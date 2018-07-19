"""
Here we collect the optimizers specifically designed for the TDC
These optimizers are to be used on lcs that have a spline already !
Typically, they are called by the runsim and runobs functions of run.py
"""

import pycs
import numpy as np
from splopt import calcknotstep


def spldiff(lcs, verbose=True, magshift=False):

	"""
	Custom spldiff optimizer for TDC

	@param lcs: list of light curves
	@param magshift: boolean, do you want to allow magnitude shifts in the pycs.spldiff.multiopt.opt_ts	function called here
	"""
	
	# We start by computing some statistics on the lcs:
	
	(lca,lcb) = lcs
	
	stats = lca.samplingstats(seasongap=30)
	sampling = stats["med"]*3.0

	if not (hasattr(lca, "vario") and hasattr(lcb, "vario")):
		lca.vario = pycs.tdc.vario.vario(lca, verbose=True)
		lcb.vario = pycs.tdc.vario.vario(lcb, verbose=True)
		print '---Vario Analysis Done---'

	knotstep = calcknotstep([lca.vario, lcb.vario])
	bokeps = np.max([sampling, knotstep/3.0])
	bokwindow = None
	
	# The stab params, quite easy :
	stabext = 300.0
	stabgap = 6.0
	stabstep = sampling
	stabmagerr = -3.0
	
	#rslc params
	pd = 10 # keep it low to avoid taking too much time running on it
	
	knots = pycs.gen.spl.seasonknots(lcs, knotstep, ingap=1)	
	if verbose:
		print "I prepared %i knots" % (len(knots))

	# Now, we call the magic function of spldiff:
	return pycs.spldiff.multiopt.opt_ts(lcs, pd=pd, knotstep=knotstep, stab=True, stabext=stabext, stabgap=stabgap, stabstep=stabstep,
					stabmagerr=stabmagerr, bokit=0, bokeps=bokeps,
					boktests=5, bokwindow=bokwindow, verbose=verbose, magshift=magshift)




def regdiff2(lcs):
	"""
	Regdiff usign scikit-learn
	Let's see if this works without any initial shift.

	@param lcs: list of light curves
	"""
	import pycs.regdiff2
	
	# We compute a theta0 parameter for the covariance function :
	for l in lcs:
		if not hasattr(l, "vario"):
			l.vario = pycs.tdc.vario.vario(l, verbose=True)
		vratio = np.clip(l.vario["vratio"], 1.0, 3.0)
		l.theta0 = 10.0**(1.3 + (vratio-0.95)*(4.4-1.3)/(2.5-0.95))

	inistep = 1.0
	radius = 100.0/inistep
	pycs.regdiff2.multiopt.opt_ts_indi(lcs, method="weights", pd=2.0, radius=radius, inistep=inistep, nit=1, verbose=True)


