"""
Here we collect the optimizers specifically designed for the TDC
These optimizers are to be used on lcs that have a spline already !
Typically, they are called by the runsim and runobs functions of run.py

"""

import pycs
import numpy as np

def spldiff(lcs, verbose=True, magshift=True):

	"""
	Custom spldiff optimizer for TDC
	here, the lcs belong the a run object and inherit its stats and knotstep
	that's why we can call them.
	
	"""
	
	# We start by computing some statistics on the lcs:
	
	(lca,lcb) = lcs
	
	stats = lca.samplingstats(seasongap=100)
	sampling = stats["med"]
	
	knotstep = lca.knotstep 
	bokeps = np.max([sampling, knotstep/3.0])
	bokwindow = None
	
	# The stab params, quite easy :
	stabext = 300.0
	stabgap = 60.0
	stabstep = sampling
	stabmagerr = -3.0
	
	#rslc params
	pd = 2 # keep it low to avoid taking too much time running on it
	
	knots = pycs.gen.spl.seasonknots(lcs, knotstep, ingap=1)	
	if verbose:
		print "I prepared %i knots" % (len(knots))
		

	# Now, we call the magic function of spldiff:
	
	
	return pycs.spldiff.multiopt.opt_ts(lcs, pd=pd, knotstep=knotstep, stab=True, stabext=stabext, stabgap=stabgap, stabstep=stabstep,
					stabmagerr=stabmagerr, bokit=0, bokeps=bokeps,
					boktests=5, bokwindow=bokwindow, verbose=verbose, magshift=magshift)




def regdiff2(lcs)
	"""
	
	"""
	import pycs.regdiff2
	
	
	
	pycs.regdiff2.multiopt.opt_ts_indi(lcs, method="weights", pd=0.25, inistep=10, nit=7, verbose=True)


