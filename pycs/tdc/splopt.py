"""
Here we collect some spline optimizers specifically designed for the TDC

"""

import pycs
import numpy as np


def calcknotstep(varios):
	"""
	Give me some outputs of vario, I try to return a good knotstep, based on the highest vratio between the curves and the sampling.
	"""
	
	vratios = np.array([out["vratio"] for out in varios])
	samplings = np.array([out["sampling"] for out in varios])
	
	# First playing around :
	"""
	# Max and min vratio :
	va = 1.0
	vb = 3.0
	# We limit the range of these values to meaningful stuff
	vratios = np.clip(vratios, va, vb)
	v = np.max(vratios)
	# Knotstep range corresponding to va and vb:
	ksa = 50.0
	ksb = 4.0#*np.min(sampling)
	# And compute the knotstep, linear scaling :	
	ks = ksa + (ksb - ksa)*(v - va)/(vb - va)
	"""
	
	# Second more serious attempt, calibrated on a bunch of curves :
	vratio = np.max(vratios)
	sampling = np.min(samplings)
	vratio = np.clip(vratio, 1.2, 5) # OK to lower this min value, but do not allow 1.0. OK to increase max value.
	ks = 14.0/(vratio - 1.0)
	ks = np.clip(ks, 2.0*sampling, 100.0)
	
	return float(ks)
	


def spl1(lcs, verbose=True):
	"""
	Custom spline optimizer for TDC
	Assumes reasonable initial time shift, but no magshift.
	
	"""
	(lca, lcb) = lcs
	
	stats = lca.samplingstats(seasongap=100)
	sampling = stats["med"]
	
	knotstep = lca.knotstep # Is the same for a and b...
	bokeps = np.max([sampling, knotstep/3.0])
	"""
	if len(lca) > 1000.0:
		knotstep = 20.
		bokeps = 10.
	else:
		knotstep = 3.0*sampling
		bokeps = knotstep/4.0
	"""
	bokwindow = None
	
	# The stab params, quite easy :
	stabext = 300.0
	stabgap = 60.0
	stabstep = sampling
	stabmagerr = -3.0
	
	knots = pycs.gen.spl.seasonknots(lcs, knotstep, ingap=1)
	#exit()
	if verbose:
		print "I prepared %i knots" % (len(knots))
	
	# We fit a spline through lca
	splinea = pycs.gen.spl.fit([lca], knots=knots,
		stab=True, stabext=stabext, stabgap=stabgap, stabstep=stabstep, stabmagerr=stabmagerr,
		bokit=0, bokeps=bokeps, boktests=5, bokwindow=bokwindow, verbose=False)
	if verbose:
		print "splinea fit done"
	
	#print "WARNING RETURN SPLINE A"
	#return splinea
	
	# We optimize the mag shift
	pycs.spl.multiopt.opt_magshift(lcs, sourcespline=splinea, verbose=False, trace=False)
	
	# And fit a spline through both curves
	spline = pycs.gen.spl.fit(lcs, knots=knots,
		stab=True, stabext=stabext, stabgap=stabgap, stabstep=stabstep, stabmagerr=stabmagerr,
		bokit=0, bokeps=bokeps, boktests=5, bokwindow=bokwindow, verbose=False)
	if verbose:
		print "spline fit done"

	# And now iteratively optimize the shits
	print "Starting opt on initial delays :"
	print pycs.gen.lc.getnicetimedelays(lcs, separator=" | ")
	for it in range(7):
		
		#pycs.spl.multiopt.opt_ts_brute(lcs, spline, movefirst=False, optml=False, r=10, step=1.0, verbose=False)
		
		pycs.spl.multiopt.opt_ts_indi(lcs, spline, method="fmin", verbose=False)
		pycs.spl.multiopt.opt_magshift(lcs, spline, verbose=False)
		
		pycs.spl.multiopt.opt_source(lcs, spline, dpmethod="extadj", bokit = 1, verbose=False)
		#print "opt_ts_brute brute done"
		#print pycs.gen.lc.getnicetimedelays(lcs, separator=" | ")
	
		print "Iteration %i done." % (it+1)
		print pycs.gen.lc.getnicetimedelays(lcs, separator=" | ")
		
	
	return spline
		
