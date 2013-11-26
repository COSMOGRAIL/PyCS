"""
Here we collect some spline optimizers specifically designed for the TDC

"""

import pycs
import numpy as np
import random



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
	(lca, lcb) = lcs  # Just nomenclature : lca is the "first one", NOT necesserily A !!!
	
	stats = lca.samplingstats(seasongap=100)
	sampling = stats["med"]
	
	if not (hasattr(lca, "vario") and hasattr(lcb, "vario")):
		lca.vario = pycs.tdc.vario.vario(lca, verbose=True)
		lcb.vario = pycs.tdc.vario.vario(lcb, verbose=True)
	
	knotstep = calcknotstep(lca.vario, lcb.vario)
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

	# And now iteratively optimize the shifts
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




def spl2(lcs, maxit=7, minchange=1.0, verbose=True):
	"""
	Custom spline optimizer for TDC
	Assumes good initial time shift, but does not care about inital magshift.
	Optimizes any ML that is present.
	
	:param maxit: maximum number of iteartions
	
	:param minchange: minimum decrease in percent of the r2. I stop if the decrease gets smaller.
	
	"""
	(lca, lcb) = lcs # Just nomenclature : lca is the "first one", NOT necesserily A !!!
	
	stats = lca.samplingstats(seasongap=100)
	sampling = stats["med"]
	
	if not (hasattr(lca, "vario") and hasattr(lcb, "vario")):
		lca.vario = pycs.tdc.vario.vario(lca, verbose=True)
		lcb.vario = pycs.tdc.vario.vario(lcb, verbose=True)
	
	knotstep = calcknotstep(lca.vario, lcb.vario)
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
	
	# We fit a spline through the curve without ML (if possible)
	nomllcs = [l for l in lcs if l.ml == None]
	if len(nomllcs) == 0: # All curves have ML, we use the first one
		spline = pycs.gen.spl.fit([lca], knots=knots,
			stab=True, stabext=stabext, stabgap=stabgap, stabstep=stabstep, stabmagerr=stabmagerr,
			bokit=0, bokeps=bokeps, boktests=5, bokwindow=bokwindow, verbose=False)
	else:
		# Spline through the fixed curve :
		spline = pycs.gen.spl.fit([nomllcs[0]], knots=knots,
			stab=True, stabext=stabext, stabgap=stabgap, stabstep=stabstep, stabmagerr=stabmagerr,
			bokit=0, bokeps=bokeps, boktests=5, bokwindow=bokwindow, verbose=False)

	if verbose:
		print "Single spline fit done"
	
	#print "WARNING RETURN SINGLE SPLINE"
	#return spline
	
	# We optimize the mag shift, moving both lca and lcb.
	pycs.spl.multiopt.opt_magshift(lcs, sourcespline=spline, verbose=False, trace=False)
	
	# If present, we do a first optimization of any ML
	pycs.spl.multiopt.opt_ml(lcs, sourcespline=spline, bokit=0, splflat=True, verbose=False)
	
	# And fit a spline through both curves
	spline = pycs.gen.spl.fit(lcs, knots=knots,
		stab=True, stabext=stabext, stabgap=stabgap, stabstep=stabstep, stabmagerr=stabmagerr,
		bokit=0, bokeps=bokeps, boktests=5, bokwindow=bokwindow, verbose=False)
	if verbose:
		print "spline fit done"

	# And now iteratively optimize the shits
	print "Starting opt on initial delays :"
	print pycs.gen.lc.getnicetimedelays(lcs, separator=" | ")
	previousr2 = spline.lastr2nostab

	for it in range(maxit):
		
		#pycs.spl.multiopt.opt_ts_brute(lcs, spline, movefirst=False, optml=False, r=10, step=1.0, verbose=False)
		
		pycs.spl.multiopt.opt_ts_indi(lcs, spline, method="brute", brutestep=0.2, bruter=10, verbose = False)
		pycs.spl.multiopt.opt_ts_indi(lcs, spline, method="fmin", verbose=False)
		pycs.spl.multiopt.opt_magshift(lcs, spline, verbose=False)
		
		pycs.spl.multiopt.opt_source(lcs, spline, dpmethod="extadj", bokit = 1, verbose=False)
		pycs.spl.multiopt.opt_ml(lcs, spline, bokit=1, splflat=True, verbose=False)
		
		#print "opt_ts_brute brute done"
		#print pycs.gen.lc.getnicetimedelays(lcs, separator=" | ")
	
		pycs.spl.multiopt.opt_source(lcs, spline, dpmethod="extadj", bokit = 0, verbose=False)
		
		r2changepercent = 100.0 * (spline.lastr2nostab - previousr2) / previousr2
		print "Iteration %i done, r2 = %8.1f (%+.2f%%)" % (it+1, spline.lastr2nostab, r2changepercent)
		print pycs.gen.lc.getnicetimedelays(lcs, separator=" | ")
		previousr2 = spline.lastr2nostab
		
		if r2changepercent < 0.0 and np.fabs(r2changepercent) < minchange:
			print "I stop, minchange reached !"
			break
	
	
	print "Timeshift stabilization and releasing of splflat :"
	for it in range(10):
		pycs.spl.multiopt.opt_ts_indi(lcs, spline, optml=True, mlsplflat=False, method="fmin", verbose = False, trace=False)
		pycs.spl.multiopt.opt_source(lcs, spline, dpmethod="extadj", bokit = 0, verbose=False, trace=False)
		r2changepercent = 100.0 * (spline.lastr2nostab - previousr2) / previousr2
		print "Iteration %i done, r2 = %8.1f (%+.2f%%)" % (it+1, spline.lastr2nostab, r2changepercent)
		print pycs.gen.lc.getnicetimedelays(lcs, separator=" | ")
		previousr2 = spline.lastr2nostab
		
	
	return spline


def splml1(lcs):
	"""
	Some rather simple spline ML model.
	Randomly added on one curve of lcs
	"""
	lctoaddto = random.choice(lcs)
	mlknotstep = 500.0
	mlbokeps = 100.0
	pycs.gen.splml.addtolc(lctoaddto, knotstep=mlknotstep, bokeps=mlbokeps)



