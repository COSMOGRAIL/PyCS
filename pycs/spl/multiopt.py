"""
Functions to optimize time shifts and microlensing between lcs, using spline fits.


By default the functions don't touch the ML and sourcespline *knots*, it's up to you to enable BOK iterations.
And check the results by eye ...

Typical order of operations : 

Put smooth microlensing on the curves
Fit a source to the curve without microlensing
opt_magshift (does not use the source)
opt_ml (does not change the source)
Fit a new source, to all the curves
Put a less smooth microlensing on them
opt_ml_source


"""
import sys, os
import pycs.gen.lc
import pycs.gen.spl
import pycs.gen.util
import pycs.gen.polyml
import numpy as np
import scipy.optimize as spopt



def opt_magshift(lcs, sourcespline=None, verbose=False, trace=False):
	"""
	If you don't give any sourcespline, this is a dirty rough magshift optimization,
	using the median mag level (without microlensing), once for all.
	We don't touch the magshift of the first curve.
	
	If you do give a sourcespline, I'll optimize the magshift of each curve one by one to match the spline.
	New : I even touch the magshift of the first curve.
	
	We use the l1-norm of the residues, not a usual r2, to avoid local minima. Important !
	This is done with the nosquare=True option to the call of r2 !
	"""
	
	if sourcespline == None:
	
		reflevel = np.median(lcs[0].getmags(noml = True))
		for l in lcs[1:]:
			level = np.median(l.getmags(noml = True))
			l.shiftmag(reflevel - level)
			if trace:
				pycs.gen.util.trace(lcs)
		if verbose:
			print "Magshift optimization done."

	else:
		
		#for l in lcs[1:]: # We don't touch the first one.
		for l in lcs:
			
			if verbose:
				print "Magshift optimization on %s ..." % (l)
			inip = l.magshift
		
			def setp(p):
				l.magshift = p
			
			def errorfct(p):	
				setp(p)
				return pycs.gen.spl.r2([l], sourcespline, nosquare=True)
			
			minout = spopt.fmin(errorfct, inip, full_output=1, xtol=0.001, disp=verbose)			
			popt = minout[0]
			if verbose:
				print "Optimal magshift: %.4f" % (popt)	
			setp(popt)
			if verbose:
				print "Magshift optimization of %s done." % (l)
		# We do not return anything
		

def opt_source(lcs, sourcespline, dpmethod="extadj", bokit = 0, bokmethod="BF", verbose = True, trace=False):
	"""
	Just the source spline, without touching the ML of the lcs.
	At each call, I update the sourcespline with the merged lcs.
	The internal knots of the sourcespline stay where they are, only the external ones are ajusted.
	"""
		
	inir2 = sourcespline.r2(nostab = True)
	if verbose:
		print "Starting source optimization ..."
		print "Initial r2 (before dp update) : %f" % (inir2)
	
	
	dp = pycs.gen.spl.merge(lcs, olddp=sourcespline.datapoints)
	sourcespline.updatedp(dp, dpmethod=dpmethod)
	
	for n in range(bokit):
		sourcespline.buildbounds(verbose=verbose)
		r2 = sourcespline.bok(bokmethod=bokmethod, verbose=verbose)
	
	if bokit == 0: # otherwise this is already done by bok.
		sourcespline.optc()		
		r2 = sourcespline.r2(nostab = True) # Important, to set sourcesplie.lastr2nostab
	
	if trace:
		pycs.gen.util.trace(lcs, [sourcespline])
	if verbose:
		print "Final r2 : %f" % (r2)
	return r2


def opt_fluxshift(lcs, sourcespline, verbose = True, trace=False):
	"""
	Optimizes the flux shift and the magshift of the lcs (not the first one)
	to get the best fit to the "sourcespline". Does not touch the microlensing, nor the spline.
	So this is a building block to be used iteratively with the other optimizers.
	Especially of the sourcespline, as we fit here even on regions not well constrained by the spline !
		
	The spline should typically well fit to the first curve.
	"""

	
	for l in lcs[1:]: # We don't touch the first one.
			
		if verbose:
			print "Fluxshift optimization on %s ..." % (l)
		
		minfs = l.getminfluxshift()
		maxfs = 100000.0
		inip = l.fluxshift
		
		
		inip = (0, 0)
		# p = (fluxshift, magshift)
		def setp(p):
			fs = p[0]*1000.0
			ms = p[1]*0.1
			
			if fs < minfs:
				l.setfluxshift(minfs, consmag = False)
			else:
				l.setfluxshift(fs, consmag = False)
			
			l.magshift = ms
			
		
		def errorfct(p):	
			setp(p)
			return pycs.gen.spl.r2([l], sourcespline)
			
		minout = spopt.fmin_powell(errorfct, inip, full_output=1, xtol=0.001, disp=verbose)			
		
		popt = minout[0]	
		setp(popt)
		if verbose:
			print "Done with %s ..." % (l)
	# We do not return anything
	


def opt_ml(lcs, sourcespline, bokit = 0, bokmethod="BF", splflat = False, verbose = True, trace=False):
	"""
	Optimizes the microlensing of the lcs (one curve after the other) so that they fit to the spline.
	I work with both polynomial and spline microlensing.
	For spline micorlensing, I can do BOK iterations to move the knots.
	
	.. note:: Does not touch the sourcespline  at all !
	
	But this it what makes the problem linear (except for the BOK iterations) for both splines and polynomial ML, and
	thus fast !
	
	Parameters for spline ML :
	
	:param bokit:
	:param bokeps:
	:param boktests:
	:param bokwindow:
	:param splflat:
	
	Parameters for poly ML :
	
	None ! We just to a linear weighted least squares on each season !
	So for poly ML, all the params above are not used at all.
	
	
	We do not return anything. Returning a r2 would make no sense, as we do not touch the sourcepline !
	
	"""
	
	if trace:
		pycs.gen.util.trace(lcs, [sourcespline])
		
	if verbose:
		print "Starting ML optimization ..."
	
	for l in lcs:
		if (l.ml != None) and (l.ml.mltype == "spline"):
			# So this is spline microlensing
			
			if verbose :
				print "Working on the spline ML of %s" % (l)
			l.ml.settargetmags(l, sourcespline)
			#l.ml.spline.maxnspace(n=1)
	
			for n in range(bokit):
				l.ml.spline.buildbounds(verbose=verbose)
				l.ml.spline.bok(bokmethod=bokmethod, verbose=verbose)
				#l.ml.spline.maxnspace(n=1)
			
			if splflat:
				l.ml.spline.optc()
				l.ml.spline.optcflat(verbose=False)
			else:
				l.ml.spline.optc()
			if trace:
				pycs.gen.util.trace(lcs, [sourcespline])


		if (l.ml != None) and (l.ml.mltype == "poly"):
			
			if verbose :
				print "Working on the poly ML of %s" % (l)
			
			# We go through the curve season by season :
			for m in l.ml.mllist:
				#print m.getparams(), m.season, m.season.indices
				
				nparams = m.nfree
				
				mlseasjds = l.jds[m.season.indices]
				mlseasjds -= np.mean(mlseasjds) # Convention for polyml, jds are "centered".
				nomlmags = l.getmags(noml = True)[m.season.indices]
				magerrs = l.magerrs[m.season.indices]
	
				absjds = l.getjds()[m.season.indices]
				targetmags = sourcespline.eval(absjds)
	
				polyparams = pycs.gen.polyml.polyfit(mlseasjds, targetmags - nomlmags, magerrs, nparams)

				m.setparams(polyparams)			

	if verbose:
		print "Done !"
	


def redistribflux(lc1, lc2, sourcespline, verbose=True, maxfrac = 0.2):
	"""
	Redistributes flux between lc1 and lc2 (assuming these curves suffer form flux sharing), so
	to minimize the r2 with respect to the sourcespline.
	I do not touch the sourcespline, but I do modify your curves in an irreversible way !
	
	:param lc1: a lightcurve
	:param lc2: another lightcurve
	:param sourcespline: the spline that the curves should try to fit to
	
	"""
	if not np.all(lc1.jds == lc2.jds):
		raise RuntimeError("I do only work on curves with identical jds !")
	
	if verbose:
		print "Starting redistrib_flux, r2 = %10.2f" % (pycs.gen.spl.r2([lc1, lc2], sourcespline))
	
	
	# The initial curves :
	lc1fluxes = lc1.getrawfluxes()
	lc2fluxes = lc2.getrawfluxes()
	
	maxamp = min(np.min(lc1fluxes), np.min(lc2fluxes)) # maximum amplitute of correction
	
	def setp(flux, i): # flux is an absolute shift in flux for point i
		lc1.mags[i] = -2.5*np.log10(lc1fluxes[i] + flux)
		lc2.mags[i] = -2.5*np.log10(lc2fluxes[i] - flux)
		
	def errorfct(flux, i):	
		setp(flux, i)
		return pycs.gen.spl.r2([lc1, lc2], sourcespline)
			
	# We can do this one point at a time ...
	#fluxes = []
	for i in range(len(lc1)):
		#print i
		out = spopt.optimize.fminbound(errorfct, -maxfrac*maxamp, maxfrac*maxamp, args=(i,), xtol=0.1, disp=True, full_output=False)			
		setp(out, i)
		#fluxes.append(out[0])
	
	#lc1.mags = -2.5*np.log10(lc1fluxes[i] + flux)
	#lc2.mags = -2.5*np.log10(lc1fluxes[i] - flux)
		
	if verbose:
		print "Done with redistrib,     r2 = %10.2f" % (pycs.gen.spl.r2([lc1, lc2], sourcespline))
	


			
# The old non-linear way to fit the poly ml :
# 		if (l.ml != None) and (l.ml.mltype in ["poly", "leg"]):
# 			
# 			inipars = l.ml.getfreeparams()
# 					
# 			def errorfct(p):
# 				# We set the ML params
# 				l.ml.setfreeparams(p)
# 				# We return a chi2 between the fixed spline, and the lightcurve :
# 				return pycs.gen.spl.r2([l], sourcespline)
# 				
# 			if verbose:
# 				print "Starting Poly ML optimization of %s ..." % (l)
# 				print "Initial pars : ", inipars
# 					
# 			minout = spopt.fmin_powell(errorfct, inipars, full_output=1, disp=verbose, maxiter=maxpowellit)
# 			popt = minout[0]
# 			if popt.shape == ():
# 				popt = np.array([popt])
# 		
# 			if verbose:
# 				print "Optimal pars : ", popt
# 			l.ml.setfreeparams(popt)
# 			
		
		


# def opt_ml_flux(lcs, sourcespline, maxit = 40, verbose=True):
# 	"""
# 	For both optimizations, we do not touch the sourcespline. Also I do not touch the knots of the ML.
# 	So this should be stable, we can iterate a lot on it.
# 	Best results when using a not very flexible ML.
# 	"""
# 	
# 	print "NOT A GOOD IDEA, THE SOURCESPLINE SHOULD UPDATE ONCE THE FLUX HAS CHANCED"
# 	
# 	it = 0
# 	while True:
# 		it +=1
# 		
# 		opt_ml(lcs, sourcespline, bokit = 0, verbose = False)
# 		opt_flux(lcs, sourcespline, verbose = False)
# 		
# 		opt_source(lcs, sourcespline, bokit = 0, verbose = False, trace=False)
# 		
# 		for l in lcs:
# 			print l
# 		
# 		if verbose:
# 			r2 = pycs.gen.spl.r2(lcs, sourcespline)
# 			print "ML + Flux optimization, iteration %i done, r2 = %.3f" % (it, r2)
# 		
# 		if it >= maxit:
# 			break
# 	
	

# def opt_ml_source(lcs, sourcespline, optflux=True, maxit = 5, stoprel=-10.0, verbose = True, trace=False):
# 	"""
# 	Iteratively optimizes the source spline coeffs and the ML, we don't touch any knots.
# 	This is something relatively unstable and potentially degenerate, especially in regions
# 	where all points have ML and if the ML has many knots, that's why we never do many iterations.
# 	
# 	I start with opt_ml, so give me a good sourcespline. It's ok to give me very wrong ML.
# 	
# 	"""
# 	
# 	inir2 = sourcespline.r2(nostab = True)
# 	if verbose:
# 		print "Starting ML+source optimization ..."
# 		print "Initial r2 : %f" % (inir2)
# 	it = 1
# 	lastr2 = inir2
# 	while True:
# 			
# 		# The microlensing...
# 		opt_ml(lcs, sourcespline, bokit = 0, verbose=False, trace=False)
# 		
# 		# The flux
# 		#if optflux:
# 		#	opt_flux(lcs, sourcespline, verbose=False, trace=False)
# 			
# 		# The sourcespline...
# 		r2 = opt_source(lcs, sourcespline, bokit = 0, verbose = False, trace=trace)
# 	
# 		# End of loop stuff :
# 		r2 = sourcespline.r2(nostab = True)
# 	
# 		if verbose:
# 			print "Iteration %02i : r2 = %f" % (it, r2)
# 
# 
# 		relimp = (lastr2 - r2)/lastr2
# 		if verbose:
# 			print "Total relative improvement : %f" % (relimp)
# 		if relimp < stoprel or it >= maxit :
# 			if verbose:
# 				print "OK, I'm done !"
# 			return r2
# 		it += 1
# 		lastr2 = r2
# 	
	



def opt_ts_powell(lcs, sourcespline, optml=True, movefirst=False, verbose = True, trace=False):
	"""
	
	If you use this again, implement mlsplflat for optml, might be important !
	
	Optimize the timeshifts of all four lightcurves AND the spline coeffs (but not the knots) so to get the best possible fit.
	
	If optml is True, I will optimize the ml for every trial delays (and then optimize the source again).
	For this I will always start from the initial settings, I don't cumulate the ML - source optimizations.
	This is required as ML - source optimizations are degenarate.
	If optml is false, this would not be needed, as the spline fitting is unique anyway.
	
	Note that even if I don't modify the sourcespline knots, I do modify the sourcespline datapoints as the time shifts move !
	
	Improvement ideas :
	see if better to not optimize the first curves time shift ?
	Don't try to optimize the timeshifts at the same time ?
	Do one after the other, with respect to the first curve ?
	Special function that takes a spline fit to the first curve as argument, to do this ?
	
	-> Yes, see opt_ts_indi below !
	
	"""
	# To start each time shift guess from the same conditions, we will internally work on copies.
	origlcs = [l.copy() for l in lcs]
	origsourcespline = sourcespline.copy()
	
		
	def errorfct(p): # p are the absolute time shifts of the n or n-1 curves
		"""
		Warning : I work only on interal copies, and do not set anything on lcs or sourcespline
		If you change me, also change the apply() below !
		"""
	
		# We make internal copies :
		mylcs = [l.copy() for l in origlcs]
		mysourcespline = origsourcespline.copy()
		
		# We set p of these internal copies : 
		# (We do an absolute update of timeshift)
		if p.shape == ():
			p = np.array([p])
		if movefirst: # this should be changed to gen.lc.gettimeshifts ...
			for (l, shift) in zip(mylcs, p):
				l.timeshift = shift
		else:
			for (l, shift) in zip(mylcs[1:], p):
				l.timeshift = shift
		# We do a first source opt :
		r2 = opt_source(mylcs, mysourcespline, verbose=False, trace=trace)
		# And then ML, and source once again :
		if optml:
			opt_ml(mylcs, mysourcespline, bokit = 0, verbose = False)
			r2 = opt_source(mylcs, mysourcespline, verbose=False, trace=False)
			
		return r2
	
	def apply(p):
		"""
		Does the same as errorfct, but on the real objects, not copies !
		"""
		if p.shape == ():
			p = np.array([p])
		if movefirst:
			for (l, shift) in zip(lcs, p):
				l.timeshift = shift
		else:
			for (l, shift) in zip(lcs[1:], p):
				l.timeshift = shift
	
		r2 = opt_source(lcs, sourcespline, verbose=False, trace=trace)
		if optml:
			opt_ml(lcs, sourcespline, bokit = 0, verbose = False)
			r2 = opt_source(lcs, sourcespline, verbose=False, trace=False)
		# These calls to opt_source will set the sourcepline.lastr2nostab !
		return r2
		
	

	if movefirst:
		inip = np.array([l.timeshift for l in lcs])
	else:
		inip = np.array([l.timeshift for l in lcs[1:]])
	
	inir2 = errorfct(inip)
	
	if verbose:
		print "Starting time shift optimization ..."
		print "Initial pars : ", inip
		print "Initial r2 : %f" % (inir2)
		
	minout = spopt.fmin_powell(errorfct, inip, full_output=1, xtol=0.001, disp=verbose)
	popt = minout[0]
	
	r2 = apply(popt) # This sets popt, and the optimal ML and source.
	
	if verbose:
		print "Final r2 : %f" % (r2)
		print "Optimal pars : ", popt
	
	return r2


def comb(*sequences):
	'''
	http://code.activestate.com/recipes/502199/
	combinations of multiple sequences so you don't have
	to write nested for loops
    
	>>> from pprint import pprint as pp
	>>> pp(comb(['Guido','Larry'], ['knows','loves'], ['Phyton','Purl']))
	[['Guido', 'knows', 'Phyton'],
	['Guido', 'knows', 'Purl'],
	['Guido', 'loves', 'Phyton'],
	['Guido', 'loves', 'Purl'],
	['Larry', 'knows', 'Phyton'],
	['Larry', 'knows', 'Purl'],
	['Larry', 'loves', 'Phyton'],
	['Larry', 'loves', 'Purl']]
	>>> 
	'''
	combinations = [[seq] for seq in sequences[0]]
	for seq in sequences[1:]:
		combinations = [comb+[item]
			for comb in combinations
			for item in seq ]
	return combinations

def opt_ts_brute(lcs, sourcespline, movefirst = True, optml=False, r=2, step=1.0, verbose = True, trace=False):
	"""
	
	If you use this again, implement mlsplflat, might be important.
	
	Given the current delays, I will explore a hypercube (r steps in each direction on each axis)
	of possible time shift combinations. And choose the shifts that gave the smallest chi2.

	For each trial shift, I optimize the sourcespline to fit all curves, and optionally also the ML of every curve.
	
	This is very slow, not really used.
	
	"""
	origlcs = [l.copy() for l in lcs]
	origsourcespline = sourcespline.copy()
	
	if movefirst:
		origshifts = np.array([l.timeshift for l in lcs])
	else:
		origshifts = np.array([l.timeshift for l in lcs[1:]])
	
	relrange = np.linspace(-1.0*r*step, 1.0*r*step, 2*r + 1)
	absrangelist = [list(relrange + origshift) for origshift in origshifts]
	
	# We want a list of combinations to explore. Its length is (2*r + 1)**len(lcs)
	combilist = comb(*absrangelist)
	if movefirst:
		assert len(combilist) == (2*r + 1)**len(lcs)
	else:
		assert len(combilist) == (2*r + 1)**(len(lcs)-1)
	

	combir2s = np.zeros(len(combilist))
	for i, absshifts in enumerate(combilist):
		if verbose:
			print "%i / %i" % (i+1 ,len(combilist))
		
		# We make internal copies :
		mylcs = [l.copy() for l in origlcs]
		mysourcespline = origsourcespline.copy()
		#mylcs = lcs
		#mysourcespline = sourcespline
		
		# We set the shifts of these internal copies : 
		
		if movefirst:
			for (l, shift) in zip(mylcs, absshifts):
				l.timeshift = shift
		else:
			for (l, shift) in zip(mylcs[1:], absshifts):
				l.timeshift = shift
		# We do a first source opt :
		r2 = opt_source(mylcs, mysourcespline, verbose=False, trace=trace)
		# And then ML, and source once again :
		if optml:
			opt_ml(mylcs, mysourcespline, bokit = 0, verbose = False)
			r2 = opt_source(mylcs, mysourcespline, verbose=False, trace=False)
		
		#print r2
		# And save the r2 : 
		combir2s[i] = r2
	
	# Find the min
	minindex = np.argmin(combir2s)
	optabshifts = combilist[minindex]
	if verbose:
		print "Best shift at index %i :" % (minindex)
		print optabshifts
	# And set the actual curves :
	if movefirst:
		for (l, shift) in zip(lcs, optabshifts):
			l.timeshift = shift
	else:
		for (l, shift) in zip(lcs[1:], optabshifts):
			l.timeshift = shift
	r2 = opt_source(lcs, sourcespline, verbose=False, trace=trace)
	if optml:
		opt_ml(lcs, sourcespline, bokit = 0, verbose = False)
		r2 = opt_source(lcs, sourcespline, verbose=False, trace=False)
	
	return r2
	


def opt_ts_indi(lcs, sourcespline, method="fmin", crit="r2", optml=False, mlsplflat=False, brutestep=1.0, bruter=5, verbose = True, trace=False):
	"""
	We shift the curves one by one so that they match to the spline, using fmin or brute force.
	A bit special : I do not touch the spline at all ! Hence I return no r2.
	Idea is to have a fast ts optimization building block.
	
	The spline should be a shape common to the joined lcs.
	No need to work on copies, as we do not change the ML or spline *iteratively*, but
	only the ML -> nothing can go wrong.
	
	:parma brutestep: step size, in days
	:param bruter: radius in number of steps
	
	"""
	
	for l in lcs:
		
		def errorfct(timeshift):
			# Set the shift :
			l.timeshift = timeshift
			# Optimize the ML for this shift :
			if optml:
				opt_ml([l], sourcespline, bokit = 0, splflat=mlsplflat, verbose = False)
			# Calculate r2 without touching the spline :
			if crit == "r2":
				error = pycs.gen.spl.r2([l], sourcespline)
				#print "Still using r2 !"
			elif crit == "tv":
				error = pycs.gen.spl.mltv([l], sourcespline)[0]
				print "Warning, using TV !"
			if verbose:
				print "%s %10.3f %10.3f" % (l.object, l.timeshift, error) 
			return error
		
		initimeshift = l.timeshift
		
		if method == "fmin":
			out = spopt.fmin(errorfct, initimeshift, xtol=0.1, ftol=0.1, maxiter=None, maxfun=100, full_output=1, disp=verbose)
			opttimeshift = float(out[0][0])
		elif method == "brute":
			
			testvals = np.linspace(-1.0*bruter*brutestep, 1.0*bruter*brutestep, 2*bruter + 1) + initimeshift
			r2vals = np.array(map(errorfct, testvals))
			minindex = np.argmin(r2vals)
			opttimeshift = float(testvals[minindex])
		
		l.timeshift = opttimeshift





	

# Version 1...
# def opt_flux(lcs, sourcespline, maxit=2, verbose = True, trace=False):
# 	"""
# 	Optimizes the flux shift, ML, and source spline.
# 	We take the first curve as a reference for the fluxshift, and tweak the fluxes of the other ones.
# 	For each trial flux shift, we optimize the ML and the sourcespline, without touchning the knots.
# 	"""
# 
# 	inir2 = sourcespline.r2(nostab = True)
# 	if verbose:
# 		print "Starting Flux+ML+source optimization ..."
# 		print "Initial r2 : %f" % (inir2)
# 	
# 	it = 1
# 	while True:
# 		for l in lcs[1:]: # We don't touch the first one.
# 			
# 			if verbose:
# 				print "Fluxshift optimization on %s ..." % (l)
# 		
# 			minfs = l.getminfluxshift()
# 			maxfs = 10000.0
# 			inifs = l.fluxshift
# 		
# 			def setfs(fs):
# 				l.setfluxshift(fs, consmag = True)
# 			
# 			def errorfct(fs):	
# 				setfs(fs)
# 				#print fs
# 				return opt_ml_source(lcs, sourcespline, maxit = 2, verbose = False, trace=trace)
# 				
# 			#def errorfctfmin(fs):	
# 			#	setfs(fs)
# 			#	print fs
# 			#	if fs < minfs or fs > maxfs:
# 			#		return 100.0 * len(l)
# 			#	return opt_ml_source(lcs, sourcespline, maxit = 2, verbose = False)
# 				
# 		
# 			#if it == 1:
# 			minout = spopt.fminbound(errorfct, minfs, maxfs, xtol=1.0, maxfun=50, full_output=1, disp=1)
# 				#inifs = minout[0]
# 				#minout = spopt.fmin(errorfctfmin, inifs, xtol=0.5, maxfun=50, full_output=1, disp=1)
# 			#else:
# 			#	minout = spopt.fmin(errorfctfmin, inifs, xtol=0.5, maxfun=50, full_output=1, disp=1)
# 			
# 			#print minout
# 			popt = minout[0]
# 			if verbose:
# 				print "Iteration %02i : r2 = %f" % (minout[3], minout[1])
# 				print "Optimal pars : ", popt
# 			setfs(popt)
# 		
# 		
# 		r2 = sourcespline.r2(nostab = True)
# 		if verbose:
# 			print "Iteration %i : r2 = %f" % (it, r2)
# 		
# 		if it >= maxit:
# 			if verbose:
# 				print "Done !"
# 			return r2
# 			
# 		it += 1




