"""
Functions to optimize timeshifts, microlensing, or all this at the same time between an arbitrary set of curves.
These are building blocks, to be assembled to make a general optimizer.

"""
import sys, os

import numpy as np
import scipy.optimize as spopt
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


import matplotlib.cm as cm
import matplotlib.colors as colors

import pycs.gen.polyml
import pycs.gen.lc
import pycs.gen.util




def opt_magshift(lcs, rawdispersionmethod, verbose=True, trace=False):
	"""
	I optimize the magnitude shifts between lightcurves. This is usually a first thing to do.
	First curve does not get shifted.
	
	You might think that this can be done by just using the median mag of a curve, but here we use the dispersion method,
	to get it right even in case of curves with non-overlapping regions.
	"""

	couplelist = [couple for couple in [[lc1, lc2] for lc1 in lcs for lc2 in lcs] if couple[0] != couple[1]]
	
	def setshifts(p): # p are the absolute magshifts.
		for (l, shift) in zip(lcs[1:], p):
			l.magshift = shift
			# We don't use the relative shiftmag() !
		if trace:
			pycs.gen.util.trace(lcs)
	
	inipars = np.array([l.magshift for l in lcs[1:]])
	if verbose:
		print "Starting magshift optimization ..."
		print "Initial pars : ", inipars

	def errorfct(p):
		setshifts(p)
		d2values = np.array([rawdispersionmethod(*couple)["d2"] for couple in couplelist])
		return np.mean(d2values)

	minout = spopt.fmin_powell(errorfct, inipars, full_output=1, disp=verbose)
	popt = minout[0]
	if popt.shape == ():
		popt = np.array([popt])
		
	if verbose:
		print "Optimal pars : ", popt
	setshifts(popt)
	# We do not return anything





def opt_ml(lcs, rawdispersionmethod, verbose=True, maxit = 3, stoprel=0.01, maxpowellit = 5, splflat=True, trace=False):
	"""
	Optimize the ml params of the lcs using the Powell method.
	Generic, give me as many lcs as you want, with either polynomial or spline microlensings.
	
	Fluxshift optimization could perhaps be included in the loop later, for now we skip this.
	
	
	If some of the lightcurves have no microlensing, I will "aim" at them for my first iteration.
	So the first iteration is special.
	"""
	
	fullcouplelist = [couple for couple in [[lc1, lc2] for lc1 in lcs for lc2 in lcs] if couple[0] != couple[1]]
	
	lastd2 = np.mean(np.array([rawdispersionmethod(*couple)["d2"] for couple in fullcouplelist]))
	loopn = 0
	
	if verbose:
		print "Starting ML optimization"
	while True:
		loopn+=1
		
		if verbose:
			print "=== Loop %i ===" % (loopn)
		for l in lcs:
			if l.ml != None:
				
				if verbose:
					print "Running on %s" % (l)
				# Here we will define a couplelist = a subset of the fullcouplelist, skipping those couples
				# not involving l. This will make it faster.
				couplelist = [couple for couple in fullcouplelist if l in couple]
				
				#for couple in couplelist:
				#	print "(%s, %s)" % (couple[0], couple[1])
				
				# Now the first iteration special case :
				if loopn == 1:
					nomllcs = [x for x in lcs if x.ml == None] # be careful not to use the variable l ...
					if len(nomllcs) == 0:
						if verbose:
							print "All curves have mircolensing, thus no special processing here."
					else:
						couplelist = [couple for couple in couplelist if couple[0] in nomllcs or couple[1] in nomllcs]
						if verbose:
							print "Aiming at curves without microlensing :"
							
							print ", ".join(["(%s, %s)" % (couple[0].object, couple[1].object) for couple in couplelist])
						
				
				#if lcs[0].ml == None and loopn == 1:
				#	# First iteration, we aim at fitting lcs[0]
				#	print "For this first loop I'll aim at fitting %s" % (lcs[0])
				#	couplelist = [[l, lcs[0]], [lcs[0], l]]
				
				#print couplelist
				#couplelist = fullcouplelist
			
				if l.ml.mltype == "spline":
				
					if splflat:
						inipars = l.ml.spline.getc(m=1)
					else:
						inipars = l.ml.spline.getc(m=0)
					
					def setp(p):
						if splflat:
							l.ml.spline.setcflat(p)
						else:
							l.ml.spline.setc(p, m=0)
						if trace:
							pycs.gen.util.trace(lcs)
							
					
					def errorfct(p):
						setp(p)
						d2values = np.array([rawdispersionmethod(*couple)["d2"] for couple in couplelist])
						return np.mean(d2values)
					
					if verbose:
						print "Starting Spline ML optimization of %s ..." % (l)
						print "Initial pars : ", inipars
					#if loopn < 3: ftol = 0.1
					#else: ftol = 0.0001
					#minout = spopt.fmin_powell(errorfct, inipars, full_output=1, disp=verbose, ftol = ftol)
					minout = spopt.fmin_powell(errorfct, inipars, full_output=1, disp=verbose, maxiter=maxpowellit)
					popt = minout[0]
					if popt.shape == ():
						popt = np.array([popt])
		
					if verbose:
						print "Optimal pars : ", popt
						
					setp(popt)
					
				if l.ml.mltype in ["poly", "leg"]:
					
					inipars = l.ml.getfreeparams()
					
					def errorfct(p):
						l.ml.setfreeparams(p)
						if trace:
							pycs.gen.util.trace(lcs)
						d2values = np.array([rawdispersionmethod(*couple)["d2"] for couple in couplelist])
						return np.mean(d2values)
					
					if verbose:
						print "Starting Poly ML optimization of %s ..." % (l)
						print "Initial pars : ", inipars
					
					minout = spopt.fmin_powell(errorfct, inipars, full_output=1, disp=verbose, maxiter=maxpowellit)
					popt = minout[0]
					if popt.shape == ():
						popt = np.array([popt])
		
					if verbose:
						print "Optimal pars : ", popt
					l.ml.setfreeparams(popt)
					
		
		d2 = np.mean(np.array([rawdispersionmethod(*couple)["d2"] for couple in fullcouplelist]))
		if verbose :
			print "Loop %i finished, dispersion = %f" % (loopn, d2)
		#if (lastd2 - d2)/lastd2 < stoprel and loopn >= 3 :
		relimp = (lastd2 - d2)/lastd2
		if verbose:
			print "Relative improvement : %f" % (relimp)
		if relimp < stoprel or loopn == maxit :
			if verbose:
				print "OK, I'm done !"
			return d2
		lastd2 = d2

	
def opt_ts_mix(lcs, rawdispersionmethod, movefirst=False, verbose=True):
	"""
	Optimize the time shifts (nothing else) with a subtle mix of optimization methods.
	
	Steps:
	* Run a simulated annealing, with a lower temperature bound and no other stop condition, i.e. a more or less
	constant number of iterations.
	* Starting from the best point found, explore a 5 x 5 cube
	* idem, but smaller
	* idem, but smaller
	* idem, but smaller
	* from this last point, simplex minimization
	
	Expect about 600 calls of the rawdispersion on the couples, i.e. 600*12 = 7200 for 4 lightcurves.
	
	.. todo:: Implement trace !
	
	"""
	
	# As we do not touch the microlensing here, we will work on copies with applied ml.
	# At the end, we just set the timeshifts of your initial curves.
	
	lcsc = [l.copy() for l in lcs]
	for l in lcsc:
		l.applyfluxshift()
		if l.ml != None:
			l.applyml() # This will make the rest a lot faster, as we avoid evaluating the polynom every time !
			
	couplelist = [couple for couple in [[lc1, lc2] for lc1 in lcsc for lc2 in lcsc] if couple[0] != couple[1]]
	
	
	#if filename != None:
	#	paramsteps = []
	#	d2steps = []

	def errorfct(timeshifts):
		if timeshifts.shape == ():
			timeshifts = np.array([timeshifts])
		pycs.gen.lc.settimeshifts(lcsc, timeshifts, includefirst = movefirst)
		d2values = np.array([rawdispersionmethod(*couple)["d2"] for couple in couplelist])
		ret = np.mean(d2values)
		#print ret
		#if filename != None:
		#	paramsteps.append(delays)
		#	d2steps.append(ret)
			
		return ret

	initshifts = pycs.gen.lc.gettimeshifts(lcsc, includefirst = movefirst)
	if verbose: 
		print "Initial time shifts : ", initshifts

	# Raw work with annealing :
	
	# Anneal no longer exists in scipy
	#res = spopt.anneal(errorfct, initshifts, schedule='fast', lower=-20, upper=20, dwell=10, learn_rate=0.2, T0 = 1.0e-5, Tf = 1.0e-7, feps=0.0, boltzmann = 0.1,full_output = 1)
	# Replacing it with 
	res = spopt.basinhopping(errorfct, initshifts, niter=100, T=1.0, stepsize=10.0, interval=50, disp=False, niter_success=10)
	# It could be that we want less iterations to be faster. For now I keep it like this, has to be investigated with a special case.
	
	ann_res = res.x # The optimal parameters
	
	if verbose:
		print "Rough optimization timeshifts : %s" % res.x
		print "Dispersion value : %7.3f" % res.fun
		#print "Final temperature : %.3g" % res[2]
		#print "Dispersion calls on couples : %i" % res[3]
		#if res[6] == 0:
		#	print "Cooled to global minimum."
		#elif res[6] == 1:
		#	print "Cooled to minimum temperature."
		#else:	
		#	print "### WARNING : %i ###" % res[6]
	

	# Medium work with iterative brute force :
		
	res = spopt.brute(errorfct, bruteranges(3,2,ann_res), full_output = 0, finish=None)
	# This would finish by default with fmin ... we do not want that.
	if verbose:
		print "Brute 1 shifts : %s" % res
		print "Brute 1 dispersion : %7.3f" % errorfct(res)
	
	res = spopt.brute(errorfct, bruteranges(1,2,res), full_output = 0, finish=None)
	if verbose:
		print "Brute 2 shifts : %s" % res
		print "Brute 2 dispersion : %7.3f" % errorfct(res)
	
	
	res = spopt.brute(errorfct, bruteranges(0.4,2,res), full_output = 0, finish=None)
	if verbose:
		print "Brute 3 shifts : %s" % res
		print "Brute 3 dispersion : %7.3f" % errorfct(res)

	res = spopt.brute(errorfct, bruteranges(0.1,2,res), full_output = 0, finish=spopt.fmin)
	# This time we let him finish with fmin.
	
	optd2 = errorfct(res)
	
	if verbose:
		print "Final shifts : %s" % res
		print "Final dispersion : %7.3f" % errorfct(res)
	
	
	#if filename != None:
	#	paramsteps = np.column_stack(paramsteps)
	#	pycs.gen.util.writepickle({"params":paramsteps, "d2": d2steps, "optparams":ann_res, "initparams":initparams}, filename)

	# Until now we worked on copies -- do not forget to shift the actual lightcurves now !
	pycs.gen.lc.settimeshifts(lcs, res, includefirst = movefirst)
	
	return optd2
	
def bruteranges(step, radius, center):
	"""
	Auxiliary function for brute force exploration.
	Prepares the "ranges" parameter to be passed to brute force optimizer
	In other words, we draw a cube ...
	radius is an int saying how many steps to go left and right of center.
	center is a list or array of the centers, it can be of any lenght.
	
	You make 2*radius + 1 steps in each direction !, so radius=2 means 5 steps thus 125 calls for 4 curves.
	"""
	
	low = - step * radius
	up = step * (radius+1)
	
	if center.shape == ():
		c = float(center)
		return [((c+low),(c+up),step)]
	else:
		return [((c+low),(c+up),step) for c in center]
	
	

		




# def opt_fluxshift(lcs, rawdispersionmethod, verbose=True):
# 	"""
# 	Optimizes the flux shift and the magshift of all your lcs.
# 	This time I'll keep it very simple : I do not touch the microlensing at all, and I don't iterate on the lcs.
# 	So this is a building block to be used iteratively with the other optimizers.
# 	Apply this function when you have already optimized a simplistic ML.
# 	"""
# 	
# 	couplelist = [couple for couple in [[lc1, lc2] for lc1 in lcs for lc2 in lcs] if couple[0] != couple[1]]
# 	
# 	for l in lcs:
# 			
# 		if verbose:
# 			print "Fluxshift optimization on %s ..." % (l)
# 		
# 		minfs = l.getminfluxshift()
# 		maxfs = 100000.0
# 		
# 		inip = (0, 0)
# 		# p = (magshift*10, fluxshift/1000) : magshift is first, for the Powell iterations.
# 		def setp(p):
# 			ms = p[0]*0.1
# 			fs = p[1]*1000.0
# 			
# 			l.magshift = ms
# 			
# 			if fs < minfs:
# 				l.setfluxshift(minfs, consmag = False)
# 			else:
# 				l.setfluxshift(fs, consmag = False)
# 	
# 		def errorfct(p):	
# 			setp(p)
# 			d2values = np.array([rawdispersionmethod(*couple)["d2"] for couple in couplelist])
# 			return np.mean(d2values)
# 
# 		minout = spopt.fmin_powell(errorfct, inip, full_output=1, xtol=0.001, disp=verbose)			
# 		
# 		popt = minout[0]	
# 		setp(popt)
# 		if verbose:
# 			print "Done with %s ..." % (l)
# 	# We do not return anything
# 


# def opt_fluxshift(lcs, rawdispersionmethod, verbose=False):
# 	"""
# 	I optimize the flux and magnitude shifts between lightcurves.
# 	I don't touch the first curve of lcs -- I need some kind of reference, otherwise I would put all your points to 0 !
# 	If the magshift is very wrong, please optimize it before calling me !
# 	Also, lightcurves should have a well fitted low order microlensing.
# 	If the microlensing is too flexible, it might already have taken some fluxshift away, so be careful...
# 	"""
# 
# 	couplelist = [couple for couple in [[lc1, lc2] for lc1 in lcs for lc2 in lcs] if couple[0] != couple[1]]
# 
# 	# First half of p are fluxshifts, then come the magshifts
# 	# This should be ideal for Powell iterations, if your magshifts are already optimized.
# 	# Otherwise, do so !
# 	
# 	n = len(lcs[1:])
# 	def setshifts(p):
# 		for (i, l) in enumerate(lcs[1:]):
# 			l.setfluxshift(p[i] * 100.0, consmag = True)
# 			l.magshift = p[n + i]
# 			#print l.fluxshift, l.magshift
# 			
# 	inipars = np.concatenate([np.array([l.fluxshift * 0.01 for l in lcs[1:]]), np.array([l.magshift for l in lcs[1:]])])
# 	
# 	if verbose:
# 		print "Starting fluxshift optimization ..."
# 		print "Initial pars : ", inipars
# 
# 	def errorfct(p):
# 		setshifts(p)
# 		d2values = np.array([rawdispersionmethod(*couple)["d2"] for couple in couplelist])
# 		r = np.mean(d2values)
# 		#print r
# 		if np.isnan(r):
# 			raise RuntimeError("Ok guys, this had to happen ... be careful with the fluxshifts.")
# 		return r
# 
# 	minout = spopt.fmin_powell(errorfct, inipars, full_output=1, disp=verbose)
# 	popt = minout[0]
# 	if popt.shape == ():
# 		popt = np.array([popt])
# 		
# 	if verbose:
# 		print "Optimal pars : ", popt
# 	setshifts(popt)
# 
# def opt_fluxshift2(lcs, rawdispersionmethod, verbose=True):
# 	"""
# 	Brute force finding of plausible fluxshifts
# 	
# 	@todo: to get this stuff to work, we would need to fit the microlensing at each stage.
# 	This is slow...
# 	"""
# 	couplelist = [couple for couple in [[lc1, lc2] for lc1 in lcs for lc2 in lcs] if couple[0] != couple[1]]
# 
# 	# Just fluxshifts here.
# 		
# 	n = len(lcs[1:])
# 	def setshifts(p):
# 		if p.shape == ():
# 			p = np.array([p])
# 		for (i, l) in enumerate(lcs[1:]):
# 			l.setfluxshift(p[i], consmag = True)
# 			
# 			
# 	inipars = np.array([l.fluxshift for l in lcs[1:]])
# 	rad = 1000.0
# 	iniranges = [(inipar-rad, inipar+rad) for inipar in inipars]
# 	
# 	if verbose:
# 		print "Starting fluxshift optimization ..."
# 		print "Initial pars : ", inipars
# 
# 	def errorfct(p):
# 		#print p
# 		setshifts(p)
# 		d2values = np.array([rawdispersionmethod(*couple)["d2"] for couple in couplelist])
# 		r = np.mean(d2values)
# 		print p, r
# 		if np.isnan(r):
# 			raise RuntimeError("Ok guys, this had to happen ... be careful with the fluxshifts.")
# 		return r
# 
# 
# 
# 	minout = spopt.brute(errorfct, iniranges, full_output=1, finish=None, Ns=5)
# 	#print minout
# 	#sys.exit()
# 
# 	#minout = spopt.anneal(errorfct, inipars, full_output=1)
# 	popt = minout[0]
# 	if popt.shape == ():
# 		popt = np.array([popt])
# 		
# 	if verbose:
# 		print "Optimal pars : ", popt
# 	setshifts(popt)
# 




	