"""
Functions to optimize shifts and microlensing between curves, minimizing a dispersion function of your choice.
Sometimes interactive plots are optional, the main goal is to return values or better change the lightcurves you pass

"""

import numpy as np
import scipy.optimize as spopt
import matplotlib.pyplot as plt

#import disps
import pycs.disp.twospec
import pycs.gen.polyml as ml

#import sys


def opt_ts_an(lc1, lc2, dispersionmethod, plot=False, verbose=True):
	"""
	A first atempt to optimize the timeshift using simulated annealing.
	We will shift lc2 in time so that it matches lc1.
	We return only the result code of the annealing.
	We do not optimize microlensing here.
	
	"""

	if plot:
		d2vals = []
		timeshifts = []
		inittimeshift = lc2.timeshift
		initd2val = dispersionmethod(lc1, lc2)['d2']
	
	def d2(timeshift):
		lc2.timeshift = timeshift
		ret = dispersionmethod(lc1, lc2)['d2']
		
		if plot:
			d2vals.append(ret)
			timeshifts.append(timeshift)
		return ret
		
	# http://www.scipy.org/doc/api_docs/SciPy.optimize.anneal.html
	res = spopt.anneal(d2, lc2.timeshift, schedule='fast', lower=-10, upper=10, dwell=6, learn_rate=0.5, T0 = 1.0e-3, Tf = 1.0e-9, full_output = 1)
	#print res
	
	# Do not forget to shift lc2 to optimal position :
	lc2.timeshift = res[0]
	
	
	if verbose :
		print "Optimal absolute shift for lc2 : %12.8f" % res[0]
		print "Dispersion value : %7.3f" % res[1]
		print "Final temperature : %.3g" % res[2]
		print "Dispersion calls : %i" % res[3]
		if res[6] == 0:
			print "Cooled to global minimum."
		elif res[6] == 1:
			print "Cooled to minimum temperature."
		else:	
			print "### WARNING : %i ###" % res[6]

	if plot:
		plt.plot(timeshifts, d2vals, color="#DDDDDD")
		plt.plot(timeshifts, d2vals, "b.")
		plt.plot([inittimeshift], [initd2val], "r.")
		plt.plot([res[0]], [res[1]], "r.")
		plt.title("Dispersion calls : %i" % len(timeshifts))
		plt.xlabel("Absolute lc2 timeshift [days]")
		plt.ylabel("Dispersion")
		plt.show()

	return res[6]
	
def opt_ts_brute(lc1, lc2, dispersionmethod, timewidth=300, timestep=0.1, verbose=True):
	"""
	"optimizes" timeshift of lc2 by calculating a full spectra ...
	
	"""
	
	# We prepare the timeshifts to explore
	timeshifts = np.arange(-(timewidth)*timestep/2.0, (timewidth+1)*timestep/2.0, timestep)
	if verbose : print "Exploring relative timeshifts from %f to %f" % (timeshifts[0], timeshifts[-1])
		
	spectrum = pycs.pelt.twospec.timeshift(lc1, lc2, timeshifts, dispersionmethod, fastret=True)
	minind = np.argmin(spectrum)
	lc2.timeshift += timeshifts[minind]
	
	if verbose : print "Minimum at absolute lc2 timeshfit : %f" %(lc2.timeshift)
	

	
	
	
	

def opt_mlts(lc1, lc2, dispersionmethod, niter = 2, startml = True, verbose = True):
	"""
	A simple try to optimize ml and timeshift : one iteration is 
	 - optimize ml
	 - then timeshift without changing ml
	
	"""
	for i in range(niter):
		print "### %s - %s delay  %f d, dispersion %f"% (lc1.object, lc2.object, lc2.timeshift - lc1.timeshift, dispersionmethod(lc1, lc2)['d2'])
		print "Iteration %i/%i" %(i+1, niter)
		
		if not (i==0 and not startml): 
			print "MICROLENSING"
			opt_ml(lc1, lc2, dispersionmethod, plot=False, verbose=verbose)
		print "TIME SHIFT"
		opt_ts_brute(lc1, lc2, dispersionmethod, verbose=verbose)
		#print "Time delay of %s with respect to %s : %8.3f days" % (lc1.object, lc2.object, lc2.timeshift - lc1.timeshift)
	print "### %s - %s delay  %f d, dispersion %f"% (lc1.object, lc2.object, lc2.timeshift - lc1.timeshift, dispersionmethod(lc1, lc2)['d2'])
		



def opt_ml(lc1, lc2, dispersionmethod, plot=False, verbose=False):
	"""
	microlensing optimization.
	Sets the eventual microlensings of lc1 and lc2 so to minimize a dispersion.
	Does not shift the curves in any other way.
	
	If plot is True, it displays the steps of the optimization.
	
	@type lc1: lightcurve
	@param lc1: This one will be passed as first argument (i.e. "reference") to the dispersionmethod
	@type lc2: lightcurve
	@param lc2: The light curve that will get interpolated by the dispersionmethod.
	
	@todo: implement a better plot ...
	
	"""

	initparams = ml.multigetfreeparams([lc1, lc2])
	if verbose:
		print "Initial params : ", initparams
		
	if len(initparams) == 0:
		raise RuntimeError, "There are no free microlensing params to optimize !"
		
	if plot:
		d2valuelist = []
	
	def d2(params):
		ml.multisetfreeparams([lc1, lc2], params)
		
		ret = dispersionmethod(lc1, lc2)["d2"]
		
		if plot:
			d2valuelist.append(ret)
		
		return ret
		
	# http://www.scipy.org/doc/api_docs/SciPy.optimize.optimize.html#fmin_bfgs
	#put disp = 0 to make them silent...

	#optparams = spopt.fmin(d2, initparams)
	#optparams = spopt.fmin_cg(d2, initparams)
	#optparams = spopt.fmin_bfgs(d2, params, epsilon =1.0e-5)
	#optparams = spopt.fmin_bfgs(d2, params, disp=1, retall=0)
	#optparams = spopt.anneal(d2, params)
	
	#optparams = spopt.fmin_ncg(d2, initparams) # this guy would need a gradient.
	#optparams = spopt.fmin_bfgs(d2, initparams) # precision loss
	
	
	# [  2.89270193e-10   4.94366701e-08  -2.63101369e-04  -4.38757115e-03]
	# [  2.95294380e-10   3.64696903e-08  -2.63882092e-04  -5.66560156e-05]
	
	
	
	optparams = spopt.fmin_powell(d2, initparams, disp=int(verbose)) # 312 , 2.592636, [  2.89270193e-10   4.94366701e-08  -2.63101369e-04  -4.38757115e-03]
	#optparams = spopt.fmin(d2, initparams) # 150, 2.614217
	#optparams = spopt.fmin(d2, initparams, xtol = 1.0e-10)
	#optparams = spopt.fmin_bfgs(d2, params, epsilon =1.0e-5)
	
	if verbose:
		print "Optimal params : ", optparams
		print repr(optparams)
	
	ml.multisetfreeparams([lc1, lc2], optparams)
	
	if plot:
		plt.semilogy(d2valuelist)
		plt.show()


#def opt_ml_2seas(lc1, lc2, dispersionmethod):
#	"""
#	This is radically different, and meant to be faster, we optimize only the ml of lc2, one season after the other.
#	
#	@todo: to implement
#	"""
#	pass
#




#def opt_mag(lclist1, lclist2, dispersionmethod, guessmagshift=0.0, magwidth = 20, magstep=0.01, plot=True):
#	"""
#	A function to specifically optimize the shift in magnitude between two light curves that should match,
#	without touching the time axis.
#	As you can have multiple objects that require the same of those magshifts, give lists of corresponding lightcurves (lists of the same length !), not juste lightcurves.
#	Example : say you have lightcurves of the same lens from Mercator and HCT, namely C{merA}, C{merB}, C{hctA}, C{hctB}. Then you should call:
#	
#		>>> magshiftopt([merA, merB], [hctA, hctB], dispersionmethod, ...)
#	
#	to find the single magshift that matches merA on hctA and merB on hctB.
#	The "figure of merit" to minimize is simply the sum of the d2 of each lightcurve pair.
#	
#	@warning: Please provide lightcurves without intrinsic magnitude shifts !
#	
#	@todo: under construction...
#	
#	@type lclist1: list of lightcurves
#	@param lclist1: These will be passed as first argument (i.e. "reference") to the dispersionmethod.
#	@type lclist2: list of lightcurves
#	@param lclist2: The light curves that will get interpolated by the dispersionmethod.
#	
#	@type dispersionmethod: function
#	@param dispersionmethod: It has to be a function that takes only two lightcurves as arguments, see the functions of specplots.
#	
#	
#	@type guessmagshift: float
#	@param guessmagshift: a first guess for the shift to apply to lc2 : steps will be centered around that value
#	@type magwidth: int
#	@param magwidth : the number of magnitude steps to explore, only used for the plot
#	@type magstep: float
#	@param magstep: the step, only used for the plot
#	
#	@type plot: boolean
#	@param plot: if True, shows a dispersion spectrum "magshift vs d2"
#	
#	@rtype: float
#	@return: the best magshift to be applied to lclist2 so that it matches lclist1
#	
#	@todo: this works, but is quite old and should be updated.
#	
#	"""
#	
#	# We check that there are no magnitude shifts applied to the lightcurves (otherwise it would 1) get messy 2) probably be a user error)
#	# This could be relaxed, if needed.
#	for lc in lclist1 + lclist2:
#		if lc.magshift != 0.0:
#			raise RuntimeError, "shiftopt.mag : lightcurve %s is already shifted in magnitude..."% str(lc)
#		if lc.ml != None:
#			raise RuntimeError, "shiftopt.mag : lightcurve %s has microlensing..."% str(lc)
#	if len(lclist1) != len(lclist2):
#		raise RuntimeError, "shiftopt.mag : lclists have different lengths !"
#		
#	print "Optimizing magnitude shift between (check order !)"
#	print "  lclist1 : %s" % " + ".join([str(l) for l in lclist1])
#	print "  lclist2 : %s" % " + ".join([str(l) for l in lclist2])
#	
#	def d2(magshift):
#		lclist2shifted = [lc.copy() for lc in lclist2] # we make a copy of the lightcurve list
#		for lc in lclist2shifted: # we shift all the lc2s
#			lc.shiftmag(magshift)
#		# we calculate the dispersion for each lc pair :
#		d2list = np.array([dispersionmethod(lc1, lc2shifted)["d2"] for (lc1, lc2shifted) in zip(lclist1, lclist2shifted)])
#		# we sum the squares of the dispersions :
#		d2value = np.sum(d2list)
#		return d2value
#		#return 4.0
#	
#	# and we simply optimize this with a "Nelder-Mead simplex algorithm" AAHHAAAA ! :
#	res = spopt.fmin(d2, guessmagshift, xtol=0.000000001, ftol=0.000000001)
#	print "Optimal magnitude shift for lclist2 : %12.8f" % res[0]
#	
#	if plot:
#		magshifts = np.arange(-(magwidth)*magstep/2.0, (magwidth+1)*magstep/2.0, magstep) + guessmagshift
#		# these shifts will be applied "as is" to the elements of lclist2 ... we do not try to cope with possible previous magnitude shifts.
#		#print magshifts
#	
#		vecd2 = np.vectorize(d2, otypes=[float])
#		# We apply it to our magshifts :
#		d2values = vecd2(magshifts)
#
#		plt.plot(magshifts, d2values, "r.")
#		
#		plt.xlabel("magshifts", size=16)
#		plt.ylabel("d2values", size=16)
#		plt.axvline(x=res, linewidth=1, color='r')
#		
#		plt.show()
#
#	
#	return res[0]


def opt_mag(lclist1, lclist2, dispersionmethod, guessmagshift=0.0, magwidth = 50, magstep=0.002, plot=False):
	"""
	A function to specifically optimize the shift in magnitude between two light curves that should match,
	without touching the time axis.
	As you can have multiple objects that require the same of those magshifts, give lists of corresponding lightcurves (lists of the same length !), not juste lightcurves.
	Example : say you have lightcurves of the same lens from Mercator and HCT, namely C{merA}, C{merB}, C{hctA}, C{hctB}. Then you should call:
	
		>>> magshiftopt([merA, merB], [hctA, hctB], dispersionmethod, ...)
	
	to find the single magshift that matches merA on hctA and merB on hctB.
	The "figure of merit" to minimize is simply the sum of the d2 dispersions of each lightcurve pair.
	
	@warning: Please provide lightcurves without intrinsic magnitude shifts !
	
	@type lclist1: list of lightcurves
	@param lclist1: These will be passed as first argument (i.e. "reference") to the dispersionmethod.
	@type lclist2: list of lightcurves
	@param lclist2: The light curves that will get interpolated by the dispersionmethod.
	
	@type dispersionmethod: function
	@param dispersionmethod: It has to be a function that takes only two lightcurves as arguments.
	
	@type guessmagshift: float
	@param guessmagshift: a first guess for the shift to apply to lc2 : steps will be centered around that value
	@type magwidth: int
	@param magwidth : the number of magnitude steps to explore, only used for the plot
	@type magstep: float
	@param magstep: the step, only used for the plot
	
	@type plot: boolean
	@param plot: if True, shows a dispersion spectrum "magshift vs d2"
	
	@rtype: float
	@return: the best magshift to be applied to lclist2 so that it matches lclist1
	
	"""
	
	# We check that there are no magnitude shifts applied to the lightcurves (otherwise it would 1) get messy 2) probably be a user error)
	# But this could be relaxed, if needed -- there is no "fundamental" reason.
	for lc in lclist1 + lclist2:
		if lc.magshift != 0.0:
			raise RuntimeError, "shiftopt.mag : lightcurve %s is already shifted in magnitude..."% str(lc)
		if lc.ml != None:
			raise RuntimeError, "shiftopt.mag : lightcurve %s has microlensing..."% str(lc)
	if len(lclist1) != len(lclist2):
		raise RuntimeError, "shiftopt.mag : lclists have different lengths !"
	
	# To be flexible we do not test that we compare only curves representing the same image : check this by yourself !
	print "Optimizing magnitude shift between (double-check respective order !)"
	print "  lclist1 : %s" % " + ".join([str(l) for l in lclist1])
	print "  lclist2 : %s" % " + ".join([str(l) for l in lclist2])
	
	def d2(magshift):
		
		for lc2 in lclist2: # we shift all the lc2s
			lc2.magshift = magshift
		# we calculate the dispersion for each lc pair :
		d2list = np.array([dispersionmethod(lc1, lc2)["d2"] for (lc1, lc2) in zip(lclist1, lclist2)])
		# we mean the squares of the dispersions :
		d2value = np.mean(d2list)
		return d2value
		# There is no need to set the magshifts back to 0 here. We apply the optimal magshift below.
	
	# http://www.scipy.org/doc/api_docs/SciPy.optimize.optimize.html#fmin_bfgs
	# and we simply optimize this with a "Nelder-Mead simplex algorithm" AAHHAAAA ! :
	res = spopt.fmin(d2, guessmagshift, xtol=0.000000001, ftol=0.000000001, disp=0)
	optmagshift = res[0]
	print "Optimal magnitude shift for lclist2 : %12.8f" % optmagshift
	
	# Now we set the shifts to the optimal value:
	for lc2 in lclist2:
		lc2.magshift = optmagshift
	
	if plot:
		magshifts = np.arange(-(magwidth)*magstep/2.0, (magwidth+1)*magstep/2.0, magstep) + guessmagshift
		# these are absolute magshifts for lclist2
	
		vecd2 = np.vectorize(d2, otypes=[float])
		# We apply it to our magshifts :
		d2values = vecd2(magshifts)

		plt.plot(magshifts, d2values, "r.")
		
		plt.xlabel("Absolute magshifts of lclist2", size=16)
		plt.ylabel("Dispersion d2", size=16)
		plt.axvline(x=optmagshift, linewidth=1, color='r')
		
		plt.show()




