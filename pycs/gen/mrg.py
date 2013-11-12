"""
Functions to match and merge lightcurves usually of different telescopes, prior to exploring time delays.
We typically adjust the magshifts and fluxshifts of curves from telescope X so that they match
the same curves from telescope Y, without touching the time axis.
Of course we need overlapping curves to do this.
"""

import numpy as np
#import matplotlib.pyplot as plt
#import scipy.optimize as spopt
#import scipy.interpolate as si

#import pycs.pelt.disps
#import pycs.pelt.multiopt

from scipy.optimize import fmin, fmin_powell


def matchtels(lcsref, lcsmatch, dispersionmethod, fluxshifts = True):
	"""
	Shifts the curves of lcsmatch so that they fit to lcsref.
	lcsmatch and lcsref are *lists* of corresponding lightcurves.
	We do not modify lscref.
	
	To help me a bit, set fluxshift and magshifts of the lcsmatch to some good 
	starting values before calling me !
	
	For now, I do :
		- one single magshift for all curves
		- one fluxshift for each curve (if fluxshifts = True)
	"""

	if len(lcsref) != len(lcsmatch):
		raise RuntimeError("Input lists must have the same length.")
	print "Matching curves : "
	for (lr, lm) in zip(lcsref, lcsmatch):
		print "   %20s -> %20s" % (lm, lr)
	
 
 	def setp(p, lcs):
 		"""
 		Sets the parameters p to the lightcurves lcs
 		All parameters have typical values around 1.0 here !
 		"""
 		
 		if fluxshifts:
 			magshift = p[0] * 0.1
 			fluxes = p[1:] * 100.0
 			for (l, fs) in zip(lcs, fluxes):
 				# We can't use the relative set methods here !
				l.fluxshift = fs
				l.magshift = magshift
 		else:
 			magshift = p[0] * 0.1
 			for l in lcs:
				l.magshift = magshift
 	
 	def inip():
 		"""
 		Returns a list of some initial values for the parameters
 		For this we use the current "settings" of these parameters of lcmatch.
 		"""
 		#if fluxshifts:
 		#	return np.ones(1 + len(lcsref))
		#else :
 		#	return [1.0]
  		
  		magpar = np.mean(np.array([l.magshift for l in lcsmatch])) * 10.0
  		if fluxshifts:
  			fluxes = [l.fluxshift for l in lcsmatch]
  			fluxpars = 0.01 * np.array(fluxes)
  			pars = np.concatenate([np.array([magpar]), fluxpars])
  			print "Initial pars : ", pars
  			return pars
  			
  		else:
  			pars = np.array([magpar])
  			print "Initial pars : ", pars
  			return pars
  			
  		
  		
	def errorfct(p):
		"""
		The error function to minimize, essentially a sum of dispersions of curve pairs.
		p is the array of parameters that will be optimized.
		p[0] = magshift, the same for all curves, ini = 0.1
		"""
	
		lcsmod = [l.copy() for l in lcsmatch]
		
		setp(p, lcsmod)
 		
 		# We calculate the error funtion :
		totdisp = 0.0
 		for (lm, lr) in zip(lcsmod, lcsref):
 			totdisp += dispersionmethod(lm, lr)["d2"]
	 	
	 	return totdisp

	# We optimize p :
	print "Starting Powell optimization ..."
	# This is very well suited, as it starts by optimizing the mag shifts, then the fluxes, then iterates.
	minout = fmin_powell(errorfct, inip(), full_output=1)
	popt = minout[0] # This might be an array, or a 0-d array. We want an array :
	#if not isinstance(popt, np.ndarray):
	#	popt = np.array([popt])
	if popt.shape == ():
		popt = np.array([popt])
	print "Optimal pars : ", popt
	setp(popt, lcsmatch)


def merge(lcss):
	"""
	Function to merge (across the first index) a list of lists of lightcurves.
	example :
	::
	
		lcss = [[EulerC2_A, EulerC2_B], [Mercator_A, Mercator_B]]
	
		# will *return* [merged_A, merged_B]
		# where merged_i is EulerC2_i.merge(Mercator_i)
	
	.. note:: I do not modify the input lightcurves !
	
	"""
	reflcs = [l.copy() for l in lcss[0]] # Important : we do not modify lcss !
	for otherlcs in lcss[1:]:
		for (reflc, otherlc) in zip(reflcs, otherlcs):
			reflc.merge(otherlc)
	return reflcs


def colourise(lcs):
	"""
	Attributes different colours to the lightcurves of the list lcs.
	
	http://www.w3schools.com/html/html_colornames.asp
	"""
	mycolours = ["red", "blue", "#009900", "purple", "brown", "magenta", "orange"]

	for (l,c) in zip(lcs, mycolours):
		l.plotcolour = c

def colorize(lcs):
	colourise(lcs)
	
