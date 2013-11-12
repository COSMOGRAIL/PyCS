# """
# Functions to optimize timeshifts, microlensing, or all this at the same time between an arbitrary set of curves.
# The top of the food chain :-)
# 
# """
# import sys
# 
# import numpy as np
# import scipy.optimize as spopt
# import matplotlib.pyplot as plt
# 
# import matplotlib.cm as cm
# import matplotlib.colors as colors
# 
# 
# #import disps
# import pycs.gen.ml as ml
# import pycs.gen.lc as lc
# import pycs.gen.util as util
# 
# 
# 
# def opt_ts(lcs, fitmethod, filename=None, verbose=True):
# 	"""
# 	Optimize the time shifts between lcs.
# 	"""
# 	
# 	# As we do not touch the microlensing here, we will work on copies with applied ml.
# 	lcsc = [l.copy() for l in lcs]
# 	for l in lcsc:
# 		if l.ml != None:
# 			l.applyml() # This will make the rest a lot faster, as we avoid evaluating the polynom every time !
# 			
# 	def chi2(delays):
# 		lc.multisettimedelays(lcsc, delays)
# 		retval = fitmethod(lcsc)["chi2n"]
# 		if verbose:
# 			print delays, retval
# 		return retval
# 	
# 	
# 	initparams = lc.multigettimedelays(lcsc)	
# 	initparams += np.random.randn(len(initparams))
# 	if verbose: print "Initial time delays : ", initparams
# 
# #	
# #	# Medium work with iterative brute force :
# #		
# #	res = spopt.brute(chi2, bruteranges(10,2,initparams), full_output = 0, finish=None)
# #	# This would finish by default with fmin ... we do not want that.
# #	if verbose:
# #		print "Brute 1 delays : %s" % res
# #		print "Brute 1 dispersion : %7.3f" % chi2(res)
# #	
# #	res = spopt.brute(chi2, bruteranges(4,2,res), full_output = 0, finish=None)
# #	if verbose:
# #		print "Brute 2 delays : %s" % res
# #		print "Brute 2 dispersion : %7.3f" % chi2(res)
# #	
# #	
# #	res = spopt.brute(chi2, bruteranges(1,2,res), full_output = 0, finish=None)
# #	if verbose:
# #		print "Brute 3 delays : %s" % res
# #		print "Brute 3 dispersion : %7.3f" % chi2(res)
# #
# #	#res = spopt.brute(d2, bruteranges(0.1,2,res), full_output = 0, finish=spopt.fmin)
# #	res = spopt.brute(chi2, bruteranges(0.1,2,res), full_output = 0, finish=None)
# #	
# 	res = spopt.fmin(chi2, initparams, disp=1, xtol=0.0001)
# 	res = spopt.fmin_powell(chi2, initparams)
# 	
# 	
# #	res = spopt.anneal(chi2, initparams, schedule='fast', lower=-10, upper=10, dwell=10, learn_rate=0.2, T0 = 1.0e-5, Tf = 1.0e-7, feps=0.0, boltzmann = 0.1,full_output = 1)
# #	if verbose:
# #		print "Annealing time delays : %s" % res[0]
# #		print "Dispersion value : %7.3f" % res[1]
# #		print "Final temperature : %.3g" % res[2]
# #		print "Dispersion calls on couples : %i" % res[3]
# #		if res[6] == 0:
# #			print "Cooled to global minimum."
# #		elif res[6] == 1:
# #			print "Cooled to minimum temperature."
# #		else:	
# #			print "### WARNING : %i ###" % res[6]
# #
# #	res = res[0]
# 	
# 	# Until now we worked on copies -- do not forget to shift the actual lightcurves now !
# 	lc.multisettimedelays(lcs, res)
# 
# 	
# 	#print res
# 	#sys.exit()
# 	
# 	optchi2 = chi2(res)
# 	#print 
# 	
# 	#if verbose:
# 	print "Final delays : \n%s" % lc.multigetnicetimedelays(lcs)
# 	
# 	
# 	print "Final chi2 : %7.3f" % optchi2
# 	
# 	
# 		
# 	return optchi2
# 	
# 	
# def opt_ml():
# 	pass
# 
# def bruteranges(step, radius, center):
# 	"""
# 	Auxiliary function for brute force exploration.
# 	Prepares the "ranges" parameter to be passed to brute force optimizer
# 	In other words, we draw a cube ...
# 	radius is an int saying how many steps to go left and right of center.
# 	center is a list or array of the centers, it can be of any lenght.
# 	
# 	You make 2*radius + 1 steps in each direction !, so radius=2 means 5 steps thus 125 calls for 4 curves.
# 	"""
# 	
# 	low = - step * radius
# 	up = step * (radius+1)
# 	return [((c+low),(c+up),step) for c in center]
# 
# 	
# 
# 
