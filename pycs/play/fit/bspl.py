# """
# B splines
# 
# 
# somewhere we should reuse coeffs here... so that a next fit looks if is has the same number of knots, and if yes, uses the previous fit as starting value.
# 
# """
# 
# import sys
# 
# from pycs.gen import *
# 
# import numpy as np
# import math
# import matplotlib.pyplot as plt
# import scipy.optimize as spopt
# import scipy.interpolate as spip
# 
# 
# 
# def fitcubbspline(x, y, yerr, t, cini=None, verbose=True):
# 	"""
# 	This is "my own" cubic B-spline fitting method, using leastsq from scipy.
# 	I know, this looks like a very naive idea from somebody who has not a clue what a spline is...
# 	But... recall that we want to 
# 	1) control the positions of the knots (probably on a fixed grid etc)
# 	2) have an irregular sampling of points (and perhaps also of knots)
# 	3) IMPORTANT : and that it may well be that we have several "points" to fit for one single JD
# 	(this last point kills sp.interpolate.splrep as well as probably all direct methods, inverse filtering etc) !
# 
# 	x y : the data
# 	t : the x-positions of the knots, WITH "prefix" and "suffix" knots !
# 	cini : initial coeffs for the knots t. If None, we start from zero.
# 	
# 
# 	We use the notation "t, c, k" from scipy.interpolate.splev etc : t the knots, c the corresponding coeffs, k the degree
# 	"""
# 
# 	
# 	k = 3 # cubic spline
# 	if cini == None:
# 		if verbose :
# 			print "Spline fit : starting from 0.0"
# 		cini = np.zeros(len(t)) # initial coeffs = 0.0
# 	else:
# 		# we check that they are compatible :
# 		if len(cini) != len(t):
# 			raise RuntimeError, "Spline fit : cini has the wrong length"
# 
# 	leastsqres = spopt.leastsq(splinefiterrfunc, cini, args=(t, k, x, y, yerr), full_output=1)
# 	# this should be faster without the full output...
# 	
# 	if verbose:
# 		print "Spline fit : %i function calls" % leastsqres[2]['nfev']
# 
# 	fittedc = leastsqres[0]
# 	tck = (t, fittedc, k)
# 
# 	return tck
# 
# 
# def splinefiterrfunc(c, t, k, x, y, yerr):
# 	"""
# 	Auxiliary function for the fit.
# 	Give me a spline (c,t,k) and some points (x, y, yerr) and I return the vector of differences.
# 	"""
# 	tck = (t, c, k)
# 	interpy = spip.splev(x, tck, der=0)
# 	return (y - interpy)/yerr
# 
# 
# 
# def knots(x, sheme = "test1"):
# 	"""
# 	Give me the x coords of some point, and I give you some knots according to a given sheme...
# 	
# 	stupid : for testing purposes
# 	
# 	test1 : I want to get the knots selected from a "fixed" absolute grid, JD = 0, n, 2n, 3n ...
# 	
# 	We add the extra coeffictients at both sides of t. See
# 	http://mathworld.wolfram.com/B-Spline.html
# 	
# 	"""
# 	# In this first step, we add the "interior" knots (i.e. containing the extremal knots, but not repeating them).
# 	
# 	if sheme == "stupid" :
# 		t = np.linspace(x[0], x[-1], 10)
# 		
# 	if sheme == "test1" :
# 		n = 10 # n = 6 -> so the grid is 0, 6, 12, ... ie equally spaced by 6
# 		first = x[0] - (x[0] % n) # this is nearly an int, but no need to convert to int().
# 		t = np.arange(first, x[-1] + n, n)
# 		
# 	# We add some extra coefficients at both sides of t
# 	prefix = np.ones(3)*t[0]
# 	suffix = np.ones(3)*t[-1]
# 	fullt = np.concatenate([prefix, t, suffix])	
# 	
# 	return fullt
# 
# 
# def cbsp(lcs, splitgap=60, oldtcks=None, verbose=True, plot=True):
# 	"""
# 	First try to get a cubic B-spline fit working, simultaneously for n lightcurves, and return a chi2 like something.
# 	Give me a list of lightcurves and I return you a value for chi2, using a specified spline fit etc
# 	
# 	oldtcks : if not "None" but a list of tcks, we will try to start the fit of the spline coeffs using these ...
# 	This typically works if the number of knots has not changed, i.e. when we optimize microlensing...
# 	
# 	"""
#  
# 	jdslist = []
# 	magslist = []
# 	magerrslist = []		
# 	for thislc in lcs:
# 		jdslist.append(thislc.getjds()[thislc.mask])
# 		magslist.append(thislc.getmags()[thislc.mask])
# 		magerrslist.append(thislc.magerrs[thislc.mask])
# 		
# 	mergedjds = np.concatenate(jdslist)
# 	mergedmags = np.concatenate(magslist)
# 	mergedmagerrs = np.concatenate(magerrslist)
# 
# 	# Now the sorting ...
# 	sortedindices = np.argsort(mergedjds)
# 	sortedjds = mergedjds[sortedindices]
# 	sortedmags = mergedmags[sortedindices]
# 	sortedmagerrs = mergedmagerrs[sortedindices]
# 
# 	# We need to find the overlapping regions ?
# 	# For a first try, let's split this as usual :
# 	first = sortedjds[:-1]
# 	second = sortedjds[1:]
# 	gapindices = np.where(second - first > splitgap)[0] + 1
# 	# make one big vector of all the indices :
# 	indices = np.arange(len(sortedjds))
# 	# split it according to the gaps :
# 	indexlist = np.hsplit(indices, gapindices)
# 	if verbose:
# 		print "We have %i splines." % len(indexlist)
# 	
# 	if oldtcks == None:
# 		# Then we do not have previous splines to start from
# 		oldtcks = [None] * len(indexlist)
# 
# 	tcks = [] # we will append here the splines from the individual spline fits, not only for plotting, also to pass them to the next call !
# 	chi2s = [] # the resulting chi2s
# 	ns = [] # the number of points for that spline
# 	
# 	for indexes, oldtck in zip(indexlist, oldtcks): # i.e. for each "season" aka "group" of points
# 		jds = sortedjds[indexes]
# 		mags = sortedmags[indexes]
# 		magerrs = sortedmagerrs[indexes]
# 
# 		t = knots(jds, sheme="test1")
# 		if (oldtck != None) and (len(t) == len(oldtck[0])) : # Then we should be able to reuse this...
# 			tck = fitcubbspline(jds, mags, magerrs, t, cini=oldtck[1], verbose=False)		
# 		else:
# 			tck = fitcubbspline(jds, mags, magerrs, t, verbose=False)
# 
# 		#if verbose:
# 		#	for (t, c) in zip(tck[0], tck[1]):
# 		#		print "t = %8.3f -> c = %8.3f" % (t, c)
# 
# 		tcks.append(tck)
# 	
# 		diffs = (mags - spip.splev(jds, tck, der=0))/magerrs
# 		chi2 = np.sum(diffs * diffs)
# 
# 		if verbose:
# 			print "chi2 : %8.3f for %i points" % (chi2, len(jds))
# 		
# 		chi2s.append(chi2)
# 		ns.append(len(jds))
# 	
# 	totchi2 = np.sum(np.array(chi2s))
# 	totn = np.sum(np.array(ns))
# 	chi2n = totchi2/float(totn)
# 	if verbose:
# 		print "tot  : %8.3f for %i points" % (totchi2, totn)
# 		print "chi2n: %8.3f" % (chi2n)
# 	
# 
# 
# 	if plot:
# 		plt.figure(figsize=(12,8))	# sets figure size
# 		axes = plt.gca()	
# 		
# 		# The points
# 		plt.errorbar(sortedjds, sortedmags, sortedmagerrs, fmt=".", color="red", ecolor="#BBBBBB")
# 		
# 		# The groups
# 		for groupindexes in indexlist:
# 			plt.axvline(sortedjds[groupindexes[0]], color = "#AAAAAA", dashes = (5,5))
# 			plt.axvline(sortedjds[groupindexes[-1]], color = "#AAAAAA", dashes = (5,5))
# 		
# 		# The spline
# 		for (tck, indexes) in zip(tcks, indexlist):
# 			xnew = np.linspace(sortedjds[indexes][0], sortedjds[indexes][-1], 1000)
# 			ynew = spip.splev(xnew,tck,der=0)
# 
# 			plt.plot(xnew, ynew, color="blue")
# 			
# 			for knot in tck[0]:
# 				plt.axvline(knot, color = "#0000AA", dashes = (2,2))
# 		
# 		
# 		# Splines may get huge, so we want to limit the axis ranges :
# 		axes.set_ylim((min(sortedmags) - 0.1, max(sortedmags) + 0.1))
# 	
# 		
# 		# Something for astronomers only : we invert the y axis direction !
# 		axes.set_ylim(axes.get_ylim()[::-1])
# 		
# 				# And we make a title for that combination of lightcurves :
# 		#plt.title("Lightcurves", fontsize=20)
# 		plt.xlabel("Days", fontsize=16)
# 		plt.ylabel("Magnitude", fontsize=16)
# 		plt.title("Spline", fontsize=16)
# 		
# 		plt.xlim([2340, 5000])
# 		
# 		plt.show()
# 
# 	
# 	return {'chi2':totchi2, 'n':totn, 'chi2n':chi2n, 'tcks':tcks}
# 
# 
# #def test(lcs, splitgap=60, usemask=True, verbose=True, plot=True):
# #	"""
# #	First try to get a cubic B-spline fit working, simultaneously for n lightcurves, and return a chi2 like something.
# #	Give me a list of lightcurves and I return you a value for chi2, using a specified spline fit etc
# #	
# #	
# #	"""
# #
# #	jdslist = []
# #	magslist = []
# #	magerrslist = []		
# #	for thislc in lcs:
# #		jdslist.append(thislc.getjds()[thislc.mask])
# #		magslist.append(thislc.getmags()[thislc.mask])
# #		magerrslist.append(thislc.magerrs[thislc.mask])
# #		
# #	mergedjds = np.concatenate(jdslist)
# #	mergedmags = np.concatenate(magslist)
# #	mergedmagerrs = np.concatenate(magerrslist)
# #
# #	# Now the sorting ...
# #	sortedindices = np.argsort(mergedjds)
# #	sortedjds = mergedjds[sortedindices]
# #	sortedmags = mergedmags[sortedindices]
# #	sortedmagerrs = mergedmagerrs[sortedindices]
# #
# #	# We need to find the overlapping regions ?
# #	# For a first try, let's split this as usual :
# #	first = sortedjds[:-1]
# #	second = sortedjds[1:]
# #	gapindices = np.where(second - first > splitgap)[0] + 1
# #	# make one big vector of all the indices :
# #	indices = np.arange(len(sortedjds))
# #	# split it according to the gaps :
# #	indexlist = np.hsplit(indices, gapindices)
# #	if verbose:
# #		print "We have %i splines." % len(indexlist)
# #	
# #
# #	tcks = [] # we will append here the splines from the individual spline fits (for plotting ...)
# #	chi2s = [] # the resulting chi2s
# #	ns = [] # the number of points for that spline
# #	
# #	for indexes in indexlist: # i.e. for each "season" aka "group" of points
# #		jds = sortedjds[indexes]
# #		mags = sortedmags[indexes]
# #		magerrs = sortedmagerrs[indexes]
# #
# #		#t = knots(jds, sheme="test1")	
# #		#tck = fitcubbspline(jds, mags, magerrs, t, verbose=False)
# #		
# #		tck = spip.splrep(jds, mags, w=(1.0/magerrs))
# #		print tck
# #		#maspline = spip.UnivariateSpline(jds, mags, w=magerrs, k=3)
# #
# #		#tck = [maspline.get_knots(), maspline.get_coeffs(), 3]
# #		#print len(maspline.get_knots())
# #
# #		#if verbose:
# #		#	for (t, c) in zip(tck[0], tck[1]):
# #		#		print "t = %8.3f -> c = %8.3f" % (t, c)
# #
# #		tcks.append(tck)
# #	
# #		diffs = (mags - spip.splev(jds, tck, der=0))/magerrs
# #		chi2 = np.sum(diffs * diffs)
# #
# #		if verbose:
# #			print "chi2 : %8.3f for %i points" % (chi2, len(jds))
# #		
# #		chi2s.append(chi2)
# #		ns.append(len(jds))
# #	
# #	totchi2 = np.sum(np.array(chi2s))
# #	totn = np.sum(np.array(ns))
# #	chi2n = totchi2/float(totn)
# #	if verbose:
# #		print "tot  : %8.3f for %i points" % (totchi2, totn)
# #		print "chi2n: %8.3f" % (chi2n)
# #	
# #
# #
# #	if plot:
# #		plt.figure(figsize=(12,8))	# sets figure size
# #		axes = plt.gca()	
# #		
# #		# The points
# #		plt.errorbar(sortedjds, sortedmags, sortedmagerrs, fmt=".", color="red", ecolor="#BBBBBB")
# #		
# #		# The groups
# #		for groupindexes in indexlist:
# #			plt.axvline(sortedjds[groupindexes[0]], color = "#AAAAAA", dashes = (5,5))
# #			plt.axvline(sortedjds[groupindexes[-1]], color = "#AAAAAA", dashes = (5,5))
# #		
# #		# The spline
# #		for (tck, indexes) in zip(tcks, indexlist):
# #			xnew = np.linspace(sortedjds[indexes][0], sortedjds[indexes][-1], 1000)
# #			ynew = spip.splev(xnew,tck,der=0)
# #
# #			plt.plot(xnew, ynew, color="blue")
# #			
# #			for knot in tck[0]:
# #				plt.axvline(knot, color = "#0000AA", dashes = (2,2))
# #		
# #		
# #		# Splines may get huge, so we want to limit the axis ranges :
# #		axes.set_ylim((min(sortedmags) - 0.1, max(sortedmags) + 0.1))
# #	
# #		
# #		# Something for astronomers only : we invert the y axis direction !
# #		axes.set_ylim(axes.get_ylim()[::-1])
# #		
# #				# And we make a title for that combination of lightcurves :
# #		#plt.title("Lightcurves", fontsize=20)
# #		plt.xlabel("Days", fontsize=16)
# #		plt.ylabel("Magnitude", fontsize=16)
# #		plt.title("Spline", fontsize=16)
# #		
# #		plt.xlim([2340, 5000])
# #		
# #		plt.show()
# #
# #	
# #	#return {'chi2':totchi2, 'n':totn, 'chi2n':chi2n}
# #
# #
# #def one(lcs, verbose=True, plot=True):
# #	"""
# #	Trying to build one big spline over all the groups...
# #	
# #	"""
# #
# #	
# #	jdslist = []
# #	magslist = []
# #	magerrslist = []		
# #	for thislc in lcs:
# #		jdslist.append(thislc.getjds()[thislc.mask])
# #		magslist.append(thislc.getmags()[thislc.mask])
# #		magerrslist.append(thislc.magerrs[thislc.mask])
# #		
# #	mergedjds = np.concatenate(jdslist)
# #	mergedmags = np.concatenate(magslist)
# #	mergedmagerrs = np.concatenate(magerrslist)
# #	
# #	# Now the sorting ...
# #	sortedindices = np.argsort(mergedjds)
# #	sortedjds = mergedjds[sortedindices]
# #	sortedmags = mergedmags[sortedindices]
# #	sortedmagerrs = mergedmagerrs[sortedindices]
# #	
# #	jds = sortedjds
# #	mags = sortedmags
# #	magerrs = sortedmagerrs
# #
# #	#t = np.linspace(jds[0], jds[-1], 10) # the knots
# #	t = np.arange(int(math.floor(jds[0])), int(math.ceil(jds[-1])), 30)
# #
# #	tck = fitcubbspline(jds, mags, magerrs, t, verbose=True)
# #
# #	if plot:
# #		plt.figure(figsize=(12,8))	# sets figure size
# #		axes = plt.gca()	
# #		
# #		# The points
# #		plt.errorbar(sortedjds, sortedmags, sortedmagerrs, fmt=".", color="red", ecolor="#BBBBBB")
# #		
# #		# The spline
# #		
# #		xnew = np.linspace(sortedjds[0], sortedjds[-1], 1000)
# #		ynew = spip.splev(xnew,tck,der=0)
# #
# #		plt.plot(xnew, ynew, color="blue")
# #			
# #		for knot in tck[0]:
# #			plt.axvline(knot, color = "#0000AA", dashes = (2,2))
# #		
# #		# Something for astronomers only : we invert the y axis direction !
# #		axes.set_ylim(axes.get_ylim()[::-1])
# #	
# #		# And we make a title for that combination of lightcurves :
# #		#plt.title("Lightcurves", fontsize=20)
# #		plt.xlabel("Days", fontsize=16)
# #		plt.ylabel("Magnitude", fontsize=16)
# #		plt.title("Spline", fontsize=16)
# #		
# #		plt.xlim([2340, 5000])
# #		
# #		plt.show()
# #
# #	
# #	diffs = (mags - spip.splev(jds, tck, der=0))/magerrs
# #	chi2 = np.sum(diffs * diffs)
# #
# #	if verbose:
# #		print "Chi2 : %8.3f" % chi2
# #
# #
# #	return chi2
# #
# #
