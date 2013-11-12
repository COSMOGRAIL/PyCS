"""
Defines dispersion functions. No optimization is done here !
We just calculate the dispersion.

* linintnp is the "default" one
* pelt95
* others should be ported into pycs from my f90 code :-(

They all share a common structure : give me 2 light curves as input, and you get a float describing the match.
By construction the dispersion is often not symmetric, i.e. if you switch lc1 and lc2 you will get a different result.
To symmetrize the techniques, use the provided symmetrize function.
I return a dictionnary, containing a field "d2" that gives the dispersion (among possible other fields).
Any microlensing is naturally taken into account, if present, simply using the getmags() method.

"""



#from pycs.gen import *

import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.optimize as spopt
import scipy.interpolate as ip



def linintnp(lc1, lc2, interpdist = 30.0, weights = True, usemask = True, plot = False):
	"""
	This is the default one, fast implementation without loops (25 times faster then maltef90).
	If usemask == True, it is a bit slower...
	If you do not change the mask, it's better to cutmask() it, then use usemask = False.
	
	"""
	if usemask :
		lc1jds = lc1.getjds()[lc1.mask] # lc1 is the "reference" curve
		lc2jds = lc2.getjds()[lc2.mask] # lc2 is the curve that will get interpolated.
		lc1mags = lc1.getmags()[lc1.mask]
		lc2mags = lc2.getmags()[lc2.mask]
		lc1magerrs = lc1.magerrs[lc1.mask]
		lc2magerrs = lc2.magerrs[lc2.mask]
	else:
		lc1jds = lc1.getjds() # lc1 is the "reference" curve
		lc2jds = lc2.getjds() # lc2 is the curve that will get interpolated.
		lc1mags = lc1.getmags()
		lc2mags = lc2.getmags()
		lc1magerrs = lc1.magerrs
		lc2magerrs = lc2.magerrs

	# --- Mask creation : which points of lc1 can be taken into account ? ---
	
	### double-sampled version ###
	
	interpdist = interpdist / 2.0
	
	# We double-sample the lc2jds:
	add_points = np.ediff1d(lc2jds)/2.0 + lc2jds[:-1]
	lc2jdsa = np.sort(np.concatenate([lc2jds, add_points])) # make this faster
	# in fact funnily this turns out to be fast... I have tried
	#lc2jdsa = np.ravel(np.vstack((lc2jds, add_points)), order="F") ... it is slower !?
	#lc2jdsa = np.ravel(np.array([lc2jds, add_points]), order="F")
	#lc2jdsa = np.array([lc2jds, add_points]).reshape(2*len(lc2jds), order="F")	
	# ... and they are not faster...
	
	# Not too wrong, but does not work as masks are not correct for small interpdists :
	
	 # BEGIN COMMENTED
	 
		
	# We want to create of mask telling us which jds of lc1 can be interpolated on lc2,
	# without doing any python loop. We will build an auxiliary array, then interpolate it.
	lc2interpcode_lr = np.cast['int8'](np.ediff1d(lc2jdsa, to_end = interpdist+1, to_begin = interpdist+1) <= interpdist) # len lc2 + 1
	lc2interpcode_sum = lc2interpcode_lr[1:] + lc2interpcode_lr[:-1] # len lc2
	# This last array contains ints :
	#	2 if both right and left distances are ok 
	#	1 if only one of them is ok
	#	0 if both are two fare away (i.e. fully isolated point)
	# We now interpolate on this (filling False for the extrapolation):
	lc1jdsmask = np.interp(lc1jds, lc2jdsa, lc2interpcode_sum, left=0, right=0) >= 1.0
	
	# we would also use np.diff and so something like :
	#diffs = np.cast['int8'](np.diff(lc2jds) <= interpdist)
	#diffr = np.append(diffs, 0)
	#diffl = np.insert(diffs, 0, 0)
	#print diffr + diffl
	
	 # END COMMENTED
	
	"""
	# old way doing two interpolations.
	lc2interpcode_r = np.ediff1d(lc2jds, to_end = interpdist+1)
	lc2interpcode_l = np.ediff1d(lc2jds, to_begin = interpdist+1)
	
	lc1jdsmask_r = np.interp(lc1jds, lc2jds, lc2interpcode_r, left=0, right=0) > interpdist
	lc1jdsmask_l = np.interp(lc1jds, lc2jds, lc2interpcode_l, left=0, right=0) > interpdist
	
	lc1jdsmask = np.logical_and(lc1jdsmask_r, lc1jdsmask_l)
	"""
			
	# --- Magnitude interpolation for all points, disregarding the mask (simpler -> faster) ---
	# These left and right values are a safety check : they should be masked, but if not, you will see them ...
	ipmagdiffs = np.interp(lc1jds, lc2jds, lc2mags, left=1.0e5, right=1.0e5) - lc1mags
	
	# --- Cutting out the mask
	okdiffs = ipmagdiffs[lc1jdsmask]
	n = len(okdiffs)
	
	if weights:
		# Then we do something similar for the errors :
		ipmagerrs = np.interp(lc1jds, lc2jds, lc2magerrs, left=1.0e-5, right=1.0e-5)
		ipsqrerrs = ipmagerrs*ipmagerrs + lc1magerrs*lc1magerrs # interpolated square errors
		oksqrerrs = ipsqrerrs[lc1jdsmask]
		
		d2 = np.sum((okdiffs*okdiffs)/oksqrerrs)/(n)
	else:
		d2 = np.sum((okdiffs*okdiffs))/(n)
	
	if plot:
		# Then we make a plot, essentially for testing purposes. Nothing important here !
		# This is not optimal, and actually quite slow.
		# Masked points are not shown at al
		
		plt.figure(figsize=(12,8))	# sets figure size
		axes = plt.gca()	
		
		# The points
		plt.plot(lc1jds, lc1mags, ".", color=lc1.plotcolour, label=lc1.object)
		plt.plot(lc2jds, lc2mags, ".", color=lc2.plotcolour, label=lc2.object)
		
		# Lines between lc2 points
		plt.plot(lc2jds, lc2mags, "-", color=lc2.plotcolour, label=lc2.object) # lines between non-masked lc2 points (= interpolation !)
		
		# Circles around non-used lc1 points
		plt.plot(lc1jds[lc1jdsmask == False], lc1mags[lc1jdsmask == False], linestyle="None", marker="o", markersize=8., markeredgecolor="black", markerfacecolor="None", color="black") # cicles around maked point
		
		# Dispersion "sticks"
		for (jd, mag, magdiff) in zip(lc1jds[lc1jdsmask], lc1mags[lc1jdsmask], okdiffs):
			plt.plot([jd, jd], [mag, mag + magdiff], linestyle=":", color="grey")

		
		# Something for astronomers only : we invert the y axis direction !
		axes.set_ylim(axes.get_ylim()[::-1])
	
		# And we make a title for that combination of lightcurves :
		#plt.title("Lightcurves", fontsize=20)
		plt.xlabel("Days", fontsize=16)
		plt.ylabel("Magnitude", fontsize=16)
		plt.title("%i interpolations" % (n), fontsize=16)
		plt.show()

	
	return {'n': n, 'd2' : d2}




def pelt95(lc1, lc2, decorlength=3.0, usemask=True, verbose=False, plot=False):
	"""
	Dispersion method as described in Pelt Kayser Refsdal 93 & 95
	For now this corresponds to D3, the third form, i.e. using only stricly neigbooring pairs,
	and no weighting.
	The fourth form uses *all* pairs within a given decorlength, non only neighboring.
	"""
	if usemask :
		lc1jds = lc1.getjds()[lc1.mask]
		lc2jds = lc2.getjds()[lc2.mask]
		lc1mags = lc1.getmags()[lc1.mask]
		lc2mags = lc2.getmags()[lc2.mask]
		lc1magerrs = lc1.magerrs[lc1.mask]
		lc2magerrs = lc2.magerrs[lc2.mask]
	else:
		lc1jds = lc1.getjds()
		lc2jds = lc2.getjds()
		lc1mags = lc1.getmags()
		lc2mags = lc2.getmags()
		lc1magerrs = lc1.magerrs
		lc2magerrs = lc2.magerrs

	
	#lc1jds = np.array([1.0, 2.0, 3.0])
	#lc1mags = np.array([5.1, 5.2, 4.9])
	#lc1magerrs = np.array([0.1, 0.1, 0.1])
	
	#lc2jds = lc1jds + 0.2
	#lc2mags = lc1mags + 1.0
	#lc2magerrs = lc1magerrs.copy()
	
	# Again, we want to avoid any loops !
	# So first we merge the lightcurves into a numpy array, introducing a last column that tracks the "origin" of each point.
	
	pts1 = np.column_stack((lc1jds, lc1mags, lc1magerrs, np.zeros(len(lc1jds))))
	pts2 = np.column_stack((lc2jds, lc2mags, lc2magerrs, np.ones(len(lc2jds))))
	pts = np.concatenate((pts1, pts2))
	
	# We sort this stuff according to the jds :
	sortedindices = np.argsort(pts[:,0])
	pts = pts[sortedindices]
	
	# So now we have a big array of mixed points, with columns :
	# jds	mags	magerrs		1/0(dpending on origin)
	
	# Now we calculate all diffs between lines k+1 and k
	diffs = pts[1:] - pts[:-1] # diffs = np.ediff1d(pts) # does not work... so we do it the hard way.
	# The trick is that we can infer good pairs from jd diffs and 1/0 diffs...
	
	# We do not care about the diffs of magnitude errors... instead we want to sum their squares
	sqrerrs = np.square(pts[1:,2]) + np.square(pts[:-1,2])
	#We can introduce this into the diffs array, in place of the error differences :
	diffs[:,2] = sqrerrs
	
	# From this table, we want to sum the square mag differences if 
	# col 0 < decorlength
	# abs(col 3) > 0.1 (meaning that the points are not from the same curve) 
	
	"""
	# The slow variant...
	goodstuff = []
	weights = []
	
	for line in diffs:
		if line[0] < decorlength and abs(line[3]) > 0.1:
			element = line[1]*line[1] / line[2]
			goodstuff.append(element)
			weights.append(2.0 * 1.0/line[2]) # This is from Pelts article.
			
	n = len(goodstuff)
	d2 = np.sum(np.array(goodstuff)) / np.sum(np.array(weights))
	#d2 = np.sum(np.array(goodstuff)) / (n-1.0)
	"""
	
	# The faster way... selection of pairs :
	pairs = diffs[np.logical_and(diffs[:,0] < decorlength, np.abs(diffs[:,3]) > 0.1)]
	# And we calculate the sum terms :
	elements = np.square(pairs[:,1]) / pairs[:,2]
	
	n = len(elements)
	#d2 = np.sum(elements) / (n-1.0)
	d2 = np.sum(elements)/(2.0 * np.sum(1.0/pairs[:,2])) # again, this is Pelt.
	# (dispersion estimator : sum a subset of half square-diffs)
	
	
	if plot:
		plt.figure(figsize=(12,8))	# sets figure size
		axes = plt.gca()
		
		# The points	
		plt.plot(lc1jds, lc1mags, ".", color=lc1.plotcolour, label=lc1.object)
		plt.plot(lc2jds, lc2mags, ".", color=lc2.plotcolour, label=lc2.object)
	
		for i,line in enumerate(diffs):
			if line[0] < decorlength and abs(line[3]) > 0.1:
				plt.plot([pts[i,0], pts[i,0] + line[0]], [pts[i,1], pts[i,1] + line[1]], linestyle="-", color="gray")
				plt.plot([pts[i,0] + line[0]/2.0, pts[i,0] + line[0]/2.0], [pts[i,1], pts[i,1] + line[1]], linestyle="-", color="black")
	
		axes.set_ylim(axes.get_ylim()[::-1])	
		plt.xlabel("Days", fontsize=16)
		plt.ylabel("Magnitude", fontsize=16)
		plt.title("%i pairs" % (n), fontsize=16)
		plt.show()

	return {"d2": d2, "n": n}		
	


def linint90(lc1, lc2, interpdist = 30.0, weights = True, verbose = False, plot = False, debug = False):
	"""
	Do not use, historical only.
	Terribly slow variant of linintnp, with plenty of loops, mimics my f90 implementation of 2007...
	
	Between two ORDERED light curves lc1 and lc2
	It does skip masked points.
	
	IT IS SLOW !!! THIS IS FOR DEMONSTRATION PURPOSES ONLY !
	
	About the maths : interpolate linearly between the points of lc2, and sum the square of the differences 
	between the points of
	lc1 and this interpolation. If you choose to weight this sum according to errors, 
	then for each point-pair (i.e. "element")
	of the sum, the weight is calculated using both the error on lc1 and the interpolated (again !) error on lc2
	
	We try to do a minimum of loops, using the fact that the light curves are ordered.
	
	lc1 will be the "reference", and we will interpolate between the points of lc2
	This algorithm absolutely NEEDS ordered light curves !!!!
	
	The way we check for the mask looks a bit un-natural, but this way should be quite fast in fact.
	Another think to note is that there is for now no special procedure for "close-by" points,
	and subtraction of very close numbers is not good...

	
	@return: a dict, containing for now :
	
		- "n" 	: the number of "pairs" that could be used
		- "d2" 	: the chi2 ( i.e. sum(squares)/(n-1) )
	
	@type lc1: lightcurve
	@param lc1: The reference light curve. It can of course be shifted and or partially masked.
	@type lc2: lightcurve
	@param lc2: The light curve that will be interpolated. Same remark.
	
	@type interpdist: float
	@param interpdist: a maximum timespan (in days) over which interpolations are done. See the code for exact details.
	If you are too conservative (low value), you might perhaps not use some isolated points at all...
	@type weights: boolean
	@param weights: True if you want to use the errorbars of the lightcurves
	@type verbose: boolean
	@param verbose: True if you want to see details about the curve matching
	@type plot: boolean
	@param plot: True if you want to see the actual interpolation as a matplotlib plot
	
	@type debug: boolean
	@param debug: Turn this on to make it even slower...
	
	"""
	
	if plot: debug = True
	
	interpjds = []		# for each interpolated point, the jd
	interpmags = []		# the corresponding mag of the reference (i.e. non-interpolated) curve
	interpmagdiffs = []	# the difference between this reference curve and the interpolation
	interpmagdifferrs = []	# some kind of error for this difference
	
	if debug: lc1notused = [] # We fill this to get a mask similar to the one used by faster implementations : this is slow, but for debugging...
	
	lc1jds = lc1.getjds()	# we make a local copy of them for fastest access.
	lc1mags = lc1.getmags()
	
	
	# To skip the masked points of lc2, the simplest is to remove them, as the algo wants consecutive points.
	# It turns out that this is really not the bottleneck, so you can safely keep all this.
	lc2masked = lc2.copy()
	if lc2masked.ml != None:
		lc2masked.applyml() # We apply the microlensing as we want to use cutmask, which would delete it.
	lc2masked.cutmask()
	lc2jds = lc2masked.getjds()
	lc2mags = lc2masked.getmags()
	
	# From here on there is no more access to the jds and mags of the lcs.
	
	
	startj = 0
	
	# loop over all points of the reference curve
	for i in range(len(lc1jds)):
		if verbose : print "### Reference point %i ###" % i
		if not lc1.mask[i] : 
			if verbose : print "masked"
			if debug: lc1notused.append(True)
			continue # we go on with the next i
		
		# we try to find lc2's closest left and right neighbors to lc1[i].
		# we look for the index j that is the left one, and j + 1 the right one.
		
		for j in range(startj, len(lc2jds)-1):
			if verbose : print "Testing j = %3i" %j
			
			if (lc2jds[j] > lc1jds[i]): # then no further j will match -> go to next i, keeping the same startj for the next i.
				if verbose : print "j = %3i is on the right side of %3i = " % (j, i)
				if debug: lc1notused.append(True)
				break
			
			
			# No, this is not that easy, thats why we already removed the masked points before.
			#if (not lc2.mask[j]) or (not lc2.mask[j+1]) : 
			#	if verbose : print "masked"
			#	continue # we go on with the next j
			
			# we check if they are left and right from our candidate :
			if (lc2jds[j] <= lc1jds[i]) and (lc2jds[j+1] >= lc1jds[i]) :
				if verbose : print "j = %3i possible : %.2f < [%.2f] < %.2f" % (j, lc2jds[j], lc1jds[i], lc2jds[j+1]) # This should be true
				# And now if they are close enough for an interpolation :
			
				#if (lc1jds[i] - lc2jds[j] < interpdist) and (lc2jds[j+1] - lc1jds[i] < interpdist) : # this was the initial, but a bit strange, way of doing it
				if (lc2jds[j+1] - lc2jds[j] < interpdist) : # I prefer this one, but it changes a bit the meaning of interpdist

					if verbose : print "\tOK ! "
					
					interpjds.append(lc1jds[i])	# we add the date of the point we found
					interpmags.append(lc1mags[i])	# handy for the plot
					
					# and now we do the actual interpolation					
					#interpmagdiffs.append( (lc2masked.mags[j] + lc2masked.mags[j+1])/2.0 - lc1.mags[i])	# this is just taking a mean... too simplistic
					interpmagdiff = lc2mags[j] + (lc2mags[j+1]-lc2mags[j])*(lc1jds[i]-lc2jds[j])/(lc2jds[j+1]-lc2jds[j]) - lc1mags[i]	# that's a linear interpolation
					interpmagdiffs.append(interpmagdiff)
					
					if weights : # then we have to interpolate the error bars, and calculate a weight for our "interpmagdiff"
						lc1error = lc1.magerrs[i] # that's the easy part
						lc2error = lc2masked.magerrs[j] + (lc2masked.magerrs[j+1]-lc2masked.magerrs[j])*(lc1jds[i]-lc2jds[j])/(lc2jds[j+1]-lc2jds[j])	# linear interpolation
						
						# these lc1error and lc2error are to be treated like standard devs :
						interpmagdifferr = math.sqrt(lc1error * lc1error + lc2error * lc2error)
						interpmagdifferrs.append(interpmagdifferr)
						
					
					# We found the neighbors, we can thus get out of the j loop. For optimization, let's remember which j was successful;
					# for the next point, we can then start to search from this j on, as our lc are ordered !
					startj = j
					if debug: lc1notused.append(False)
					break # we get out of the j loop -> next i
					
				else :
					if verbose : print "\tNo, they are too far away."
					# In this case, there is no point looking further : the next j will be on the right of lc1[i]. So same as above :
					startj = j
					if debug: lc1notused.append(True)
					break # we get out of the j loop, 
					
	
			if verbose and j == len(lc2jds)-2: print "I did not find convenient neighbours."
			if debug and j == len(lc2jds)-2: lc1notused.append(True)
	
	n = len(interpmagdiffs)			
	#print "Ok, %i points interpolated." % n
	
	if debug:
		lc1notused = np.array(lc1notused)
		#print n, len(lc1notused)
	
	if plot:
		# Then we make a plot, essentially for testing purposes. Nothing important here !
		# This is not optimal, and actually quite slow.
		
		plt.figure(figsize=(12,8))	# sets figure size
		axes = plt.gca()
		
		# The points	
		plt.plot(lc1.getjds(), lc1.getmags(), ".", color=lc1.plotcolour, label=lc1.object)
		plt.plot(lc2.getjds(), lc2.getmags(), ".", color=lc2.plotcolour, label=lc2.object)
		#plt.errorbar(lc1.getjds(), lc1.getmags(), lc1.magerrs, fmt=".", color=lc1.plotcolour, ecolor="#BBBBBB", label=lc1.object) # the acutal points, with errorbars
		#plt.errorbar(lc2.getjds(), lc2.getmags(), lc2.magerrs, fmt=".", color=lc2.plotcolour, ecolor="#BBBBBB", label=lc2.object)
		
		# Lines between all lc2points
		plt.plot(lc2jds, lc2mags, "-", color=lc2masked.plotcolour, label=lc2masked.object)
		
		# Circles around masked points (masked in the lightcurves)
		#plt.plot(lc1.getjds()[lc1.mask == False], lc1.getmags()[lc1.mask == False], linestyle="None", marker="o", markersize=8., markeredgecolor="black", markerfacecolor="None", color="black") # cicles around maked point
		#plt.plot(lc2.getjds()[lc2.mask == False], lc2.getmags()[lc2.mask == False], linestyle="None", marker="o", markersize=8., markeredgecolor="black", markerfacecolor="None", color="black") # cicles around maked point

		#Circles around points of lc1 that cannot be used :
		if debug: plt.plot(lc1.getjds()[lc1notused], lc1.getmags()[lc1notused], linestyle="None", marker="o", markersize=8., markeredgecolor="black", markerfacecolor="None", color="black") # cicles around maked point
		
		# Dispersion "sticks"
		for (interpjd, interpmag, interpmagdiff) in zip(interpjds, interpmags, interpmagdiffs):
			plt.plot([interpjd, interpjd], [interpmag, interpmag + interpmagdiff], linestyle=":", color="grey") # lines between the original point and the interpolation.
				
		
#		Not up to date, should be rewritten
#		if lc1.showlabels:
#			for i, label in enumerate(lc1.labels):
#				if label != "":
#					axes.annotate(label, (lc1.jds[i], mlmagslc1[i]), xytext=(7, -6), textcoords='offset points',size=12, color = lc1.plotcolour)
#		if lc2.showlabels:
#			for i, label in enumerate(lc2.labels):
#				if label != "":
#					axes.annotate(label, (lc2.jds[i], mlmagslc2[i]), xytext=(7, -6), textcoords='offset points',size=12, color = lc2.plotcolour)	
#		
#		if weights:
#			for (interpjd, interpmag, interpmagdifferr) in zip(interpjds, interpmags, interpmagdifferrs):
#				yerrlabelpos = interpmag # + 0.5*interpmagdiff
#				axes.annotate("%.3f" % interpmagdifferr, (interpjd, yerrlabelpos), xytext=(7, 0), textcoords='offset points', size=12, color = "black")
#				
				
		# Something for astronomers only : we invert the y axis direction !
		axes.set_ylim(axes.get_ylim()[::-1])
	
		# And we make a title for that combination of lightcurves :
		#plt.title("Lightcurves", fontsize=20)
		plt.xlabel("Days", fontsize=16)
		plt.ylabel("Magnitude", fontsize=16)
		plt.title("%i interpolations" % (n), fontsize=16)
		plt.show()

	
	
	# Ok, and now we calculate a dispersion from our interpmagdiffs
	
	if n > 1 :
		if not weights: # simple case, trivial calculation :
			d2 = np.sum(np.array(interpmagdiffs)*np.array(interpmagdiffs))/(n)
		else:
			#interpmagdiffweights = (array(interpmagdifferrs).mean()/array(interpmagdifferrs)) # we normalize with the mean -> no effect if all errors are identical
			# Not sure if this is optimal ... Do we want something close to 1.0 like a chi2 or not ?
			# let's try to not normalize this :
			interpmagdiffweights = (1.0/np.array(interpmagdifferrs))
			
			d2 = np.sum((np.array(interpmagdiffs)*np.array(interpmagdiffs))*(interpmagdiffweights*interpmagdiffweights))/(n)
	
		return {'n': n, 'd2' : d2}
	else :
		return {'n': 0, 'd2' : 0.0}





def symmetrize(lc1, lc2, dispersionmethod):
	"""
	Calls your dispersion method on (A,B) and on (B,A) and returns the "average".
	"""
	res1 = dispersionmethod(lc1, lc2)
	res2 = dispersionmethod(lc2, lc1)
	#return {'n': 0.5*(res1['n']+res2['n']), 'd2' : 0.5*(res1['d2']+res2['d2'])}
	#return {'n': (res1['n']+res2['n']), 'd2' : 0.5*(res1['d2']+res2['d2'])}
	n = res1['n']+res2['n']
	d2 = float(res1['n']*res1['d2'] + res2['n']*res2['d2'])/float(n)
	return {'n': n, 'd2' : d2}
	
	
	


