from pycs.gen import *

import numpy as np
#import math
import matplotlib.pyplot as plt
#import scipy.optimize as spopt
#import scipy.interpolate as ip



def chi2(lcs, plot=False):
	"""
	lcs : a list of lightcurve objects
	"""
	
	# We merge the lcs into three numpy arrays
	jdslist = []
	magslist = []
	magerrslist = []		
	for thislc in lcs:
		jdslist.append(thislc.getjds()[thislc.mask])
		magslist.append(thislc.getmags()[thislc.mask])
		magerrslist.append(thislc.magerrs[thislc.mask])
		
	mergedjds = np.concatenate(jdslist)
	mergedmags = np.concatenate(magslist)
	mergedmagerrs = np.concatenate(magerrslist)

	# Now the sorting ...
	sortedindices = np.argsort(mergedjds)
	jds = mergedjds[sortedindices]
	mags = mergedmags[sortedindices]
	errs = mergedmagerrs[sortedindices]
	
	plt.plot(jds, mags, "b.")
	
	# Something for astronomers only : we invert the y axis direction !
	axes = plt.gca()
	axes.set_ylim(axes.get_ylim()[::-1])

	plt.xlabel("Days")
	plt.ylabel("Magnitude")
	plt.show()	
	
	
	