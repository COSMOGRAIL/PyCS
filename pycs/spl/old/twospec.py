"""

"""

import numpy as np
import matplotlib.pyplot as plt
import pycs.spl.multiopt

def calcr2s(lc1, lc2, reltimeshifts, spline, trace=False):
	"""
	I calcuate the r2 for an array of relative time shifts.
	To be compated to calcd2 of pycs.pelt.twospec !
	"""
	
	lc2abstimeshifts = reltimeshifts + lc2.timeshift

	def r2(lc2abstimeshift):
		
		# We work with copies at every trial time delay, to always start from the same position.
		mylc1 = lc1.copy()
		mylc2 = lc2.copy()
		mylc2.timeshift = lc2abstimeshift
		myspline = spline.copy()
		
		return pycs.spl.multiopt.opt_source([mylc1, mylc2], myspline, verbose=False, trace=trace)

	# We vectorize this before applying it to our abstimeshifts
	vecr2 = np.vectorize(r2, otypes=[np.ndarray])
	r2s = vecr2(lc2abstimeshifts)
		
	return r2s
	
	 

def plot(lc1, lc2, timewidth=30, timestep=1.0, optml=False):
	"""
	Interactively plot a timeshift spectrum.
	"""
	
	# We prepare the timeshifts to explore
	timeshifts = np.arange(-(timewidth)*timestep/2.0, (timewidth+1)*timestep/2.0, timestep)
	print "Exploring timeshifts from %f to %f" % (timeshifts[0], timeshifts[-1])
		
	# This corresponds to the following actual time delays :
	timedelays = timeshifts + lc2.timeshift - lc1.timeshift
	print "... that is time delays from %f to %f"% (timedelays[0], timedelays[-1])
	
	# We calculate the spectra :
	spectrum = timeshift(lc1, lc2, timeshifts)
	
	# And show all this in a plot :	
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	
	ax1.plot(timedelays, spectrum, "r.")
	ax1.set_ylabel("chi2")
	ax1.set_xlabel("Time delay of %s with respect to %s [days]" % (lc1.object, lc2.object))
	#for tl in ax1.get_yticklabels():
    	#	tl.set_color('r')

	#ax2 = ax1.twinx()
	#ax2.plot(timedelays, spectrum["n"], "b.")
	#for tl in ax2.get_yticklabels():
    	#	tl.set_color('b')
	#ax2.set_ylabel("Number of terms")	# this causes a crash, probably a bug...
	
	#plt.figtext(0.13, 0.9, "Magnitude shift of %s with respect to %s : %6.4f" % (lc2.object, lc1.object, lc2.magshift - lc1.magshift))
	
	plt.figtext(0.14, 0.96, "%s" % (str(lc1)))
	plt.figtext(0.14, 0.93, "versus")
	plt.figtext(0.14, 0.90, "%s" % (str(lc2)))
	
	plt.show()

