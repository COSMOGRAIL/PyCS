"""
Functions to graphically explore the magnitude- and time-shift parameter space for any dispersion function of your choice.
"""


from numpy import *
import matplotlib.pyplot as plt
#import pycs.gen.ml
import pycs.pelt.dispersions as dispersions

import time	
import pycs.gen.lc as lc





def timeshiftspectrum(lc1, lc2, dispersionmethod, timewidth=30, timestep=1.0, mlopt=False, savefig=False):
	"""
	The most important plot... dispersion versus time shfit.
	Shows the simple 2D plot obtained by varying only the time shift (what I previously did in fortran).
	Magnitude shift is kept "as is" -- but you can choose to optimize ml of the two curves, this is one form of mag shift ...
	The plot is centered around the present time shifts, so shift your curves by hand to match more or less.
	
	In the code, we will shift lc2 around, but this has no importance : there is a perfect symmetry.
	So if your dispersionmethod is symmetric versus curve exchange, the timeshiftspectrum will be the same.
	
	@type mlopt: boolean
	@param mlopt: Do you want to optimize the microlensing at each step ? 
	"""

	
	lc1 = lc1.copy() # for the microlensing, just in case...
	lc2 = lc2.copy()
	
	timeshifts = arange(-(timewidth)*timestep/2.0, (timewidth+1)*timestep/2.0, timestep)
	print "Exploring timeshifts from %f to %f" % (timeshifts[0], timeshifts[-1])
	

	# these are the corresponding absolute shifts, the important time delay !
	timedelays = timeshifts + lc2.timeshift - lc1.timeshift
	print "... that is time delays from %f to %f"% (timedelays[0], timedelays[-1])
	
	
	def d2(timeshift):
		lc2shifted = lc2.copy()
		#lc2.resettimeshift()	# NO ! this would kill any previous time shift ! And in fact it does not speed up to reuse the paramerters.
		lc2shifted.shifttime(timeshift)
		#lc2.shifttime(timeshift)
		if mlopt:
			#n = n+1
			#print "%4i / %4i" % (n, timewidth)
			print "Additional shift : %f days." % timeshift	# Additional, as there might already be another one.
			dispersions.optimizeml(lc1, lc2shifted, dispersionmethod)

		res = dispersionmethod(lc1, lc2shifted)
		
		if savefig:
			filename = "%.0f.png" % (time.time()*100)
			#print filename
			lc.display([lc1, lc2shifted], showerrorbars=False, legendloc="upper left", silent=True, filename=filename, figsize=(9,6), message = "d2 = %12.6f, n = %4i" % (res["d2"], res["n"]))
						
		return ([res["d2"], res["n"]])
	
	#vecd2 = vectorize(d2, otypes=[float]) # this was for the d2-only output
	vecd2 = vectorize(d2, otypes=[ndarray])

	# We apply it to our shifts :
	#(d2values, nvalues) = vecd2(timeshifts)
	
	bla = vecd2(timeshifts)
	d2values = array([el[0] for el in bla])
	nvalues = array([el[1] for el in bla])
	#print d2values
	
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	
	ax1.plot(timedelays, d2values, "r.")
	ax1.set_ylabel("Dispersion d2", size=16, color='r')
	ax1.set_xlabel("Time delay of %s with respect to %s [days]" % (lc1.object, lc2.object), size=16)
	for tl in ax1.get_yticklabels():
    		tl.set_color('r')

	ax2 = ax1.twinx()
	ax2.plot(timedelays, nvalues, "b.")
	for tl in ax2.get_yticklabels():
    		tl.set_color('b')
	#ax2.set_ylabel("Number of terms")	# this causes a crash, probably a bug...
	
	
	#plt.figtext(0.13, 0.9, "Magnitude shift of %s with respect to %s : %6.4f" % (lc2.object, lc1.object, lc2.magshift - lc1.magshift))
	
	plt.figtext(0.14, 0.96, "%s" % (str(lc1)))
	plt.figtext(0.14, 0.93, "versus (mlopt = %s)" % (str(mlopt)))
	plt.figtext(0.14, 0.90, "%s" % (str(lc2)))
	
	#plt.imshow(Z, interpolation='nearest', extent = extent, aspect='auto')
	#plt.colorbar()
	
	#plt.title("%s versus %s, mlopt = %s" % (str(lc1), str(lc2), str(mlopt)), size=10)
	#plt.xlabel("Time delay of %s with respect to %s [days]" % (lc1.object, lc2.object), size=16)
	#plt.ylabel("Dispersion d2", size=16)
	
	plt.show()








def dispsurface(lc1, lc2, dispersionmethod, timewidth=20, timestep=1.0, magwidth=20, magstep=0.01):
	"""
	Shows the dispersion surface (sqrt of d2) obtained by varying both the time shift and magnitude shift.
	
	How to specify the region to explore ? Before passing lightcurves to this function, shift them so that they match more or less.
	Then, the parameters like timewidth and timestep indicate the region/step to explore around this central position. Nevertheless, 
	the axes labels on the graphics show I{absolute} shifts, i.e. between the original, unshifted light curves !
	
	I'm not sure if this is the way to go, but it seems rather convienient for now.
	
	@todo: Send to ds9, if it is possible to pass the "x" and "y" axis meanings.
	
	@return: None
		
	@type lc1: lightcurve
	@param lc1: A reference lightcurve that gets passed to the disperionmethod
	@type lc2: lightcurve
	@param lc2: the light curve that will be shifted, then interpolated by the dispersionmethod.
	This means that shifts in the output graph refer to shifts to apply to lc2 so that it matches lc1.
	
	@type dispersionmethod: function
	@param dispersionmethod: It has to be a function that takes only two lightcurves as arguments,
	lc1 and lc2, and returns a dictionnary containting a dispersion with key "d2" (that way it is quite generic).
	You can prepare such a function this using lambda, like here :
	
		>>> mydisp = lambda lc1, lc2 : dispersions.maltef90(lc1, lc2, interpdist = 30.0, weights=True)
		>>> specplots.dispsurface(mercator_a, mercator_b, mydisp)

	@type timewidth: int
	@param timewidth: number of samples you want along the time axis
	@type timestep: float
	@param timestep: step alont the time axis
	@type magwidth: int
	@param magwidth: number of samples you want along the magnitude axis
	@type magstep: float
	@param magstep: step alont the magnitude axis

		
	"""
	
	raise RuntimeError, "I need to be updated and tested."
	
	# these are the shifts to apply to lc2 :
	#timeshifts = arange(-timewidth*timestep/2.0, timewidth*timestep/2.0 + timestep, timestep)
	#magshifts = arange(-magwidth*magstep/2.0, magwidth*magstep/2.0 + magstep, magstep)
	
	timeshifts = arange(-(timewidth)*timestep/2.0, (timewidth+1)*timestep/2.0, timestep)
	magshifts = arange(-(magwidth)*magstep/2.0, (magwidth+1)*magstep/2.0, magstep)
	
	# these are the corresponding absolute shifts, for the axes labels to show up nicely :
	extent = (-timewidth*timestep/2.0 - timestep/2.0 + lc2.timeshift - lc1.timeshift, 
		timewidth*timestep/2.0 + timestep/2.0 + lc2.timeshift - lc1.timeshift,
		magwidth*magstep/2.0 + magstep/2.0 + lc2.magshift - lc1.magshift,
		-magwidth*magstep/2.0 - magstep/2.0 + lc2.magshift - lc1.magshift)
	

	X,Y = meshgrid(timeshifts, magshifts)
	
	def z(timeshift, magshift):
		lc2shifted = lc2.copy()
		lc2shifted.shifttime(timeshift)
		lc2shifted.shiftmag(magshift)
		#lc2shifted.jds =+ timeshift
		#lc2shifted.mags =+ magshift
		res = dispersionmethod(lc1, lc2shifted)
		return sqrt(res["d2"])
		#return res["n"]	# if you want to plot the number of used points...
	
	
	# Now we pimp our z function, so that it accepts numpy arrays as input, and also returns a numpy array (it will be applied "element-wise" to input arrays)
	vecz = vectorize(z, otypes=[float])
	
	# We apply it to our grid :
	Z = vecz(X,Y)
	
	#print Z.shape
	#print timeshifts
	#print magshifts
	
	#plt.pcolormesh(X, Y, Z)
	#plt.colorbar()
	
	plt.imshow(Z, interpolation='nearest', extent = extent, aspect='auto')
	plt.colorbar()
	
	plt.title("Dispersion d", size=18)
	plt.xlabel("Time delay of %s with respect to %s [days]" % (lc1.object, lc2.object), size=16)
	plt.ylabel("Magnitude shift of %s with respect to %s" % (lc2.object, lc1.object), size=16) # yes, this is right
	
	plt.show()


def dispsurfacetofits(lc1, lc2, dispersionmethod, filename="disp.fits", timewidth=20, timestep=1.0, magwidth=20, magstep=0.01):
	"""Very similar to dispsurface, but saves to a FITS file"""
	
	raise RuntimeError, "I need to be updated and tested."
	
	import pyfits
	
	timeshifts = arange(-(timewidth)*timestep/2.0, (timewidth+1)*timestep/2.0, timestep)
	magshifts = arange(-(magwidth)*magstep/2.0, (magwidth+1)*magstep/2.0, magstep)
	

	X,Y = meshgrid(timeshifts, magshifts)
	
	def z(timeshift, magshift):
		lc2shifted = lc2.copy()
		lc2shifted.shifttime(timeshift)
		lc2shifted.shiftmag(magshift)
		#lc2shifted.jds =+ timeshift
		#lc2shifted.mags =+ magshift
		res = dispersionmethod(lc1, lc2shifted)
		return sqrt(res["d2"])
		#return res["n"]	# if you want to plot the number of used points...
	
	# Now we pimp our z function, so that it accepts numpy arrays as input, and also returns a numpy array (it will be applied "element-wise" to input arrays)
	vecz = vectorize(z, otypes=[float])
	# We apply it to our grid :
	Z = vecz(X,Y)
	
	hdu = pyfits.PrimaryHDU(Z.transpose()) # The transpose is to give a coherent orientation with DS9.

	hdu.header.update("WCSNAME", "Chi square")	# Name for this WCS (you could make multiple WCSs... -> see article)

	hdu.header.update("CUNIT1", "Day")	# unit of axis (max 3 characters)
	hdu.header.update("CRVAL1", 0.0)	# value at reference position
	hdu.header.update("CRPIX1", 50.0)	# position of reference (in pixels, as float)
	hdu.header.update("CDELT1", 0.1)	# step per pixel
	hdu.header.update("CTYPE1", "Time    ")	# This has to be 8 characters long !


	hdu.header.update("CUNIT2", "Mag")
	hdu.header.update("CRVAL2", 0.0)
	hdu.header.update("CRPIX2", 50.0)
	hdu.header.update("CDELT2", 0.01)
	hdu.header.update("CTYPE2", "Magshift")

	if os.path.isfile(filename):
		os.remove(filename)
	hdu.writeto(filename)

	


