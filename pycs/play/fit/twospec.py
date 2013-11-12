# """
# Functions to test ranges of time-delays
# 
# 
# """
# 
# import numpy as np
# import matplotlib.pyplot as plt
# #import pycs.gen.ml
# #import pycs.pelt.dispersions as dispersions
# 
# #import time	
# import pycs.gen.util as util
# import pycs.gen.lc as lc
# #from pycs.pelt import twoopt
# 
# 
# 
# def timeshift(lc1, lc2, timeshifts, fitmethod):
# 	"""
# 	The under-the-hood function to calculate dispersion spectra.
# 	Do not use this directly.
# 	
# 	lc1 and lc2 are two lightcurve objects. lc2 will be shifted by timeshifts (an array) in a relative way to it's present shift.
# 		
# 	We return an array of the dispersion values.
# 	"""
# 
# 	lc2origtimeshift = lc2.timeshift # to restore this later ...
# 
# 	abstimeshifts = timeshifts + lc2origtimeshift
# 
# 	# We define how to calculate a dispersion
# 	def chi2(abstimeshift):
# 		
# 		lc2.timeshift = abstimeshift	 # The choice of shifting lc2 or lc1 should make no fundamental difference
# 		#if optml:
# 		#	twoopt.opt_ml(lc1, lc2, dispersionmethod)
# 			
# 		#if fastret:
# 		return fitmethod([lc1, lc2])["chi2n"]
# 		#else:
# 		#	res = dispersionmethod(lc1, lc2)
# 		#	return ([res["d2"], res["n"]]) # so this is tuple
# 
# 	# We vectorize this before applying it to our abstimeshifts
# 	vecd2 = np.vectorize(chi2, otypes=[np.ndarray])
# 	
# 	#if fastret:
# 	retstuff = vecd2(abstimeshifts)
# 				
# 	#else:
# 	#	bla = vecd2(abstimeshifts)
# 	#	d2values = np.array([el[0] for el in bla])
# 	#	nvalues = np.array([el[1] for el in bla])
# 	#	retstuff = {"d2":d2values, "n":nvalues}
# 	
# 	# We reset lc2 at least to the same timeshift (but note that microlensing might have changed).
# 	lc2.timeshift = lc2origtimeshift
# 	
# 	return retstuff
# 	
# 	 
# 
# def plot(lc1, lc2, fitmethod, timewidth=30, timestep=1.0, optml=False):
# 	"""
# 	Interactively plot a timeshift spectrum.
# 	"""
# 	
# 	# We prepare the timeshifts to explore
# 	timeshifts = np.arange(-(timewidth)*timestep/2.0, (timewidth+1)*timestep/2.0, timestep)
# 	print "Exploring timeshifts from %f to %f" % (timeshifts[0], timeshifts[-1])
# 		
# 	# This corresponds to the following actual time delays :
# 	timedelays = timeshifts + lc2.timeshift - lc1.timeshift
# 	print "... that is time delays from %f to %f"% (timedelays[0], timedelays[-1])
# 	
# 	# We calculate the spectra :
# 	spectrum = timeshift(lc1, lc2, timeshifts, fitmethod)
# 	
# 	# And show all this in a plot :	
# 	fig = plt.figure()
# 	ax1 = fig.add_subplot(111)
# 	
# 	ax1.plot(timedelays, spectrum, "r.")
# 	ax1.set_ylabel("Chi2", size=16, color='r')
# 	ax1.set_xlabel("Time delay of %s with respect to %s [days]" % (lc1.object, lc2.object), size=16)
# 	#for tl in ax1.get_yticklabels():
#     	#	tl.set_color('r')
# 
# 	#ax2 = ax1.twinx()
# 	#ax2.plot(timedelays, spectrum["n"], "b.")
# 	#for tl in ax2.get_yticklabels():
#     	#	tl.set_color('b')
# 	#ax2.set_ylabel("Number of terms")	# this causes a crash, probably a bug...
# 	
# 	#plt.figtext(0.13, 0.9, "Magnitude shift of %s with respect to %s : %6.4f" % (lc2.object, lc1.object, lc2.magshift - lc1.magshift))
# 	
# 	plt.figtext(0.14, 0.96, "%s" % (str(lc1)))
# 	plt.figtext(0.14, 0.93, "versus (optml = %s)" % (str(optml)))
# 	plt.figtext(0.14, 0.90, "%s" % (str(lc2)))
# 	
# 	plt.show()
# 
# 
# 	
# #def boot(lc1, lc2, dispersionmethod, n=20, timewidth=30, timestep=1.0, optml=False, filename="boot.pkl"):
# #	"""
# #	"Stacks" dispersion spectra obtained by repeatedly bootstrapping the lightcurves.
# #	
# #	It turns out that this is not usefull.
# #	
# #	optml : if false, ml is not optimized at all, if true, it is done for every bootstrap (slow !)
# #	It's a good idea to optimize the ml before calling this technique, then you can
# #	reasonnably choose optml = False.
# #
# #	"""
# #
# #	lc1c = lc1.copy() # We keep the originals unchanged, as we will apply the microlensing if optml=False
# #	lc2c = lc2.copy() # Calm down, this is only called once.
# #	
# #	# The timeshifts to explore.
# #	timeshifts = np.arange(-(timewidth)*timestep/2.0, (timewidth+1)*timestep/2.0, timestep)
# #	print "Exploring timeshifts from %f to %f" % (timeshifts[0], timeshifts[-1])
# #		
# #	# This corresponds to the following actual time delays :
# #	timedelays = timeshifts + lc2c.timeshift - lc1c.timeshift
# #	print "... that is time delays from %f to %f"% (timedelays[0], timedelays[-1])
# #		
# #	if optml == False: # then we apply ml :
# #		if lc2c.ml != None:
# #			lc2c.applyml()	# This will make the rest a lot faster, as we avoid evaluating the polynom every time !
# #					# (recall that it deletes the ml object, thus no more evaluation)
# #			print "Applied ml to a copy of lc2"
# #		if lc1c.ml != None:
# #			lc1c.applyml()
# #			print "Applied ml to a copy of lc1"
# #	
# #	# At this stage we backup the mags, to restore them after each bootstrapping.
# #	lc1corigmags = lc1c.mags.copy()
# #	lc2corigmags = lc2c.mags.copy()
# #	
# #	d2array = []
# #
# #	for b in range(n):
# #		
# #		
# #		print "Bootstrap %4i / %4i" % (b+1, n)
# #		lc2c.bootstrap()
# #		lc1c.bootstrap()
# #
# #		d2row = timeshift(lc1c, lc2c, timeshifts, dispersionmethod, fastret=True, optml=optml)
# #		d2array.append(d2row) 
# #		
# #		# Crucial ... we revert to the original mags :
# #		lc1c.mags = lc1corigmags.copy()
# #		lc2c.mags = lc2corigmags.copy()
# #	
# #	
# #	d2array = np.column_stack(d2array)
# #	d2means = np.mean(d2array, axis=1)
# #	d2vars = np.var(d2array, axis=1)
# #	
# #	# bug ? don't ask, sqrt or std do not work as they should here ...
# #	d2stds = np.zeros(len(d2vars))
# #	for i in range(len(d2vars)):
# #		d2stds[i] = np.sqrt(d2vars[i])
# #	
# #	# Ok, now let's do this once with the non-bootstrapped curves :
# #	
# #	d2origs = timeshift(lc1c, lc2c, timeshifts, dispersionmethod, fastret=True, optml=optml)
# #	
# #	util.writepickle({'lc1':lc1, 'lc2':lc2, 'n':n, 'optml':optml, 'timedelays':timedelays, 'd2means':d2means, 'd2stds':d2stds, 'd2origs':d2origs}, filename)
# #	
# #
# #def bootplot(filename="boot.pkl"):
# #	
# #	filedict = util.readpickle(filename)
# #	
# #	lc1 = filedict['lc1']
# #	lc2 = filedict['lc2']
# #	optml = filedict['optml']
# #	
# #	print "Plotting %i bootstrap dispersions." % filedict['n']
# #	
# #	fig = plt.figure()
# #	ax = fig.add_subplot(111)
# #	
# #	plt.plot(filedict['timedelays'], filedict['d2origs'], "b.")
# #	
# #	plt.errorbar(filedict['timedelays'], filedict['d2means'], filedict['d2stds'], fmt=".", color="red", ecolor="#BBBBBB")
# #	
# #	plt.figtext(0.14, 0.96, "%s" % (str(lc1)))
# #	plt.figtext(0.14, 0.93, "versus (optml = %s)" % (str(optml)))
# #	plt.figtext(0.14, 0.90, "%s" % (str(lc2)))
# #		
# #	ax.set_ylabel("Dispersion d2", size=16)
# #	ax.set_xlabel("Time delay of %s with respect to %s [days]" % (lc1.object, lc2.object), size=16)
# #	
# #	plt.show()
# #	
# 
# 
# #
# #def time(lc1, lc2, dispersionmethod, timewidth=30, timestep=1.0, mlopt=False, savefig=False):
# #	"""
# #	The most important plot... dispersion versus time shfit.
# #	Shows the simple 2D plot obtained by varying only the time shift (what I previously did in fortran).
# #	Magnitude shift is kept "as is" -- but you can choose to optimize ml of the two curves, this is one form of mag shift ...
# #	The plot is centered around the present time shifts, so shift your curves by hand to match more or less.
# #	
# #	In the code, we will shift lc2 around, but this has no importance : there is a perfect symmetry.
# #	So if your dispersionmethod is symmetric versus curve exchange, the timeshiftspectrum will be the same.
# #	
# #	@type mlopt: boolean
# #	@param mlopt: Do you want to optimize the microlensing at each step ? 
# #	"""
# #
# #	
# #	lc1 = lc1.copy() # for the microlensing, just in case...
# #	lc2 = lc2.copy()
# #	
# #	timeshifts = arange(-(timewidth)*timestep/2.0, (timewidth+1)*timestep/2.0, timestep)
# #	print "Exploring timeshifts from %f to %f" % (timeshifts[0], timeshifts[-1])
# #	
# #
# #	# these are the corresponding absolute shifts, the important time delay !
# #	timedelays = timeshifts + lc2.timeshift - lc1.timeshift
# #	print "... that is time delays from %f to %f"% (timedelays[0], timedelays[-1])
# #	
# #	
# #	def d2(timeshift):
# #		lc2shifted = lc2.copy()
# #		#lc2.resettimeshift()	# NO ! this would kill any previous time shift ! And in fact it does not speed up to reuse the paramerters.
# #		lc2shifted.shifttime(timeshift)
# #		#lc2.shifttime(timeshift)
# #		if mlopt:
# #			#n = n+1
# #			#print "%4i / %4i" % (n, timewidth)
# #			print "Additional shift : %f days." % timeshift	# Additional, as there might already be another one.
# #			dispersions.optimizeml(lc1, lc2shifted, dispersionmethod)
# #
# #		res = dispersionmethod(lc1, lc2shifted)
# #		
# #		if savefig:
# #			filename = "%.0f.png" % (time.time()*100)
# #			#print filename
# #			lc.display([lc1, lc2shifted], showerrorbars=False, legendloc="upper left", silent=True, filename=filename, figsize=(9,6), message = "d2 = %12.6f, n = %4i" % (res["d2"], res["n"]))
# #						
# #		return ([res["d2"], res["n"]])
# #	
# #	#vecd2 = vectorize(d2, otypes=[float]) # this was for the d2-only output
# #	vecd2 = vectorize(d2, otypes=[ndarray])
# #
# #	# We apply it to our shifts :
# #	#(d2values, nvalues) = vecd2(timeshifts)
# #	
# #	bla = vecd2(timeshifts)
# #	d2values = array([el[0] for el in bla])
# #	nvalues = array([el[1] for el in bla])
# #	#print d2values
# #	
# #	fig = plt.figure()
# #	ax1 = fig.add_subplot(111)
# #	
# #	ax1.plot(timedelays, d2values, "r.")
# #	ax1.set_ylabel("Dispersion d2", size=16, color='r')
# #	ax1.set_xlabel("Time delay of %s with respect to %s [days]" % (lc1.object, lc2.object), size=16)
# #	for tl in ax1.get_yticklabels():
# #    		tl.set_color('r')
# #
# #	ax2 = ax1.twinx()
# #	ax2.plot(timedelays, nvalues, "b.")
# #	for tl in ax2.get_yticklabels():
# #    		tl.set_color('b')
# #	#ax2.set_ylabel("Number of terms")	# this causes a crash, probably a bug...
# #	
# #	
# #	#plt.figtext(0.13, 0.9, "Magnitude shift of %s with respect to %s : %6.4f" % (lc2.object, lc1.object, lc2.magshift - lc1.magshift))
# #	
# #	plt.figtext(0.14, 0.96, "%s" % (str(lc1)))
# #	plt.figtext(0.14, 0.93, "versus (mlopt = %s)" % (str(mlopt)))
# #	plt.figtext(0.14, 0.90, "%s" % (str(lc2)))
# #	
# #	#plt.imshow(Z, interpolation='nearest', extent = extent, aspect='auto')
# #	#plt.colorbar()
# #	
# #	#plt.title("%s versus %s, mlopt = %s" % (str(lc1), str(lc2), str(mlopt)), size=10)
# #	#plt.xlabel("Time delay of %s with respect to %s [days]" % (lc1.object, lc2.object), size=16)
# #	#plt.ylabel("Dispersion d2", size=16)
# #	
# #	plt.show()
# #	#
# #
# 
# 
# 
