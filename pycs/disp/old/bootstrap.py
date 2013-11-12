#"""
#
#We make dispersion curves like in specplots, but this time mean them over n bootstrapped evaluations
#
#
#Ok all this stuff needs to be "merged" with the specplots stuff...
#
#"""
#
#import sys
#
#from pycs.gen import *
#import dispersions
#
##from pycs.pelt import specplots
#
#
#
#import numpy as np
#import matplotlib.pyplot as plt
##import scipy.optimize as spopt
#
#
#
#def boot(lc1, lc2, dispersionmethod, n=20, timewidth=30, timestep=1.0, mlopt=False, filename="boot.pkl"):
#	"""
#	
#	"Stacks" dispersion spectra obtained by repeatedly bootstrapping the lightcurves. 
#	
#	mlopt : if false, ml is not optimized at all, if true, it is done for every bootstrap (slow !)
#	It's a good idea to optimize the ml before calling this technique, then you can
#	reasonnably choose mlopt = False.
#
#	"""
#
#	lc1c = lc1.copy() # We keep the originals unchanged, as we might apply the microlensing.
#	lc2c = lc2.copy() # Calm down, this is only called once.
#	
#	# Ok, we prepare the timeshifts to explore.
#	timeshifts = np.arange(-(timewidth)*timestep/2.0, (timewidth+1)*timestep/2.0, timestep)
#	abstimeshifts = lc2c.timeshift + timeshifts # these are the absolute values to which we will set lc2.
#	print "Exploring timeshifts from %f to %f" % (timeshifts[0], timeshifts[-1])
#		
#	# This corresponds to the following actual time delays :
#	timedelays = timeshifts + lc2c.timeshift - lc1c.timeshift
#	print "... that is time delays from %f to %f"% (timedelays[0], timedelays[-1])
#		
#	if mlopt == False: # then we apply ml :
#		if lc2c.ml != None:
#			lc2c.applyml()	# This will make the rest a lot faster, as we avoid evaluating the polynom every time !
#					# (recall that it deletes the ml object, thus no more evaluation)
#		if lc1c.ml != None:
#			lc1c.applyml()
#	
#	
#	# We define a dispersion, depending only on the abstimeshift of lc2 :
#	def d2(abstimeshift):
#		
#		lc2c.timeshift = abstimeshift	 	# The choice of shifting lc2 or lc1 makes no difference at all
#		if mlopt:
#			dispersions.optimizeml(lc1c, lc2c, dispersionmethod)
#				
#		return dispersionmethod(lc1c, lc2c)["d2"]
#
#	vecd2 = np.vectorize(d2, otypes=[np.ndarray])
#	
#	# Now the n-loop:
#
#	d2array = []
#	
#	lc1corigmags = lc1c.mags.copy()
#	lc2corigmags = lc2c.mags.copy()
#	
#	for b in range(n):
#		
#		print "Bootstrap %4i / %4i" % (b+1, n)
#		lc2c.bootstrap()
#		lc1c.bootstrap()
#		
#		#lc.display([lc1c, lc2c])
#		
#		d2row = vecd2(abstimeshifts)
#		d2array.append(d2row) 
#		
#		# Crucial ... we revert to the original mags :
#		lc2c.mags = lc1corigmags.copy()
#		lc1c.mags = lc2corigmags.copy()
#	
#	# and we revert lc2c to its original timeshift :
#	lc2c.timeshift = lc2.timeshift
#	
#	d2array = np.column_stack(d2array)
#	d2means = np.mean(d2array, axis=1)
#	d2vars = np.var(d2array, axis=1)
#	
#	# bug ? don't ask, sqrt or std do not work as they should here ...
#	d2stds = np.zeros(len(d2vars))
#	for i in range(len(d2vars)):
#		d2stds[i] = np.sqrt(d2vars[i])
#	
#	# Ok, now let's do this once with the non-bootstrapped curves :
#	
#	d2origs = vecd2(timeshifts)
#	
#	div.writepickle({'lc1':lc1, 'lc2':lc2, 'n':n, 'mlopt':mlopt, 'timedelays':timedelays, 'd2means':d2means, 'd2stds':d2stds, 'd2origs':d2origs}, filename)
#	
#
#
#
#def plot(filename="boot1.pkl"):
#	
#	filedict = div.readpickle(filename)
#	
#	lc1 = filedict['lc1']
#	lc2 = filedict['lc2']
#	mlopt = filedict['mlopt']
#	
#	print "Plotting %i bootstrap dispersions." % filedict['n']
#	
#	fig = plt.figure()
#	ax = fig.add_subplot(111)
#	
#	plt.plot(filedict['timedelays'], filedict['d2origs'], "b.")
#	
#	plt.errorbar(filedict['timedelays'], filedict['d2means'], filedict['d2stds'], fmt=".", color="red", ecolor="#BBBBBB")
#	
#	plt.figtext(0.14, 0.96, "%s" % (str(lc1)))
#	plt.figtext(0.14, 0.93, "versus (mlopt = %s)" % (str(mlopt)))
#	plt.figtext(0.14, 0.90, "%s" % (str(lc2)))
#		
#	ax.set_ylabel("Dispersion d2", size=16)
#	ax.set_xlabel("Time delay of %s with respect to %s [days]" % (lc1.object, lc2.object), size=16)
#	
#	plt.show()
#	
#	
#	
#
#
#	