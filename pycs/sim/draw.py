"""
This module generates fake curves by "drawing" them from common sourcesplines, adding noise, tweaking microlensing.
"""

import sys, os, time
import numpy as np

import pycs.gen.lc
import pycs.gen.util
import pycs.sim.src

import scipy.ndimage.filters


def sample(lc, spline):
	"""
	You give me a lightcurve, and I will set the lc.mags according to the spline.
	I DO modify the lc object in place !
	I do NOT modify the spline, of course.
	
	If your lightcurve has mircolensing, magshifts, or fluxshift, and of course also time shifts,
	I will set its mags so that it matches the spline with all of these applied !
	
	If you think about it, this is WONDERFUL !!!
	To do this, I have to do something like a reverse getmags().
	
	Note that I do not add any random noise here, this is sampling only.
	"""
	
	#lc.telescopename = "sim" 
	#lc.object += "sim"
	lc.telescopename += "sim" 
	lc.commentlist.append("Magnitudes simulated by source %s" % (spline))
	
	if lc.fluxshift != 0.0:
		if (lc.ml != None):
			lc.mags = spline.eval(lc.getjds()) - lc.ml.calcmlmags(lc) - lc.magshift
			lc.mags -= lc.calcfluxshiftmags(inverse=True) # Don't put this on the above line ... uses lc.mags !
		else:
			lc.mags = spline.eval(lc.getjds()) - lc.magshift
			lc.mags -= lc.calcfluxshiftmags(inverse=True)
	
	else:
		if (lc.ml != None):
			lc.mags = spline.eval(lc.getjds()) - lc.ml.calcmlmags(lc) - lc.magshift
		else:
			lc.mags = spline.eval(lc.getjds()) - lc.magshift	



def saveresiduals(lcs, spline):
	"""
	You give me some optimized lcs, and the fitting spline.
	This function saves the residuals of the optimized lightcurves (as an attribute), in terms of non-fluxshifted magnitudes.
	This means that you will have to ADD these residuals (in terms of magnitude) *before any fluxshifts*
	to the spline mags in order to recover the lc.
	This is made so that you can reuse this particular shotnoise once fake curves are drawn.
	
	We have to keep this function separate from draw(), as you might want to shift the curves after saving the residuals...
	
	.. warning:: Call this function if your spline matches to the lcs !
	
	"""
	
	
	for l in lcs:
		rawmags = l.mags.copy()
		samplelc = l.copy()
		sample(samplelc, spline)
		l.residuals = rawmags - samplelc.mags
		
		# Naive residuals (don't take into account the fluxshift problem) :
		#l.residuals = l.getmags() - spline.eval(jds = l.getjds()) 
		# So this is wrong, in case there is a fluxshift
		# Nice try, but this was wrong as well :
		"""
		# We want to compensate for the fluxshift, that "deformed" these residuals :
		# Easier to think in terms of flux...
		# The following is a bit too explicit, but hopefully clear.
		shiftedlcfluxs = 10**(-0.4 * lcmags)
		shiftedsplinefluxs = 10**(-0.4 * splinemags)
		nonshiftedlcfluxs = shiftedlcfluxs - l.fluxshift
		#nonshiftedsplinefluxs = shiftedsplinefluxs - l.fluxshift
		nonshiftedsplinefluxs = shiftedsplinefluxs
		nonshiftedlcmags = -2.5 * np.log10(nonshiftedlcfluxs)
		nonshiftedsplinemags = -2.5 * np.log10(nonshiftedsplinefluxs)
		
		l.residuals = nonshiftedlcmags - nonshiftedsplinemags
		"""
	
def transfershifts(lcs, reflcs, transferml=True):
	"""
	I will put (copies) of all the shifts and ML from the reflcs onto the lcs.
	I will check that you didn't mix up any curves.
	"""
	if len(lcs) != len(reflcs):
		raise RuntimeError("Not the same lengths !")
	
	for (l, refl) in zip(lcs, reflcs):
		if l.object != refl.object:
			raise RuntimeError("Seems that you mixed up the curves !")
		
		l.timeshift = refl.timeshift
		l.fluxshift = refl.fluxshift
		l.magshift = refl.magshift
		if transferml == True:
			if refl.ml == None:
				l.rmml()
			else:
				l.ml = refl.ml.copy()
		else:
			l.rmml() # Seems safe to do this, it cannot harm at least. Leaving the ML could be unwanted.
	
	
	

def draw(lcs, spline, shotnoise=None, shotnoisefrac=1.0, tweakml=None, scaletweakresi=True, tweakspl=None, keepshifts=True, keeptweakedml=False, trace=False, tracedir="draw"):
	"""
	Wrapper to produce one set of fake lightcurves that are similar to your lcs.
	Give me some lightcurves, that are optimized in fluxshift, magshift, microlensing, timeshift, and the resulting spline fit
	representing the intrinsic variations.
	I will tweak the source spline and the microlensing etc, and return to you a list of fake lightcurves.
	These lightcurves are just blank data points, "as observed". I will build them from scratch.
	But note that by default the drawn light curve will be equiped with the same shifts and ML as your lcs.
	
	.. note:: I do NOT modify lcs or spline !
	
	.. note:: It is perfectly ok to shift your input curves in time after the optimization,
		to make me build fake lightcurves with other delays !
		So if you want me to draw curves that look different from yours, it is perfectly ok to give me lcs that do not match to the spline !
		Same logic applies for all other shifts and microlensing.
	
	:param lcs: reference lightcurves to "immitate". I will use their epochs, their time/mag/flux shifts, their microlensing, their errorbars.
	
	:param spline: reference spline from which I will draw my magnitudes.
	
	:param shotnoise: Select among [None, "magerrs", "res", "mcres", "sigma"] 
		
		It tells what kind of noise to add to the fake mags.
		This noise will, as required, be added to the "observed" data (i.e. not fluxshifted).
	
	:param shotnoisefrac: a multiplier of the shotnoise added to the curve. Set to 0.5 and I'll add only "half" of the noise ...
	
	:param tweakml: either a function, or a list of functions, that takes a list of lightcurves and tweaks their ML in place.
		I will use this tweaked ML to draw them.
		If you give a single function, I will apply it to all the curves
		If you give a list of functions, I will apply them to the respective curves of your lcs (i.e., give the functions in the order corresponding to your lcs !).
	
	:param scaletweakresi: scales the "residuals" obtained by tweakml according to saved residuals
	
	:param tweakspl: a function that takes a spline and returns a tweaked spline. I will use this on
		the sourcespline you pass me, before drawing from it.
	
	:param keepshifts: by default I will set the time/flux/mag/ML shifts from your lcs also to the fake curves.
	
	
	:param keeptweakedml: if keepshifts is True, and keeptweakedml is True, I will keep the tweaked ML, not the input ML, on the 
		output curves.
	
	
	.. note:: I will tweak the ML only of those curves that have spline ML. You probably want me to tweak the ML of all your curves !
		So be sure to add some spline ML to all your curves before calling me !
	
	
	"""
	
	if isinstance(tweakml, list) :
		assert len(tweakml) == len(lcs)
	
	# We build a tweaked copy of the "intrinsic" spline. It will be common to all curves.
	if tweakspl != None:
		tweakedspline = tweakspl(spline)
	else:
		tweakedspline = spline # No need to make a copy
	
	
	# For the lcs, I'll not modify them in place, but return new ones.
	# Nevertheless, I work on a copy, as I'll call sample().
	copylcs = [l.copy() for l in lcs]
	fakelcs = []

	
	for i, l in enumerate(copylcs):

		# Before tweaking the microlensing, we make a copy (might be None)
		
		if l.ml != None:
			origml = l.ml.copy()
		else:
			origml = None
			
		# We now tweak the microlensing, before calling sample(). Cleaner this way, so that we can keep it.
		# Hmm, I guess we could make this afterwards as well (then apply it)
		
		if tweakml != None:
			# Then it should be either just a function, or a list of functions :
			if isinstance(tweakml, list) :
				thislctweakml = tweakml[i]
			else:
				thislctweakml = tweakml
			# We run the tweaking :
			if l.ml != None:
				thislctweakml([l]) # modifies in place
				tweakedml = l.ml.copy()
			else:
				print "WARNING: curve %s has no ML to tweak !" % (str(l))
				tweakedml = None
		else:
			tweakedml = origml	

		# And we sample the curve from the source, changing l in place.
		sample(l, tweakedspline)
		
		# Rescale the residuals ?

		# seems tempting here to do an "autoscaling" so that the sigma of the drawn residuals automatically matches the observed ones.
		# This is not possible, as we first need to redo the spline curve shifting to make a fair comparision !
		
		if (scaletweakresi == True) and  (l.ml != None): # Otherwise no need to do anything 
		
			# We sample from the same spline, but without tweaked ml :		
			lorigml = l.copy()
			lorigml.ml = origml
			# Sampling :
			sample(lorigml, tweakedspline)
			
			# Only the l.mags are different between l and lorigml.
			# Now we can compute the "residuals" that are due to the tweaked ml :
				
			tweakresis = l.mags - lorigml.mags
			#print np.std(tweakresis)


			# We scale them using previously saved residuals :
			if not hasattr(l, 'residuals'):
				raise RuntimeError("Save the residuals first !")
			
			scalings = np.fabs(l.residuals)
			# We normalise this a first time, to get values around 1.0 :
			scalings /= np.mean(scalings)
			#print np.mean(scalings)
			
			smoothscalings = scipy.ndimage.filters.median_filter(scalings, size=7, mode='constant', cval=1.0)
			# We normalise these again :
			smoothscalings /= np.mean(smoothscalings)
			#print np.mean(smoothscalings)
			
			"""
			import matplotlib.pyplot as plt
			plt.plot(l.jds, scalings)
			plt.plot(l.jds, smoothscalings, color="red")
			plt.show()
			exit()
			"""
			
			#print np.std(tweakresis*smoothscalings)
						
			l.mags = lorigml.mags + tweakresis*smoothscalings
			
		
		# And add purely uncorrelated noise on top of all this.
		# In future we should think about the fluxshifts here,
		# and not blindly add this according to the errorbars.
		
		if shotnoise == "magerrs":
			l.montecarlomags(f = shotnoisefrac)
		elif shotnoise == "none" or shotnoise == None:
			pass
		elif shotnoise == "res":
			if not hasattr(l, 'residuals'):
				raise RuntimeError("Save the residuals first !")
			l.mags += shotnoisefrac * l.residuals.copy()
			l.commentlist.append("Added previously saved residuals !")
		elif shotnoise == "mcres":
			# We use the residuals as 1 sigma of a gaussian to draw the shotnoise from.
			if not hasattr(l, 'residuals'):
				raise RuntimeError("Save the residuals first !")
			l.mags += l.residuals * shotnoisefrac * np.random.randn(len(l))
			l.commentlist.append("Monte Carlo with previously saved residuals as sigma !")
		elif shotnoise == "sigma":
			# We use the std of the residuals as 1-sigma amplitude of white noise to add.
			if not hasattr(l, 'residuals'):
				raise RuntimeError("Save the residuals first !")
			sigma = np.std(l.residuals)
			l.mags += sigma * shotnoisefrac * np.random.randn(len(l))
			l.commentlist.append("White noise with std of residuals as sigma !")
		else:
			raise RuntimeError("Couldn't understand your shotnoise.")
		
		
		
		
		# Usually we won't want to return the tweaked ML.
		# So let's build a brand new lightcurve.
		# So to be sure that this new one will "forget" all timeshifts, ML, fluxshits etc if we want so.
		
		jds = l.jds.copy() # We do not use getjds (we simulate "real data" taken at the same absolute epochs)
		mags = l.mags.copy() # Idem, we do not use getmags (as we might keep magshift, ML, etc)
		magerrs = l.magerrs.copy()
				
		#telescopename = "Fake(%s)" % (sl.telescopename)
		telescopename = l.telescopename
		object = l.object
		
		fakel = pycs.gen.lc.factory(jds, mags, magerrs, telescopename=telescopename, object=object, verbose=False)
		
		fakel.plotcolour = l.plotcolour
		#fakel.ploterrorbars = False
		
		# In any case we want to keep the true shifts of each curve stored under another attribute :
		fakel.truetimeshift = l.timeshift
	
		# And we might even want to keep the simulation shifts that make the points fit to the spline :
		if keepshifts:
			fakel.timeshift = l.timeshift
			fakel.fluxshift = l.fluxshift
			fakel.magshift = l.magshift
			if keeptweakedml == False:
				fakel.ml = origml # Yes, it's that simple ! fakel and l have the same jds, after all.
			else:
				fakel.ml = tweakedml
		
		fakelcs.append(fakel)
	
	# For trace, I save the tweaked ML and the tweaked spline.
	if trace:
		pycs.gen.util.trace(lclist=copylcs, splist=[tweakedspline], tracedir = tracedir)

	
	return fakelcs
	

def multidraw(lcs, spline=None, optfctnots=None, onlycopy=False, n=20, npkl=5, simset="draw", simdir=None, shotnoise=None, shotnoisefrac=1.0, truetsr=8.0, tweakml=None, scaletweakresi=True, tweakspl=None, shuffle=True, verbose=True, trace=False, destpath ='./'):
	"""
	Even higher wrapper to produce mock + tweaked lightcurves, and save them into a directory (as pickle files),
	in preparation for analysing them with :py:func:`pycs.sim.run.multirun`
	
	The curves I return are "like observed" : they have no shifts, no ML. Just datapoints.
	
	
	
	:param lcs: The starting-point lightcurves from which I will draw the mock data. They must match (!), i.e.,
		have appropriate microlensing and be optimized somehow.
	:type lcs: list
	
	:param spline: The source spline used to draw the new curves. Is not used (-> None) if onlycopy=True, or if optfct is specified.
	:type spline: spline
	
	:param optfctnots: A function to fit a spline and the ML + all shifts execpt for timeshifts. 
		It is called after setting the true time shifts. Put at None if you want to use always the same current ML and spline.
	:type optfctnots: function
	
	:param onlycopy: If True, I will simply save copies of the input lcs, not drawing anything.
	
	:param n: number of lightcurve-lists to simulate per pickle
	:param npkl: number of pickle files
	
	:param simset: give a name to your simulation !
	
	:param simdir: where should I put these simulations ?
	
	:param shotnoise: Select among None, "magerrs", "mcres", "res".
		See definitions in :py:func:`pycs.sim.draw.draw`. 
	:type shotnoise: string
	
	:param truetsr: radius of exploration for unifomly slectrec true time shifts (passed to multidraw)
		Not used if draw = False, of course.
		
	:param shuffle: Shuffle the curves before feeding them into optfctnots ?
		If you use this, be sure to pycs.gen.lc.objsort the curves before !
	
	.. todo:: Possibility to add shot noise even with draw = False.
		Like in the good old days. -> Separate addshotnoise from multidraw ?


	.. note:: I will tweak the ML only of those curves that have spline ML. You probably want me to tweak the ML of all your curves !
		So be sure to add some spline ML to all your curves before calling me !
	

	"""
	
	if simdir == None:
		destdir = destpath+"sims_%s" % (simset)
	else:
		destdir = destpath+simdir
	if verbose:
		print "Now thowing dice into %s ..." % destdir
	
	# We prepare the destination directory
	if not os.path.isdir(destdir):
		os.mkdir(destdir)
	else:
		if verbose:
			print "The directory exists, I'll add my new curves."
	
	# Print out some info
	if verbose:
		print "Input shifts :"
		print pycs.gen.lc.getnicetimeshifts(lcs, separator=" | ")
		print "Input delays :"
		print pycs.gen.lc.getnicetimedelays(lcs, separator=" | ")
	
	origshifts = np.array([l.timeshift for l in lcs]) # the mean shifts for the simulations
	
	
	# In case we trace :
	if trace:
		# We put the "real curves" once in each tracedir :
		
		if spline == None:
			splist = []
		else:
			splist = [spline]
		pycs.gen.util.trace(lclist=lcs, splist=splist, tracedir = "trace_sims_%s_tweak" % (simset))
		
		rawlcs = [l.copy() for l in lcs]
		for l in rawlcs: # to get the unshifted curve :
			l.resetshifts()
		pycs.gen.util.trace(lclist=rawlcs, splist=[], tracedir = "trace_sims_%s_draw" % (simset))
	
	
			
	for i in range(npkl):
		
		pklfilepath = os.path.join(destdir, "%i_%.5f.pkl" % (i+1, float(time.time())))
		
		if onlycopy == True:
			if verbose:
				print "Preparing %i identical copies for pkl %i/%i ..." % (n, (i+1), npkl)
			simlcslist = [[l.copy() for l in lcs] for ni in range(n)]
			# We remove any microlensing or shifts :
			for simlcs in simlcslist:
				for l in simlcs:
					l.resetshifts()
	
		else:
			if verbose:
				print "Drawing %i simulations for pkl %i/%i ..." % (n, (i+1), npkl)
					
			simlcslist = []
			
			# The absolute shifts to set the the lcs :
			sampleshifts = [np.random.uniform(low = -truetsr, high = truetsr, size=(len(lcs))) + origshifts for i in range(n)]
		
			for (i, shifts) in enumerate(sampleshifts):
				
				# So this loop is run for every simulated data set.
				# We work on copies of lcs, as we will change them !
				lcscopies = [l.copy() for l in lcs]
				
				# We set the absolute true shifts :
				for (l, shift) in zip(lcscopies, shifts):
					l.timeshift = shift

				if optfctnots == None: # Then we just call draw on these time-shifted curves, using the provided spline and ML etc.
					simlcs = draw(lcscopies, spline, shotnoise=shotnoise, shotnoisefrac=shotnoisefrac, tweakml=tweakml, scaletweakresi=scaletweakresi, tweakspl=tweakspl, keepshifts=False, keeptweakedml=False, trace=trace, tracedir="trace_sims_%s_tweak" % (simset))
				
				else:
					# We fit the custom ML + all shifts except time and get a new spline
					# We start on lcscopies
					
					
					if shuffle:
						pycs.gen.lc.shuffle(lcscopies)
					
					indispline = optfctnots(lcscopies) # Sets ML, mag and flux shifts, but does not change the time shifts.
					
					# Very important : we sort the lightcurves after this optimization !
					# So that the tweakml will correspond.
					if shuffle:
						pycs.gen.lc.objsort(lcscopies, verbose=False)
						
					saveresiduals(lcscopies, indispline) # in case mcres

					simlcs = draw(lcscopies, indispline, shotnoise=shotnoise, shotnoisefrac=shotnoisefrac, tweakml=tweakml, scaletweakresi=scaletweakresi, tweakspl=tweakspl, keepshifts=False, keeptweakedml=False, trace=trace, tracedir="trace_sims_%s_tweak" % (simset))
						
				simlcslist.append(simlcs)
	
		
		# We save the simlcslist into the pkl file
		pycs.gen.util.writepickle(simlcslist, pklfilepath, verbose=verbose)
	
		
		# The trace with tweak is already done by multidraw. We add a trace of the drawn curves :
		if trace:
			for simlcs in simlcslist:
				pycs.gen.util.trace(lclist=simlcs, splist=[], tracedir = "trace_sims_%s_draw" % (simset))
		
		
	
	
	
def shareflux(lc1, lc2, frac=0.01):
	"""
	I add "noise" to lc1 and lc2 by randomly sharing flux between the two sources.
	
	:param frac: The stddev of the gaussian "noise" in flux, with respect to the minimum flux in the curves.
	
	"""
	if not np.all(lc1.jds == lc2.jds):
		raise RuntimeError("I do only work on curves with identical jds !")
	
	#lc1fs = lc1.getrawfluxes()
	#lc2fs = lc2.getrawfluxes()

	minshift = np.fabs(max(lc1.getminfluxshift(), lc2.getminfluxshift()))
	shifts = frac * minshift * np.random.randn(len(lc1))
	shifts = np.clip(shifts, -minshift+1.0, minshift-1.0) # To garantee that we won't get negative fluxes
	
	
	lc1.addfluxes(shifts)
	lc2.addfluxes(-shifts)






	
	
