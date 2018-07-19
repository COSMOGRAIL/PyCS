"""
Module containing all we need to manipulate lightcurves. Could ideally be used by all our curve-shifting algorithms.
	
"""
import sys
import os
import copy as pythoncopy
import cPickle as pickle
import numpy as np
import util

# For the sort and shuffle stuff :
import random
import operator


class lightcurve:
	"""
	A class to handle light curves of lensed QSO images.
	
	The actual information is stored into independent "vectors" (1D numpy arrays) : one for time, one for magnitude...
	It would seem nicer to have some kind of iterable lightcurve object, using perhaps python lists/dicts,
	but this would be slower to manipulate with numpy I guess.
	
	An exception to this is the properties attribute : a list of dicts, one for each data point, usually all containing the same fields.
	
	Note that a lightcurve MUST always be chronologically ordered, this allows faster algorithms.
	See methods sort() and validate().
	Importation functions should automatically sort().
	
	
	There is only one remaining philosophy about handling shifts of lightcurves :
	
		- use shifttime() and cie, this is fast and you do not change the actual data (recommended)
		
	"""

	# Documentation is quite explicit and intensive here, there is no need for the same level of detail for other
	# functions...
	def __init__(self, telescopename = "Telescope", object = "Object", plotcolour = "red"):
		"""
		The plain constructor, binding all the instance variables.
		It gives a small default lightcurve with 5 points, labels, and a mask to play with.
				
		@todo: nothing to do, that's just a test of the "todo" field :-)
		
		@type	telescopename: string
		@param	telescopename: The name of the telescope that took this curve
		@type	object: string
		@param	object: The name of the object
		@type	plotcolour: "colour"
		@param	plotcolour: A colour to use when plotting this lightcurve. In principle this is a string,
		you can use matplotlib names like "red", "blue" etc, or "#FF0000" ...
		
		@note: Usually you will use I{factory-functions} like "importfromvirginie", and not explicitely call
		this constructor.

		"""
		
		# Some various simple attributes :
		
		self.telescopename = telescopename 
		"""@type: string
		@ivar: Choose a name for the telescope
		"""
		
		self.object = object
		"""@type: string
		@ivar: A name for the object"""
		
		self.plotcolour = plotcolour
		"""@type: "colour"
		@ivar: "red", "blue", "#00FF00"..."""
		
		
		# The data itself : here we give some default values for testing purposes :
		
		self.jds = np.array([1.1, 2.0, 3.2, 4.1, 4.9])	# 
		"""@type: 1D float array
		@ivar: The "dates", as floats, in HJD, MHJD, ... these jds should always be in 
		chronological order please (and of course always correspond to the other arrays)."""
		
		self.mags = np.array([-10.0, -10.1, -10.2, -10.3, -10.2])
		"""@type: 1D float array
		@ivar: Magnitudes"""
		
		self.magerrs = np.array([0.1, 0.2, 0.1, 0.1, 0.1])
		"""@type: 1D float array
		@ivar: Errors of magnitudes, as floats"""
		
		# Aside of these, we keep one more array of mags, called origmags, that is handy for stuff like bootstrapping etc.
		# It turns out all this infrastructure was to complicated... It's easier to do this for your own, if you want to
		# be able to come back.
		#self.origmags = self.mags.copy()
		#"""@type: 1D float array
		#@ivar: Kind of a backup-copy of the original mags : allows you to "apply" microlensing
		#or add bootstrapping etc to the primary mags, and still be able to restore original mags.
		#"""
		
		# Other attributes :
		
		self.mask = np.array([True, True, False, True, True]) 
		"""@type: 1D boolean array
		@ivar: A boolean mask, allows to "switch off" some points.
		True means "included", False means "skip this one". This is nice, it directly allows you to
		use the mask as indices; example (suppose C{mylc} is a lightcurve) :
		
			>>> mylc.mask[34] = False
			>>> mylc.mags[mylc.mask].std()
			
		masks the 35th point and gives you the standard deviation of the unmasked mags. Masked points
		show up on graphs as circled in black, and are naturally skipped by dispersion methods etc.
		
		"""
		
		
		# Can we make it even cooler ? Yes we can :
		
		self.labels = ["0", "1", "2 = Huge problem", "3", "4"]
		"""@type: list of strings
		@ivar: Each point has a string label. Plots can show these labels.
		Setting these is simple :
		
			>>> mylc.labels[42] = "This point is funny"
		
		We can use them for all kind of things (debugging...).
		"""
		
		
		self.commentlist = ["Test lightcurve"]
		"""@type: list of strings
		@ivar: This is a list of string comments for each full lightcurve, and NOT one comment for each point.
		You can "append" comments at any time; they can be shown by various functions, and will trace the 
		history of the lightcurve.
		"""
		
		self.properties = [{"fwhm":1.0}, {"fwhm":1.0}, {"fwhm":1.0}, {"fwhm":1.0}, {"fwhm":1.0}]
		"""
		:ivar properties: This is a list of dicts (exactly one per data point) in which you can put whatever you want to.
			You can use them to carry image names, airmasses, sky levels, telescope names etc.
			Of course this slows down a bit things like merging etc.
		:type properties: list of dicts of whatever
		
		"""
		
		
		# Other useful settings for plots
		
		self.ploterrorbars = True
		"""@type: boolean
		@ivar: A flag to show or hide errorbars in plots of this lightcurve."""
		
		self.showlabels = False
		"""@type: boolean
		@ivar: A flag to show or hide labels in plots of this lightcurve."""

		#self.hidemasked = False # Not yet implemented... not sure if this is needed
		
		
		# Values of the two obvious possible shifts.
		
		self.timeshift = 0.0
		"""@type: float
		@ivar: A float giving the shift in time (same unit as jds) that was applied to the lightcurve. 
		This is updated at every shift, so that a lightcurve is always "aware" of its eventual shift. 
		"""
		
		self.magshift = 0.0
		"""@type: float
		@ivar: A float giving the shift in magnitude (similar to timeshift)."""
		
		self.fluxshift = 0.0
		"""@type: float
		@ivar: A float giving the shift in flux (!). As we work with magnitudes, a shift in flux is a bit special.
		Note by the way that a flux shift does not commute with a magnitude shift !
		"""
		
		# And now the microlensing : by default there is none :
		
		self.ml = None
		"""@type: microlensings object
		@ivar: The optional microlensings object allows to take
		microlensing into account during the optimizations and time-delay quests."""
		
			
			
	# I explicitly define how str(mylightcurve) will look. This allows a nice "print mylightcurve" for instance !
	def __str__(self):
		"""
		We explicitly define how str(mylightcurve) will look. "[Euler/A]", or "[Euler/A](0.00,0.10)"
		in case the lightcurve is shifted (here by 0.1 mag).
		For more info, look at longstr().
		This will show up in plot legends etc. Plus, at any time, it allows you to :
		
			>>> print mylc
			
		Microlensing "constrains" are shown like |12|, meaning in this case
		two seasons : first one with 1 parameter, second one with two.
		"""
		retstr = "[%s/%s]" % (self.telescopename, self.object)
		
		if self.timeshift != 0.0 or self.magshift != 0.0 or self.fluxshift != 0.0:
			retstr += "(%.3f,%.3f,%.0f)" % (self.timeshift, self.magshift, self.fluxshift)
		
		if self.ml != None:
			retstr += "%s" % (self.ml)
			#if self.hideml == True:
			#	retstr += "X|"
			
		return retstr
	
	def __len__(self):
		"""
		We define what len(mylightcurve) should return : the number of points.
		I do not exclude masked points -- this is the full length.
		"""
		return len(self.jds)
	
	
	def samplingstats(self, seasongap=60.0):
		"""
		Calculates some statistics about the temporal sampling of the lightcurve object.
		I do not take into account the mask !
		
		:param seasongap: Minimal interval in days to consider a gap to be a season gap.
		
		I return a dict with the following keys :
		
		* nseas : number of seasons
		* meansg : mean length (in days, of course) of the season gaps
		* minsg : min season gap
		* maxsg : max season gap
		* stdsg : standard deviation of the season gaps
		* med : the median sampling interval inside seasons
		* mean : the mean sampling interval inside seasons
		* min : ...
		* max
		* std
		
		"""
		
		self.validate()
		
		stats = {"len":len(self)}
		
		first = self.jds[:-1]
		second = self.jds[1:]
		gaps = second - first
		
		seasongaps = gaps[gaps >= seasongap]
		intervals = gaps[gaps < seasongap] # The "normal" gaps
		
		if len(seasongaps) == 0:
			seasongaps = [0]
			stats["nseas"] = 1 # Number of seasons
		else:
			stats["nseas"] = int(len(seasongaps) + 1) # Number of seasons
		stats["meansg"] = np.mean(seasongaps)
		stats["minsg"] = np.min(seasongaps)
		stats["maxsg"] = np.max(seasongaps)
		stats["stdsg"] = np.std(seasongaps)
		
		stats["med"] = np.median(intervals)
		stats["mean"] = np.mean(intervals)
		stats["max"] = np.max(intervals)
		stats["min"] = np.min(intervals)
		stats["std"] = np.std(intervals)
		
		return stats
		
		
		
		
	
	def commonproperties(self, notonlycommon = False):
		"""
		Returns a list of strings containing the property field names that are common to
		all data points.
		If notonlycommon is True, then I return all fields, not only the ones common to
		all my data points.
		"""
		propkeys = [set(props.keys()) for props in self.properties]
		if notonlycommon:
			return sorted(list(reduce(lambda x,y: x.union(y), propkeys)))
		else:
			return sorted(list(reduce(lambda x,y: x.intersection(y), propkeys)))
	
	def longinfo(self, title="None"): # note that the name info() is reserved for python
		"""
		Returns a multi-line string (small paragraph) with information
		about the current state of a lightcurve.
		"""

		if title != "None":
			titlestr = "\n\t%s" % (title.upper())
		else:
			titlestr = ""

		samplingstats = self.samplingstats()
		string_list = ["- "*30, titlestr, "\n\t", str(self),
		"\n", "%i points (total), %i of which are masked" % (len(self), self.mask.tolist().count(False)),
		"\n", "%i seasons (gap: >60), gap length : %.1f +/- %.1f days" % (samplingstats["nseas"], samplingstats["meansg"], samplingstats["stdsg"]),
		"\n", "Sampling : median %.1f, mean %.1f, max %.1f, min %.2f days" % (samplingstats["med"], samplingstats["mean"], samplingstats["max"], samplingstats["min"]),
		
		"\nShifts : (%.5f,%.5f,%.2f) [days, mag, flux]" % (self.timeshift, self.magshift, self.fluxshift),
		"\nColour : ", self.plotcolour,
		"\nCommon properties : ", ", ".join(self.commonproperties()),
		"\n   All properties : ", ", ".join(self.commonproperties(notonlycommon = True)),
		"\nComments :"]
		for comment in self.commentlist :
			string_list.append("\n   %s" %comment)
		if self.ml != None:
			#if self.hideml == True:
			#	mlstr = "\nHIDDEN Microlensing %s" % (str(self.ml))
			#else:
			mlstr = "\nMicrolensing : %s" %  (str(self.ml))
			#mlstr = mlstr + self.ml.longinfo().replace("\n", "\n   ")
			string_list.append(mlstr)
		string_list.extend(["\n", "- "*30])
		return ''.join(string_list)
	
	
	def printinfo(self, title="None"):
		""" Prints the above longinfo about the lightcurve."""
		print self.longinfo(title)


	# Elementary lightcurve methods related to shifts
	
	
	
			
	def shifttime(self, days):
		"""
		Update the value of timeshift, relative to the "current" shift.
		The jds array is not shifted; we just keep track of a timeshift float. So this is fast.
		"""
		self.timeshift += float(days)
		
	def shiftmag(self, deltamag):
		"""
		Shifts the curve in mag, with respect to the "current" shift. Same remark as for shifttime().
		@todo: think about changing the errors as well... -- no, the error in magnitude is a relative error.
		"""
		self.magshift += float(deltamag)



	
	def shiftflux(self, flux):
		"""
		Idem, but a shift in flux.
		@todo: Hmm, here we should change the errors as well -- but flux shifts are usuall very small, no ?
		"""
		self.fluxshift += float(flux)
	
	def setfluxshift(self, flux, consmag = False):
		"""
		An absolute setter for the fluxshift, that directly tries to correct
		for the equivalent magnitude shift with respect to the mean magnitude of the curve.
		So be careful : I change magshift as well !
		"""
		if consmag == True:
		
			oldfsmags = self.calcfluxshiftmags()
			self.fluxshift = float(flux)
			newfsmags = self.calcfluxshiftmags()
			magcor = - np.mean(newfsmags - oldfsmags)
			#print magcor
			self.shiftmag(magcor) # to compensate for the fluxshift
		else:
			self.fluxshift = float(flux)
		
	

	def getjds(self):
		"""
		A getter method returning the jds + timeshift. This one does the actual shift.
		Because of the addition, we return a copy.
		"""
		return self.jds + self.timeshift
	
	def getminfluxshift(self):
		"""
		Returns the minimal flux shift (i.e. largest negative fluxshift) so that the flux of all points in this curve
		remain positive.
		-> returns - the minimum flux of the curve.
		This does not depend on timeshifts, magshifts, or microlensing.
		When optimizing fluxshifts, you should check by yourself that you respect this minimum value.
		"""
		return - np.min(self.getrawfluxes())
	
	def getrawfluxes(self):
		"""
		Returns the raw fluxes corresponding to self.mags
		"""
		
		return 10.0**(self.mags/-2.5)
		
	
	def calcfluxshiftmags(self, inverse=False):
		"""
		Returns an array of mags that you have to add to the lc.mags if you want to take into account for the fluxshift.
		Be careful with the order of things ! Always take into account this fluxshift *before* the magshift / microlensing.
		
		Inverse = True is a special feature, it returns the mags to subtract from fluxshiftedmags to get back the original mags.
		Perfectly clear, no ? Just changes signs in the fomula.
		You could also use - calcfluxshiftmags(self, inverse=True) and reverse the sign of the fluxshift, but this is more convenient.
		"""
		if self.fluxshift != 0.0:
			# We need to include a check :
			if inverse:
				shifts = 2.5*np.log10( (-self.fluxshift * np.ones(len(self)) / (10.0**(self.mags/-2.5)))+ 1.0)
			else:
				shifts = -2.5*np.log10( (self.fluxshift * np.ones(len(self)) / (10.0**(self.mags/-2.5)))+ 1.0)
			if np.all(np.isnan(shifts) == False) == False: # If there is a nan in this...
				print "Ouch, negative flux !"
				return np.zeros(len(self))
			else:
				return shifts
			
		else:
			return np.zeros(len(self))
		
	
	def addfluxes(self, fluxes):
		"""
		Give me fluxes as a numpy array as long as self.mags, and I'll add those to self.mags.
		I work directly on self.mags, thus I do not take into account magshifts or ml.
		As this changes the actual self.mags, it cannot be undone !
		"""	
		# We check that we won't go into negative fluxes :
		if not np.all(fluxes > self.getminfluxshift()):
			raise RuntimeError("That would give negative fluxes ...")
		
		newfluxes = self.getrawfluxes() + fluxes
		self.mags = -2.5*np.log10(newfluxes)
		self.commentlist.append("Added some fluxes")
		
	
	
	def getmags(self, noml = False):
		"""
		A getter method returning magnitudes. Now this is non trivial, as magnitudes are influenced by :
			- fluxshift
			- magshift
			- microlensing (= like a variable magshift)
			
		We "add" this stuff in this order, i.e. fluxshift always comes first.
		We always return a copy of the array, to prevent you from changing the actual lightcurve.
		
		
		lc.mags "+" fluxshift + magshift + microlensings if ml is present,
		and just lc.mags "+" fluxshift + magshift if not.
		You can overwrite this : if noml = True, I don't include the ml even if a microlensing is present !
		
		A bit ugly all these ifs, but I think it makes it faster (the fluxshift calculation is slow).
		
		"""
		if self.fluxshift != 0.0:
			if (self.ml != None) and (noml == False):
				return self.mags + self.calcfluxshiftmags() + self.magshift + self.ml.calcmlmags(self)
			else:
				return self.mags + self.calcfluxshiftmags() + self.magshift	
		else:
			if (self.ml != None) and (noml == False):
				return self.mags + self.magshift + self.ml.calcmlmags(self)
			else:
				return self.mags + self.magshift	


#	def getflux(self):
# 		"""
# 		For debugging ...
# 		No, see getrawfluxes !
# 		"""
# 		return 10.0**(self.getmags()/-2.5)
# 		
	
	def getmagerrs(self):
		"""
		Hmm, how does fluxshift and stuff influence the magerrs ?
		Copy is important here as we don't do any calculation.
		"""
		return self.magerrs.copy()
	
	
	
	# Now, stuff for microlensing and manipulating the magnitudes :
	
	def addml(self, microlensing):
		"""
		Adds a microlensing object (by add we mean it replaces one if it was already present) -- we are not just setting
		the parameters here ! Note that we COPY the microlensing object, so that the one added is independent from yours.
		"""
		
		if self.ml != None:
			print "I replace an existing mircolensing."
			
		
		self.ml = microlensing.copy() # this copy is important if you append the "same" new ml object to different lcs.
		
		self.ml.checkcompatibility(self)
		self.commentlist.append('Setup of microlensing "%s"'% self.ml)
		
		
	def rmml(self):
		"""Simply deletes the microlensing object"""
		if self.ml != None : self.commentlist.append('Microlensing "%s" removed.'% self.ml)
		self.ml = None
	
	def resetml(self):
		"""
		Resets the microlensing to "zero".
		* for polynomial ML the coefficients are put to zero
		* for spline ML we reset the coeffs to zero
		"""
		if self.ml != None:
			self.ml.reset()
		
		#self.commentlist.append("ML Reset")
		
		
	
	def resetshifts(self, keeptimeshift=False):
		"""
		Removes all shifts and microlensing, getting back to the observed data.
		Does keep the mask.
		"""
		
		self.rmml()
		self.fluxshift = 0.0
		self.magshift = 0.0
		#self.commentlist = ["Reset"]
		#self.commentlist.append("###### Full Reset #######")
		
		if keeptimeshift == False:
			self.timeshift = 0.0
		
	
	def applyfluxshift(self):	
		"""
		It adds the fluxshift-float to the present flux, then puts this fluxshift-float to 0. So that "nothing" changes as seen from
		the outside.
		This is used for instance when you want to merge different lightcurves with different fluxshifts.
		Needless to say, use this carefully, only if it really makes sense for what you want to do.
		
		Note that we do not touch microlensing here, it remains in place and does not change its meaning in any way. 
		"""
		
		self.mags += self.calcfluxshiftmags()
		self.commentlist.append("CAUTION : fluxshift of %f APPLIED" % (self.fluxshift))
		self.fluxshift = 0.0 # as this is now applied.
	
	def applymagshift(self):	
		"""
		It adds the magshift-float to the present mags, then puts this magshift-float to 0. So that "nothing" changes as seen from
		the outside.
		This is used for instance when you want to merge different lightcurves with different magshifts.
		Needless to say, use this carefully, only if it really makes sense for what you want to do.
		
		Note that we do not touch microlensing here, it remains in place and does not change its meaning in any way. 
		"""
		
		if self.fluxshift != 0.0:
			raise RuntimeError, "Apply the fluxshift before applying the magshift !"
		
		self.mags += self.magshift
		self.commentlist.append("CAUTION : magshift of %f APPLIED" % (self.magshift))
		self.magshift = 0.0 # as this is now applied.


	def applyml(self):
		""" 
		We "add" the microlensing to the actual mags, then remove the microlensing object.
		The advantage is that getmags() gets faster, as it does not have to evaluate the microlensing anymore.
		It also allows you to do some other tricks like tweaking a curve, shifting seasons etc if you want.
		So as a conclusion : use this when you know what you are doing.
		
		We do not touch magshift, it remains perfectly valid afterwards.
		But, if there is a *fluxshift*, we will have to make sure that it was "applied" before !
		"""
		
		if self.ml == None:
			raise RuntimeError, "Hey, there is no ml associated to this lightcurve !"
		
		if self.fluxshift != 0.0:
			raise RuntimeError, "Apply the fluxshift before applying the ML !"
			# This is really important. One possibility would be that this function applies the fluxshift ?
		
		self.mags += self.ml.calcmlmags(self)
		self.commentlist.append('CAUTION : microlensing %s APPLIED' % self.ml)
		self.rmml() # very important, we get rid of the present stuff. This prevents you from applying it twice etc.
	

	# And various other features in random order.

	def setindexlabels(self):
		"""First point gets "0", second point gets "1", etc; handy to identify troublemakers on a plot !"""
		
		self.labels = [str(i) for i in range(len(self.jds))]
		#self.showlabels = True
	
	def setjdlabels(self):
		"""Points are labeled by their mjd rounded to .1 days. Use these labels to write a skiplist."""
		
		self.labels = ["%.1f" % (jd) for jd in self.jds]
		#self.showlabels = True
	
	def maskskiplist(self, filepath, searchrange=0.2, accept_multiple_matches=False, verbose=True):
		"""
		I mask points according to a skiplist. I do not modify the mask for points that are not on the skiplist,
		i.e. I do not reset the mask in any way.
		The structure of the skiplist is one "point" per line, where the first element in the line
		is the MJD (for instance as identified on a plot using setjdlabels() !).
		For each point from the skiplist, I will search for the corresponding point in the lightcurve
		within searchrange days. I will warn you in case of anything fishy.
		
		
		:param filepath: file to read
		:type filepath: string or path
		:param searchrange: range in which I'll look for points (in days).
		:type searchrange: float	
		
		
		"""
		skippoints = util.readidlist(filepath, verbose=verbose)
		
		if verbose:
			if self.hasmask():
				print "Note : %i epochs are already masked." % (np.sum(self.mask == False))

		for skippoint in skippoints:
			skipjd = float(skippoint[0])
			indices = np.argwhere(np.fabs(self.jds - skipjd) <= searchrange)
			#print indices
			if len(indices) == 0:
				if verbose:
					print "Warning, epoch %s from skiplist not found in %s !" % (skippoint[0], str(self))
			elif len(indices) > 1:
				if verbose:
					print "Warning, multiple matches for epoch %s from skiplist !" % (skippoint[0])
				if accept_multiple_matches:
					if verbose:
						print "I mask all of them..."
					for index in indices:
						if self.mask[index] == False:
							if verbose:
								print "Epoch %s is already masked." % (skippoint[0])
						else:
							self.mask[index] = False
			elif len(indices) == 1:
				index = indices[0]
				if self.mask[index] == False:
					if verbose:
						print "Epoch %s is already masked." % (skippoint[0])
				else:
					self.mask[index] = False
		
		if verbose:
			print "Done with maskskiplist, %i epochs are now masked." % (np.sum(self.mask == False))
	
	def maskinfo(self):
		"""
		Returns a description of the masked points and available properties of them.
		Note that the output format can be directly used as a skiplist.
		"""
		
		cps = self.commonproperties()
		lines = []
		maskindices = np.argwhere(self.mask == False)
		for maskindex in maskindices:
			comment = ", ".join(["%s : %s" % (cp, self.properties[maskindex][cp]) for cp in cps])
			txt = "%.1f    %s" % (self.jds[maskindex], comment)
			lines.append(txt)
			
		txt = "\n".join(lines)
		txt = "# %i Masked points of %s :\n" % (np.sum(self.mask == False), str(self)) + txt
		return txt
		
		
	def clearlabels(self):
		"""Sets label = "" for all points.
		You could use this before setting the labels of just some of the points "by hand", for instance."""
		
		self.labels = [""] * len(self)
	
	def clearproperties(self):
		"""Removes all properties"""
		
		self.properties = [{} for i in range(len(self))]
		
	def copy(self):
		"""
		Returns a "deep copy" of the lightcurve object. Try to avoid this within loops etc ... it is slow !
		Typically if you want to optmize time and mag shifts, think about using local backups of lc.getmags() and lc.getjds() etc ... 
		
		We use the copy module, imported as "pythoncopy" to avoid confusion with this method.
		
		@rtype: lightcurve
		@return: A copy of the lightcurve.
		
		""" 
		return pythoncopy.deepcopy(self)

	def cutmask(self):
		"""
		Erases all masked points from the lightcurve.
		
		This is one way of handling the mask, but perhaps not the best/fastest one,
		depending on what you want do !
		If you write a curve shifting algorithm, probably your code will be more elegant if you teach him to simply skip
		masked points, instead of making a copy and removing them using this method !
		
		.. warning:: Cutmask will change the meaning of seasons, that's why we are forced to delete any present microlensing. If you want to "keep" it, you'll have to apply it first.
		"""
		
		self.jds = self.jds[self.mask]
		self.mags = self.mags[self.mask]
		self.magerrs = self.magerrs[self.mask]
		self.labels = [self.labels[i] for i in range(len(self.mask)) if self.mask[i]]
		# This is a bit harder, as labels is a plain python list and not a numpy array.
		# Be careful not to use self.jds in the range() at this point ...
		# Same is true for the properties :
		self.properties = [self.properties[i] for i in range(len(self.mask)) if self.mask[i]]
		# finally we change the mask itself :
		self.mask = self.magerrs >= 0.0	# This should be True for all !
		
		if self.ml != None:
			print "WARNING : cutmask() just removed your microlensing !"
			self.rmml()
		# Note that removing microlensing is important, as the parameters change
		# their meanings, and seasons might change as well etc.
	
	def hasmask(self):
		"""
		Returns True if some points are masked (i.e., if there is a False), False otherwise.
		"""
		return not np.alltrue(self.mask)
	
	def clearmask(self):
		"""
		Sets all the mask to True.
		"""
		self.mask = np.ones(len(self), dtype=np.bool)
	
	def validate(self, verbose=False):
		"""
		Checks the "health" and coherency of a lightcurve.
		Are there as many mags as jds, etc ? No return value, but raises RuntimeError in case
		of failure.
		
		.. note:: Here we also check that the curve is "ordered", i.e. the jds are monotonously increasing.
			At this stage it is OK to have several points with the *same* julian date. 
			In fact such points are a problem for splines, but this is addressed once you use splines.
		
		"""
		
		ndates = len(self.jds)
		if len(self.mags) != ndates or len(self.magerrs) != ndates or len(self.mask) != ndates or len(self.labels) != ndates or len(self.properties) != ndates :
			raise RuntimeError, "Incoherence in the length of your lightcurve !"
	
		# I postulate that lightcurves shorter then 2 points make no sense (the next test, and seasons() etc would crash...)
		if ndates < 2:
			raise RuntimeError, "Your lightcurve is too short... give me at least 2 points."
	
		# some algorithms have been optimized so that they NEED ordered lightcurves, i.e. values in self.jds must be increasing.
		# here we check if this is the case for a given lightcurve.
		first = self.jds[:-1]
            	second = self.jds[1:]
		if not np.alltrue(np.less_equal(first,second)): # Note the less_equal : ok to have points with same JD.
			raise RuntimeError, "Your lightcurve is not sorted !"
		
		# we check if the errors are positive:
		if not np.alltrue(self.magerrs > 0.0):
			raise RuntimeError, "Your lightcurve has negative or null errors !"
		# note that we use this sometimes in the code to generate masks etc... so don't just
		# remove this check.
		
		# the associated microlensing object if present
		if self.ml != None:
			self.ml.checkcompatibility(self)
			
		if verbose:
			print "%s : validation done !" % (self)

	def sort(self):
		"""
		We sort the lightcurve points according to jds.
		Of course all the arrays and lists have to be sorted.
		
		The meaning of seasons is changed, that's why we have to delete any microlensing.
		"""
		# Now look how nice this is (so proud...) :
		sortedindices = np.argsort(self.jds)
		self.jds = self.jds[sortedindices]
		self.mags = self.mags[sortedindices]
		self.magerrs = self.magerrs[sortedindices]
		self.mask = self.mask[sortedindices]
		self.labels = [self.labels[i] for i in sortedindices]	# trick as labels is not a numpy array
		self.properties = [self.properties[i] for i in sortedindices]	# trick as properties is not a numpy array
		
		if self.ml != None:
			print "WARNING : sort() just removed your microlensing !"
			self.rmml()
			
			
	def montecarlomags(self, f = 1.0, seed=None) :
		"""
		We add gaussian noise to all mags according to their errors.
		We do not care about any microlensing of shifts here, but directly tweak the self.mags.
		The order of the points is not touched.
		
		I modify all points, even masked ones.
		
		:param seed: if None, the clock is used, if not, the given seed.
		:type seed: int or None
		
		:param f: a multiplier
		
		"""
	
		self.commentlist.append("Monte Carlo on mags !") # to avoid confusions.
		rs = np.random.RandomState(seed) # we create a random state object, to control the seed.	
		self.mags += rs.standard_normal(self.mags.shape) * f * self.magerrs # and here we do the actual bootstrapping !

	def montecarlojds(self, amplitude=0.5, seed=None, keepml=True) :
		"""
		We add a UNIFORM random variable from -amplitude to +amplitude to each jd.
		This is not trivial, as we have to take care about the order.
		
		We need sorting, so normally the microlensing would be removed.
		But if you know that, using only a small amplitude, the seasons will remain of the same length after the bootstrapping,
		you can use keepml=True (default), and we will "keep" the microlensing. Note that the microlensing functions are defined
		with respect to the jds : so as the jds are bootstrapped, the mircolensing is slightly affected, but this is fine as typically
		you just want to keep a rough estimate of the microlensing, and run an ML optimization anyway !
		Hence this obscure option.
		
		I modify all points, even masked ones.
		"""
		
		self.commentlist.append("Monte Carlo on jds !") # to avoid confusions.
		rs = np.random.RandomState(seed) # we create a random state object, to control the seed.
		self.jds += (rs.uniform(low=0.0, high=2.0*amplitude, size=self.jds.shape) - amplitude)
		# uniform distribution. Yes, this is a bit strange, but low cannot be negative.
		
		# And now everything is fine but the curve might not be sorted, so :
		if keepml == False or self.ml == None:
			self.sort() # bye bye microlensing
		else:
			ml_i_wanna_keep = self.ml.copy()
			self.rmml()
			self.sort()
			self.addml(ml_i_wanna_keep)
			# Isn't this wonderful ?!
	
	def pseudobootstrap(self):
		"""
		We randomly mask some points, but without duplicating them.
		Real bootstap would be to draw N points from a lightcurve of N points with replacements.
		But we do not replace here, instead we do as if we would replace, but then skip any "duplicates".
		
		I respect mask : masked points will stay masked, I do as if they were not there.
		"""
		
		indices = np.arange(len(self))[self.mask] # the indices we want to bootstrap : only unmasked ones.
		
		if indices.size == 1:
			#print "WARNING : not enough points to bootstrap !"
			#self.commentlist.append("Bootstrap error !")
			raise RuntimeError, "Not enough points to bootstrap !"
		
		draws = np.random.randint(0, high=indices.size, size=indices.size) # indexes of the indices that are drawn
		# We remove the duplicates from these draws:
		#uniquedraws = np.array(sorted(list(set(list(draws)))))
		
		# Faster and easier :
		newmask = np.zeros(len(self), dtype=np.bool) # an array of False, as long as the full curve
		newmask[indices[draws]] = True # drawn points are set True
		self.mask = newmask
	
		self.commentlist.append("Pseudobootstraped !")
		
		
		
	def bootstrap(self):
		"""
		Real bootstap = we draw the points with replacements.
		You will get points with identical jds in the curve !
		
		I respect mask : masked points will stay masked, I do as if they were not there.
		
		Hmm, a bit slow... would need to rearrange all the properties etc
		
		In fact no, should be ok, do similar then above
		
		"""
		pass
		
	
	
	
	def merge(self, otherlc):
		"""
		Merges a lightcurve into the current one. Warning : it takes the other lightcurve as it comes, with all shifts
		or eventual microlensings -- "applied" or just specified !
		In fact it behaves as if every shift or microlensing was applied before the merging.
				
		
		It's up to you to check that you don't merge lightcurves with timeshifts, which would probably be nonsense.
		We delete an eventual current microlensing object -> You will have to define new seasons and parameters if you want.
		
		:param otherlc: A lightcurve to merge into the current one.
		:type otherlc: lightcurve
		
		The origin of each point of the otherlc is appended to the labels.
		Settings and colour etc from the current curve are not changed.
		After the merging process, the lightcurve is sorted and validated.
		
		"""
		# Let's warn the user if there are timeshifts in the curves :
		if self.timeshift != 0.0 or otherlc.timeshift != 0.0 :
			print "WARNING : you ask me to merge time-shifted lightcurves !"
		
		# and microlensing :
		if self.ml != None or otherlc.ml != None :
			print "WARNING : I am merging lightcurves with possible microlensing !"
			
		# for the magnitude shift, for otherlc it is quite common, but not for self :
		if self.magshift != 0.0 or self.fluxshift != 0.0:
			print "WARNING : merging into a lightcurve with magshift or fluxshift : everything gets applied ! "
		
		# Frist we just concatenate all the values into new numpy arrays :
		
		concjds = np.concatenate([self.getjds(), otherlc.getjds()])
		concmags = np.concatenate([self.getmags(), otherlc.getmags()])
		# To calculate the ML, we need the untouched input "self" curve.
		# Hence this new variable (and not directly using self.mags etc) !
		concmagerrs = np.concatenate([self.magerrs, otherlc.magerrs])
		concmask = np.concatenate([self.mask, otherlc.mask])
		
		
		self.jds = concjds
		self.mags = concmags
		self.magerrs = concmagerrs
		self.mask = concmask
		self.resetshifts() # As we used getjds() and getmags() above ...
		
		
		# We change the new labels so that they tell us from which lightcurve the point was from.
		newlabels = [label + "(from %s)" % str(otherlc) for label in otherlc.labels] # This way we make a copy
		self.labels.extend(newlabels)
		
		# And now the properties :
		self.properties.extend(otherlc.properties)
		
		self.commentlist.append("Merged with %s" % str(otherlc) )
		#self.commentlist.extend(otherlc.commentlist)
		self.telescopename = self.telescopename + "+%s" % otherlc.telescopename
		
		# The very essential :
		self.sort()
		
		# Just to be sure that everything went fine :
		self.validate()
			

	def rdbexport(self, filename=None, separator="\t", writeheader=True, rdbunderline=True, properties=None):
		"""
		Writes the lightcurve into an "rdb" file, that is tab-separeted columns and a header.
		Note that any shifts/ML are taken into account. So it's a bit like if you would apply the
		shifts before writing the file.
		
		Includes mask column only if required (if there is a mask)
		
		
		:param filename: where to write the file
		:type filename: string or path
		:param separator: how to separate the collumns
		:type separator: string
		:param writeheader: include rdb header ?
		:type writeheader: boolean	
		:param properties: properties of the lightcurves to be include in the file.
		:type properties: list of strings, e.g. ["fwhm", "ellipticity"]
			
		"""
		import csv
		
		self.validate() # Good idea to keep this here, as the code below is so ugly ...
	
		if filename == None:
			filename = "%s_%s.rdb" % (self.telescopename, self.object)
		
		# We include a "mask" column only if mask is not True for all points
		if False in self.mask:
			colnames = ["mhjd", "mag", "magerr", "mask"]
			data = [self.getjds(), self.getmags(), self.magerrs, self.mask]
		
		else:
			colnames = ["mhjd", "mag", "magerr"]
			data = [self.getjds(), self.getmags(), self.magerrs]
		
		
		# Now we do some special formatting for the cols mhjd, mag, and magerr
		data[0] = map(lambda mhjd: "%.8f" % (mhjd), data[0]) # formatting of mhjd
		data[1] = map(lambda mhjd: "%.8f" % (mhjd), data[1]) # formatting of mhjd
		data[2] = map(lambda mhjd: "%.8f" % (mhjd), data[2]) # formatting of mhjd
		
		data = map(list, list(zip(*data))) #  list to make it mutable
		
		# We add further columns
		if properties == None:
			properties = []
		colnames.extend(properties)
		for i in range(len(self.jds)):
			for property in properties:
				data[i].append(self.properties[i][property])

		underline = ["="*n for n in map(len, colnames)]	
		
		outfile = open(filename, "wb") # b needed for csv
		writer = csv.writer(outfile, delimiter=separator)
	
		if writeheader :
			writer.writerow(colnames)
			if rdbunderline:
				writer.writerow(underline)
		writer.writerows(data)
		
		outfile.close()
		print "Wrote %s into %s." % (str(self), filename)


# Factory functions : give me some stuff and I make a lightcurve object from this.

def factory(jds, mags, magerrs=None, telescopename="Unknown", object="Unknown", verbose=False):
	"""Returns a valid lightcurve object from the provided arrays.
	The numpy arrays jds and mags are mandatory. If you do not specify a third array containing the magerrs,
	we will calculate them "automatically" (all the same value), to avoid having 0.0 errors.
	
	@type	jds: 1D numpy array
	@param	jds: julian dates
	@type	mags: 1D numpy array
	@param	mags: magnitudes
	@type	magerrs: 1D numpy array
	@param	magerrs: optional magnitude errors
	
	@todo: improve it and use this in file importing functions
	
	"""
	# Make a brand new lightcurve object :
	newlc = lightcurve()
	
	# Of couse we can/should check a lot of things, but let's be naive :
	
	newlc.jds = np.asarray(jds)
	newlc.mags = np.asarray(mags)
	
	if magerrs is None:
		newlc.magerrs = np.zeros(len(newlc.jds)) + 0.1
	else:
		newlc.magerrs = np.asarray(magerrs)

	if len(newlc.jds) != len(newlc.mags) or len(newlc.jds) != len(newlc.magerrs):
		raise RuntimeError, "lightcurve factory called with arrays of incoherent lengths"

	newlc.mask = newlc.magerrs >= 0.0	# This should be true for all !

	newlc.properties = [{}] * len(newlc.jds)

	newlc.telescopename = telescopename
	newlc.object = object

	newlc.setindexlabels()
	newlc.commentlist = []

	newlc.sort() # not sure if this is needed / should be there

	newlc.validate()

	if verbose: print "New lightcurve %s with %i points" % (str(newlc), len(newlc.jds))

	return newlc


def flexibleimport(filepath, jdcol=1, magcol=2, errcol=3, startline=1, flagcol=None, propertycols=None, telescopename="Unknown", object="Unknown", plotcolour="red", verbose = True, absmagerrs=False):
	"""
	The general form of file reading. We read only one lightcurve object.
	To comment a line in the input file, start the line with # like in python.

	:param jdcol: The column number containing the MHJDs. First column is number 1, not 0 !
	:param magcol: The column number containing the magnitudes.
	:param errcol: The column number containing the magnitude errorbars.
	:param flagcol: (default : None) The colum containing a mask (if available). This one should contain False (= will be masked) or True (not masked).

	:param propertycols: (default : None) is a dict : ``propertycols = {"fwhm":7, "skylevel":8}`` means that col 7 should be read in as property "fwhm" and 8 as "skylevel".
	:type propertycols: dictionary


	"""
	if verbose : print "Reading \"%s\"..." % (os.path.basename(filepath))
	rdbfile = open(filepath, "r")
	rdbfilelines = rdbfile.readlines()[startline-1:] # we directly "skip" the first lines of eventual headers ...
	rdbfile.close()

	jds = []
	mags = []
	magerrs = []
	flags = []
	properties = []

	elsperline = None

	for i, line in enumerate(rdbfilelines) :


		if line[0] == "#":
			continue

		if len(line.strip()) < 5:
			print "Skipping empty line %i : %s" % (i+startline, repr(line))
			continue

		elements = line.split() # line is a string, elements is a list of strings

		# We check the consistency of the number of elements...
		if elsperline != None:
			if len(elements) != elsperline:
				raise RuntimeError, "Parsing error in line %i, check columns : \n%s" % (i+startline, repr(line))
		elsperline = len(elements)

		jds.append(float(elements[jdcol-1]))
		mags.append(float(elements[magcol-1]))
		magerrs.append(float(elements[errcol-1]))

		if flagcol != None :
			strflag = str(elements[flagcol-1])
			if strflag == "True":
				flags.append(True)
			elif strflag == "False":
				flags.append(False)
			else:
				print "Flag error in line %i : %s" % (i+startline, repr(line))
				flags.append(True)
		else :
			flags.append(True)

		propdict = {} # an empty dict for now
		if propertycols != None :
			for (propname, propcol) in propertycols.items():
				propdict[propname] = str(elements[propcol-1])
		properties.append(propdict)

	if absmagerrs:
		magerrs = np.abs(np.array(magerrs))

	# Make a brand new lightcurve object :
	newlc = factory(np.array(jds), np.array(mags), magerrs=np.array(magerrs), telescopename=telescopename, object=object)
	newlc.properties = properties[:]
	newlc.mask = np.array(flags[:])
	newlc.plotcolour = plotcolour
	nbmask = np.sum(newlc.mask==False)
	commentcols = "(%i, %i, %i)" % (jdcol, magcol, errcol)
	newlc.commentlist.insert(0, "Imported from %s, columns %s" % (os.path.basename(filepath), commentcols))
	if verbose: print "%s with %i points imported (%i of them masked)." % (str(newlc), len(newlc.jds), nbmask)
	return newlc




def rdbimport(filepath, object="Unknown", magcolname="mag", magerrcolname="magerr", telescopename="Unknown", plotcolour="red", mhjdcolname="mhjd", flagcolname = None, propertycolnames = "lcmanip", verbose = True, absmagerrs=False):
	"""
	The relaxed way to import lightcurves, especially those from cosmouline, provided they come as rdb files.
	Don't care about column indices, just give the column headers that you want to read.

	Propertycolnames is a list of column names that I will add as properties.
	Possible settings :
	"lcmanip" : (default) I will import the standard cols from lcmanip / cosmouline, if those are available.
	None : I will not import any properties
	["property1", "property2"] : just give a list of properties to import.

	The default column names are those used by the method :py:meth:`pycs.gen.lc.lightcurve.rdbexport` : "mhjd", "mag", "magerr", "mask".

	We use flexibleimport under the hood.

	"""

	if verbose : print "Checking header of \"%s\"..." % (os.path.basename(filepath))
	rdbfile = open(filepath, "r")
	rdbfilelines = rdbfile.readlines()
	rdbfile.close()
	headers = rdbfilelines[0].split()
	underlines = rdbfilelines[1].split()
	if map(len, headers) != map(len, underlines):
		raise RuntimeError, "Error in parsing headers"
	#headerindices = np.array(range(len(headers))) + 1 # +1 as we use human convention in flexibleimport

	# We build the default property names :
	if propertycolnames == "lcmanip": # Then we put it the default set, but only if available.
		propertycolnames = set(["telescope", "fwhm", "ellipticity", "airmass", "relskylevel", "normcoeff", "nbimg"]).intersection(set(headers))

	# We check if the headers you want are available :
	checknames = [mhjdcolname, magcolname, magerrcolname]
	if flagcolname : checknames.append(flagcolname)
	if propertycolnames : checknames.extend(propertycolnames)
	for name in checknames:
		if name not in headers:
			raise RuntimeError, 'I cannot find a column named "%s" in your file !' % (name)


	jdcol = headers.index(mhjdcolname) + 1 # +1 as we use human convention in flexibleimport
	magcol = headers.index(magcolname) + 1
	errcol = headers.index(magerrcolname) + 1

	if flagcolname != None :
		flagcol = headers.index(flagcolname) + 1
	else:
		flagcol = None


	if propertycolnames != None :
		propertycols = {}
		for propertycolname in propertycolnames:
			propertycols[propertycolname] = headers.index(propertycolname) + 1
	else:
		propertycols = None

	newlc = flexibleimport(filepath, jdcol=jdcol, magcol=magcol, errcol=errcol, startline=3, flagcol=flagcol, propertycols=propertycols, telescopename=telescopename, object=object, verbose=verbose, absmagerrs=absmagerrs)
	newlc.plotcolour = plotcolour
	return newlc


def multidisplay(setlist=[],
	title=None, style=None, showlegend=True, showlogo=False, logopos="left", showdates=False, showdelays=False, nicefont=False, text=None,
	jdrange=None, magrange=None, initfigsize=(12,8), plotsize=(0.08, 0.96, 0.09, 0.95), showgrid=False,
	markersize=6, showerrorbars=True, errorbarcolour = "#BBBBBB", capsize=3, knotsize=0.015,
	legendloc = "best", showspldp=False, colourprop=None, hidecolourbar=False, transparent=False,
	collapseref=False, jdmintickstep=100, magmintickstep=0.2, filename="screen", verbose=False):

	#todo: this is a mess ! clean
	"""
	Function that uses matplotlib to plot a **list** of lightcurves/splines/GPRs, either on screen or into a file.
	It uses lightcurve attributes such as ``lc.plotcolour`` and ``lc.showlabels``, and displays masked points
	as black circles. It's certainly a key function of pycs.
	You can also put tuples (lightcurve, listofseasons) in the lclist, and the seasons will be drawn.
	This function is intended both for interactive exploration and for producing publication plots.

	:param lclist: A list of lightcurve objects [lc1, lc2, ...] you want to plot.
	:type lclist: list

	:param splist: A list of spline or rslc (e.g., GPR) objects to display on top of the data points.
	:type splist: list

	:param title: Adds a title to the plot, center top of the axes, usually used for lens names.
		To nicely write a lens name, remember to use raw strings and LaTeX mathrm, e.g. :
		``title = r"$\mathrm{RX\,J1131-1231}$"``
	:type title: string

	:param style: A shortcut to produce specific kinds of stylings for the plots.
		Available styles:

			* ``homepagepdf`` : for cosmograil homepage, ok also with long magnitude labels (like -13.2)

	:type style: string

	:param showlegend: Automatic legend (too technical/ugly for publication plots, uses str(lightcurve))
	:type showlegend: boolean

	:param showlogo: Adds an EPFL logo + www.cosmograil.org on the plot.
	:type showlogo: boolean

	:param logopos: Where to put it, "left" or "right"
	:type logopos: string

	:param showdates: If True, the upper x axis will show years and 12 month minor ticks.
	:type showdates: boolean

	:param showdelays: If True, the relative delays between the curves are written on the plot.
	:type showdelays: boolean

	:param nicefont: Sets default to serif fonts (terrible implementation, but works)
	:type nicefont: boolean

	:param text:
		Generic text that you want to display on top of the plot, in the form : [line1, line2, line3 ...]
		where line_i is (x, y, text, kwargs) where kwargs is e.g. {"fontsize":18} and x and y are relative positions (from 0 to 1).
	:type text: list

	:param jdrange: Range of jds to plot, e.g. (53000, 54000).
	:type jdrange: tuple

	:param magrange: Range of magnitudes to plot, like ``magrange = [-11, -13]``.
		If you give only a float, like ``magrange=1.0``, I'll plot +/- this number around the mean curve level
		Default is None -> automatic mag range.
	:type magrange: tuple

	:param figsize: Figure size (width, height) in inches, for interactive display or savefig.
	:type figsize: tuple

	:param plotsize: Position of the axes in the figure (left, right, bottom, top), default is (0.065, 0.96, 0.09, 0.95).
	:type plotsize: tuple

	:param showgrid: Show grid, that is vertical lines only, one for every year.
	:type showgrid: boolean

	:param markersize: Size of the data points, default is 6
	:type markersize: float

	:param showerrorbars: If False, the ploterrorbar settings of the lightcurves are disregared and no error bars are shown.
	:type showerrorbars: boolean

	:param errorbarcolour: Color for the error bars
	:type errorbarcolour: string

	:param capsize: Size of the error bar "ticks"
	:type capsize: float

	:param knotsize: For splines, the length of the knot ticks, in magnitudes.
	:type knotsize: float

	:param legendloc: Position of the legend. It is passed to matplotlib legend. It can be useful to specify this if you want to make animations etc.

			* string	int
			* upper right	1
			* upper left	2
			* lower left	3
			* lower right	4
			* right	5
			* center left	6
			* center right	7
			* lower center	8
			* upper center	9
			* center	10

	:type legendloc: string or int

	:param showspldp: Show the acutal data points of spline objects (for debugging etc)
	:type showspldp: boolean

	:param colourprop: If this is set I will make a scatter plot with points coloured according to the given property, disregarding the lightcurve attribute plotcolour.
		Format : (property_name, display_name, minval, maxval), where display_name is a "nice" version of property_name, like "FWHM [arcsec]" instead of "seeing".
		Note that this property will be used in terms of a float. So you cannot use this for properties that are not floats.
	:type colourprop: tuple

	:param hidecolourbar: Set to True to hide the colourbar for the colourprop
	:type hidecolourbar: boolean

	:param transparent: Set transparency for the plot, if saved using filename
	:type transparent: boolean

	:param collapseref: Plot one single dashed line as the reference for the microlensing.
		Use this if you would otherwise get ugly overplotted dashed lines nearly at the same level ...
		This option is a bit ugly, as it does not correct the actual microlensing curves for the collapse.
	:type collapseref: boolean

	:param jdmintickstep: Minor tick step for jd axis
	:type jdmintickstep: float

	:param magmintickstep: Minor tick step for mag axis
	:type magmintickstep: float

	:param filename: If this is not "screen", I will save the plot to a file instead of displaying it. Try e.g. "test.png" or "test.pdf". Success depends on your matplotlib backend.
	:type filename: string

	:param verbose: Set to True if you want me to print some details while I process the curves
	:type verbose: boolean


	"""

	import matplotlib as mpl
	import matplotlib.pyplot as plt
	import matplotlib.font_manager as fm
	from matplotlib.ticker import MultipleLocator, FormatStrFormatter
	import matplotlib.dates
	import matplotlib.lines

	if style == None:
		pass
	elif style=="homepagepdf":
		figsize=(10,5)
		plotsize=(0.09, 0.97, 0.10, 0.95)
		showlogo=True
		nicefont=False
		showdelays=False
		showlegend=False
		showdates=True
		errorbarcolour="#777777"
		markersize=3.0
		capsize=0
		jdmintickstep=50
		magmintickstep=0.2
		showgrid=True
	else:
		raise RuntimeError("I do not know the style %s" % (style))


	# warning discarded : user must be careful (bad...)
	#if not (isinstance(lclist, list) or isinstance(lclist, tuple)):
		#raise TypeError, "Hey, give me a LIST of lightcurves !"

	if colourprop != None:
		(colourpropname, colournicename, colourminval, colourmaxval) = colourprop

	labelfontsize = 14
	if nicefont:
		#mpl.rcParams['font.size'] = 20
		mpl.rcParams['font.family'] = 'serif'
		#labelfontsize = 20
	else:
		labelfontsize = 14


	figsize = [initfigsize[0]]
	if len(setlist) != 1:
		figsize.append(15)
	else:
		figsize.append(initfigsize[1])
	fig = plt.figure(figsize=figsize)	# sets figure size
	fig.subplots_adjust(left = plotsize[0], right=plotsize[1], bottom=plotsize[2], top=plotsize[3])

	axes = plt.gca()

	if verbose : print "Plotting %i lightcurves and %i splines ..." % (len(lclist), len(splist))

	reflevels = [] # only used for collapseref



	for ind,set in enumerate(setlist):

		lclist = set[0]
		if len(set) != 3:
			if len(set) == 1:
				splist = []
				comment = ''
			else:
				if len(set) == 2:
					if isinstance(set[1],str):
						splist  = []
						comment = set[1]
					else:
						splist  = set[1]
						comment = ''
				else:
					raise TypeError, 'Give me a good setlist ! (i.e. [lcs,[spline],somecaption])'
		else:
			splist  = set[1]
			comment = set[2]


		#initialize the subplot we will draw on:
		plt.subplot(len(setlist),1,ind+1)

		# The lightcurves :
		for curve in lclist :

			if type(curve).__name__ == 'tuple': # then we have both a lightcurve and a season to plot

				actualcurve = curve[0]
				curveseasons = curve[1]

				if not isinstance(curveseasons, list):
					raise TypeError, "lc.display wants LISTs of seasons, not individual seasons !"
				for curveseason in curveseasons:
					# the x lims :
					(x1, x2) = curveseason.getjdlims(actualcurve)
					# for y, lets take the median of that season
					y = np.median(actualcurve.getmags()[curveseason.indices])

					# we make this robust even with old versions of matplotlib, so no fancy arrows here.
					plt.plot([x1, x2], [y, y], color = actualcurve.plotcolour, dashes = (1,1))
					axes.annotate(str(curveseason), ((x1 + x2)/2.0, y), xytext=(-50, -15), textcoords='offset points', size=10, color = actualcurve.plotcolour)
					#plt.axvline(seasonjdlims[0], color = curve[0].plotcolour, dashes = (5,5))
					#plt.axvline(seasonjdlims[1], color = curve[0].plotcolour, dashes = (5,5))

				curve = curve[0] # for the rest of this loop, curve is now only the lightcurve.


			if verbose: print "#   %s -> %s\n\t%s" % (curve, str(curve.plotcolour), "\n\t".join(curve.commentlist))
			#if verbose and (curve.ml != None):
			#	print curve.ml.longinfo()

			tmpjds = curve.getjds()
			tmpmags = curve.getmags() # to avoid calculating the microlensing each time we need it

			if colourprop != None:
				scattervalues = np.array([float(propertydict[colourpropname]) for propertydict in curve.properties])
				plt.scatter(tmpjds, tmpmags, s=markersize, c=scattervalues, vmin=colourminval, vmax=colourmaxval, edgecolors="None")
			else:

				if curve.ploterrorbars and showerrorbars:
					plt.errorbar(tmpjds, tmpmags, curve.magerrs, fmt=".", markersize = markersize, markeredgecolor=curve.plotcolour, color=curve.plotcolour, ecolor=errorbarcolour, capsize=capsize, label=str(curve), elinewidth=0.5)
					#plt.errorbar(tmpjds, tmpmags, curve.magerrs, linestyle="-", marker=".", color=curve.plotcolour, ecolor="#BBBBBB", label=str(curve))
					#
				else :
					plt.plot(tmpjds, tmpmags, marker=".", markersize = markersize, linestyle="None", markeredgecolor=curve.plotcolour, color=curve.plotcolour, label=str(curve))

			# We plot little round circles around masked points.
			plt.plot(tmpjds[curve.mask == False], tmpmags[curve.mask == False], linestyle="None", marker="o", markersize=8., markeredgecolor="black", markerfacecolor="None", color="black")

			# And now we want to graphically display the microlensing in a nice way. This costs some cpu but anyway
			# for a display it's fine.
			#if curve.ml != None and curve.hideml == False:
			if curve.ml != None:
				if curve.ml.mltype in ["poly", "leg"]:
					for sfct in curve.ml.mllist:
						smoothml = sfct.smooth(curve)
						if not collapseref:
							plt.plot(smoothml['jds'], smoothml['refmags'], color=curve.plotcolour, dashes = ((3,3))) # the new ref
						else:
							reflevels.append(np.mean(smoothml['refmags']))
						plt.plot(smoothml['jds'], smoothml['refmags'] + smoothml['ml'], color=curve.plotcolour)
				if curve.ml.mltype == "spline":
					smoothml = curve.ml.smooth(curve)

					if not collapseref:
						plt.plot(smoothml['jds'], np.zeros(smoothml["n"]) + smoothml['refmag'], color=curve.plotcolour, dashes = ((3,3))) # the new ref
					else:
						reflevels.append(smoothml['refmag'])

					plt.plot(smoothml['jds'], smoothml['refmag'] + smoothml['ml'], color=curve.plotcolour)
					# We want to overplot the knots
					#plt.plot(smoothml['knotjds'], smoothml['knotmags'] + smoothml["refmag"], color=curve.plotcolour, linestyle="none", marker=".")
					if getattr(curve.ml.spline, "showknots", True) == True:
						plt.errorbar(smoothml['knotjds'], smoothml['knotmags'] + smoothml["refmag"], knotsize*np.ones(len(smoothml['knotjds'])), capsize=0, ecolor=curve.plotcolour, linestyle="none", marker="", elinewidth=1.5)

			# Labels if wanted :
			if curve.showlabels:
				for i, label in enumerate(curve.labels):
					if label != "":
						#axes.annotate(label, (curve.jds[i], curve.mags[i]))
						if len(label) > 4: # Probably jd labels, we write vertically :
							axes.annotate(label, (tmpjds[i], tmpmags[i]), xytext=(-3, -70), textcoords='offset points',size=12, color = curve.plotcolour, rotation = 90)
						else:	# horizontal writing
							axes.annotate(label, (tmpjds[i], tmpmags[i]), xytext=(7, -6), textcoords='offset points',size=12, color = curve.plotcolour)

		if collapseref and len(reflevels) != 0:
			print "WARNING : collapsing the refs %s" % (reflevels)
			plt.axhline(np.mean(np.array(reflevels)), color="gray", dashes = ((3,3))) # the new ref


		# The supplementary objects
		if len(splist) != 0:
			for stuff in splist:

				# We do some stupid type checking. But I like this as it does not require
				# to import spline and gpr etc.
				if hasattr(stuff, "knottype"): # Then it's a spline
					spline = stuff
					if verbose: print "#   %s -> %s" % (str(spline), str(spline.plotcolour))

					npts = (spline.datapoints.jds[-1] - spline.datapoints.jds[0])*2.0
					xs = np.linspace(spline.datapoints.jds[0], spline.datapoints.jds[-1], npts)
					ys = spline.eval(jds = xs)
					plt.plot(xs, ys, "-", color=spline.plotcolour, zorder=+20, label=str(spline))
					# For the knots, we might not want to show them (by default we do show them) :
					if getattr(spline, "showknots", True) == True:
						#plt.errorbar(spline.getinttex(), spline.eval(jds = spline.getinttex()), 0.015*np.ones(len(spline.getinttex())), capsize=0, ecolor=spline.plotcolour, linestyle="none", marker="", elinewidth=1.5, zorder=40, barsabove=True)

						ax = plt.gca()
						knotxs = spline.getinttex()
						knotys = spline.eval(knotxs)
						for (knotx, knoty) in zip(knotxs, knotys):
							l = matplotlib.lines.Line2D([knotx,knotx],[knoty-knotsize,knoty+knotsize], zorder=30, linewidth=1.5, color=spline.plotcolour)
 							ax.add_line(l)

					if showspldp: # The datapoints of the spline (usually not shown)
						plt.plot(spline.datapoints.jds, spline.datapoints.mags, marker = ",", linestyle="None", color=spline.plotcolour, zorder=-20)

#				if hasattr(stuff, "regfct"): # Then it's a GPR
#
# 					gpr = stuff
# 					npts = (gpr.jds[-1] - gpr.jds[0])*2.0
# 					xs = np.linspace(gpr.jds[0], gpr.jds[-1], npts)
# 					(ys, yerrs) = gpr.regfct(xs)
# 					#print "regfct evaluated"
# 					plt.plot(xs, ys, "-", color=gpr.plotcolour, zorder=+20, label=str(gpr))
# 					xf = np.concatenate((xs, xs[::-1]))
#       	  			yf = np.concatenate((ys+yerrs, (ys-yerrs)[::-1]))
#	         			plt.fill(xf, yf, facecolor = gpr.plotcolour, alpha=0.1, edgecolor = (1,1,1))

 				if hasattr(stuff, "pd"): # Then it's a rslc

					rs = stuff
					#plt.plot(rs.getjds(), rs.mags, "-.", color=rs.plotcolour)
					plt.plot(rs.getjds(), rs.mags, "-", color=rs.plotcolour)
					xf = np.concatenate((rs.getjds(), rs.getjds()[::-1]))
        				yf = np.concatenate((rs.mags+rs.magerrs, (rs.mags-rs.magerrs)[::-1]))
        				plt.fill(xf, yf, facecolor = rs.plotcolour, alpha=0.2, edgecolor = (1,1,1), label=str(rs))


		# Astronomers like minor tick marks :
		minorxLocator = MultipleLocator(jdmintickstep)
		axes.xaxis.set_minor_locator(minorxLocator)
		#minorLocator = MultipleLocator(1) # so to have a tick every day
		#axes.xaxis.set_minor_locator(minorLocator)


		# Something for astronomers only : we invert the y axis direction !
		axes.set_ylim(axes.get_ylim()[::-1])


		if colourprop != None and hidecolourbar == False:
			cbar = plt.colorbar(orientation='vertical', shrink=1.0, fraction=0.065, pad=0.025)
			cbar.set_label(colournicename)

		# And we make custom title :

		if title == "None" or title == None or title == "none":
			#plt.title("Lightcurves", fontsize=18)
			pass
		else:
			#plt.title(title, fontsize=18)
			axes.annotate(title, xy=(0.5, 1.0), xycoords='axes fraction', xytext=(0, -10),
				textcoords='offset points', ha='center', va='top', fontsize=30)

		if jdrange != None:
			plt.xlim(jdrange[0], jdrange[1])


		plt.xlabel("HJD - 2400000.5 [day]", fontsize = labelfontsize)
		plt.ylabel("Magnitude (relative)", fontsize = labelfontsize)

		if showdelays:
			txt = getnicetimedelays(lclist, separator="\n")

			#axes.annotate(txt, xy=(0.0, 1.0), xycoords='axes fraction', xytext=(6, -6),
				#textcoords='offset points', ha='left', va='top')


			#stupid way of doing things...but it works and I'm not in the mood of doing it cleaner right now...
			if len(setlist)==1:
				ytxtcoord = [0.99]
			if len(setlist)==2:
				ytxtcoord = [0.99,0.445]
			if len(setlist)==3:
				ytxtcoord = [0.99,0.63,0.28]
			if len(setlist)==4:
				ytxtcoord = [0.99,0.73,0.47,0.21]
			if len(setlist)==5:
				ytxtcoord = [0.995,0.787,0.582,0.374, 0.167]
			if len(setlist)==6:
				ytxtcoord = [0.995,0.826, 0.655,  0.48, 0.31, 0.14]

			# you want more subplots ? Feel free to add ytxtcoord...
			'''
			plt.text(0.01, ytxtcoord[ind], txt,
				horizontalalignment='left', verticalalignment='top',
				transform = axes.transAxes)

			plt.text(0.4, ytxtcoord[ind], comment,
				horizontalalignment='left', verticalalignment='top',
				transform = axes.transAxes)
			legendloc = 1
			if verbose:
				print "Delays between plotted curves :"
				print txt
			'''


		if showlegend and (len(lclist) > 0 or len(splist) > 0):
			plt.legend(loc = legendloc, numpoints = 1, prop = fm.FontProperties(size = 12))

		if magrange != None:
			if type(magrange) == float or type(magrange) == int:
				# We find the mean mag of the stuff to plot :
				allmags = []
				for l in lclist:
					allmags.extend(l.getmags())
				meanlevel = np.mean(np.array(allmags))
				plt.ylim(meanlevel+magrange, meanlevel-magrange)
			else:
				plt.ylim(magrange[0], magrange[1])

		if showdates: # Be careful when you change something here, it could mess up the axes.
			# Especially watch out when you change the plot range.
			# This showdates stuff should come at the very end
			minjd = axes.get_xlim()[0]
			maxjd = axes.get_xlim()[1]
			#axes.set_xlim(minjd, maxjd)
			yearx = axes.twiny()
			yearxmin = util.datetimefromjd(minjd + 2400000.5)
			yearxmax = util.datetimefromjd(maxjd + 2400000.5)
			yearx.set_xlim(yearxmin, yearxmax)
			yearx.xaxis.set_minor_locator(matplotlib.dates.MonthLocator())
			yearx.xaxis.set_major_locator(matplotlib.dates.YearLocator())
			yearx.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))
			yearx.xaxis.tick_top()
			#yearx.set_xlabel("Date")

		minoryLocator = MultipleLocator(magmintickstep)
		axes.yaxis.set_minor_locator(minoryLocator)


		if showgrid:
			plt.grid(zorder=2000)

		if text != None:
			for line in text:
				axes.text(line[0], line[1], line[2], transform=axes.transAxes, fontsize=18)

		if showlogo:

			# The EPFL logo :
			from matplotlib._png import read_png
			from matplotlib.offsetbox import OffsetImage, AnnotationBbox
			logodir = os.path.dirname(__file__)
			im = read_png(os.path.join(logodir, "epfl.png"))
			imagebox = OffsetImage(im, zoom=0.13, interpolation="sinc", resample = True)

			if logopos == "left":
				ab = AnnotationBbox(imagebox,  xy=(0.0, 0.0), xycoords='axes pixels', xybox = (52, 30),
					boxcoords="offset points",
					pad=0.0, frameon=False
				)
				axes.add_artist(ab)
				axes.annotate("COSMOGRAIL.org", xy=(0.0, 0.0), xycoords='axes fraction', fontsize=16, xytext=(105, 7),
					textcoords='offset points', ha='left', va='bottom', color="gray")

			if logopos == "right":
				ab = AnnotationBbox(imagebox,  xy=(1.0, 0.0), xycoords='axes fraction', xybox = (-200, 30),
					boxcoords="offset points",
					pad=0.0, frameon=False
				)
				axes.add_artist(ab)
				axes.annotate("COSMOGRAIL.org", xy=(1.0, 0.0), xycoords='axes fraction', fontsize=16, xytext=(-10, 7),
					textcoords='offset points', ha='right', va='bottom', color="gray")

			# Alternative possibility (just to keep the idea) :
			"""
			try:
				import Image
			except ImportError:
				print "Couldn't import PIL ! Therefore I won't be able to display the cosmograil logo."
			else:
				im = Image.open('epfl.png')
				height = im.size[1]
				print height
				im = np.array(im).astype(np.float) / 255
				fig.figimage(im, 0, fig.bbox.ymax - height)
				# With newer (1.0) versions of matplotlib, you can
				# use the "zorder" kwarg to make the image overlay
				# the plot, rather than hide behind it... (e.g. zorder=10)
				fig.figimage(im, 0, fig.bbox.ymax - height)
			"""

	if filename == "screen":
		plt.show()
	else:
		plt.savefig(filename, transparent=transparent)
		print "Plot written to %s" % filename
		plt.close() # this seems important so that the plot is not displayed when a next plt.show() is called.


def display(lclist=[], splist=[],
	title=None, titlexpos=None, style=None, showlegend=True, showlogo=False, logopos="left", showdates=False, showdelays=False, nicefont=False, text=None, keeponlygrid=False,
	jdrange=None, magrange=None, figsize=(12,8), plotsize=(0.08, 0.96, 0.09, 0.95), showgrid=False,
	markersize=6, showerrorbars=True, showdatapoints=True, errorbarcolour = "#BBBBBB", capsize=3, knotsize=0.015,
	legendloc = "best", showspldp=False, colourprop=None, hidecolourbar=False, transparent=False,
	collapseref=False, hidecollapseref=False, jdmintickstep=100, magmintickstep=0.2, filename="screen", showinsert=None, insertname=None, verbose=False, ax=None):
	"""
	Function that uses matplotlib to plot a **list** of lightcurves/splines/GPRs, either on screen or into a file.
	It uses lightcurve attributes such as ``lc.plotcolour`` and ``lc.showlabels``, and displays masked points
	as black circles. It's certainly a key function of pycs.
	You can also put tuples (lightcurve, listofseasons) in the lclist, and the seasons will be drawn.
	This function is intended both for interactive exploration and for producing publication plots.

	:param lclist: A list of lightcurve objects [lc1, lc2, ...] you want to plot.
	:type lclist: list

	:param splist: A list of spline or rslc (e.g., GPR) objects to display on top of the data points.
	:type splist: list

	:param title: Adds a title to the plot, center top of the axes, usually used for lens names.
		To nicely write a lens name, remember to use raw strings and LaTeX mathrm, e.g. :
		``title = r"$\mathrm{RX\,J1131-1231}$"``
	:type title: string

	:param titlexpos: Specify where you want your to anchor the center of your title in the x axis. the value is in the x-axis fraction.  Default is center. (x = 0.5)
	:type titlexpos: float

	:param style: A shortcut to produce specific kinds of stylings for the plots.
		Available styles:

			* ``homepagepdf`` : for cosmograil homepage, ok also with long magnitude labels (like -13.2)
			
	:type style: string

	:param showlegend: Automatic legend (too technical/ugly for publication plots, uses str(lightcurve))
	:type showlegend: boolean
	
	:param showlogo: Adds an EPFL logo + www.cosmograil.org on the plot.
	:type showlogo: boolean
	
	:param logopos: Where to put it, "left" or "right"
	:type logopos: string

	:param showdates: If True, the upper x axis will show years and 12 month minor ticks.
	:type showdates: boolean

	:param showdelays: If True, the relative delays between the curves are written on the plot.
	:type showdelays: boolean
	
	:param nicefont: Sets default to serif fonts (terrible implementation, but works)
	:type nicefont: boolean
	
	:param text:
		Generic text that you want to display on top of the plot, in the form : [line1, line2, line3 ...]
		where line_i is (x, y, text, kwargs) where kwargs is e.g. {"fontsize":18} and x and y are relative positions (from 0 to 1). 
	:type text: list
	
	:param jdrange: Range of jds to plot, e.g. (53000, 54000).
	:type jdrange: tuple
	
	:param magrange: Range of magnitudes to plot, like ``magrange = [-11, -13]``.
		If you give only a float, like ``magrange=1.0``, I'll plot +/- this number around the mean curve level
		Default is None -> automatic mag range.
	:type magrange: tuple
	
	:param figsize: Figure size (width, height) in inches, for interactive display or savefig.
	:type figsize: tuple
	
	:param plotsize: Position of the axes in the figure (left, right, bottom, top), default is (0.065, 0.96, 0.09, 0.95).
	:type plotsize: tuple
	
	:param showgrid: Show grid, that is vertical lines only, one for every year.
	:type showgrid: boolean

	:param markersize: Size of the data points, default is 6
	:type markersize: float
	
	:param showerrorbars: If False, the ploterrorbar settings of the lightcurves are disregared and no error bars are shown.
	:type showerrorbars: boolean

	:param showdatapoints: If False, no data points are shown. Useful if you want e.g. to plot only the microlensing
	:type showerrorbars: boolean

	:param keeponlygrid: If True, keeps the yearly grid from showdates but do not display the dates above the plot.
	:type keeponlygrid: boolean

	:param showinsert: If True, display the insertname image in the top-right corner of the main image
	:type showinsert: boolean

	:param insertname: path to the image you want to insert
	:type instername: string

	:param errorbarcolour: Color for the error bars
	:type errorbarcolour: string
	
	:param capsize: Size of the error bar "ticks"
	:type capsize: float

	:param knotsize: For splines, the length of the knot ticks, in magnitudes.
	:type knotsize: float
	
	:param legendloc: Position of the legend. It is passed to matplotlib legend. It can be useful to specify this if you want to make animations etc.
		
			* string	int
			* upper right	1
			* upper left	2
			* lower left	3
			* lower right	4
			* right	5


			* center left	6
			* center right	7
			* lower center	8
			* upper center	9
			* center	10
	
	:type legendloc: string or int

	:param showspldp: Show the acutal data points of spline objects (for debugging etc)
	:type showspldp: boolean
	
	:param colourprop: If this is set I will make a scatter plot with points coloured according to the given property, disregarding the lightcurve attribute plotcolour.
		Format : (property_name, display_name, minval, maxval), where display_name is a "nice" version of property_name, like "FWHM [arcsec]" instead of "seeing".
		Note that this property will be used in terms of a float. So you cannot use this for properties that are not floats.
	:type colourprop: tuple
	
	:param hidecolourbar: Set to True to hide the colourbar for the colourprop
	:type hidecolourbar: boolean

	:param transparent: Set transparency for the plot, if saved using filename
	:type transparent: boolean

	:param collapseref: Plot one single dashed line as the reference for the microlensing.
		Use this if you would otherwise get ugly overplotted dashed lines nearly at the same level ...
		This option is a bit ugly, as it does not correct the actual microlensing curves for the collapse.
	:type collapseref: boolean
	
	:param jdmintickstep: Minor tick step for jd axis
	:type jdmintickstep: float
	
	:param magmintickstep: Minor tick step for mag axis
	:type magmintickstep: float
	
	:param filename: If this is not "screen", I will save the plot to a file instead of displaying it. Try e.g. "test.png" or "test.pdf". Success depends on your matplotlib backend.
	:type filename: string

	:param ax: if not None, I will return what I plot in the given matplotlib axe you provide me with instead of plotting it.
	:type ax: matplotlib axes

	:param verbose: Set to True if you want me to print some details while I process the curves
	:type verbose: boolean


	"""
	
	import matplotlib as mpl
	import matplotlib.pyplot as plt
	import matplotlib.font_manager as fm
	from matplotlib.ticker import MultipleLocator, FormatStrFormatter
	import matplotlib.dates
	import matplotlib.lines
	
	if style == None:
		pass
	elif style in ["homepagepdf", "homepagepdfnologo"]:
		figsize=(10,5)
		plotsize=(0.09, 0.97, 0.10, 0.95)
		showlogo=True
		if style == "homepagepdfnologo":
			showlogo=False
		nicefont=False
		showdelays=False
		showlegend=False
		showdates=True
		errorbarcolour="#777777"
		markersize=5.0
		capsize=0
		jdmintickstep=50
		magmintickstep=0.2
		showgrid=True
		transparent=False

	elif "2m2" in style:
		figsize=(10,5)
		plotsize=(0.09, 0.97, 0.10, 0.95)
		showlogo = True
		if "nologo" in style:
			showlogo=False
		logopos="right"
		nicefont=True
		showdelays=False
		if "showdelays" in style:
			showdelays = True
		showlegend=False
		#showdates=False
		errorbarcolour="#777777"
		markersize=7.0
		capsize=0
		jdmintickstep=50
		if "largeticks" in style:
			jdtickstep=100
		magmintickstep=0.2
		#showgrid=False
		transparent=False

		
		
	elif style=="posterpdf":
		figsize=(10,5.5)
		plotsize=(0.09, 0.97, 0.10, 0.95)
		showlogo=False
		nicefont=True
		showdelays=False
		showlegend=False
		showdates=True
		errorbarcolour="#777777"
		markersize=4.0
		capsize=0
		jdmintickstep=50
		magmintickstep=0.2
		showgrid=False
		transparent=True
		title=None

	elif style=="internal":
		figsize=(10,5.5)
		plotsize=(0.09, 0.97, 0.10, 0.95)
		showlogo=False
		nicefont=True
		showdelays=False
		showdates=True
		errorbarcolour="#777777"
		markersize=5.0
		capsize=0
		jdmintickstep=50
		magmintickstep=0.2
		showgrid=True
		transparent=False



	else:
		raise RuntimeError("I do not know the style %s" % (style))
	
	
	
	if not (isinstance(lclist, list) or isinstance(lclist, tuple)):
		raise TypeError, "Hey, give me a LIST of lightcurves !"

	if colourprop != None:
		(colourpropname, colournicename, colourminval, colourmaxval) = colourprop
	
	labelfontsize = 14
	if nicefont:
		#mpl.rcParams['font.size'] = 20
		mpl.rcParams['font.family'] = 'serif'
		#labelfontsize = 20
	else:
		labelfontsize = 14

	if ax == None:
		fig = plt.figure(figsize=figsize)	# sets figure size
		fig.subplots_adjust(left = plotsize[0], right=plotsize[1], bottom=plotsize[2], top=plotsize[3])
		axes = plt.gca()

	else:
		ihaveax = True
		axes = ax
	
	if verbose : print "Plotting %i lightcurves and %i splines ..." % (len(lclist), len(splist))

	reflevels = [] # only used for collapseref

	# The lightcurves :
	for curve in lclist :
		if showdatapoints:
			if type(curve).__name__ == 'tuple': # then we have both a lightcurve and a season to plot

				actualcurve = curve[0]
				curveseasons = curve[1]

				if not isinstance(curveseasons, list):
					raise TypeError, "lc.display wants LISTs of seasons, not individual seasons !"
				for curveseason in curveseasons:
					# the x lims :
					(x1, x2) = curveseason.getjdlims(actualcurve)
					# for y, lets take the median of that season
					y = np.median(actualcurve.getmags()[curveseason.indices])

					# we make this robust even with old versions of matplotlib, so no fancy arrows here.
					axes.plot([x1, x2], [y, y], color = actualcurve.plotcolour, dashes = (1,1))
					axes.annotate(str(curveseason), ((x1 + x2)/2.0, y), xytext=(-50, -15), textcoords='offset points', size=10, color = actualcurve.plotcolour)
					#plt.axvline(seasonjdlims[0], color = curve[0].plotcolour, dashes = (5,5))
					#plt.axvline(seasonjdlims[1], color = curve[0].plotcolour, dashes = (5,5))

				curve = curve[0] # for the rest of this loop, curve is now only the lightcurve.


			if verbose: print "#   %s -> %s\n\t%s" % (curve, str(curve.plotcolour), "\n\t".join(curve.commentlist))
			#if verbose and (curve.ml != None):
			#	print curve.ml.longinfo()

			tmpjds = curve.getjds()
			tmpmags = curve.getmags() # to avoid calculating the microlensing each time we need it

			if colourprop != None:
				scattervalues = np.array([float(propertydict[colourpropname]) for propertydict in curve.properties])
				axes.scatter(tmpjds, tmpmags, s=markersize, c=scattervalues, vmin=colourminval, vmax=colourmaxval, edgecolors="None")
			else:

				if curve.ploterrorbars and showerrorbars:
					axes.errorbar(tmpjds, tmpmags, curve.magerrs, fmt=".", markersize = markersize, markeredgecolor=curve.plotcolour, color=curve.plotcolour, ecolor=errorbarcolour, capsize=capsize, label=str(curve), elinewidth=0.5)
					#plt.errorbar(tmpjds, tmpmags, curve.magerrs, linestyle="-", marker=".", color=curve.plotcolour, ecolor="#BBBBBB", label=str(curve))
					#
				else :
					axes.plot(tmpjds, tmpmags, marker=".", markersize = markersize, linestyle="None", markeredgecolor=curve.plotcolour, color=curve.plotcolour, label=str(curve))

			# We plot little round circles around masked points.
			axes.plot(tmpjds[curve.mask == False], tmpmags[curve.mask == False], linestyle="None", marker="o", markersize=8., markeredgecolor="black", markerfacecolor="None", color="black")

		# And now we want to graphically display the microlensing in a nice way. This costs some cpu but anyway
		# for a display it's fine.
		#if curve.ml != None and curve.hideml == False:
		if curve.ml != None:
			if curve.ml.mltype in ["poly", "leg"]:
				for sfct in curve.ml.mllist:
					smoothml = sfct.smooth(curve)
					if not collapseref:
						axes.plot(smoothml['jds'], smoothml['refmags'], color=curve.plotcolour, dashes = ((3,3))) # the new ref
					else:
						reflevels.append(np.mean(smoothml['refmags']))
					axes.plot(smoothml['jds'], smoothml['refmags'] + smoothml['ml'], color=curve.plotcolour)
			if curve.ml.mltype == "spline":
				smoothml = curve.ml.smooth(curve)

				if not collapseref:
					axes.plot(smoothml['jds'], np.zeros(smoothml["n"]) + smoothml['refmag'], color=curve.plotcolour, dashes = ((3,3))) # the new ref
				else:
					reflevels.append(smoothml['refmag'])

				if hasattr(curve, "hideml"):
					if not curve.hideml:
						axes.plot(smoothml['jds'], smoothml['refmag'] + smoothml['ml'], color=curve.plotcolour)
				else:
					axes.plot(smoothml['jds'], smoothml['refmag'] + smoothml['ml'], color=curve.plotcolour)
				# We want to overplot the knots
				if hasattr(curve, "hideml"):
					if not curve.hideml:
						if getattr(curve.ml.spline, "showknots", True) == True:
							axes.errorbar(smoothml['knotjds'], smoothml['knotmags'] + smoothml["refmag"], knotsize*np.ones(len(smoothml['knotjds'])), capsize=0, ecolor=curve.plotcolour, linestyle="none", marker="", elinewidth=1.5)
				else:
					if getattr(curve.ml.spline, "showknots", True) == True:
						axes.errorbar(smoothml['knotjds'], smoothml['knotmags'] + smoothml["refmag"], knotsize*np.ones(len(smoothml['knotjds'])), capsize=0, ecolor=curve.plotcolour, linestyle="none", marker="", elinewidth=1.5)


		# Labels if wanted :
		if curve.showlabels:
			for i, label in enumerate(curve.labels):
				if label != "":
					#axes.annotate(label, (curve.jds[i], curve.mags[i]))
					if len(label) > 4: # Probably jd labels, we write vertically :
						axes.annotate(label, (tmpjds[i], tmpmags[i]), xytext=(-3, -70), textcoords='offset points',size=12, color = curve.plotcolour, rotation = 90)
					else:	# horizontal writing
						axes.annotate(label, (tmpjds[i], tmpmags[i]), xytext=(7, -6), textcoords='offset points',size=12, color = curve.plotcolour)

	if collapseref and len(reflevels) != 0:
		print "WARNING : collapsing the refs %s" % (reflevels)
		if not hidecollapseref:
			axes.axhline(np.mean(np.array(reflevels)), color="gray", dashes = ((3,3))) # the new ref

	
	# The supplementary objects
	if len(splist) != 0:
		for stuff in splist:
			
			# We do some stupid type checking. But I like this as it does not require
			# to import spline and gpr etc.
			if hasattr(stuff, "knottype"): # Then it's a spline
				spline = stuff
				if verbose: print "#   %s -> %s" % (str(spline), str(spline.plotcolour))
				
				npts = (spline.datapoints.jds[-1] - spline.datapoints.jds[0])*2.0
				xs = np.linspace(spline.datapoints.jds[0], spline.datapoints.jds[-1], npts)
				ys = spline.eval(jds = xs)
				axes.plot(xs, ys, "-", color=spline.plotcolour, zorder=+20, label=str(spline))
				# For the knots, we might not want to show them (by default we do show them) :
				if getattr(spline, "showknots", True) == True:
					if ax != None:
						ax.errorbar(spline.getinttex(), spline.eval(jds = spline.getinttex()), 0.015*np.ones(len(spline.getinttex())), capsize=0, ecolor=spline.plotcolour, linestyle="none", marker="", elinewidth=1.5, zorder=40, barsabove=True)
						knotxs = spline.getinttex()
						knotys = spline.eval(knotxs)
						for (knotx, knoty) in zip(knotxs, knotys):
							l = matplotlib.lines.Line2D([knotx,knotx],[knoty-knotsize,knoty+knotsize], zorder=30, linewidth=1.5, color=spline.plotcolour)
							ax.add_line(l)
					else:
						axes = plt.gca()
						knotxs = spline.getinttex()
						knotys = spline.eval(knotxs)
						for (knotx, knoty) in zip(knotxs, knotys):
							l = matplotlib.lines.Line2D([knotx,knotx],[knoty-knotsize,knoty+knotsize], zorder=30, linewidth=1.5, color=spline.plotcolour)
							axes.add_line(l)
				
				if showspldp: # The datapoints of the spline (usually not shown)
					axes.plot(spline.datapoints.jds, spline.datapoints.mags, marker = ",", linestyle="None", color=spline.plotcolour, zorder=-20)
			
#			if hasattr(stuff, "regfct"): # Then it's a GPR
# 				
# 				gpr = stuff
# 				npts = (gpr.jds[-1] - gpr.jds[0])*2.0
# 				xs = np.linspace(gpr.jds[0], gpr.jds[-1], npts)
# 				(ys, yerrs) = gpr.regfct(xs)
# 				#print "regfct evaluated"
# 				plt.plot(xs, ys, "-", color=gpr.plotcolour, zorder=+20, label=str(gpr))
# 				xf = np.concatenate((xs, xs[::-1]))
#         			yf = np.concatenate((ys+yerrs, (ys-yerrs)[::-1]))
#         			plt.fill(xf, yf, facecolor = gpr.plotcolour, alpha=0.1, edgecolor = (1,1,1))

 			if hasattr(stuff, "pd"): # Then it's a rslc
				
				rs = stuff
				#plt.plot(rs.getjds(), rs.mags, "-.", color=rs.plotcolour)
				axes.plot(rs.getjds(), rs.mags, "-", color=rs.plotcolour)
				xf = np.concatenate((rs.getjds(), rs.getjds()[::-1]))
        			yf = np.concatenate((rs.mags+rs.magerrs, (rs.mags-rs.magerrs)[::-1]))
        			plt.fill(xf, yf, facecolor = rs.plotcolour, alpha=0.2, edgecolor = (1,1,1), label=str(rs))
				
	
	# Astronomers like minor tick marks :
	minorxLocator = MultipleLocator(jdmintickstep)
	axes.xaxis.set_minor_locator(minorxLocator)

	if style:
		if "largeticks" in style:
			majorxLocator = MultipleLocator(jdtickstep)
			axes.xaxis.set_major_locator(majorxLocator)
	#minorLocator = MultipleLocator(1) # so to have a tick every day
	#axes.xaxis.set_minor_locator(minorLocator)
	
	
	# Something for astronomers only : we invert the y axis direction !
	axes.set_ylim(axes.get_ylim()[::-1])
	
	
	if colourprop != None and hidecolourbar == False:
		cbar = plt.colorbar(orientation='vertical', shrink=1.0, fraction=0.065, pad=0.025)
		cbar.set_label(colournicename) 

	# And we make custom title :
	
	if title == "None" or title == None or title == "none":
		#plt.title("Lightcurves", fontsize=18)
		pass
	else:
		#plt.title(title, fontsize=18)
		if titlexpos == None:
			axes.annotate(title, xy=(0.5, 1.0), xycoords='axes fraction', xytext=(0, -4),
			textcoords='offset points', ha='center', va='top', fontsize=25)
		else:
			axes.annotate(title, xy=(titlexpos, 1.0), xycoords='axes fraction', xytext=(0, -4),
			textcoords='offset points', ha='center', va='top', fontsize=25)
	if jdrange != None:
		axes.set_xlim(jdrange[0], jdrange[1])
	
	
	axes.set_xlabel("HJD - 2400000.5 [day]", fontsize = labelfontsize)
	axes.set_ylabel("Magnitude (relative)", fontsize = labelfontsize)
	
	if showdelays:
		txt = getnicetimedelays(lclist, separator="\n")
		axes.annotate(txt, xy=(0.77, 0.77), xycoords='axes fraction', xytext=(6, -6),
			textcoords='offset points', ha='left', va='top')
		#plt.text(0.01, 0.99, txt,
		#	horizontalalignment='left', verticalalignment='top',
		#	transform = axes.transAxes)
		legendloc = 1
		if verbose:
			print "Delays between plotted curves :"
			print txt
	
	if showlegend and (len(lclist) > 0 or len(splist) > 0):
		axes.legend(loc = legendloc, numpoints = 1, prop = fm.FontProperties(size = 12))
	
	if magrange != None:
		if type(magrange) == float or type(magrange) == int:
			# We find the mean mag of the stuff to plot :
			allmags = []		
			for l in lclist:
				allmags.extend(l.getmags())
			meanlevel = np.mean(np.array(allmags))
			axes.set_ylim(meanlevel+magrange, meanlevel-magrange)
		else:
			axes.set_ylim(magrange[0], magrange[1])
	
	if showdates: # Be careful when you change something here, it could mess up the axes.
		# Especially watch out when you change the plot range.
		# This showdates stuff should come at the very end
		minjd = axes.get_xlim()[0]
		maxjd = axes.get_xlim()[1]
		#axes.set_xlim(minjd, maxjd)
		yearx = axes.twiny()
		yearxmin = util.datetimefromjd(minjd + 2400000.5)
		yearxmax = util.datetimefromjd(maxjd + 2400000.5)
		yearx.set_xlim(yearxmin, yearxmax)
		yearx.xaxis.set_minor_locator(matplotlib.dates.MonthLocator())
		yearx.xaxis.set_major_locator(matplotlib.dates.YearLocator())
		yearx.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%Y'))
		yearx.xaxis.tick_top()
		if keeponlygrid:
			yearx.set_xticklabels([])
		#yearx.set_xlabel("Date")
	
	minoryLocator = MultipleLocator(magmintickstep)
	axes.yaxis.set_minor_locator(minoryLocator)
	
	
	if showgrid:
		plt.grid(zorder=20)
	
	if text != None:
		for line in text:
			axes.text(line[0], line[1], line[2], transform=axes.transAxes, **line[3])

	if showinsert:
		assert insertname != None
		from matplotlib._png import read_png
		from matplotlib.offsetbox import OffsetImage, AnnotationBbox
		im = read_png(insertname)
		imagebox = OffsetImage(im, zoom=0.5, interpolation="sinc", resample = True)
		ab = AnnotationBbox(imagebox,  xy=(1.0, 1.0), xycoords='axes fraction', xybox = (-75, -75),
				boxcoords="offset points",
				pad=0.0, frameon=False
			)
		axes.add_artist(ab)

	if showlogo:
	
		# The EPFL logo :	
		from matplotlib._png import read_png
		from matplotlib.offsetbox import OffsetImage, AnnotationBbox
		logodir = os.path.dirname(__file__)
		im = read_png(os.path.join(logodir, "epfl.png"))
		imagebox = OffsetImage(im, zoom=0.13, interpolation="sinc", resample = True)
	
		if logopos == "left":
			ab = AnnotationBbox(imagebox,  xy=(0.0, 0.0), xycoords='axes pixels', xybox = (52, 30),
				boxcoords="offset points",
				pad=0.0, frameon=False
			)
			axes.add_artist(ab)
			axes.annotate("COSMOGRAIL.org", xy=(0.0, 0.0), xycoords='axes fraction', fontsize=16, xytext=(105, 7),
				textcoords='offset points', ha='left', va='bottom', color="gray")
				
			## name lightcurves:
			if 0:	
				axes.annotate("A", xy=(0.0,0.0), xycoords='axes fraction', fontsize=25 , xytext=(20,260),
					textcoords='offset points', ha='center', va='bottom', color="red")
				axes.annotate("B", xy=(0.0,0.0), xycoords='axes fraction', fontsize=25 , xytext=(20,150),
					textcoords='offset points', ha='center', va='bottom', color="green")					
				axes.annotate("C", xy=(0.0,0.0), xycoords='axes fraction', fontsize=25 , xytext=(20,125),
					textcoords='offset points', ha='center', va='bottom', color="blue")						
				axes.annotate("D", xy=(0.0,0.0), xycoords='axes fraction', fontsize=25 , xytext=(20,80),
					textcoords='offset points', ha='center', va='bottom', color="purple")
					
							
		if logopos == "right":
			ab = AnnotationBbox(imagebox,  xy=(1.0, 0.0), xycoords='axes fraction', xybox = (-200, 30),
				boxcoords="offset points",
				pad=0.0, frameon=False
			)
			axes.add_artist(ab)
			axes.annotate("COSMOGRAIL.org", xy=(1.0, 0.0), xycoords='axes fraction', fontsize=16, xytext=(-10, 7),
				textcoords='offset points', ha='right', va='bottom', color="gray")
				
		if logopos == "center":
			ab = AnnotationBbox(imagebox,  xy=(0.55, 0.0), xycoords='axes fraction', xybox = (-80, 30),
				boxcoords="offset points",
				pad=0.0, frameon=False
			)
			axes.add_artist(ab)
			axes.annotate("COSMOGRAIL.org", xy=(0.55, 0.0), xycoords='axes fraction', fontsize=16, xytext=(40, 7),
				textcoords='offset points', ha='center', va='bottom', color="gray")				
				
		
		# Alternative possibility (just to keep the idea) :
		"""
		try:
			import Image
		except ImportError:
			print "Couldn't import PIL ! Therefore I won't be able to display the cosmograil logo."
		else:
			im = Image.open('epfl.png')
			height = im.size[1]
			print height
			im = np.array(im).astype(np.float) / 255
			fig.figimage(im, 0, fig.bbox.ymax - height)
			# With newer (1.0) versions of matplotlib, you can 
			# use the "zorder" kwarg to make the image overlay
			# the plot, rather than hide behind it... (e.g. zorder=10)
			fig.figimage(im, 0, fig.bbox.ymax - height)
		"""

	if ax != None:
		return

	if filename == "screen":
		plt.show()
	else:

		plt.savefig(filename, transparent=transparent)
		#if verbose:
		print "Plot written to %s" % filename
		plt.close() # this seems important so that the plot is not displayed when a next plt.show() is called.





def displayrange(lcs, margin=0.05):
	"""
	returns a plausible range of mags and hjds to plot, so that you can keep this fixed in your plots
	"""
	mags = []
	jds = []
	for l in lcs:
		mags.extend(l.getmags())
		jds.extend(l.getjds())
	
	magrange = np.max(mags) - np.min(mags)
	jdrange = np.max(jds) - np.min(jds)
	return ((np.min(jds)-margin*jdrange, np.max(jds)+margin*jdrange), (np.max(mags)+margin*magrange, np.min(mags)-margin*magrange))


# 
# def multigettimedelays(lcs):
# 	"""
# 	Old way ... not used for spline stuff anymore
# 	Returns the timedelays of the n-1 last curves with respect to the first one,
# 	so lc2-lc1, lc3-lc1, ...
# 	We call these numbers time delays of lc1 with respect to lcx. 
# 	So this function "translates" time shifts into time delays.
# 	"""
# 	print "Deprecation warning :-)"
# 	return np.array([l.timeshift - lcs[0].timeshift for l in lcs[1:]]) # this is a copy
# 	
# 	
# def multisettimedelays(lcs, tds):
# 	"""
# 	Old way ... not used for spline stuff anymore
# 	The inverse function of multigettimedelays : you give n-1 time delays, and I set the n-1 last timeshifts of
# 	the lightcurve list so to that they would give these time delays.
# 	The timeshift of the first curve is not changed.
# 	"""
# 	
# 	print "Deprecation warning :-)"
# 	if len(lcs) == 2: 
# 		# i.e. 2 lcs and 1 td
# 		# because 0-d arrays cannot be indexed, we need a special solution :
# 		lcs[1].timeshift = float(tds) + lcs[0].timeshift
# 
# 	else:
# 		for i in range(1, len(lcs)):
# 			#lcs[i].magshift == 0.0 # not sure why I once wrote this here ...
# 			lcs[i].timeshift = tds[i-1] + lcs[0].timeshift
# 	

def getnicetimedelays(lcs, separator="\n", sorted = False):
	"""
	Returns a formatted piece of text of the time delays between a list of lc objects.
	
	:param separator: try ``" | "`` to get all the delays in one line.
	:type separator: string
	
	:param sorted: If True, I will sort my output according to l.object.
		But of course I will **not** modify the order of your lcs !
		By default this is False : if the curves are B, A, C and D in this order,
		I would return BA, BC, BD, AC, AD, CD
	
	:type sorted: boolean
	
	.. warning:: This is the function that **defines** the concept of "delay" (and its sign) used by
		pycs. **The delay AB corresponds to timeshift(B) - timeshift(A)**
		Other places where this is hard-coded : 
		:py:func:`pycs.sim.plot.hists`
	
	"""
	n = len(lcs)
	if sorted == True:
		worklcs = objsort(lcs, ret=True, verbose = False)
	else:
		worklcs = lcs
	couples = [(worklcs[i], worklcs[j]) for i in range(n) for j in range(n) if i < j]
	return separator.join(["%s%s = %+7.2f" % ( lc1.object, lc2.object, lc2.timeshift - lc1.timeshift) for (lc1, lc2) in couples])

def getnicetimeshifts(lcs, separator="\n"):
	"""
	I return the timeshifts as a text string, use mainly for checks and debugging.
	"""
	return separator.join(["%s  %+7.2f" % (l.object, l.timeshift) for l in lcs])


def gettimeshifts(lcs, includefirst=True):
	"""
	I simply return the **absolute** timeshifts of the input curves as a numpy array.
	This is used by time shift optimizers, and usually not by the user.
	
	:param includefirst: If False, I skip the first of the shifts.
		For timeshift optimizers, it's indeed often ok not to move the first curve.
	:type includefirst: boolean
	
	"""
	if includefirst == True:
		return np.array([l.timeshift for l in lcs])
	if includefirst == False:
		return np.array([l.timeshift for l in lcs[1:]])

def settimeshifts(lcs, shifts, includefirst=False):
	"""
	I set the timeshifts like those returned from :py:func:`gettimeshifts`.
	An do a little check.
	"""
	
	if includefirst == True:
		assert len(lcs) == len(shifts)
		for (l, shift) in zip(lcs, shifts):
				l.timeshift = shift
	
	elif includefirst == False:
		#print lcs, shifts
		assert len(lcs)-1 == len(shifts)
		for (l, shift) in zip(lcs[1:], shifts):
				l.timeshift = shift

	
def shuffle(lcs):
	"""
	I scramble the order of the objects in your list, in place.
	"""
	random.shuffle(lcs) # That was easier than I thought it would be when starting to wright this function ...
	
	
def objsort(lcs, ret=False, verbose=True):
	"""
	I sort the lightcurve objects in your list (in place) according to their object names.
	
	:param ret: If True, I do not sort them in place, but return a new sorted list (containing the same objects).	
	:type ret: boolean
	
	"""

	# Maybe we start with some checks :
	diff_objects = set([l.object for l in lcs])
	if len(diff_objects) != len(lcs):
		raise RuntimeError("Cannot sort these objects : %s" % ", ".join([l.object for l in lcs]))
	
	# The actual sorting ...
	if ret == False:
		lcs.sort(key = operator.attrgetter('object'))
		if verbose:
			print "Sorted lcs, order : %s" % ", ".join([l.object for l in lcs])
		return
	if ret == True:
		return sorted(lcs, key = operator.attrgetter('object'))
	
	


	
	
	

