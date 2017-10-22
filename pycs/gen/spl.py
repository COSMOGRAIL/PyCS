"""
Module defining the Spline class, something easy to wrap around SciPy splines.
Includes BOK algorithms (Mollinari et al)

Some rules of splrep (k = 3)
	- do not put more then 2 knots between data points.
	- splrep wants inner knots only, do not give extremal knots, even only "once".

"""

import numpy as np
import sys
import pycs.gen.util
import copy as pythoncopy
import matplotlib.pyplot as plt
import scipy.optimize as spopt
import scipy.interpolate as si


class DataPoints():
	"""
	An ultralight version of a lightcurve, made for fast computations.
	Can be "merged" from a list of lightcurves, see factory function below.
	
	A Spline object has such a DataPoints object as attribute.
	
	ATTENTION
	Datapoints are expected to be ALWAYS SORTED BY JDS, and no two datapoints have the same jd !
	See the splitup option of the constructor.
	Note that this is not the case for lightcurves ! Hence the existence of datapoints.
	Should be enforced in every function that builds datapoints.
	
	ABOUT STAB POINTS
	With scipy splines, we always get the last knots at the extrema of data points.
	So to get knots "outside" of the real datapoints, we have to insert fake points.
	And while we are at it, these fake points can also be used to stabilize the spline in 
	gaps.
	The mask is used to differentiate between actual data points and "stabilization points"
	that are inserted to make the spline behave well at the extrema and in season gaps.
	It is modified by the two addgappts and addextpts.
	
	The info about stabpoints is written into the object,
	so that they can be reconstrucuted from any new jds and mags.
	"""

	def __init__(self, jds, mags, magerrs, splitup=True, deltat=0.000001, sort=True, stab=False, 
		stabext=300.0, stabgap = 30.0, stabstep = 5.0, stabmagerr = -2.0, stabrampsize = 0, stabrampfact = 1.0):
		"""
		Constructor
		Always leave splitup and sort on True ! Only if you know that you are already
		sorted you can skip them.
		You cannot specify a mask, I do this myself. (could be done in principle).
		
		stab : do you want stabilization points ?
		Don't forget to run splitup, sort, and addstab again if you change the data !
		"""
		
		self.jds = jds
		self.mags = mags
		self.magerrs = magerrs
		
		self.stab = stab
		self.stabext = stabext
		self.stabgap = stabgap
		self.stabstep = stabstep
		self.stabmagerr = stabmagerr
		
		self.stabrampsize = stabrampsize
		self.stabrampfact = stabrampfact
			
		self.mask = np.ones(len(self.jds), dtype=np.bool) # an array of True
		
		self.deltat = deltat
		if splitup:
			self.splitup()	
		elif sort: # If we do the splitup, we sort anyway.
			self.sort()
			
		self.putstab()
			
	
# 	def update(self, jds, mags, magerrs):
# 		"""
#		NOT NEEDED ANYMORE, JUST CALL MERGE AND GIVE AN OLDDP. SAFER.
#
# 		Give me some new datapoints (no stabs) (already splitup and sorted, by definition), I'll update myself.
# 		In fact everything might move !
# 		"""
# 		if newdatapoints.stab = True:
# 			raise RuntimeError("Give me points without stab !")
# 		self.jds = newdatapoints.jds
# 		self.mags = newdatapoints.mags
# 		self.magerrs = newdatapoints.magerrs
# 		self.mask = np.ones(len(self.jds), dtype=np.bool)
# 		self.addstab() # runs only if stab = True
	
	def splitup(self):
		"""
		TO WRITE !!!
		We avoid that two points get the same jds...
		Note that this might change the order of the jds,
		but only of very close ones, so one day it would be ok to leave the mags as they are.
		
		"""
		self.jds += self.deltat * np.random.randn(len(self.jds))
		self.sort()
	
	def sort(self):
		"""
		Absolutely mandatory, called in the constructor.
		"""
		sortedindices = np.argsort(self.jds)
		self.jds = self.jds[sortedindices]
		self.mags = self.mags[sortedindices]
		self.magerrs = self.magerrs[sortedindices]
		self.mask = self.mask[sortedindices]
		self.validate()
	
	def validate(self):
		"""
		We check that the datapoint jds are increasing strictly :
		"""
		first = self.jds[:-1]
            	second = self.jds[1:]
		if not np.alltrue(np.less(first,second)): # Not less_equal ! Strictly increasing !
			raise RuntimeError, "These datapoints don't have strcitly increasing jds !"
		
	def rmstab(self):
		"""
		Deletes all stabilization points
		"""
		self.jds = self.jds[self.mask]
		self.mags = self.mags[self.mask]
		self.magerrs = self.magerrs[self.mask]
		self.mask = np.ones(len(self.jds), dtype=np.bool)
	
	def putstab(self):
		"""
		Runs only if stab is True.
		I will :
		add datapoints (new jds, new mags, new magerrs)
		modify the mask = False for all those new datapoints.
		"""
		if self.stab == True:
		
			# We start by deleting any previous stab stuff :
			
			self.rmstab()
			self.addgappts()
			self.addextpts()
		else:
			pass
	
	
	def calcstabmagerr(self):
		"""
		Computes the mag err of the stabilisation points.
		"""
		if self.stabmagerr >= 0.0:
			return self.stabmagerr
		else:
			return - self.stabmagerr * np.median(self.magerrs)
			
		
	
	def addgappts(self):
		"""
		We add stabilization points with low weights into the season gaps
		to avoid those big excursions of the splines.
		This is done by a linear interpolation across the gaps.
		"""
		
		absstabmagerr = self.calcstabmagerr()
		
		gaps = self.jds[1:] - self.jds[:-1] # has a length of len(self.jds) - 1
		gapindices = np.arange(len(self.jds) - 1)[gaps > self.stabgap] # indices of those gaps that are larger than stabgap
		
		for n in range(len(gapindices)):
			i = gapindices[n]
			a = self.jds[i]
			b = self.jds[i+1]
			
			newgapjds = np.linspace(a, b, float(b-a)/float(self.stabstep))[1:-1]
			newgapindices = i + 1 + np.zeros(len(newgapjds))
			newgapmags = np.interp(newgapjds, [a, b], [self.mags[i], self.mags[i+1]])
			newgapmagerrs = absstabmagerr * np.ones(newgapmags.shape)
			newgapmask = np.zeros(len(newgapjds), dtype=np.bool)
			
			self.jds = np.insert(self.jds, newgapindices, newgapjds)
			self.mags = np.insert(self.mags, newgapindices, newgapmags)
			self.magerrs = np.insert(self.magerrs, newgapindices, newgapmagerrs)
			self.mask = np.insert(self.mask, newgapindices, newgapmask)
			
			gapindices += newgapjds.size # yes, as we inserted some points the indices change.
	
		# If you change this structure, be sure to check SplineML.settargetmags as well !
	
		self.validate()
	
	def addextpts(self):
		"""
		We add stabilization points at both extrema of the lightcurves
		This is done by "repeating" the extremal points, and a ramp in the magerrs
		"""
		
		absstabmagerr = self.calcstabmagerr()
		
		extjds = np.arange(self.jds[0], self.jds[0] - self.stabext, -1*self.stabstep)[::-1][:-1]
		extmags = self.mags[0] * np.ones(extjds.shape)
		extmagerrs = absstabmagerr * np.ones(extjds.shape)
		for i in range(1, self.stabrampsize+1):
			extmagerrs[-i] += (self.stabrampsize +1 -i) * absstabmagerr * self.stabrampfact
		extindices = np.zeros(extjds.shape)
		mask = np.zeros(len(extjds), dtype=np.bool)
		self.jds = np.insert(self.jds, extindices, extjds)
		self.mags = np.insert(self.mags, extindices, extmags)
		self.magerrs = np.insert(self.magerrs, extindices, extmagerrs)
		self.mask = np.insert(self.mask, extindices, mask)
		
		# And the same at the other end :
		
		extjds = np.arange(self.jds[-1], self.jds[-1] + self.stabext, self.stabstep)[1:]
		extmags = self.mags[-1] * np.ones(extjds.shape)
		extmagerrs = absstabmagerr * np.ones(extjds.shape)
		for i in range(0, self.stabrampsize):
			extmagerrs[i] += (self.stabrampsize -i) * absstabmagerr * self.stabrampfact
		extindices = len(self.jds) + np.zeros(extjds.shape)
		mask = np.zeros(len(extjds), dtype=np.bool)
		self.jds = np.insert(self.jds, extindices, extjds)
		self.mags = np.insert(self.mags, extindices, extmags)
		self.magerrs = np.insert(self.magerrs, extindices, extmagerrs)
		self.mask = np.insert(self.mask, extindices, mask)
	
		self.validate()
		
	def getmaskbounds(self):
		"""
		Returns the upper and lower bounds of the regions containing stabilization points.
		This is used when placing knots, so to put fewer knots in these regions.
		Crazy stuff...
		"""
	
		maskindices = np.where(self.mask == False)[0]
		#print maskindices
		
		if len(maskindices) < 3:
			print "Hmm, not much masked here ..."
			return (np.array([]), np.array([]))
		else:
			lcuts = maskindices[np.where(maskindices[1:] - maskindices[:-1] > 1)[0] + 1]
			lcuts = np.insert(lcuts, 0, maskindices[0])
			ucuts = maskindices[np.where(maskindices[1:] - maskindices[:-1] > 1)[0]]
			ucuts = np.insert(ucuts, len(ucuts), maskindices[-1])
			return (lcuts, ucuts)

	def ntrue(self):
		"""
		Returns the number of real datapoints (skipping stabilization points)
		"""
		return np.sum(self.mask)



def merge(lcs, olddp=None, splitup=True, deltat=0.000001, sort=True, stab=False,
	stabext=300.0, stabgap = 30.0, stabstep = 5.0, stabmagerr = 2.0, stabrampsize = 0, stabrampfact = 1.0):
	"""
	Factory function for DataPoints objects, starting from lightcurves.
	Takes a list of lightcurves and quickly concatenate the jds, mags, and magerrs.
	
	Instead of specifying all the stab point parameters, you can give me an old datapoints object,
	and I will reuse its settings... This is useful if you want to "update" the data points.
	
	If overlap is True, I will keep only points that are "covered" by all four lightcurves !
	This is useful when you want to build a first source spline, and your microlensing is messy at the borders.
	NOT YET IMPLEMENTED ...
	"""

	jds = np.concatenate([l.getjds() for l in lcs])
	mags = np.concatenate([l.getmags() for l in lcs])
	magerrs = np.concatenate([l.getmagerrs() for l in lcs])

	if olddp == None:
		return DataPoints(jds, mags, magerrs, splitup=splitup, deltat=deltat, sort=sort,
			stab=stab, stabext=stabext, stabgap=stabgap, stabstep=stabstep, stabmagerr=stabmagerr,
			stabrampsize=stabrampsize, stabrampfact=stabrampfact)
	else:
		return DataPoints(jds, mags, magerrs, splitup=splitup, sort=sort,
			deltat=olddp.deltat,
			stab=olddp.stab, stabext=olddp.stabext, stabgap=olddp.stabgap, stabstep=olddp.stabstep, stabmagerr=olddp.stabmagerr,
			stabrampsize=olddp.stabrampsize, stabrampfact=olddp.stabrampfact)


class Spline():
	"""
	A class to represent a spline, that is essentially a set of knots and coefficients.
	As finding knots and coefficients requires access to some data points, these are included
	in the form of a DataPoints object.
	
	Abount knots :
	Spline.t are all the knots, including extremas with multiplicity.
	But splrep wants internal knots only ! By internal we mean : not even the data extremas !
	Spline.getintt() returns only these internal knots.
	"""
	
	def __init__(self, datapoints, t = None, c = None, k = 3, bokeps = 2.0, boktests = 5, bokwindow = None, plotcolour="black"):
		"""
		t : all the knots (not only internal ones !)
		c : corresponding coeffs
		k : degree : default = cubic splines k=3 -> "order = 4" ???
		whatever ... 3 means that you can differentiate twice at the knots.
		
		"""
		
		#self.origdatapoints = datapoints
		self.datapoints = datapoints
		
		# At this point we know that your datapoint jds are monotonously increasing. This is tested
		# by validate() of datapoints.
		
		self.t = t # the array of knots
		self.c = c # the coeffs	
		self.k = k
		
		self.bokeps = bokeps
		self.boktests = boktests
		self.bokwindow = bokwindow
		
		self.knottype = "none"
		self.plotcolour = plotcolour
		self.showknots = True
		
		# Bounds, for BOK
		self.lims = None
		self.l = None
		self.u = None
		
		# We want to keep trace of the r2 of a spline.
		self.lastr2nostab = 0.0 # without stab points (the real thing)
		self.lastr2stab = 0.0 # with stab points (usually not so interesting)
		
		# If you did not give me a t&c, I'll make some default ones for you :
		try:
			if (self.t is None):
				self.uniknots(2) # This also puts self.c to 0s
		except:
			if (len(self.t) == 0):
				self.uniknots(2) # This also puts self.c to 0s
		
		
	def __str__(self):
		"""
		Returns a string with:
		* degree
		* knot placement
		* number of intervals
		
		"""
		#return "Spline of degree %i, %i knots (%i inner knots), and %i intervals." % (self.k, len(self.t), len(self.getintt()), self.getnint())
		
		if len(self.knottype) > 6: # That's a string
			knottext = "%il%ib" % (self.knottype.count("l"), self.knottype.count("b"))
		else:
			knottext = self.knottype
		return "~%i/%s/%i~" % (self.k, knottext, self.getnint())
	
	
	def copy(self):
		"""
		Returns a "deep copy" of the spline.
		""" 
		return pythoncopy.deepcopy(self)
	
	
	def shifttime(self, timeshift):
		"""
		Hard-shifts your spline along the time axis.
		By "hard-shift", I mean that unlike for a lightcurve, the spline will not know that it was shifted !
		It's up to you to be sure that you want to move it.
		We shift both the datapoints and the knots.
		"""
		
		self.t += timeshift
		self.datapoints.jds += timeshift
		
		

	def updatedp(self, newdatapoints, dpmethod="stretch"):
		"""
		
		Replaces the datapoints of the spline, and makes sure that the knots
		stay compatible.
		
		If you tweaked your datapoints, I will have to tweak my knots to make sure
		that my external knots fit. Hence this method !
		Due to the splitup, this is needed even if you just tweaked the mags !
		And anyway in this case I have to rebuild the stab points.
		
		.. warning :: IT'S UP TO YOU TO CHECK THAT YOU DON'T REPLACE DATATOINTS WITH DIFFERENT STAB SETTINGS
			Anyway it would work, just look ugly !

		Replaces the datapoints (jds, mags, and magerrs) touching the knots and coeffs as less as possible.
		Note that we also have to deal with stab points here !
		
		This is made for instance for time shifts that only very slightly change the datapoints, and you don't want to
		optimize the knots all the time from scratch again.
		The current knots are "streched" (keeping their relative spacings) accross the new datapoints.
		
		Options for "dpmethod" :
		- "stretch" : changes all the knots
		- "extadj" : does not touch the internal knots, but adjusts the external ones only, to
		fit the new datapoints. Probably the method to use when optimizing time shifts.
		- "leave" : does not touch the knots -> ok to evaluate the spline,
		but you will not be able to fit it anymore, as the external knots don't correspond to datapoints.	
		
		
		.. todo:: In principle, why don't we just update the real datapoints here, and leave the stab as
			they are ?
		
		"""
		
		if dpmethod == "stretch":
		
			oldmin = self.datapoints.jds[0] # This includes potential stab points
			oldmax = self.datapoints.jds[-1]
			
			newmin = newdatapoints.jds[0] # Idem
			newmax = newdatapoints.jds[-1]
			
			oldknots = self.getinttex()
			#print oldknots
		
			# we will stretch the oldknots by a factor a :
			a = (newmax - newmin)/(oldmax - oldmin)
					
			newknots = newmin + a*(oldknots-oldmin)
			
			# We set the new datapoints:
			self.datapoints = newdatapoints
			self.setinttex(newknots)
		
		elif dpmethod == "extadj" :
			
			intknots = self.getintt()
			self.datapoints = newdatapoints
			
			# Ok, now the newdatapoints might be narrower or wider than the knots, we have to deal with this.
			# If they are wider, it's easy : setint will put move the external knot on the external datapoint.
			# If they are narrower, it's trickier : we have to remove some extra knots, so to really just keep the "internal" ones.
			# to feed into setintt.
			
			#if True: # works as well, but maybe faster to test first :
			if (self.datapoints.jds[0] >= intknots[0]) or (self.datapoints.jds[-1] <= intknots[-1]):
			
				keepmask = np.ones(intknots.shape, dtype=np.bool)
				for i in range(len(intknots)): # Starting from the left ...
					if intknots[i] <= self.datapoints.jds[0]:
						keepmask[i] = False
					else:
						break
				for i in range(len(intknots))[::-1]: # And now the right ...
					if intknots[i] >= self.datapoints.jds[-1]:
						keepmask[i] = False
					else:
						break
				
				#nkick = np.sum(keepmask == False)
				#if nkick != 0:
				#	print "I'll kick %i knots !" % (nkick)
				
				# And finally, we apply the mask .
				intknots = intknots[keepmask]
			
			self.setintt(intknots) # This automatically adjusts the extremal knots.
			

		elif dpmethod == "leave" :
			
			knots = self.getinttex()
			self.datapoints = newdatapoints
			
			# We quickly check the boundaries
			if ( knots[0] >= self.datapoints.jds[0] ) or ( knots[-1] <= self.datapoints.jds[-1] ):
				raise RuntimeError("Your newdatapoints are to wide for the current knots !")
			
		else:
			raise RuntimeError("Don't know this updatedp method !")

			
		# We reset any bounds just to be sure.
		self.lims = None
		self.l = None
		self.u = None


	def uniknots(self, nint, n=True):
		"""
		Uniform distribution of internal knots across the datapoints (including any stab points).
		We don't make a difference between stab and real points.
		
		:param nint: The number of intervals, or the step
		:param n:
			If True, nint is the number of intervals (== piecewise polynoms) you want.
			If False : nint is a step in days you want between the knots (approximately).
		:type n: boolean
		
		.. note:: I also put all coeffs back to 0.0 !
		
		"""

		#intt = np.linspace(self.datapoints.jds[0], self.datapoints.jds[-1], step+1)[1:-1] # we remove the extremas	
		a = self.datapoints.jds[0]
		b = self.datapoints.jds[-1]
		if n:
			intt = np.linspace(a, b, nint + 1)[1:-1]
		else:
			intt = np.linspace(a, b, float(b-a)/float(nint))[1:-1]
			
		if len(intt) == 0:
			raise RuntimeError("I am uniknots, and I have only 0 (zero) internal knots ! Increase this number !")
			
		self.setintt(intt)
		self.knottype = "u"
		
		# Important : we put some 0 coeffs to go with the new knots
		self.resetc()


	def resetc(self):
		"""
		Sets all coeffs to 0.0 -- if you want to start again your fit, keeping the knot positions.
		"""
		self.c = np.zeros(len(self.t))
		
		
	def reset(self):
		"""
		Calls uniknots, i.e. resets both coeffs and knot positions, keeping the same number of knots.
		"""
		self.uniknots(self.getnint() ,n=True)
		

	
	def buildbounds(self, verbose = True):
		"""
		Build bounds for bok.
		By default I will make those bounds as wide as possible, still respecting epsilon.
		The parameter epsilon is the minimum distance two knots can have.
		If you give me a window size, I will not make the bounds as wide as possible, but only put them
		0.5*window days around the current knots (still respecting all this epsilon stuff of course).
		
		I look where your current knots are, and for each knots I build the bounds so that
		epsilon distance is respected between adjacent upper and lower bounds.
		But, there might already be knots only epsilon apart.
		So I'm a bit tricky, not so straightforward as my predecessors.
		
		Knots at the extrema are not allowed to move.
		
		Requires existing knots, puts lims in between them, and builds the bounds.
		@todo: Optimize me using numpy ! This is experiemental code for now.
		"""
		
		if verbose:
			print "Building BOK bounds (bokeps = %.3f, bokwindow = %s) ..." % (self.bokeps, self.bokwindow)
		
		knots = self.getinttex() # Including extremal knots (once).
		n = len(knots)
		
		# We start by checking the knot spacing
		knotspacings = knots[1:] - knots[:-1]
		if not np.alltrue(knotspacings > 0.0):
			raise RuntimeError("Ouch, your knots are not sorted !")
		minspace = np.min(knotspacings)
		if verbose :
			print "Minimal knot spacing : %.3f" % (minspace)
		if minspace < self.bokeps - 0.00001: # Rounding errors, we decrease epsilon a bit...
			# If this does still happens, then it was not just a rounding error ...
			# Yes it still happens, due to updatedp stretch ...
			raise RuntimeError("Knot spacing min = %f, epsilon = %f" % (minspace, self.bokeps))
		
		
		# Loop through the knots. 
		lowers = [knots[0]] # First knot is not allowed to move
		uppers = [knots[0]]
		for i in range(1, n-1): # Internal knots
			tk = knots[i] # this knot
			pk = knots[i-1] # previous knot
			nk = knots[i+1] # next knot
			
			# First we build the wide bounds :
			guessl = 0.5*(pk + tk) + 0.5*self.bokeps
			if guessl >= tk:
				guessl = tk
			
			guessu = 0.5*(nk + tk) - 0.5*self.bokeps
			if guessu <= tk:
				guessu = tk
			
			# Now we see if the use wants a narrower window within those bounds :
			if self.bokwindow != None:	
				if tk - 0.5*self.bokwindow >= guessl:
					guessl = tk - 0.5*self.bokwindow
				if tk + 0.5*self.bokwindow <= guessu:
					guessu = tk + 0.5*self.bokwindow
	
			lowers.append(guessl)
			uppers.append(guessu)
		
		# And now this last knot, doesn't move, like the first one:
		lowers.append(knots[-1])
		uppers.append(knots[-1])
		self.l = np.array(lowers)
		self.u = np.array(uppers)
		self.knottype += "l"
		if verbose:
			print "Buildbounds done."
		
	
	def bok(self, bokmethod="BF", verbose=True, trace=False):
		"""
		We optimize the positions of knots by some various techniques.
		We use fixed bounds for the exploration, run buildbounds (with low epsilon) first.
		This means that I will not move my bounds.
		
		For each knot, i will try ntestpos linearly spaced positions within its bounds.
		In this version, the bounds are included : I might put a knot on a bound !
		The way the bounds are placed by buildbounds ensures that in any case the minimal
		distance of epsilon is respected.
		
		Using this sheme, it is now possible to iteratively call mybok and buildbounds in a loop
		and still respect epsilon at any time.
		
		
		bokmethods :
			- MCBF : Monte Carlo brute force with ntestpos trial positions for each knot
			- BF : brute force, deterministic. Call me twice
			- fminind : fminbound on one knot after the other.
			- fmin :global fminbound
		
		Exit is automatic, if result does not improve anymore...
		"""
		intknots = self.getintt() # only internal, the ones we will move
		nintknots = len(intknots)
		weights = 1.0/self.datapoints.magerrs
		
		def score(intknots, index, value):
			modifknots = intknots.copy()
			modifknots[index] = value
			return si.splrep(self.datapoints.jds, self.datapoints.mags, w=weights, xb=None, xe=None, k=self.k, task=-1, s=None, t=modifknots, full_output=1, per=0, quiet=1)[1]
		
		iniscore = score(intknots, 0, intknots[0])
		lastchange = 1
		lastscore = iniscore
		iterations = 0
		if verbose:
			print "Starting BOK-%s on %i intknots (boktests = %i)" % (bokmethod, nintknots, self.boktests)
		
		if bokmethod == "MCBF":
		
			while True:
				if lastchange >= 2*nintknots: # somewhat arbitrary, but why not.
					break
				i = np.random.randint(0, nintknots) # (inclusive, exclusive)
			
				testknots = np.linspace(self.l[i+1], self.u[i+1], self.boktests)
				# +1, as u and l include extremal knots...
				# So we include the extremas in our range to test.
			
				testscores = np.array([score(intknots, i, testknot) for testknot in testknots])
				bestone = np.argmin(testscores)
			
				bestscore = testscores[bestone]
				if bestscore < lastscore:
					lastchange = 0
				intknots[i] = testknots[bestone] # WE UPDATE the intknots array !
				lastscore = bestscore
				lastchange += 1
				iterations += 1
				
				if trace:
					self.optc()
					pycs.gen.util.trace([], [self])
		
		if bokmethod == "BF":
			
			intknotindices = range(nintknots) # We could potentially change the order, just to see if that makes sense.
			# No, it doesn't really help
			#mid = int(len(intknotindices)/2.0)
			#intknotindices = np.concatenate([intknotindices[mid:], intknotindices[:mid][::-1]])
			
			for i in intknotindices:
				testknots = np.linspace(self.l[i+1], self.u[i+1], self.boktests)
				# +1, as u and l include extremal knots...
				# So we include the extremas in our range to test.
			
				testscores = np.array([score(intknots, i, testknot) for testknot in testknots])
				bestone = np.argmin(testscores)
			
				bestscore = testscores[bestone]
				intknots[i] = testknots[bestone] # WE UPDATE the intknots array !
				iterations += 1
				
				if trace:
					self.optc()
					pycs.gen.util.trace([], [self])
		
		if bokmethod == "fminind":
			intknotindices = range(nintknots)
			for i in intknotindices:
				
				def target(value):
					return score(intknots, i, value)
				#inival = intknots[i]
				#bounds = (self.l[i+1], self.u[i+1])
				
				out = spopt.fminbound(target, self.l[i+1], self.u[i+1], xtol=0.01, maxfun=100, full_output=1, disp=1)
				#print out
				optval = out[0]
				bestscore = out[1]
				
				intknots[i] = optval # WE UPDATE the intknots array !
				iterations += 1
				
				if trace:
					self.optc()
					pycs.gen.util.trace([], [self])
					
		if bokmethod == "fmin":
				
			def target(modifknots):
				
				#iterations += 1
				#if trace:
				#	self.optc()
				#	pycs.gen.util.trace([], [self])
				
				return si.splrep(self.datapoints.jds, self.datapoints.mags, w=weights, xb=None, xe=None, k=self.k, task=-1, s=None, t=modifknots, full_output=1, per=0, quiet=1)[1]
		
			bounds = [(a, b) for (a, b) in zip(self.l[1:-1], self.u[1:-1])]
			
			out  = spopt.fmin_l_bfgs_b(target, intknots, approx_grad=True, bounds=bounds, m=10, factr=1e7, pgtol=1.e-05, epsilon=1e-04, iprint=-1, maxfun=15000)
			#out = spopt.fminbound(target, self.l[1:-1], self.u[1:-1], xtol=0.01, maxfun=1000, full_output=1, disp=3)
			#print out
			intknots = out[0]
			bestscore = out[1]
		
		
		
		# relative improvement :
		relimp = (iniscore - bestscore)/iniscore
		self.knottype += "b"
		self.setintt(intknots)
		
		#pycs.gen.lc.display([],[self])
		#self.display()
		self.optc() # Yes, not yet done !
		finalr2 = self.r2(nostab=True)
		if verbose:
			print "r2 = %f (without stab poins)" % finalr2
			print "Done in %i iterations, relative improvement = %f" % (iterations, relimp)
			# We count all datapoints here, as score returns the full chi2 including stab pts.
		
		return finalr2



		
	# Some stuff about knots :
	
	def getintt(self):
		"""
		Returns the internal knots (i.e., not even the datapoints extrema)
		This is what you need to feed into splrep !
		There are nint - 1 such knots
		"""
		return self.t[(self.k+1):-(self.k+1)].copy() # We cut the outer knots.
	
	def getinttex(self):
		"""
		Same as above, but we include the extremal points "once".
		"""
		return self.t[(self.k):-(self.k)].copy()
	
	def knotstats(self):
		"""
		Returns a string describing the knot spacing
		"""
		knots = self.getinttex()
		spacings = knots[1:] - knots[:-1]
		return " ".join(["%.1f" % (spacing) for spacing in sorted(spacings)])

	
	def setintt(self, intt):
		"""
		Give me some internal knots (not even containing the datapoints extrema),
		and I build the correct total knot vector t for you.
		I add the extremas, with appropriate multiplicity.
		
		@TODO: check consistency of intt with datapoints !
		"""
		
		# Ok a quick test for consisency :
		
		if len(intt) == 0:
			raise RuntimeError("Your list of internal knots is empty !")
		
		if not self.datapoints.jds[0] < intt[0]:
			raise RuntimeError("Ouch.")
		if not self.datapoints.jds[-1] > intt[-1]:
			raise RuntimeError("Ouch.")
		#assert self.datapoints.jds[0] < intt[0] # should we put <= here ?
		#assert self.datapoints.jds[-1] > intt[-1]
		
		pro = self.datapoints.jds[0] * np.ones(self.k+1)
		post = self.datapoints.jds[-1] * np.ones(self.k+1)
		
		self.t = np.concatenate((pro, intt, post))
		
	
	
	def setinttex(self, inttex):
		"""
		Including extremal knots
		"""
		#pro = self.datapoints.jds[0] * np.ones(self.k)
		#post = self.datapoints.jds[-1] * np.ones(self.k)
		pro = inttex[0] * np.ones(self.k)
		post = inttex[-1] * np.ones(self.k)
		
		self.t = np.concatenate((pro, inttex, post))
	
	
	def getnint(self):
		"""
		Returns the number of intervals
		"""
		return(len(self.t) - 2* (self.k + 1) + 1)
	
	# Similar stuff about coeffs :
		
	def getc(self, m=0):
		"""
		Returns all active coefficients of the spline, the ones it makes sense to play with.
		The length of this guy is number of intervals - 2 !
		"""
		return self.c[m:-(self.k + 1 + m)].copy()
		
		
	def setc(self, c, m=0):
		"""
		Puts the coeffs from getc back into place.
		"""
		self.c[m:-(self.k + 1 + m)] = c
	
	def getco(self, m=0):
		"""
		Same as getc, but reorders the coeffs in a way more suited for nonlinear optimization
		"""
		c = self.getc(m=m)
		mid = int(len(c)/2.0)
		return np.concatenate([c[mid:], c[:mid][::-1]])
	
	def setco(self, c, m=0):
		"""
		The inverse of getco.
		"""
		mid = int(len(c)/2.0)
		self.setc(np.concatenate([c[mid+1:][::-1], c[:mid+1]]), m=m)
		
	
	def setcflat(self, c):
		"""
		Give me coeffs like those from getc(m=1), I will set the coeffs so that the spline extremas
		are flat (i.e. slope = 0).
		"""
		
		self.setc(c, m=1)
		self.c[0] = self.c[1]
		self.c[-(self.k + 2)] = self.c[-(self.k + 3)]
	
		
	def setcoflat(self, c):
		"""
		idem, but for reordered coeffs.
		"""
		mid = int(len(c)/2.0)
		self.setcflat(np.concatenate([c[mid:][::-1], c[:mid]]))
		
	
	
	
	def r2(self, nostab=True, nosquare=False):
		"""
		Evaluates the spline, compares it with the data points and returns a weighted sum of residuals r2.
		
		If nostab = False, stab points are included 
		This is precisely the same r2 as is used by splrep for the fit, and thus the same value as 
		returned by optc !
		
		This method can set lastr2nostab, so be sure to end any optimization with it.
		
		If nostab = True, we don't count the stab points
		"""
		
		if nostab == True :
			splinemags = self.eval(nostab = True, jds = None)
			errs = self.datapoints.mags[self.datapoints.mask] - splinemags
			werrs = errs/self.datapoints.magerrs[self.datapoints.mask]
			if nosquare:
				r2 = np.sum(np.fabs(werrs))
			else:
				r2 = np.sum(werrs * werrs)
			self.lastr2nostab = r2
		else :
			splinemags = self.eval(nostab = False, jds = None)
			errs = self.datapoints.mags - splinemags
			werrs = errs/self.datapoints.magerrs
			if nosquare:
				r2 = np.sum(np.fabs(werrs))
			else:
				r2 = np.sum(werrs * werrs)
			self.lastr2stab = r2
			
		return r2
		
		#if red:
		#	return chi2/len(self.datapoints.jds)
	

	def tv(self):
		"""
		Returns the total variation of the spline. Simple !
		http://en.wikipedia.org/wiki/Total_variation
		
		"""
		
		# Method 1 : linear approximation
		
		ptd = 5 # point density in days ... this is enough !
		
		a = self.t[0]
		b = self.t[-1]
		x = np.linspace(a, b, int((b-a) * ptd))
		y = self.eval(jds = x)
		tv1 = np.sum(np.fabs(y[1:] - y[:-1]))
		#print "TV1 : %f" % (tv1)
		
		return tv1	

		# Method 2 : integrating the absolute value of the derivative ... hmm, splint does not integrate derivatives ..		
		#si.splev(jds, (self.t, self.c, self.k))
		
	
	def optc(self):
		"""
		Optimize the coeffs, don't touch the knots
		This is the fast guy, one reason to use splines :-)
		Returns the chi2 in case you want it (including stabilization points) !
		
		Sets lastr2stab, but not lastr2nostab !
		
		"""

		out = si.splrep(self.datapoints.jds, self.datapoints.mags, w=1.0/self.datapoints.magerrs, xb=None, xe=None, k=self.k, task=-1, s=None, t=self.getintt(), full_output=1, per=0, quiet=1)
		# We check if it worked :
		if not out[2] <= 0: 
			raise RuntimeError("Problem with spline representation, message = %s" % (out[3]))
		
		self.c = out[0][1] # save the coeffs
			
		#import matplotlib.pyplot as plt
		#plt.plot(self.datapoints.jds, self.datapoints.magerrs)
		#plt.show()
		
		self.lastr2stab = out[1]
		return out[1]
	
	def optcflat(self, verbose = False):
		"""
		Optimizes only the "border coeffs" so to get zero slope at the extrema
		Run optc() first ...
		This has to be done with an iterative optimizer
		"""
		
		full = self.getc(m=1)
		inip = self.getc(m=1)[[0, 1, -2, -1]] # 4 coeffs
		
		def setp(p):
			full[[0, 1, -2, -1]] = p
			self.setcflat(full)
		
		if verbose:
			print "Starting flat coeff optimization ..."
			print "Initial pars : ", inip

		def errorfct(p):
			setp(p)
			return self.r2(nostab=False) # To get the same as optc would return !

		minout = spopt.fmin_powell(errorfct, inip, full_output=1, disp=verbose)
		popt = minout[0]
		if popt.shape == ():
			popt = np.array([popt])
		
		if verbose:
			print "Optimal pars : ", popt
		setp(popt)
		return self.r2(nostab=False) # We include the stab points, like optc does.
		# This last line also updates self.lastr2 ...
	

	
	def eval(self, jds = None, nostab = True):
		"""
		Evaluates the spline at jds, and returns the corresponding mags-like vector.
		By default, we exclude the stabilization points !
		If jds is not None, we use them instead of our own jds (in this case excludestab makes no sense)
		"""
		if jds is None:
			if nostab:
				jds = self.datapoints.jds[self.datapoints.mask]
			else:
				jds = self.datapoints.jds
		else:
			# A minimal check for non-extrapolation condition should go here !
			pass
		
		fitmags = si.splev(jds, (self.t, self.c, self.k))
		# By default ext=0 : we do return extrapolated values
		return fitmags
	
	
	
	def display(self, showbounds = True, showdatapoints = True, showerrorbars=True, figsize=(16,8)):
		"""
		A display of the spline object, with knots, jds, stab points, etc.
		For debugging and checks.
		"""
	
		fig = plt.figure(figsize=figsize)
		
		if showdatapoints:
			if showerrorbars:
				mask = self.datapoints.mask
				plt.errorbar(self.datapoints.jds[mask], self.datapoints.mags[mask], yerr=self.datapoints.magerrs[mask], linestyle="None", color="blue")
				if not np.alltrue(mask):
					mask = mask == False
					plt.errorbar(self.datapoints.jds[mask], self.datapoints.mags[mask], yerr=self.datapoints.magerrs[mask], linestyle="None", color="gray")

			else:
				plt.plot(self.datapoints.jds, self.datapoints.mags, "b,")
				
	
		if (self.t != None) :
			
			if getattr(self, "showknots", True) == True:
				for knot in self.t:
					plt.axvline(knot, color="gray")
			
			# We draw the spline :
			xs = np.linspace(self.datapoints.jds[0], self.datapoints.jds[-1], 1000)
			ys = self.eval(jds = xs)
			plt.plot(xs, ys, "b-")
		
		if showbounds :
			if (self.l != None) and (self.u != None) :
				for l in self.l:
					plt.axvline(l, color="blue", dashes=(4, 4))
				for u in self.u:
					plt.axvline(u, color="red", dashes=(5, 5))
		
	
		axes = plt.gca()
		axes.set_ylim(axes.get_ylim()[::-1])
	
		plt.show()
		


# Some functions to interact directly with lightcurves :


def fit(lcs, knotstep=20.0, n=None, knots=None, stab=True,
	stabext=300.0, stabgap=20.0, stabstep=5.0, stabmagerr=-2.0, stabrampsize=0, stabrampfact=1.0,
	bokit=1, bokeps=2.0, boktests=5, bokwindow=None, k=3, verbose=True):
	"""
	The highlevel function to make a spline fit.
	
	lcs : a list of lightcurves (I will fit the spline through the merged curves)
	
	Specify either 
	knotstep : spacing of knots
	or
	n : how many knots to place
	or 
	knots : give me actual initial knot locations, for instance prepared by seasonknots.
	
	stab : do you want to insert stabilization points ?
	stabext : number of days to the left and right to fill with stabilization points
	stabgap : interval of days considered as a gap to fill with stab points.
	stabstep : step of stab points
	stabmagerr : if negative, absolte mag err of stab points. If positive, the error bar will be stabmagerr times the median error bar of the data points.
	
	
	bokit : number of BOK iterations (put to 0 to not move knots)
	bokeps : epsilon of BOK
	boktests : number of test positions for each knot
	
	
	"""
	dp = merge(lcs, stab=stab, stabext=stabext, stabgap=stabgap, stabstep=stabstep, stabmagerr=stabmagerr, stabrampsize=stabrampsize, stabrampfact=stabrampfact)
	
	s = Spline(dp, k=k, bokeps=bokeps, boktests=boktests, bokwindow=bokwindow)
	
	if knots==None:
		if n == None:
			s.uniknots(nint = knotstep, n = False)
		else :
			s.uniknots(nint = n, n = True)
	else:
		s.setintt(knots)
	#if stab:
	#	s.unistabknots(stabknotn,n=True)
	
	for n in range(bokit):
		s.buildbounds(verbose=verbose)
		s.bok(bokmethod="BF", verbose=verbose)
		
	s.optc()
	s.r2(nostab=True) # This is to set s.lastr2nostab
	
	return s


def seasonknots(lcs, knotstep, ingap, seasongap=60.0):
	"""
	A little helper to get some knot locations inside of seasons only
	
	knotstep is for inside seasons
	ingap is the number of knots inside gaps.
	
	"""
	knots = []
	
	#knotstep = 10
	dp = merge(lcs, splitup=True, deltat=0.000001, sort=True, stab=False)
		
	gaps = dp.jds[1:] - dp.jds[:-1]
	gapindices = list(np.arange(len(dp.jds)-1)[gaps > seasongap])
	
	# knots inside of seasons :
	a = dp.jds[0]
	for gapi in gapindices:
		b = dp.jds[gapi]
		#print (a, b)
		knots.append(np.linspace(a, b, float(b - a)/float(knotstep)))
		a = dp.jds[gapi+1]
	b = dp.jds[-1]
	knots.append(np.linspace(a, b, float(b - a)/float(knotstep)))
	
	
	# knots inside of gaps	
	for gapi in gapindices:
		a = dp.jds[gapi]
		b = dp.jds[gapi+1]
		knots.append(np.linspace(a, b, ingap+2)[1:-1])
	
	
	knots = np.concatenate(knots)
	knots.sort()	
	return knots
	
	#print gapindices
	
	
	"""
	for n in range(len(gapindices)):
		i = gapindices[n]
		a = self.jds[i]
		b = self.jds[i+1]
			
			newgapjds = np.linspace(a, b, float(b-a)/float(self.stabstep))[1:-1]
			newgapindices = i + 1 + np.zeros(len(newgapjds))
			newgapmags = np.interp(newgapjds, [a, b], [self.mags[i], self.mags[i+1]])
			newgapmagerrs = absstabmagerr * np.ones(newgapmags.shape)
			newgapmask = np.zeros(len(newgapjds), dtype=np.bool)
			
			self.jds = np.insert(self.jds, newgapindices, newgapjds)
		
	
	knotstep

	"""



def r2(lcs, spline, nosquare=False):
	"""
	I do not modify the spline (not even its datapoints) !
	Just evaluate the quality of the match, returning an r2 (without any stab points, of course).

	This is used if you want to optimize something on the lightcurves without touching the spline.
	
	Of course, I do not touch lastr2nostab or lastr2stab of the spline ! So this has really nothing
	to do with source spline optimization !
	
	
	"""
	
	myspline = spline.copy()
	newdp = pycs.gen.spl.merge(lcs, stab=False) # Indeed we do not care about stabilization points here.
	myspline.updatedp(newdp, dpmethod="leave")
	return myspline.r2(nostab=True, nosquare=nosquare)



def mltv(lcs, spline, weight=True):
	"""
	Calculates the TV norm of the difference between a lightcurve (disregarding any microlensing !) and the spline.
	I return the sum over the curves in lcs.
	
	Also returns a abs(chi) like distance between the lcs without ML and the spline
	
	If weight is True, we weight the terms in sums according to their error bars.
	
	
	Idea : weight the total variation somehow by the error bars ! Not sure if needed, the spline is already weighted.
	"""
	
	#import matplotlib.pyplot as plt
	tv = 0.0
	dist = 0.0
	for l in lcs:
		
		# We have a spline, and a lightcurve
		lmags = l.getmags(noml = True) # We get the mags without ML (but with mag and fluxshift !)
		ljds = l.getjds() # Inluding any time shifts.
		
		# Evaluating the spline at those jds :
		splinemags = spline.eval(ljds)
		
		# The residues :
		res = lmags - splinemags
		
		#plt.plot(ljds, res, "r.")
		#plt.show()
		
		if weight == False:
			tv += np.sum(np.fabs(res[1:] - res[:-1]))
			dist += np.sum(np.fabs(res))
		
		else:
			magerrs = l.getmagerrs()
			a = res[1:]
			aerrs = magerrs[1:]
			b = res[:-1]
			berrs = magerrs[:-1]
			
			vari = np.fabs(a - b)
			varierrs = np.sqrt(aerrs * aerrs + berrs * berrs)
			
			tv += np.sum(vari/varierrs)
			dist += np.sum(np.fabs(res) / np.fabs(magerrs))
		
			
	return (tv, dist)
	


def optcmltv(lcs, spline, verbose=True):
	"""
	I will optimize the coefficients of the spline so to minimize the mltv.
	I do not use the microlensing of the lcs at all !
	
	Simple powell optimization, slow. A pity.
	
	Add BOK and time shifts in there and it might be bingo !
	
	Would be more efficient if we add knots on the fly
	"""
	
	
	inic = spline.getc(m=2)
	
	def setc(c):
		spline.setc(c, m=2)
			
	def errorfct(c):
		setc(c)
		(tv, dist) = mltv(lcs, spline, weight=False)
		print "put weight"
		return tv + 0.1*spline.tv()

	minout = spopt.fmin_powell(errorfct, inic, full_output=1, disp=verbose)
	copt = minout[0]

	# We find a common shift to all coeffs so that the level matches
	
	meanc = np.mean(spline.getc(m=2))
	meanmag = np.mean(np.concatenate([l.getmags(noml = True) for l in lcs]))
	
	setc(copt)
	spline.c += meanmag - meanc	
	
	
	
	
	
