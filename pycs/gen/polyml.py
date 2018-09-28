"""

One season -- one microlensing
This potentially allows to optimize one season after the other, thus reducing the number of simultaneous parameters,
but of course you can still choose to optimize evrything at once (when using dispersion).
When using a spline fit, optimizing a polynom corresponds to weighted to linear least squares -- fast !



Idea for the future : "mutable" objects should be used here, with pointers in tables to pass the positions of each ml parameter very quickly.
The present implementation is not that elegant... we make to many copies.

We still need a more sophisticated factory function, so that you can set masks on individual params.
This is not possible for now, only because of the lack of such a factory.
(edit : not sure if really needed, I never use these...)
Let's forget about these "free params" and masks.


"""

# code to test legendre poly :
# import sys
# import numpy as np
# 
# import matplotlib.pyplot as plt
# import scipy.special as ss
# 
# 
# xs = np.linspace(-1, 1, 1000)
# 
# 
# for n in range(2, 2):
# 	poly = ss.legendre(n,0)
# 	ys = poly(xs)
# 	print n, ys[0]
# 	plt.plot(xs, ys)
# 
# 
# plt.xlim(-1, 1)
# plt.ylim(-1,1)
# plt.show()



import lc
import sea
import numpy as np
import copy as pythoncopy
import scipy.special as ss
import pycs.gen.sea



def polyfit(jds, mags, magerrs, nparams):
	"""
	I define my own polyfit, as numpy polyfit does not accept weights until numpy version 1.5
	This is linear least squares -- can only be used when optimizing polyml so to fit a "model" like a spline.
	
	What we want to do :
	fitp = np.polyfit(jds, mags, deg = nparams-1)
	
	:param jds: x values
	:param mags: y values
	:param magerrs: y errors
	:param nparams: number of parameters of polynom to fit.
	
	
	Code to test that the "heteroskedasticity correction" works :
	
	::
		
		import matplotlib.pyplot as plt

		x = np.linspace(0, 10, 20)
		ystd = 1.0 * np.ones(len(x))
		yerrs = ystd * np.random.randn(len(x))
		yerrs[-3] = 40.0
		ystd[-3] = 40.0 
		print yerrs

		y = 5.0 + 2.0*x - 1.0*x*x + yerrs

		finex = np.linspace(0, 10, 100)

		plt.errorbar(x, y, ystd)
	
		polyp = polyfit(x, y, ystd, 3)
		#polyp = np.polyfit(x, y, 3)

		finey = np.polyval(polyp, finex)

		plt.plot(finex, finey)
		plt.show()
		

	
	"""
	
	pows = np.arange(nparams)

	# We do "generalized least squares" aka "heteroskedasticity correction"
	a = np.column_stack([(jds ** pow)/magerrs for pow in pows[::-1]]) 
	fitp = np.linalg.lstsq(a, mags/magerrs)[0]
	
	# Not using weights, this would be :
	#a = np.column_stack([(jds ** pow) for pow in pows[::-1]]) 
	#fitp = np.linalg.lstsq(a, mags)[0]
	
	# Those fitp are directly the polynom coeffs as returned by np.polyfit
	return fitp





class seasonfct:
	"""
	A seasonfct object is one season + a certain polynom (simple poly or legendre for now -- from outside you do not see what it is) + some stuff
	to play with that.
	By itself, it makes no sense. It needs to be linked to a lightcurve : this is done by adding it to lc.ml, wich is a
	microlensing object, i.e. essentially a list of seasonfct objects.
	
	The idea is that the number of params and the mask remain fixed for one given seasonfct object.
	This makes the code simpler and faster.
	So once the instance is made, only the values of the parameters can be updated.
	If you want to change something else, simply make a new object.
	
	One idea of the relatively simple get and set methods is that one day we might want to use something else
	then polynoms, and at that time all this might be handy.
	
	for the mltype argument, see discussion of ml.factory function.
	
	"""

	def __init__(self, season, mltype="poly", params = np.array([0.0]), mask = np.array([True])):
		"""
		Also see the factory functions below ... they will be much easier to use for your task.
		Default values shown above : 1 parameter, i.e. a simple additive constant for the given season, which is free to be optimized.
		"""
		
		# Some basic checks :
		#if not isinstance(season, sea.season):
		#	raise RuntimeError, "microlensing : I want a season object."
		if len(params) != len(mask):
			raise RuntimeError, "I want params and mask to be of the same length !"
		
		
		self.season = season 
		"""@type: season
		@ivar: The season for which this microlensing if for.
		"""
		
		if mltype not in ["poly", "leg"]:
			raise RuntimeError, "Sorry I do not know this microlensing mltype."
		self.mltype = mltype
		"""@type: string
		@ivar: The type of microlensing
		"""
		
		self.params = np.asarray(params) 
		"""@type: ndarray
		@ivar: A 1D numpy array of params (floats)
		"""
		
		self.mask = np.asarray(mask)
		"""@type: ndarray
		@ivar: boolean array : True if the parameter is "free".
		"""
		
		self.nfree = np.bincount(self.mask)[-1]
		"""@type: int
		@ivar: How many free parameters
		"""
		# This remains fixed !!!
		
	def copy(self):
		return pythoncopy.deepcopy(self)


	def __str__(self):
		"""
		This is just the number of params of the polynom (i.e. degree + 1)
		(Not the number of free params !)
		"""
		return "%i" % (len(self.params))
		
	def longinfo(self):
		"""
		Returns a longer description of the object.
		"""
		if self.mltype == "poly":
		
			strlist=["Plain poly with %s params on %s :" % (str(self), str(self.season))]
			for (i, p, m) in zip(range(len(self.params)-1, -1, -1), self.params, self.mask):
				if m:
					paramstr = "\n deg %2i : %g (free)" % (i, p)
				else:
					paramstr = "\n deg %2i : %g" % (i, p)
				strlist.append(paramstr)
			return "".join(strlist)
			
		
		if self.mltype == "leg":
		
			strlist=["Legendre poly with %s params on %s :" % (str(self), str(self.season))]
			for (n, p, m) in zip(range(len(self.params)), self.params, self.mask):
				if m:
					paramstr = "\n deg %2i : %g (free)" % (n, p)
				else:
					paramstr = "\n deg %2i : %g" % (n, p)
				strlist.append(paramstr)
			return "".join(strlist)
		
		
	def printinfo(self):
		"""
		Prints the longer description of the object.
		"""
		print self.longinfo()
		
	
	#def resetparams(self):
	#	self.params = zeros(len(self.params))
		
	def setparams(self, p):
		"""
		Remember that you are not allowed to change the number of parameters !
		"""
		if len(p) == len(self.params):
			self.params = p # no further testing ...
		else :
			raise RuntimeError, "Wrong number of parameters !"
		
	def getparams(self):
		return self.params 	# perhaps put a copy here ? Let's see how we use this.
					# it would be nice to return pointers one day...
		
	def setfreeparams(self, p):
		"""
		Here we should really do some testing ... be careful when using this !!!
		"""
		if len(p) == self.nfree:
			self.params[self.mask] = p # no further testing ...
		else :
			raise RuntimeError, "Wrong number of free parameters !"
		
	
	def getfreeparams(self):
		return self.params[self.mask] # put a copy here ?
	
	def validate(self):
		"""
		@todo: code this in case of ideas...
		"""
		pass
		
	def checkcompatibility(self, lightcurve):
		"""
		Not sure if needed.
		"""
		self.season.checkcompatibility(lightcurve)

	def calcmlmags(self, lightcurve):
		"""
		Returns a "lc.mags"-like array made using the ml-parameters.
		It has the same size as lc.mags, and contains the microlensing to be added to them.
		The lightcurve object is not changed !
		
		For normal use, call getmags() from the lightcurve.
		
		Idea : think about only returning the seasons mags to speed it up ? Not sure if reasonable, as no seasons defined outside ?
		"""
		jds = lightcurve.jds[self.season.indices] # Is this already a copy ? It seems so. So no need for an explicit copy().
		# We do not need to apply shifts (i.e. getjds()), as anyway we "center" the jds.
		
		
		# Old method :
		if self.mltype == "poly":
		
			refjd = np.mean(jds)
			jds -= refjd # This is apparently safe, it does not shifts the lightcurves jds.
		
			allmags = np.zeros(len(lightcurve.jds))
			allmags[self.season.indices] = np.polyval(self.params, jds) # probably faster then +=
			return allmags
		
		
		# Legendre polynomials :
		if self.mltype == "leg":
		
			rjd = (np.max(jds) - np.min(jds))/2.0
			cjd = (np.max(jds) + np.min(jds))/2.0
			jds = (jds - cjd)/rjd
		
			allmags = np.zeros(len(lightcurve.jds))
		
			for (n, p) in enumerate(self.params):
				allmags[self.season.indices] += p * ss.legendre(n)(jds)
		
			return allmags
		
	
		

	def smooth(self, lightcurve):
		"""
		Only for plotting purposes : returns jds, mlmagshifts, and refmags with a tight and regular sampling,
		over the range given by the season.
		Note that this time we are interested in the acutal shifted jds, so that it looks right when plotted !
		We return arrays that can directly be plotted to illustrate the microlensing.
		
		TODO : return refmag, not refmags !
		"""

		jds = lightcurve.getjds()[self.season.indices]

		# Old method :
		if self.mltype == "poly":
		
			refjd = np.mean(jds)
			jds -= refjd
		
			refmag = np.median(lightcurve.getmags()) 	# So the reference magnitude is evaluated from the entire lightcurve.
							# Independent on seasons.
		
			smoothtime = np.linspace(jds[0], jds[-1], 50)
			smoothml =  np.polyval(self.params, smoothtime)
			smoothtime += refjd # important, to get the time back at the right place.
			refmags = np.zeros(50) + refmag
			return {"jds":smoothtime, "ml":smoothml, "refmags":refmags}
		
		
		# Legendre polynomials :
		
		if self.mltype == "leg":
		
			rjd = (np.max(jds) - np.min(jds))/2.0
			cjd = (np.max(jds) + np.min(jds))/2.0
			jds = (jds - cjd)/rjd
		
			refmag = np.median(lightcurve.getmags()) 	# So the reference magnitude is evaluated from the entire lightcurve.
							# Independent on seasons.
		
			smoothtime = np.linspace(-1, 1, 300)
			smoothml = np.zeros(len(smoothtime))
			for (n, p) in enumerate(self.params):
				smoothml += p * ss.legendre(n)(smoothtime)
			
			smoothtime = smoothtime * rjd + cjd # important, to get the time back at the right place.
			refmags = np.zeros(300) + refmag
			return {"jds":smoothtime, "ml":smoothml, "refmags":refmags}
		
		


class microlensing:
	"""
	Contains a list of seasonfct objects OF ONE SAME LIGHTCURVE OBJECT, and some methods to treat them.
	
	You probably do not want your seasonfct seasons to overlap. But no problem if they do.
	
	Again : do not change the contents of such a microlensing object otherwise then retrieving and updating the params
	through the provided methods !
	Do not change the list of seasonfct objects in any other way !
	1) build the seasonfct objects
	2) put them into a microlensing object
	3) do not touch the seasonfct objects anymore if you do not know what you are doing.
	
	
	Instead, make a brand new object if you need, using factory functions provided below or to be written.

	"""


	def __init__(self, mllist):
		self.mllist = mllist 
		"""@type: list
		@ivar: The list of seasonfct objects to be applied.
		"""
		
		self.nfree = np.sum(np.array([sfct.nfree for sfct in self.mllist]))
		"""@type: int
		@ivar: The number of free parameters.
		"""
		
		self.mltype = self.mllist[0].mltype # ok this is a bit mean...
		
	
	def copy(self):
		return pythoncopy.deepcopy(self)

	def __str__(self):
		return "".join(["|%s/" % self.mltype] + ["%s" % (m) for m in self.mllist] + ["|"])
		
	def longinfo(self):
		return "\n".join(["%s" % (m.longinfo()) for m in self.mllist])
	
	def printinfo(self):
		print self.longinfo()

	def getfreeparams(self):
		"""
		Straightforward.
		"""
		return np.concatenate([sfct.getfreeparams() for sfct in self.mllist])
		# I tested this, the concatenate makes a copy, otherwise 
		# we would be mutable, that would be nice to solve these unelegant setfreeparams problem.
	
	#def getfreeshape(self):
	#	return getfreeshape(self.mllist)
		
	def setfreeparams(self, p):
		"""
		This method distributes the params on the microlensing objects as fast as I could ... it's a bit delicate and unelegant.
		As we want to be fast -- do not mess with this and give a p with exactly the right size.
		"""
	
		if len(p) == self.nfree:
			startindex = 0
			for sfct in self.mllist:
				stopindex = startindex + sfct.nfree
				if len(p) == 1:
					sfct.setfreeparams(p)
				else:
					sfct.setfreeparams(p[startindex : stopindex])
				startindex += sfct.nfree		
			#if startindex != len(params):
			#	raise RuntimeError, "You provided too many params !"
		else :
			raise RuntimeError, "Wrong number of free parameters !"
	
	def reset(self):
		"""
		Puts all coefs back to 0.0
		Allows to start from scratch without having to rebuild a new ML.
		"""
		self.setfreeparams(0.0 * self.getfreeparams()) # A bit sad, I agree :)
		

	def checkcompatibility(self, lightcurve):
		for sfct in self.mllist:
			sfct.checkcompatibility(lightcurve)
		
	def calcmlmags(self, lightcurve):
		"""
		Returns one a "lc.mags"-like array made using the parameters of all seasonfct objects.
		This array has the same size as lc.mags, and contains the microlensing to be added to lc.mags.
		
		Idea : perhaps we can make this faster by not calling calcmags of the seasonfct objects ?
		"""

		allmags = np.zeros(len(lightcurve.jds))
		for microlensing in self.mllist:
			allmags += microlensing.calcmlmags(lightcurve)
		return allmags


	def stats(self, lightcurve):
		"""
		Calculates some statistics on the microlensing deformation evaluated at the same sampling
		than the lightcurve.
		The idea is to get the flux ratios, in case of calm microlensing.
		
		"""
		
		mlmags = self.calcmlmags(lightcurve)
		
		mlmean = np.mean(mlmags)
		mlstd = np.std(mlmags)
		
		return {"mean":mlmean, "std":mlstd}



def multigetfreeparams(lclist):
	"""
	Give me a list of lightcurves (with or without ml !), I give you a single flat array of parameters of all MLs concatenated.
	Note that the order of lightcurves in lclist is important ! You will have to give the same order for multisetfreeparams() !
	
	For now a bit simple ... but it's not that slow.
	"""
	
	params = np.array([]) # we start with an empty array
	for curve in lclist :
		if curve.ml != None:
			#print "Microlensing for %s" % str(curve)
			params = np.append(params, curve.ml.getfreeparams())
		else:
			#print "No microlensing for %s" % str(curve)
			pass
	
	#print "Joined free params : ", params
	if len(params) == 0:
		print "WARNING : there are no free ml params !"
	return params 



def multisetfreeparams(lclist, params):
	"""
	Be careful to respect the order of lcs in lclist ... otherwise everything gets messed up.
	Again this seems a bit slower then it really is -- prefer simplicity.
	"""

	startindex = 0
	for curve in lclist:
		if curve.ml != None:
			stopindex = startindex + curve.ml.nfree
			if stopindex > len(params):
				raise RuntimeError, "Something is fishy with your params..."
			
			if len(params) == 1: # special solution needed if only 1 parameter :-(
				curve.ml.setfreeparams(params)
			else:
				curve.ml.setfreeparams(params[startindex : startindex + curve.ml.nfree])
			startindex += curve.ml.nfree
			
	if startindex != len(params):
		raise RuntimeError, "You provided too many params !"
			


# What we still need is a method that collects and sets all the free params of a list of lightcurves.
	

def factory(seasons, nparams, mltype="poly"):
	"""
	A factory function to create a microlensings object filled by seasonfct objects.
	seasons is a list of season objects
	nparams is an array or list of ints. "default" = one constant per season.
	All parameters will be set to 0.0, and free to be optimized.
	
	mltype = "poly" : simple polynomial microlensing, very stupid but fast, ok for degree <= 3
	default type.
	
	mltype = "leg" : legendre polynomials, very clever but slow :-) These are fine for degree <= 25
	Above deg 25, some numerical errors set in. Could perhaps be rewritten to make this faster.
	Or implement in C !!!!
	
	"""
	#if nparams == "default":
	#	nparams = ones(len(seasons), dtype="int")
	
	if len(nparams) != len(seasons):
		raise RuntimeError, "Give as many nparams as they are seasons !"
	
	mllist = []
	for (season, n) in zip(seasons, nparams):
		if n != 0:
			p = np.zeros(n, dtype="float")
			mask = p > -1.0
			sfct = seasonfct(season, mltype, p, mask)
			mllist.append(sfct)

	return microlensing(mllist)


def addtolc(l, seasons=None, nparams=1, autoseasonsgap = 60.0):
	"""
	Adds polynomial ML to the lightcurve l.
	Top level function, to make it really easy ! We just wrap the factory above.
	
	If seasons = None, I will make some autoseasons. You can specify autoseasonsgap.
	Else, seasons should be a list of seasons (see factory function above)
	
	If nparams is an int, each season will get a polynom of nparams parameters.
	1 = contant, 2 = slope, ...
	Else nparams should be a list of ints (see factory function above)
	
	If you want one single polynom for the full curve, just set autoseasonsgap to 10000.0 ...
	
	::
	
		
		pycs.gen.polyml.addtolc(l, nparams=1) # One constant on each season.
		
	
	
	"""
	
	if seasons == None:
		seasons = pycs.gen.sea.autofactory(l, seasongap=autoseasonsgap)
	if type(nparams) == int:
		nparams = [nparams] * len(seasons)
	
	m = pycs.gen.polyml.factory(seasons, nparams, mltype="poly")
	
	l.addml(m)
	


#class microlensing:
#	"""
#	A microlensing object is a representation of microlensing for a specific lightcurve.
#	
#	Note that parameters defining the ml are defined to be relative to the mean julian date in each season.
#	This way, simply shifting a curve does not intrinsically change the aspect of the microlensing for a given set of parameters.
#	
#	Microlensing is represented by a 0th, 1st, or 2nd order polynom added to each season.
#	You have to specify the seasons.
#		
#	You have to specify what kind of microlensing you want to optimize for each season, i.e. what kind of constrains to have.
#	You do so by specifying a string containing as many chars as there are seasons :
#		- "2" for 2nd degree polynom
#		- "1" for 1st degree polynom
#		- "0" for 0 degree polynom = simple magnitude shift.
#		- "-" for no microlensing
#	"""
#	
#	def __init__(self, seasons, constrains="default"):
#		"""
#		seasons in the same format as given by lc.autoseasosn(), that is a list of arrays of indices.
#		The default value for constains will be like "-000" i.e. first season fixed and following seasons represented with a constant.
#		
#		"""
#		
#		sea.validateseasons(seasons)
#		
#		self.seasons = seasons
#		"""@type: list of season objects
#		@ivar: Like given for instance by the lc.autoseasons method.
#		"""
#		
#		if constrains=="default":
#			self.constrains = "-" + "0"*(len(seasons)-1)
#			"""@type: string
#			@ivar: A string of same length then the seasons, encoding how each season's microlensing should be treated"""
#		else :
#			self.constrains = constrains
#		if len(self.constrains) != len(self.seasons):
#			raise RuntimeError, "Microlensing : seasons and constrains have different lengths !"
#		
#		
#		self.params = zeros((len(seasons), 3)) # 3 as max order is 2 for now -> a, b c
#		"""@type: 2D array
#		@ivar: A table containing the actual parameters of the microlensing. Don't access them directly from outside,
#		use the getparams and setparams methods ! params[i] can be seen as a 1D array with elements [c, b, a] in
#		a x^2  +  b x  +  c
#		"""
#		
#		self.pmask = [] #  we fill this below
#		"""@type: 2D array
#		@ivar: array of same size as params, tells which one are free (True = free)."""
#		for code in self.constrains :
#			if code == "2":
#				self.pmask.append([True, True, True])
#			elif code == "1":
#				self.pmask.append([True, True, False])
#			elif code == "0":
#				self.pmask.append([True, False, False])
#			elif code == "-":
#				self.pmask.append([False, False, False])
#			else:
#				raise RuntimeError, "Microlensing : illegal char %s in constrains '%s' !" % (code, self.constrains)
#	
#		# we make an array out of this :
#		self.pmask = array(self.pmask)
#
#		#self.pshape = self.params.shape
#		#"""@type: tuple of ints
#		#@ivar: shape of the params table """
#
#		self.nbrfree = len(self.params[self.pmask].ravel())
#		"""@type: int
#		@ivar: gives the number of free params (i.e. length of the flatparams from setparams and getparams),
#		according to the constrains string"""
#	
#	def __str__(self):
#		"""Representation of the microlensing in form of a string of constrains"""
#		return "|%s|" % (self.constrains)
#	
#	def resetparams(self):
#		"""Simply puts the whole parameter array back to 0.0"""
#		
#		self.params = zeros((len(self.seasons), 3))
#		
#	
#	def setparams(self, flatparams):
#		"""The goal of this set and get methods are to handle flat parameter arrays, making it easier to optmize those from outside this class.
#		
#		@type flatparams: array
#		@param flatparams: A flat array of floats, containing values for the free parameters.
#		"""
#		
#		if len(flatparams) != self.nbrfree :
#			raise RuntimeError, "Microlensing : wrong length of flatparams, need %i." % self.nbrfree
#		
#		self.params[self.pmask] = flatparams # yes, this works. we update only the free parameters, shape stays
#		
#	
#	
#	def getparams(self):	
#		"""Returns all the needed free parameters as one flat 1D array
#		
#		@rtype: 1D array
#		@return: A flat array of floats containing the free parameters that should get optimized.
#		"""
#		return self.params[self.pmask]
#	
#	
#	def printinfo(self):
#		"""Prints a short info about the status of this object, for debugging purposes."""
#		
#		print "Microlensing info :"
#		for i, season in enumerate(self.seasons):
#			print "Season %s, degree %s : (%s)" % (i, season, self.constrains[i], ", ".join(["%f" % param for param in self.params[i]]))
#		print "This gives : (%s)" % ", ".join(["%f" % param for param in self.getparams()])	
#
#
#	def checkcompatibility(self, lightcurve):
#		"""This is a bit ridiculous... but anyway. We check that seasons are not longer then they should.
#		At this stage we do not check that seasons validate... this should be done before (constructor).
#		"""
#	
#		for season in self.seasons:
#			season.checkcompatibility(lightcurve)
#		
#		
#
#	def calcmags(self, lightcurve):
#		"""Returns a "mags"-like vector made using the current parameters microlensing applied to the mags of a lightcurve.
#		The lightcurve object is not changed !
#		We do not take into accound any non-applied magshift. For normal use, call getmags() from the lightcurve.
#		"""
#		
#		magscopy = lightcurve.mags.copy() 	# this way we do not modify the lightcurve object
#							# We WANT to directly access the mags attribute here !
#		for paramline, season in zip(self.params, self.seasons):
#			jds = lightcurve.jds[season.indices]
#			refjd = mean(jds)
#			magshifts = paramline[0] + paramline[1]*(jds - refjd) + paramline[2]*((jds - refjd)**2)
#			magscopy[season.indices] += magshifts
#			
#		return magscopy
#
#
#def getparams(lclist):
#	"""returns the microlensing params of a list of lightcurves, as provided by ml.getparams, concatenated
#	into one 1D array.
#	Note that the order of lightcurves in lclist is important ! You will have to give the same order for
#	setmlparams() !
#	
#	@type lclist: list of lightcurves
#	@param lclist: list of lightcurves you want to have the ml params from.
#	
#	@rtype: 1D array
#	@return: 1D array of concatenated ml params from these lightcurves.
#	"""
#	
#	params = array([]) # we start with an empty array
#	for curve in lclist :
#		if curve.ml != None:
#			#print "Microlensing for %s" % str(curve)
#			params = append(params, curve.ml.getparams())
#		else:
#			#print "No microlensing for %s" % str(curve)
#			pass
#	
#	#print "Joined free params : ", params
#	if len(params) == 0:
#		print "WARNING : there are no free ml params !"
#	return params 
#	
#	
#def setparams(params, lclist):
#	"""
#	Be careful to respect the order of lcs in lclist, if not everything gets messed up.
#	
#	@type params: 1D array
#	@param params: An array of params, similar as the one got from getmlparams().
#	
#	@type lclist: list of lightcurves
#	@param lclist: List of lightcurves to which you want to have the ml params set.
#	"""
#	
#	startindex = 0
#	for curve in lclist:
#		if curve.ml != None:
#		
#			stopindex = startindex + curve.ml.nbrfree
#			if stopindex > len(params):
#				raise RuntimeError, "ml.setparams : something is fishy with your params..."
#			
#			#print startindex, startindex + curve.ml.nbrfree
#			curve.ml.setparams(params[startindex : startindex + curve.ml.nbrfree])
#			#print "Set parmams for %s : %s" %(str(curve), str(params[startindex : startindex + curve.ml.nbrfree])) 
#			startindex += curve.ml.nbrfree
#	
#	if startindex != len(params):
#		raise RuntimeError, "ml.setparams : you provided too many params !"
#			
#
#		
#class seasonslopeml:
#	"""
#	Describes a microlensing modelled by "magshift = m * days + h" for each season of a lightcurve.
#	Note that parameters defining the slopes are defined to be relative to the mean julian date in each season.
#	This means that shifted curves will not have different parameters regarding microlensing !
#	
#	
#	ms is an array of slopes (in mags/timeunit)
#	hs is an array of steps (in mags)
#	"""
#	
#	def __init__(self, seasons):
#	
#		if len(seasons) == 0:
#			raise RuntimeError, "Wrong format for seasons !"
#	
#		self.nbpoints = sum(array(map(len, seasons)))
#		self.seasons = seasons
#		self.params = array([0.0] * 2*len(seasons))
#		self.distributeparams()
#		self.description = "Flat"
#			
##	def validate(self):
##		"""We check if ms, hs etc are coherent. Also we check if seasons is meaningfull."""
##		pass
#	
#	def distributeparams(self):
#		self.ms = self.params[:len(self.seasons)]
#		self.hs = self.params[-len(self.seasons):]
#				
#	def describe(self):
#		self.distributeparams()
#		self.description = ",".join(["(%f,%f)"%(m,h) for (s,m,h) in zip(arange(len(self.ms))+1,self.ms,self.hs)])
#	
#	def printinfo(self):
#		self.describe()
#		string_list = ["- "*30, "\n\t", "Number of seasons : %i"%len(self.seasons), "\n", self.description, "\n", "- "*30]
#		print ''.join(string_list)
#
#	
##	def iscompatible(self, lightcurve):
##		""" """
##		pass
#	
#	def applyto(self, lightcurve):
#		"""Applies this to the mags of a lightcurve"""
#		
#		if self.nbpoints != len(lightcurve.jds):
#			raise RuntimeError, "%s not compatible with %s"%(self.description, str(lightcurve))
#		
#		self.distributeparams()	
#		for h, m, season in zip(self.hs, self.ms, self.seasons):
#			jds = lightcurve.jds[season]
#			refjd = mean(jds)
#			magshifts = (jds - refjd)*m + h
#			lightcurve.mags[season] += magshifts
#				
#
#		self.describe()
#		lightcurve.commentlist.append("ML:" + self.description)
#
#	def copy(self):
#		"""Returns a "deep copy" of the object.""" 
#		return pythoncopy.deepcopy(self)

#def optimizeml(lc1, lc2, dispersionmethod):
#	"""
#	
#	@type lc1: lightcurve
#	@param lc1: This one will be passed as first argument (i.e. "reference") to the dispersionmethod.
#	@type lc2: lightcurve
#	@param lc2: The light curve that will get interpolated by the dispersionmethod.
#	
#	@todo: this has to be rewritten and tested !
#	"""
#
#	params = lc2.ml.getparams()
#	print "Initial params : ", params
#	
#	def d2(params):
#	
#		lc2ml = lc2.copy()
#		ml2.params = params
#		ml2.applyto(lc2ml) 
#		d2val = dispersionmethod(lc1, lc2ml)["d2"]
#		return d2val
#			
#	#res = spopt.fmin(d2, params, xtol=0.001, ftol=0.0001)
#	res = spopt.fmin_bfgs(d2, params)
#	print "Optimal params : ", res
#	
#	ml2.params = res
#	ml2.distributeparams
#	return ml2


#class seasonslopemlcurve(lightcurve):
#	"""
#	A microlensing affected curve in which microlensing consists of individual slopes for each season.
#	It inherits a lightcurve, and adds some new attributes and methods.
#	
#	Note that parameters defining the microlensing are defined to be relative to the first point in each season.
#	This means that shifted curves will not have different parameters regarding microlensing !
#	"""
#
#	# The constructor (will never be used, but serves to introduce new attributes) :
#	
#	def __init__(self, lc):
#		
#		# We call the parent constructor :
#		#lightcurve.__init__(self)
#		
#		
#		self.nbpoints = len(self.jds)
#		self.seasons = [arange(nbpoints)] # only one season for this demo
#		"""@type: List of 1D integer arrays
#		@ivar: Indexes of the seasons, as a list of 1D integer arrays, as returned for instance by autoseasons method."""
#		
#		self.hs = [0.3]		# mag shift of seasons. must have same length then seasons.
#		self.ms = [0.0001]	# slope of seasons
#		
#		self.magshifts = zeros(nbpoints)
#		self.description = "zeros"
#		self.intendedfor = str(lightcurve)
#




#class mlcurve():
#
#	def __init__(self, nbpoints):
#		
#		self.magshifts = zeros(nbpoints)
#		self.description = "zeros"
#		#self.intendedfor = str(lightcurve)
#
#
#	def applyto(self, lightcurve):
#		lightcurve.mags += self.magshifts
#		lightcurve.commentlist.append("ML : " + self.description)
#
#
##def mlfunction(lightcurve, 
#
#
## lightcurve object : add origmags and origjds ? Easier when handling ML no, add them to child object
#
#class seasonslopemlcurve(lightcurve):
#	"""
#	A microlensing affected curve. It is a lightcurve, but with bonus attributes.
#	keeps track of the original mags (redefine shiftmag and shifttime, but warning ?)
#	and commentlist
#	adds parameters, numer of seasons
#	a method applyml
#	"""
#	pass
#
#
#class seasonslopes():
#	""" note that parameters are absolute. If you shift a curve, they will not be the same. Do we want this ?
#	we could make them relative (to the first point of each season)
#	"""
#
#	def __init__(self, lightcurve, seasons):
#	
#		self.lc = lightcurve.copy() # we keep a copy
#		self.origmags = pythoncopy.copy(lightcurve.mags
#		self.origcommentlist = lightcurve.commentlist
#		self.ms = [0.0] * len(seasons)
#		self.hs = [0.0] * len(seasons)
#		self.description = "flat"
#		
#		nbpoints = sum(array(map(len, seasons)))
#		if nbpoints != len(self.lc.jds):
#			raise RuntimeError, "seasonslopes : params must have same length then seasons !"		
#
#	def getparams(self):
#		pass
#	
#	def setparams(self):
#		pass
#	
#	def validate(self):
#		pass
#	
#	def describe(self):
#		
#	
#	
#	def apply(self):
#		"""
#		just change self.lc, do not return anything
#		"""
#	
#"""
#other idea : make an object lightcurveml that inherits from lightcurve
#"""
#
#
#def seasonsteps(seasons, steps):
#	"""
#	factory function for mlcurve
#	seasons is a list of arrays of indices
#	steps are steps to apply to these seasons
#	returns a mlcurve object.
#	"""
#	
#	if len(seasons) != len(steps):
#		raise RuntimeError, "seasonsteps : steps must have same lenght then seasons !"
#	
#	nbpoints = sum(array(map(len, seasons)))
#	
#	mlc = mlcurve(nbpoints)
#	for step, season in zip(steps, seasons):
#		mlc.magshifts[season] = step
#	
#	mlc.description = "%i seasonsteps" % len(steps)
#		
#	return mlc
#
#
#def seasonslopes(lightcurve, seasons, m, h)
#	"""
#	m is a list of slopes (in mags/timeunit)
#	h is a list of steps (in mags)
#	"""
#	
#	if len(seasons) != len(m) or len(seasons) != len(h) :
#		raise RuntimeError, "seasonslopes : params must have same length then seasons !"
#	
#	mstring = ",".join(["%f"%item for item in m])
#	hstring = ",".join(["%f"%item for item in h])
#	mlc.description = "%i seasonslopes (m=%s, h=%s)" % (len(seasons), mstring, hstring)
#	return mlc
