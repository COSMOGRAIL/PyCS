"""
What we call a season is in fact a way of specifying any "connected" range of points in a lightcurve.
Seasons are needed for microlensing, polynomial methods, etc.
Functions that want or return "seasons" usually handle lists of these season objects.
And as such lists of seasons are so common, there are a few functions below that directly handle lists of seasons.


"""


import lc
from numpy import *
#import copy as pythoncopy


class season:
	"""
	A season is essentially a numpy integer array containing the indexes of the lightcurve to consider.
	I made this class to allow for a bit more explicit error checking, and display capabilities etc.
	"""
	
	def __init__(self, indices, name=""):
		"""
		Prefer to use the factory functions below.
		
		@type	indices: int array
		@param	indices: The included indices of the season
		"""
	
		#self.indices = arange(minindex, maxindex+1)
		self.indices = indices
		"""@type: array of ints
		@ivar: An array of indices, that codes a season. You can use it like lc.jds[season.indices] to get the jds of that season.
		"""
		
		self.name = name
		
		#self.origin = "None"
		
		self.validate()
	
	def __str__(self):
		return "Season %s [%i, %i]" % (self.name, self.indices[0], self.indices[-1])
	
	def getjdlims(self, lc):
		"""
		Returns a tuple of two jds giving the extremities of the season.
		"""
		self.checkcompatibility(lc)
		tempjds = lc.getjds()
		return(tempjds[self.indices[0]], tempjds[self.indices[-1]])
		
	
	def validate(self):
		"""
		This checks some inherent properties of a season, that must be true independent of any lightcurve.
		"""
		
		if type(self.indices).__name__ != "ndarray":
			raise RuntimeError, "A seasons indices must be an array !"
		if self.indices.dtype != dtype('int'):
			raise RuntimeError, "A seasons indices must be an array of ints !"

		if len(self.indices) > 1:
			first = self.indices[:-1]
            		second = self.indices[1:]
			if not alltrue(less_equal(first,second)): # this is very fast, better then a loop
				raise RuntimeError, "A season must be sorted !"
			if self.indices[-1] - self.indices[0] + 1 != len(self.indices):
				print "WARNING : your season has holes !?"
	
	def checkcompatibility(self, lc):
		"""
		This checks if the season is compatible with the given lightcurve.
		"""
		
		if (max(self.indices) + 1) > len(lc.jds):
			raise RuntimeError, "Season not compatible, index %i too long for %s" % (max(self.indices), lightcurve)
		

# Factory functions to make such seasons, directly in form of lists :

def autofactory(l, seasongap=60, tpe='seas'):
	"""
	A factory function that automatically create season objects. Formerly known as the lightcurve method "autoseasons".
	Returns a LIST of season objects, describing the "shape" of the seasons of lc.
	In this case, these seasons are determined simply by looking for gaps larger then seasonsgap in the lightcurve.
	
		mode option added 11/2013:
		
		mode='seasons' (default option) ,  do create the season objects descibed above
		mode='interseasons' , return the extremity days between seasons, useful for cutting seasons in a synthetic lightcurve.
		
		
	Small example : suppose you have a lightcurve called C{mylc}, then :
	
		>>> myseasons = sea.autofactory(mylc)
		>>> mylc.mags[myseasons[0].indices] += 0.2
	
	... increases the mags of the first season by 0.2 !
	Or perhaps you want to mask the second season :
	
		>>> mylc.mask[myseasons[1].indices] = False

	@rtype: list
	@return: list of season objects.
	
	@todo: rewrite the example above.
		
	"""
	if not isinstance(l, lc.lightcurve):
		raise RuntimeError, "Please give me a lightcurve, not something else."
	
	# we find the indices just before season gaps 
	tempjds = l.getjds()
	first = tempjds[:-1]
	second = tempjds[1:]
	gapindices = where(second - first > seasongap)[0] + 1
            		
	# make one big vector of all the indices :
	indices = arange(len(tempjds))
	# split it according to the gaps :
	indexlist = hsplit(indices, gapindices)
	seasons = [season(ind, "a%i" % (i+1)) for (i, ind) in enumerate(indexlist)]		
	
	#seasonlengths = [len(season) for season in returnlist]
	#print "Lengths of seasons :", seasonlengths
	validateseasons(seasons)
	
	if tpe == 'seas':
		return seasons
	
	  
	if tpe == 'interseasons':
		
		"""
		Here, we play with the seasons, and return a list of 2d tuples (jdsmin, jdsmax)
		which determine the extremity days of the inter-season. It is useful e.g. to play with
		days out of observations in a synthetic lightcurve (as done in the spldiff method)		 
		"""
	    	
		# by the way, this is not optimized at all, but as the operation is nearly costless we don't care...(for now)
		          	
		beg=[]
		end=[]
		begint=[]
		endint=[]
		intsea=[]
		
		for sea in seasons:
			beg.append(sea.getjdlims(l)[0])
			end.append(sea.getjdlims(l)[1])
		
		for ind in arange(len(seasons)-1):
			begint.append(end[ind])  
			endint.append(beg[ind+1])
		
		for ind in arange(len(seasons)-1):
			intsea.append((begint[ind],endint[ind]))
		
		return intsea
			
	

	


def manfactory(lc, jdranges):
	"""
	A factory function that returns a LIST of season objects.
	As input, you give a lightcurve object as well as a list of JD ranges.
	Note that if the lightcurve is shifted, shifted jds are used !
	
	Example:
		>>> seas1 = sea.manfactory(lc1, [[2800, 3500], [3500, 4260]])

	@rtype: list
	@return: list of season objects.
	
	"""
	
	returnlist = []
	indices = arange(len(lc.jds))
	tempjds = lc.getjds()
	for i, jdrange in enumerate(jdranges):
		seasindices = indices[logical_and(tempjds > jdrange[0], tempjds < jdrange[1])]
		if len(seasindices) == 0:
			raise RuntimeError, "Empty seasons, check your ranges !"
		newseas = season(seasindices, "m%i" % (i+1))
		
		returnlist.append(newseas)
		
	validateseasons(returnlist)
	return returnlist
	
	


def validateseasons(seasons):
	"""
	Checks if your seasons (i.e. a list of season objects) can really be considered as seasons ...
	This is typically called by the microlensing constructor.
	Seasons can well be "incomplete", i.e. you can give only one season for your 4 year lightcurve, but we do not allow for overlapping seasons,
	as this would -up microlensing.
	"""
	
	# First we check the general structure.
	if type(seasons).__name__ != "list":
		raise RuntimeError, "Your seasons are not even a list..."
	if not len(seasons) >= 1: 
		raise RuntimeError, "Seasons are emtpy. No seasons = no sense, so I prefer to stop here."
	
	for season in seasons:
		season.validate()
		
	# Now some higher level properties :
	
	allindices = concatenate([s.indices for s in seasons])
	if len(allindices) != len(list(set(list(allindices)))):
		raise RuntimeError, "No overlapping seasons please !"
	

def printinfo(seasons):
	for s in seasons:
		print str(s)



def easycut(lclist, keep=(1), seasongap=60, mask=False, verbose=False):
	"""
	Wrapper function to easily cut out some seasons from a list of lightcurves.
	Seasons are determined for each lightcurve individually
	
	keep : a tuple of "indices" (human convention, starting at 1) of seasons to keep.
	
	mask : if True, I will not cut the seasons, but only mask them out.
	"""
	
	for l in lclist:
		seasons = autofactory(l, seasongap=seasongap)
		if verbose:
			print "%i seaons in %s" % (len(seasons), str(l))
		
		# We check user provided indices ...
		for humanindex in keep:
			if humanindex -1 not in range(len(seasons)):
				raise RuntimeError("There is no season %i !" % (humanindex))
		
		
		seasonstomask = [season for (i, season) in enumerate(seasons) if i+1 not in keep] 
		
		if l.hasmask():
			print "Warning : %s has already a mask, I might get cut !" % (str(l))
		
		for season in seasonstomask:
			l.mask[season.indices] = False
		
		if mask == True:
			continue
		else :
			l.cutmask()

		
		
		
		
		
	
