"""
Microlensing represented by splines.
Rewamped in July 2011
"""


import numpy as np
import pycs.gen.spl
import copy as pythoncopy


class SplineML():
	"""
	A lightcurve can have such a microlensing object as attribute "ml".
	For now we have only one spline per lightcurve (i.e., no seasons).
	
	The splml.spline.datapoints are relative, first point is always at 0.0
	That's why we store a jdoffset EDIT : No, we don't store it anymore.
	The spline just starts at 0.0, nobody here cares about this offset... use the timeshift of my lightcurve,
	if you want to know where I am located on an absolute time axis.
	
	We dont use stabililzation points for microlensing splines (not needed I hope)
	Hence, the datapoints.jds are directly the jds - jdoffset from the associated lightcurve
	(They will not change)
	EDIT : we now use stab points, but still, the jds remain fixed.
	
	The datapoint mags = what has to be added to the lightcurve so that it matches the sourcespline (or something else).
	This is the same convention set as for ML by polynoms.
	
	"""
	
	def __init__(self, spline):
		"""
		
		EDIT : I do not store jdoffset anymore, makes no sense.
		jdoffset = Offset to add to spline.datapoints so to get back jds as they were at the moment of the factory call.
		Hmm, this is a bit stange, I think nobody needs to know this offset.
		If you want to plot the ML for instance, you should not use it, but take the current lc.timeshift, as the curve
		might have been shifted in time, and this jdoffset does not get updated !
		"""
	
		self.spline = spline
		self.mltype = "spline"
		
		
		# Question : why do we need this ?
		# settargetmags used them ! No, not anymore !
		#self.stab=stab
		#self.stabgap = stabgap
		#self.stabstep = stabstep
		#self.staberrs = staberrs
		
		# See remark above
		#self.jdoffset = jdoffset
		
	
	def __str__(self):
		return "|" + str(self.spline) + "|"
	
	#def longinfo(self):
	#	return "Bla"
	
	#def printinfo(self):
	#	print self.longinfo()
	
		
	def copy(self):
		return pythoncopy.deepcopy(self)
	
	def checkcompatibility(self, lightcurve):
		"""
		It should be checked here that the self.datapoints are comptabile with the lightcurve
		(i.e. same "relative" jds)
		"""
		# Stupid example :
		#if self.spline.datapoints.ntrue() != len(lightcurve):
		#	raise RuntimeError("Ouch, this lightcurve is not compatible with its SplineML datapoints !")
		# We should go further :
		if not np.all(np.fabs(self.spline.datapoints.jds[self.spline.datapoints.mask] - (lightcurve.jds - lightcurve.jds[0])) < 0.001):
			raise RuntimeError("Ouch, this lightcurve is no longer compatible with its SplineML datapoints !")

	def settargetmags(self, lightcurve, sourcespline):
		"""
		We update the self.spline.datapoints.mags so that a fit will results in a microlensing approximation.
		
		Sourcespline is a spline object that is the reference where you want your lightcurve to go.
		
		Only the mags are updated !
		Indeed for this microlensing spline I don't care if your lightcurve was shifted in time.
		But of course a shift in time with restpect to the sourcespline will give you differnent mags !
		
		I hope you did not remove points or other tweaks.
		
		@todo: treat magerrs with a getter method !
		
		EDIT : Let's try to add some stabilization points in the datapoints...
		The constructor of the SplineML now defines how to calculate these stabilisation points.
		These settings are permanent for one SplineML, hence like without any stabpoints, all we need to do
		is to update the mags.
		No need to update the stab point jds, or the datapoints.mask etc.
		But, to update the mags, we need to calculate the mags for these stab points.
		
		"""
	
		self.checkcompatibility(lightcurve)

		#if self.stab == False:
		#	#Then its very easy.
		#	self.spline.datapoints.mags = sourcespline.eval(jds = lightcurve.getjds()) - lightcurve.getmags(noml = True) # Note that we do not want to include the ML here of course !
		#	self.spline.datapoints.magerrs = lightcurve.magerrs # just in case they have changed, see todo above.
		#
		#	import matplotlib.pyplot as plt
		#	plt.plot(self.spline.datapoints.jds, self.spline.datapoints.mags, "r.")
		#	plt.show()
		#
		# But the following works in any case :
		
		jdoffset = lightcurve.getjds()[0] # That's the only stuff we'll use from the lightcurve.
			
		# Evaluating the sourcespline at the real jds :
		sourcerealjds = self.spline.datapoints.jds[self.spline.datapoints.mask] + jdoffset
		sourcerealmags = sourcespline.eval(sourcerealjds)
			
		# Getting the real mags :
		realpointmags = lightcurve.getmags(noml = True)  # Note that we do not want to include the ML here of course !
		realtargetmags = sourcerealmags - realpointmags
			
		# And now we need to interpolate the values for the stab points.
		# Is it better to interpolate the real mags, or directly the correction (the targetmags) ?
		# Perhaps better to interpolate the correction, as we do not want the ML to follow the intrinsic spline
		# anyway. We assume the ML is a simple as possible.
			
		sourcestabjds = self.spline.datapoints.jds[self.spline.datapoints.mask == False] + jdoffset
			
		# No we have :
		#sourcerealjds
		#realtargetmags
		#sourcestabjds
		# And want to know :
		#stabtargetmags	
		stabtargetmags = np.interp(sourcestabjds, sourcerealjds, realtargetmags, left=100.0, right=100.0)
		# left and right should never be needed ...
			
		#import matplotlib.pyplot as plt
		#plt.plot(sourcerealjds, realtargetmags, "r.")
		#plt.plot(sourcestabjds, stabtargetmags, "b.")
		#plt.show()
			
		# The update itself :
		self.spline.datapoints.mags[self.spline.datapoints.mask] = realtargetmags
		self.spline.datapoints.mags[self.spline.datapoints.mask == False] = stabtargetmags
			
		# The magerrs of the stabpoints do not change, we update just the real magerrs, in case.
		self.spline.datapoints.magerrs[self.spline.datapoints.mask] = lightcurve.magerrs
			
	
	def replacespline(self, newspline):
		"""
		If you want to change the spline of this SplineML object
		I should maybe check that you didn't change stab attribute stuff here ?
		"""
		olddp = self.spline.datapoints
		self.spline = newspline.copy()
		# And now we want to update the splines datapoints, by principle. Even if typically you will not fit this spline
		# anymore, we need at least to be able to evaluate it !
		self.spline.updatedp(olddp, dpmethod="extadj")
		#print olddp.jds[0], self.spline.datapoints.jds[0]
		#print olddp.jds[-1], self.spline.datapoints.jds[-1]
		
	
	def reset(self):
		"""
		Puts all coeffs back to zero, and redistributes the knots uniformly.
		The number of knots does not change of course !
		"""
		self.spline.reset()

	def calcmlmags(self, lightcurve=None):
		"""
		Required by lc (for getmags, applyml, etc...)
		Returns a mags-like vector containing the mags to be added to the lightcurve.
		
		We don't need the lightcurve object !
		Yes, as we do not care about time shifts, returning only mags.
		I leave it in the arguments as it is present in the calcmlmags of the polynomial splines.
		"""
		
		return self.spline.eval(nostab = True) # Of course, we do not want to evaluate the stab points !
		
	
	def smooth(self, lightcurve, n=1000):
		"""
		Returns a dict of points to plot when displaying this microlensing.
		Just for plots.
		
		Warning : we suppose here that lightcurve is the same as was used to build the ml !
		
		Here we directly return stuff that can be plotted, so we *do* care about timeshifts.
		
		n is the number of points you want. the more the smoother.
		"""
		
		#jds = lightcurve.getjds()
		
		#smoothtime = np.linspace(jds[0], jds[-1], n) - self.jdoffset
		#smoothml = self.spline.eval(smoothtime)
		#smoothtime += self.jdoffset
		smoothtime = np.linspace(self.spline.datapoints.jds[0], self.spline.datapoints.jds[-1], n)
		smoothml = self.spline.eval(jds = smoothtime)
		refmag = np.median(lightcurve.getmags())
		
		# We also return the knots, in the form of points upon the smooth curve
		# These can then be plotted with errorbars, for instance.
		
		#knotjds = self.spline.getintt()
		knotjds = self.spline.getinttex()
		knotmags = self.spline.eval(jds = knotjds)
		
		# Setting the correct offset, so that the microlensing is shown along with the lightcurve :
		jdref = lightcurve.getjds()[0]
		#jdref = self.jdoffset # NO, this is wrong : it might well be that the lightcurve was shifted in time since the microlensing was created !
		smoothtime += jdref
		knotjds += jdref
		
		return {"n":n, "jds":smoothtime, "ml":smoothml, "refmag":refmag, "knotjds":knotjds, "knotmags":knotmags}
		
	

	# YES, there is place for setparams and getparams here !!!
	# This allows to use spline ML with dispersion techniques.
	# Hmm, I'm directly accessing the spline, easier.


def addtolc(lc, targetlc=None, n = 5, knotstep=None, stab=True, stabgap=30.0, stabstep=3.0, stabmagerr=1.0,
	bokeps = 10.0, boktests = 10, bokwindow = None):
	"""
	Adds a SplineML to the lightcurve.
	SplineML splines have NO external stabilization points (stabext = 0.0) !
	We ignore any mask of the lc, cut it before putting this ML.
	This is just about putting things in place. Of course we cannot optimize anything here !

	If targetlc is not None, then pass me another light curve that already have a microlensing spline. I will "copy" that microlensing to your lc and adjust the knots coefficients accordingly.

	The stab stuff inserts stab points into gaps only.


	We use uniknots, n is the number of uniform intervals you want.
	Or specify knotstep : if this is specified, I don't use n.
	
	::
	
		pycs.gen.splml.addtolc(l, knotstep=200)

	"""
	lcjds = lc.getjds()
	jdoffset = lcjds[0] # We can do this, as jds is sorted.
	lcjds -= jdoffset # so the first true datapoint of a ML spline is at 0.0

	if targetlc == None:
		dpmags = np.zeros(len(lcjds))
		dpmagerrs = np.ones(len(lcjds))

	else:
		targetml = targetlc.ml
		assert targetml.mltype == "spline"
		dpmags = targetml.spline.eval(lcjds)
		dpmags -= np.mean(dpmags)
		dpmagerrs = np.ones(len(lcjds)) * np.median(targetml.spline.datapoints.magerrs) # This is a bit arbitrary...but seems we don't really care

	dp = pycs.gen.spl.DataPoints(lcjds, dpmags, dpmagerrs, splitup=True, sort=True,
		stab=stab, stabext=0.0, stabgap=stabgap, stabstep=stabstep, stabmagerr=stabmagerr)
	
	s = pycs.gen.spl.Spline(dp, bokeps=bokeps, boktests=boktests, bokwindow=bokwindow)
	#s.display()
	
	if knotstep == None:
		s.uniknots(n, n=True)
	else:
		s.uniknots(knotstep, n=False)
	#s.maxnspace(n=1)
	
	# And we add the SplineML to our curve :
	#lc.addml(SplineML(s, stab=stab, stabgap=stabgap, stabstep=stabstep, stabmagerr=stabmagerr))
	lc.addml(SplineML(s))

	# Update the knots coeffs if you copied the ml from targetlc
	if targetlc != None:
		lc.ml.spline.optc()

	
# 	if knottype == "uni":
# 		spline.uniknots(nint)
# 	elif knottype == "equi":
# 		spline.equiknots(nint)
# 	else:
# 		raise RuntimeError("I don't knot this knowtype ;-)")
# 	
# 	return SplineML(spline, jdoffset)


	
# 	datapoints = pycs.gen.spl.DataPoints(jds-jdoffset, np.zeros(len(jds)), np.ones(len(jds)), splitup=False, sort=False)
# 	spline = pycs.gen.spl.Spline(datapoints)
	

	
# def factory(lc, nint=5, knottype="uni"):
# 	"""
# 	Takes a lightcurve, returns a SplineML that you can add to your lightcurve.
# 	
# 	For all this spline microlensing stuff, we use the shifted jds as seen from the outside,
# 	but we allways rescale them so that the first point is on 0.0 in the inside.
# 	This way the microlensing can be conserved if you shift a curve.
# 	But if the lc gets shifted, do not use this offset to get back the correct jds !
# 	
# 	"""
# 	
# 	jds = lc.getjds()
# 	jdoffset = jds[0] # We can do this, as jds is sorted.
# 	
# 	datapoints = pycs.gen.spl.DataPoints(jds-jdoffset, np.zeros(len(jds)), np.ones(len(jds)), splitup=False, sort=False)
# 	spline = pycs.gen.spl.Spline(datapoints)
# 	
# 	if knottype == "uni":
# 		spline.uniknots(nint)
# 	elif knottype == "equi":
# 		spline.equiknots(nint)
# 	else:
# 		raise RuntimeError("I don't knot this knowtype ;-)")
# 	
# 	return SplineML(spline, jdoffset)
	

	
	
	
