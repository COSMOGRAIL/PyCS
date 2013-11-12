"""
Defines a class that represents a gaussian process regression (GPR)

"""
import numpy as np
import pymcgp
import pycs.gen.spl

class GPR():
	"""
	Class representing a gaussian process regression.
	evaluated on a fine grid.
	So it is possible to shift and subtract or add gaussian processes etc,
	doing simple calculations on this grid.
	
	"""

	def __init__(self, lcs, plotcolour=None):
		"""
		lcs is a list of lightcurves, will be merged
		"""
		
		
		#spline = pycs.gen.spl.fit(lcs, knotstep = 15.0, stab=True, stabext=300.0, stabgap=20.0, stabstep=5.0, staberr=0.1, bokit = 2, bokeps = 2.0, boktests = 5, bokwindow=None,  k=3, verbose=True)
		
		#self.spline = spline
		
		dp = pycs.gen.spl.merge(lcs, splitup=False, sort=True, stab=False)
		self.jds = dp.jds
		self.mags = dp.mags
		self.magerrs = dp.magerrs
	
	
		if plotcolour == None:
			if len(lcs) == 1:
				self.plotcolour = lcs[0].plotcolour
			else:
				self.plotcolour = "black"
	
		
		self.pd = 5 # point density, n per day.
		self.pad = 200 # padding, in days.
		
		
		self.regression()
		
	
	def regression(self):
	
		# Prior on mean function :
		mean_mag = np.mean(self.mags)
		def meanprior(query):
			return (0.0 * query + mean_mag)
		
		# Ok we try this, but with a spline
		#def meanprior(jds):
		#	inshape = jds.shape
		#	return spline.eval(jds=jds.reshape(inshape[0])).reshape(inshape)
	
		self.regfct = pymcgp.regression(self.jds, self.mags, self.magerrs, meanprior)
		
		#npts = (gpr.jds[-1] - gpr.jds[0])*2.0
		#		xs = np.linspace(gpr.jds[0], gpr.jds[-1], npts)
		#		(ys, yerrs) = gpr.regfct(xs)