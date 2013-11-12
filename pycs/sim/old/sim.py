"""
Defines some functions that return simulated lightcurves, for testing purposes
"""

import lc
from numpy import *

def sinuslike(seed=None, object="sim", phase=0.0, amp=0.5, period=450.0, length=1000.0, jdsampling=8.0, jdsamplingstd=4.0, mag=-12.0, magstdmean=0.07, magstdstd=0.03, seasonlength=250.0, seasonlengthstd=90.0):
	"""
	First simple lightcurve simulation to explore curve shifting methods.
	All timelike parameters (phase, period, length, jdsampling, jdsamplingstd, seasongap, etc) are expressed in days.
	
	@type seed: int or None
	@param seed: if None, the clock is used, if not, the given seed.
	
	@rtype: lightcurve
	@return: a simulated ligthcurve.
	"""
	simlc = lc.lightcurve()	# we start with a default one
	
	# create the random seed
	rs = random.RandomState(seed) 
	simlc.jds = arange(0.0, length, jdsampling) # make the base jd array 
	simlc.jds = simlc.jds + jdsamplingstd*rs.randn(len(simlc.jds)) # randomize the sampling
	simlc.jds.sort() # and ensure that we are chronological
	
	seasonmask = (simlc.jds + seasonlengthstd*rs.rand(len(simlc.jds))) % 365.0 < seasonlength # simulate the seasons
	#seasonmask = simlc.jds % 365.0 < seasonlength 
	simlc.jds = simlc.jds[seasonmask] # cut out the good points
	
	simlc.mags = mag + amp * sin(2.0*pi/period*(simlc.jds + phase)) # create the magnitudes
	simlc.magerrs = magstdmean*ones(len(simlc.jds)) + magstdstd*rs.randn(len(simlc.jds))
	simlc.magerrs = abs(simlc.magerrs)
	
	simlc.mask = simlc.magerrs >= 0.0 # That's true for all of them
	
	simlc.clearlabels()

	simlc.object = object
	simlc.telescopename = "sim.sinuslike"
	simlc = simlc.bootstrap(seed=seed) # we pass the seed to bootstrap method
	simlc.commentlist = ["phase=%f"%phase]
	simlc.ploterrorbars = True
	
	simlc.setorigs()
	
	simlc.validate()
	return simlc
