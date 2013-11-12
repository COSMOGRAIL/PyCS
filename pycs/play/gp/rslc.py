"""
Defines a class that represents a regularly sampled lightcurve

"""
import sys
import numpy as np
import pymcgp
import pycs.gen.spl
import copy as pythoncopy
import scipy.optimize as spopt


import matplotlib.pyplot as plt

class rslc():
	"""
	A regularly sampled lightcurve, typically obtained by regression.
	To make such a rslc from a usual lightcurve object, look at the factory function below.
	One idea is that we want to be able to add and subtract those, propagating errors.
	There is no "microlensing" or similar stuff -- only time shifts.
	"""

	def __init__(self, jds, mags, magerrs, pad, pd, timeshift=0.0, name="Name", plotcolour="black"):
		
		self.jds = jds
		self.mags = mags
		self.magerrs = magerrs
		self.plotcolour = plotcolour
		self.name = name
		self.timeshift = timeshift
		self.pad = pad
		self.pd = pd
	
	def __str__(self):
		retstr = "[RS:%s]" % (self.name)
		if self.timeshift != 0.0:
			retstr += "(%.3f)" % (self.timeshift)
		return retstr

	
	def shifttime(self, timeshift):
		self.timeshift += timeshift
	
	def copy(self):
		return pythoncopy.deepcopy(self)

	def getjds(self):
		return self.jds + self.timeshift


	def mask(self, maxmagerr = 0.1, target = 20.0):
		self.magerrs[self.magerrs > maxmagerr] = target 
		
		
	def wtv(self, method = "weights"):
		"""
		Return some weighted total variation.
		Usuall called on a "difference" lightcurve.
		"""
		
		#return np.sum(np.fabs(self.mags[1:] - self.mags[:-1]))
		
		#mask = self.magerrs < maxmagerr
		
		
		if method == "weights":
			dys = self.mags[1:] - self.mags[:-1]
			dyws = 1.0 / (0.5*(self.magerrs[1:] + self.magerrs[:-1]))
		
			out = np.sum(np.fabs(dys) * dyws) / np.sum(dyws)
		
		if method == "simple":
			out = np.sum(np.fabs(self.mags[1:] - self.mags[:-1]))
		
		return out
		
		

def factory(l, pad=300, pd=2, plotcolour=None):
	"""
	Give me a lightcurve, I return a regularly sampled light curve, by performing some regression.
	
	:param pad: the padding, in days
	:param pd: the point density, in points per days.
	
	The points live on a regular grid in julian days, 0.0, 0.1, 0.2, 0.3 ...
	"""
	
	if plotcolour == None:
		plotcolour = l.plotcolour
	
	name = l.object
	
	jds = l.jds.copy()
	timeshift = l.timeshift
	
	mags = l.getmags(noml=True)
	magerrs = l.getmagerrs()
	
	minjd = np.round(jds[0] - pad)
	maxjd = np.round(jds[-1] + pad)
		
	npts = int(maxjd - minjd)*pd
		
	rsjds = np.linspace(minjd, maxjd, npts) # rs for regularly sampled
		
	# The regression itself
	
	mean_mag = np.mean(mags)
	def meanprior(query):
		return (0.0 * query + mean_mag)
		
	regfct = pymcgp.regression(jds, mags, magerrs, meanprior)
	(rsmags, rsmagerrs) = regfct(rsjds)
	
	return rslc(rsjds, rsmags, rsmagerrs, pad, pd, timeshift=timeshift, name=name, plotcolour=plotcolour)


def subtract(rs1, rs2):
	"""
	I subtract rs2 from rs1.
	This means I keep the jds and timeshift of rs1, and only change the mags and magerrs,
	interpolating rs2.
	
	The difference has no time shift.
	
	"""
	
	newjds = rs1.getjds()
	
	newmags = rs1.mags.copy()
	newmagerrs = rs1.magerrs.copy()
	newpad = rs1.pad
	newpd = rs1.pd
	newname = "%s(%+.1f)-%s(%+.1f)" % (rs1.name, rs1.timeshift, rs2.name, rs2.timeshift) 
	
	# We interpolate rs2 at the positions of rs1
	newrs2mags = np.interp(rs1.getjds(), rs2.getjds(), rs2.mags, left=np.nan, right=np.nan)
	newrs2magerrs = np.interp(rs1.getjds(), rs2.getjds(), rs2.magerrs, left=np.nan, right=np.nan)

	# These arrays contain NaN at one of there extremas.
	newmags -= newrs2mags
	newmagerrs = np.sqrt(rs1.magerrs*rs1.magerrs + newrs2magerrs*newrs2magerrs)
	
	# The NaN are now propagated in newmags and newmagerrs
	# We cut them :
	
	nanmask = np.isnan(newmags)
	#nnan = np.sum(nanmask)
	#print nnan/newpd
	
	newjds = newjds[nanmask == False]
	newmags = newmags[nanmask == False]
	newmagerrs = newmagerrs[nanmask == False]
	
	return rslc(newjds, newmags, newmagerrs, newpad, newpd, timeshift=0.0, name=newname, plotcolour="black")

	
def wtvdiff(rs1, rs2, method):
	"""
	Returns the wtv of the difference between 2 curves
	"""
	out = subtract(rs1, rs2).wtv(method)
	#print out
	return float(out)

def opt_ts(rslcs, method="weights", verbose=True):
	"""
	I optimize the timeshifts between the rslcs to minimize the wtv between them.
	Note that even if the wtvdiff is only about two curves, we cannot split this into optimizing
	AB AC AD in a row, as this would never calculate BC, and BC is not contained into AB + AC.
	"""
	rslcsc = [rs.copy() for rs in rslcs]
	
	# No need for reverse combis, as wtvdiff is symmetric.
	#couplelist = [couple for couple in [[rs1, rs2] for rs1 in rslcsc for rs2 in rslcsc] if couple[0] != couple[1]]
	
	indexes = np.arange(len(rslcsc))
	indlist = [c for c in [[i1, i2] for i1 in indexes for i2 in indexes] if c[1] > c[0]]
	couplelist = [[rslcsc[i1], rslcsc[i2]] for (i1, i2) in indlist]
	
	inishifts = [rs.timeshift for rs in rslcsc[1:]]

	def errorfct(timeshifts):
		if timeshifts.shape == ():
			timeshifts = np.array([timeshifts])
		for (rs, timeshift) in zip(rslcsc[1:], timeshifts):
			rs.timeshift = timeshift
		
		tvs = np.array([wtvdiff(rs1, rs2, method=method) for (rs1, rs2) in couplelist])
		ret = np.sum(tvs)
		print timeshifts, ret
		return ret
	
	if verbose:
		print "Starting time shift optimization ..."
		print "Initial pars : ", inishifts
		
	minout = spopt.fmin_powell(errorfct, inishifts, xtol=0.001, full_output=1, disp=verbose)
	#minout = spopt.fmin_bfgs(errorfct, inishifts, maxiter=None, full_output=1, disp=verbose, retall=0, callback=None)
	
	popt = minout[0]
	
	r2 = errorfct(popt) # This sets popt, and the optimal ML and source.
	
	if verbose:
		print "Optimal pars : ", popt
	
	
	
	
	
	