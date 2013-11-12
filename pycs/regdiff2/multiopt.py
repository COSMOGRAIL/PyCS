"""
Optimization using regdiff

Would be great to have an opt_flux in here !

"""
import sys
import numpy as np

import rslc
import pycs.gen.lc



def opt_ts(lcs, method="weights", pd=5, inistep=5, iniradius=5, nit=5, finish=True, verbose=True):
	"""
	Give me lightcurves (with more or less good initial time shifts)
	I run a regression on them, optimize regdiff, and set their delays to the optimal values.
	This one is meant to measure accurate delays.
	
	:param pd: the point density, in points per days.
	"""
	
	if verbose:
		print "Starting regdiff opt_ts, initial time delays :"
		print "%s" % (pycs.gen.lc.getnicetimedelays(lcs, separator=" | "))	
	
	
	rss = [pycs.regdiff2.rslc.factory(l, pad=100, pd=pd) for l in lcs]
	# The time shifts are transfered to these rss, any microlensing is disregarded

	if verbose:
		print "Regressions done."
	
	minwtv = pycs.regdiff2.rslc.opt_ts(rss, method=method, inistep=inistep, iniradius=iniradius, nit=nit, finish=True, verbose=True)
	
	for (l, r) in zip(lcs, rss):
		l.timeshift = r.timeshift
		l.commentlist.append("Timeshift optimized with regdiff.")

	
	if verbose:
		print "Optimization done ! Optimal time delays :"
		print "%s" % (pycs.gen.lc.getnicetimedelays(lcs, separator=" | "))	

	return minwtv



def opt_ts_indi(lcs, method="weights", pd=0.5, inistep=3, nit=5, verbose=True):
	"""
	I optimize the time shifts one by one with respect to the first curve only.
	I do not shift the first curve.
	This is a function to robustly get first guess time delays.
	"""
	
	jdss = [l.getjds() for l in lcs]
	lengths = map(lambda x: np.max(x) - np.min(x), jdss)
	minlength = np.min(lengths)
	radius = minlength / (2.0 * (float(inistep)))
	
	if verbose:
		print "I will search for time shifts in a range of %.1f days (%i steps)." % (minlength, 2*radius+1)
		
	if verbose:
		print "Computing GPRs..."
	rss = [pycs.regdiff2.rslc.factory(l, pd=pd, pad=0.0) for l in lcs]
	# The time shifts are transfered to these rss, any microlensing is disregarded
	if verbose:
		print "Regressions done."
	
	for rs in rss[1:]:
		minwtv = pycs.regdiff2.rslc.opt_ts([rss[0], rs], method=method, inistep=inistep, iniradius=radius, nit=nit, finish=False, verbose=True)
	
	for (l, r) in zip(lcs, rss):
		l.timeshift = r.timeshift
		l.commentlist.append("Timeshift optimized with regdiff.")

	
	if verbose:
		print "Optimization done ! Optimal time delays :"
		print "%s" % (pycs.gen.lc.getnicetimedelays(lcs, separator=" | "))	



	
	
	
	
	
	