"""
Optimization using regdiff

Would be great to have an opt_flux in here !

"""
import sys
import numpy as np

import rslc
import pycs.gen.lc



def opt_ts(lcs, method="weights", pd=5, pow=1.5, amp=2.0, scale=200.0, errscale=5.0, verbose=True):
	"""
	Give me lightcurves (with more or less good initial time shifts)
	I run a regression on them, optimize regdiff, and set their delays to the optimal values.
	
	:param pd: the point density, in points per days.
	
	The parameters pow, amp, scale, errscale are passed to the GPR, see its doc (or explore their effect on the GPR before runnign this...)
	"""
	
	if verbose:
		print "Starting regdiff opt_ts, initial time delays :"
		print "%s" % (pycs.gen.lc.getnicetimedelays(lcs, separator=" | "))	
	
	
	rss = [pycs.regdiff.rslc.factory(l, pd=pd, pow=pow, amp=amp, scale=scale, errscale=errscale) for l in lcs]
	# The time shifts are transfered to these rss, any microlensing is disregarded

	if verbose:
		print "Regressions done."
	
	minwtv = pycs.regdiff.rslc.opt_ts(rss, method=method, verbose=True)
	
	for (l, r) in zip(lcs, rss):
		l.timeshift = r.timeshift
		l.commentlist.append("Timeshift optimized with regdiff.")

	
	if verbose:
		print "Optimization done ! Optimal time delays :"
		print "%s" % (pycs.gen.lc.getnicetimedelays(lcs, separator=" | "))	

	return minwtv
	
	
	
	
	
	
	