"""
Optimization using regdiff

Would be great to have an opt_flux in here !

"""
import sys
import numpy as np

import rslc
import pycs.gen.lc
import pycs.spl.multiopt


def opt_ts(lcs, method="weights", pd=5,plotcolour=None,knotstep=20.0, n=None, stab=True,stabext=300.0, stabgap=20.0, stabstep=5.0,
		stabmagerr=-2.0, stabrampsize=0, stabrampfact=1.0, bokit=1, bokeps=2.0, boktests=5, bokwindow=None, k=3, verbose=True, magshift=True):
	"""
	Give me lightcurves (with more or less good initial time shifts)
	I run a spline regression on them, optimize the difference, and set their delays to the optimal values.
	I also optimise the magnitude shifts for display purpose, but it does not affect the time shifts at all !
	
	:param pd: the point density of the regularly sampled lightcurves, in points per days.
		"""
	
	if verbose:
		print "Starting spldiff opt_ts, initial time delays :"
		print "%s" % (pycs.gen.lc.getnicetimedelays(lcs, separator=" | "))	
	
	
	
	rss = [pycs.spldiff.rslc.factory(l, pd=pd,plotcolour=plotcolour,knotstep=knotstep, n=n, stab=stab, stabext=stabext, stabgap=stabgap, stabstep=stabstep,
						stabmagerr=stabmagerr, stabrampsize=stabrampsize, stabrampfact=stabrampfact, bokit=bokit, bokeps=bokeps, boktests=boktests,
						bokwindow=bokwindow, k=k, verbose=verbose)[0] for l in lcs]
	if magshift:						
		spline = pycs.spldiff.rslc.factory(lcs[0], pd=pd,plotcolour=plotcolour,knotstep=knotstep, n=n, stab=stab, stabext=stabext, stabgap=stabgap, stabstep=stabstep,
						stabmagerr=stabmagerr, stabrampsize=stabrampsize, stabrampfact=stabrampfact, bokit=bokit, bokeps=bokeps, boktests=boktests,
						bokwindow=bokwindow, k=k, verbose=verbose)[1]				
						
	# The time shifts are transfered to these rss, any microlensing is disregarded

	if verbose:
		print "Regressions done."
	

	minwtv = pycs.spldiff.rslc.opt_ts(rss, method=method, verbose=verbose)

	
	for (l, r) in zip(lcs, rss):
		l.timeshift = r.timeshift
		l.commentlist.append("Timeshift optimized with spldiff.")
		
	if magshift:
		if verbose:
			print "Now, I will optimise the shifts in magnitude"
			print "I will shift your lightcurves so that they match a similar spline"
		pycs.spl.multiopt.opt_magshift(lcs,spline)		

	
	if verbose:
		print "Optimization done ! Optimal time delays :"
		print "%s" % (pycs.gen.lc.getnicetimedelays(lcs, separator=" | "))	

	

	return minwtv
	
	
	
	
	
	
	
