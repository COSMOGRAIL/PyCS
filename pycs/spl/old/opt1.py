"""
First generation optimization functions, made here for RXJ1131.
A priori these are not up to date anmyore.

Something to remember about the philosophy of these optimizations : the delays are all we want.
For given delays, at any time, it is easy to build a microlensing or even find fluxshifts.
So it seems ok in a first step to iterate on this, keeping only the delays after each iteration.

The functions must return a final optimal spline. This spline object also contains its final r2, saved into the pkl, and used for plots etc.

"""

import sys
sys.path.append("../")			

import pycs.gen.util
import pycs.gen.lc
import pycs.gen.splml
import pycs.gen.spl
import pycs.spl.multiopt



def rough(lcs, movefirst=True, mln=20, knotstep=20.0, iterations=3, verbose=False):
	"""
	rough = no BOK
	The first l in lcs is used as a reference for some fits.
	I assume that the timeshifts are already quite ok.
	
	At each new iteration, the time shifts are the only thing kept from the previous iteration.
	The reste is redone from scratch.
	Hence you can reduce the number of iterations if you are confident about your initial delays.
	
	There is no risk here to to extra iterations, no divergence or so.
	"""
	
	for it in range(iterations):
		
		pycs.spl.multiopt.opt_magshift(lcs)
		
		# We add some rough microlensing to all but the first curve
		lcs[0].rmml()
		for l in lcs[1:]:
			l.rmml()
			pycs.gen.splml.addtolc(l, mln)
		
		# Spline on the first curve
		roughsinglespline = pycs.gen.spl.fit([lcs[0]], knotstep = knotstep, stab=True, bokit = 0, verbose=verbose)
		#print pycs.gen.spl.r2([lcs[0]], roughsinglespline)
		
		pycs.spl.multiopt.opt_ml(lcs, roughsinglespline, bokit=0, verbose=verbose)
		
		# Rough spline, all
		roughspline = pycs.gen.spl.fit(lcs, knotstep = knotstep, stab=True, bokit = 0, verbose=verbose)
		
		# Quick time shift optimization
		r2 = pycs.spl.multiopt.opt_ts_brute(lcs, roughspline, optml=True, r=2, step=3.0, movefirst=movefirst, verbose = verbose)
		r2 = pycs.spl.multiopt.opt_ts_powell(lcs, roughspline, optml=True, movefirst=movefirst, verbose = verbose)

		print "Done with rough iteration %i/%i, r2 = %.2f" % (it+1, iterations, r2)
		
		if verbose:
			print pycs.gen.lc.multigetnicetimedelays(lcs)
		
	return roughspline
	
	

def fine(lcs, movefirst=True, mln=20, knotstep=15.0, mlbokit=2, splbokit=2, iterations=3, verbose=False):
	"""
	Similar to rough, but with some BOK.
	Also we directly start with a spline fit on all curve, so you should already have a good ML -> run
	rough first !
	"""
	
	for it in range(iterations):
		
		# We build a better spline, with BOK to get better knots
		finespline = pycs.gen.spl.fit(lcs, knotstep = knotstep, stab=True, bokit = splbokit, verbose=verbose)

		lcs[0].rmml()
		for l in lcs[1:]:
			l.rmml()
			pycs.gen.splml.addtolc(l, mln)
		
		# BOK on the ML, starting from scratch, to get better knots
		pycs.spl.multiopt.opt_ml(lcs, finespline, bokit = mlbokit, verbose=verbose)
	
		# And do a opt_ts
		r2 = pycs.spl.multiopt.opt_ts_brute(lcs, finespline, optml=True, r=2, step=1.0, movefirst=movefirst, verbose = verbose)
		r2 = pycs.spl.multiopt.opt_ts_powell(lcs, finespline, optml=True, movefirst=movefirst, verbose = verbose)
	
		print "Done with fine iteration %i/%i, r2 = %.2f" % (it+1, iterations, r2)	
		if verbose:
			print pycs.gen.lc.multigetnicetimedelays(lcs)
	
	return finespline
	
	

def combi(lcs, movefirst=True, mln=20, knotstep=15, mlbokit=2, splbokit=2, rit = 2, fit = 2, verbose=False):


	roughspline = rough(lcs, movefirst=movefirst, mln=mln, knotstep=knotstep, iterations=rit, verbose=verbose)
	finespline = fine(lcs, movefirst=movefirst, mln=mln, knotstep=knotstep, mlbokit=mlbokit, splbokit=splbokit, iterations=fit, verbose=verbose)
	return finespline


