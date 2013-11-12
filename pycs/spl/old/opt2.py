"""
Second generation optimizers, with much more iterative optimizations of the delays themselves.

The functions must return a final optimal spline. This spline object also contains its final r2, saved into the pkl, and used for plots etc.

"""

import sys
sys.path.append("../")			

import pycs.gen.util
import pycs.gen.lc
import pycs.gen.splml
import pycs.gen.spl
import pycs.spl.multiopt



def indiopt(lcs, roughit=3, fineit=10, mlstep=100, splstep=20, verbose=True):
	"""
	Progressively increases the complexity of ML
	Fast !
	September 2011
	"""
	
	if verbose:
		print "Starting indiopt on initial delays :"
		print pycs.gen.lc.multigetnicetimedelays(lcs, separator=" | ")
	
	# ====== Part 1 : rough ========
	pycs.spl.multiopt.opt_magshift(lcs)
	
	# Putting in place a first common spline and a rough ML
	spline = pycs.gen.spl.fit([lcs[0]], knotstep = 60, stab=True, bokit = 0, verbose=False)
	lcs[0].rmml()
	for l in lcs[1:]:
		l.rmml()
		pycs.gen.splml.addtolc(l, knotstep=360)
	pycs.spl.multiopt.opt_ml(lcs, spline, bokit=0, splflat=True, verbose=False)
	spline = pycs.gen.spl.fit(lcs, knotstep = 30, stab=True, bokit = 0, verbose=False)
	
	# Iteratively refining it (without any bok) :
	if verbose:
		print "Rough itetations :"
	for it in range(roughit):
		pycs.spl.multiopt.opt_ml(lcs, spline, bokit=0, splflat=True, verbose=False)
		spline = pycs.gen.spl.fit(lcs, knotstep = 30, stab=True, bokit=0, boktests=5, bokeps=5.0, verbose=False)
		pycs.spl.multiopt.opt_ts_indi(lcs, spline, optml=False, method="brute", brutestep=2.0, bruter=10, verbose = False, trace=False)
		pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 0, verbose=False, trace=False)
		thisr2 = spline.lastr2nostab
		if verbose:
			print "%s    (Iteration %2i, r2 = %8.1f)" % (pycs.gen.lc.multigetnicetimedelays(lcs, separator=" | "), it+1, thisr2)	
	
	
	# ====== Part 2 : fine ========
	# We start from scratch, new ML and spline :
	if verbose:
		print "Setting up fine spline :"
	
	lcs[0].rmml()
	for l in lcs[1:]:
		l.rmml()
		pycs.gen.splml.addtolc(l, knotstep=mlstep)
	pycs.spl.multiopt.opt_ml(lcs, spline, bokit=1, bokwindow = 30, boktests=10, splflat=True, verbose=False)
	spline = pycs.gen.spl.fit(lcs, knotstep = splstep, stab=True, bokit=1, boktests=5, bokeps=6.0, verbose=False)
	pycs.spl.multiopt.opt_ts_indi(lcs, spline, optml=True, method="brute", brutestep=0.5, bruter=10, verbose = False, trace=False)
	pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 0, verbose=False, trace=False)
	pycs.spl.multiopt.opt_ts_indi(lcs, spline, optml=True, method="fmin", verbose = False, trace=False)
	pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 0, verbose=False, trace=False)
	if verbose:
		print "%s" % (pycs.gen.lc.multigetnicetimedelays(lcs, separator=" | "))	
	
	# ====== Part 3 : refinement ========
	if verbose:
		print "Fine iterations :"
		
	for it in range(fineit):
		
		pycs.spl.multiopt.opt_ts_indi(lcs, spline, optml=False, method="brute", brutestep=0.5, bruter=2, verbose = False, trace=False)
		pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 0, verbose=False, trace=False)
		pycs.spl.multiopt.opt_ts_indi(lcs, spline, optml=False, method="fmin", verbose = False, trace=False)
		pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 0, verbose=False, trace=False)
		pycs.spl.multiopt.opt_ml(lcs, spline, bokit=1, bokwindow=10.0, boktests=5, splflat=False, verbose=False)
		pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 1, bokwindow=5.0, boktests=5, bokeps=6.0, verbose=False, trace=False)
		
		thisr2 = spline.lastr2nostab
		if verbose:
			print "%s    (Iteration %2i, r2 = %8.1f)" % (pycs.gen.lc.multigetnicetimedelays(lcs, separator=" | "), it+1, thisr2)	

	# ====== Part 4 : stabilizing time shifts ========
	if verbose:
		print "Timeshift stabilization :"

	for it in range(5):
		
		pycs.spl.multiopt.opt_ts_indi(lcs, spline, optml=False, method="fmin", verbose = False, trace=False)
		pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 0, verbose=False, trace=False)
		thisr2 = spline.lastr2nostab
		if verbose:
			print "%s    (Iteration %2i, r2 = %8.1f)" % (pycs.gen.lc.multigetnicetimedelays(lcs, separator=" | "), it+1, thisr2)	


	return spline
	
	
	"""
	# ====== Part 2 : iterative ========
	# At each iteration we throw everything over board and start from scratch !

	if verbose:
		print "Fine iterations :"
		
	for it in range(fineit):
		#lastr2 = spline.lastr2nostab
		
		# Putting in place a finer ML
		for l in lcs:
			l.rmml()
			pycs.gen.splml.addtolc(l, knotstep=180)
		pycs.spl.multiopt.opt_ml(lcs, spline, bokit=1, splflat=True, verbose=False)
		
		# Putting in place a new spline
		spline = pycs.gen.spl.fit(lcs, knotstep = 20, stab=True, bokit=1, boktests=5, bokeps=5.0, verbose=False)
		
		# Finer time shift optimization
		pycs.spl.multiopt.opt_ts_indi(lcs, spline, optml=False, method="brute", brutestep=0.5, bruter=10, verbose = False, trace=False)
		pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 0, verbose=False, trace=False)
		pycs.spl.multiopt.opt_ts_indi(lcs, spline, optml=False, method="fmin", verbose = False, trace=False)
		pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 0, verbose=False, trace=False)
		
		thisr2 = spline.lastr2nostab
		if verbose:
			print "%s    (Iteration %i, r2 = %8.1f)" % (pycs.gen.lc.multigetnicetimedelays(lcs, separator=" | "), it+1, thisr2)	
	"""

