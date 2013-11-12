"""
Third generation
The functions must return a final optimal spline. This spline object also contains its final r2, saved into the pkl, and used for plots etc.


.. todo:: [Solved ! See opt4 ... ] Think about ways to keep the ML somehow, not very nice to redefine it inside of these functions.
	Evenutally split up this indiopt into 3 calls of one same function that keeps the ML ?
	This would allow to have different MLs for each curve ! (important !)
	But on the other hand we run this on drawn curves (that can be different -> different ML)
	Ok we would need to copy the input lcs onto the drawn curves.
	
	The input of multirun needs a good ML to draw !
	So at the end better if the optfct defines its own ML.

"""

import sys
sys.path.append("../")			

import pycs.gen.util
import pycs.gen.lc
import pycs.gen.splml
import pycs.gen.spl
import pycs.spl.multiopt



def indiopt(lcs, spline=None, rmlstep = 200, medit=10, fineit=10, mlstep=120, splstep=20, verbose=True, trace=False, tracedir = "opt3"):
	"""
	Progressively increases the complexity of ML
	September 2011
	
	If spline =! None, I keep the current ML and spline, skipping the rough part.
	
	"""
	
	if verbose:
		print "Starting indiopt on initial delays :"
		print pycs.gen.lc.getnicetimedelays(lcs, separator=" | ")
	
	
	# ====== Part 1 : rough ========
	# No BOK, splflat = True, no iterative tweaking of the source
	
	if spline == None:
		if verbose:
			print "Setting up rough splines ..."
		
		pycs.spl.multiopt.opt_magshift(lcs)
	
		# Putting in place a first common spline and a rough ML
		spline = pycs.gen.spl.fit([lcs[0]], knotstep = 60, stab=True, bokit = 0, verbose=False)
		lcs[0].rmml()
		for l in lcs[1:]:
			l.rmml()
			pycs.gen.splml.addtolc(l, knotstep=rmlstep)
		pycs.spl.multiopt.opt_ml(lcs, spline, bokit=0, splflat=True, verbose=False)
		spline = pycs.gen.spl.fit(lcs, knotstep = 30, stab=True, bokit = 0, verbose=False)
	
		# Iteratively refining it (without any bok, starting from new splines) :
		if verbose:
			print "Rough iterations :"
		for it in range(3):
			if trace:
				pycs.gen.util.trace(lcs, [spline], tracedir = tracedir)
			pycs.spl.multiopt.opt_ml(lcs, spline, bokit=0, splflat=True, verbose=False)
			spline = pycs.gen.spl.fit(lcs, knotstep = 30, stab=True, bokit=0, boktests=5, bokeps=5.0, verbose=False)
			pycs.spl.multiopt.opt_ts_indi(lcs, spline, optml=False, method="brute", brutestep=2.0, bruter=10, verbose = False, trace=False)
			pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 0, verbose=False, trace=False)
			if verbose:
				print "%s    (Iteration %2i, r2 = %8.1f)" % (pycs.gen.lc.getnicetimedelays(lcs, separator=" | "), it+1, spline.lastr2nostab)	
			
	else:
		if verbose:
			print "Reusing given spline, skipping rough estimation."
	
	
	# ====== Part 2 : med ========
	
	# We start from scratch, new ML and spline.
	# splflat = True for all ML optimizations.
	
	if medit != 0:
		if verbose:
			print "Setting up med splines :"
		
		lcs[0].rmml()
		for l in lcs[1:]:
			l.rmml()
			pycs.gen.splml.addtolc(l, knotstep=rmlstep)
		pycs.spl.multiopt.opt_ml(lcs, spline, bokit=2, bokwindow = 30, boktests=10, splflat=True, verbose=False)
		spline = pycs.gen.spl.fit(lcs, knotstep = 20, stab=True, bokit=2, boktests=5, bokeps=6.0, verbose=False)
		pycs.spl.multiopt.opt_ts_indi(lcs, spline, optml=True, method="brute", brutestep=0.5, bruter=10, verbose = False, trace=False)
		pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 0, verbose=False, trace=False)
		pycs.spl.multiopt.opt_ts_indi(lcs, spline, optml=True, method="fmin", verbose = False, trace=False)
		pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 0, verbose=False, trace=False)
		if verbose:
			print "%s" % (pycs.gen.lc.getnicetimedelays(lcs, separator=" | "))	
		
	
		if verbose:
			print "Med iterations :"
			
		for it in range(medit):
			if trace:
				pycs.gen.util.trace(lcs, [spline], tracedir = tracedir)
			pycs.spl.multiopt.opt_ts_indi(lcs, spline, optml=False, method="brute", brutestep=0.5, bruter=2, verbose = False, trace=False)
			pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 0, verbose=False, trace=False)
			pycs.spl.multiopt.opt_ts_indi(lcs, spline, optml=False, method="fmin", verbose = False, trace=False)
			pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 0, verbose=False, trace=False)
			pycs.spl.multiopt.opt_ml(lcs, spline, bokit=1, bokwindow=10.0, boktests=5, splflat=True, verbose=False)
			pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 1, bokwindow=5.0, boktests=5, bokeps=6.0, verbose=False, trace=False)
			if verbose:
				print "%s    (Iteration %2i, r2 = %8.1f)" % (pycs.gen.lc.getnicetimedelays(lcs, separator=" | "), it+1, spline.lastr2nostab)	
	
	else:
		if verbose:
			print "Skipping med splines."
	
	
	# ====== Part 3 : fine ========
	# We start from scratch, new ML and spline :
	
	
	if fineit != 0:
		if verbose:
			print "Setting up fine splines :"
	
		for l in lcs:
			l.rmml()
			pycs.gen.splml.addtolc(l, knotstep=mlstep)
			
		pycs.spl.multiopt.opt_ml(lcs, spline, bokit=2, bokwindow = 30, boktests=10, splflat=True, verbose=False)
		spline = pycs.gen.spl.fit(lcs, knotstep = splstep, stab=True, bokit=2, boktests=5, bokeps=6.0, verbose=False)
		if verbose:
			print "%s" % (pycs.gen.lc.getnicetimedelays(lcs, separator=" | "))	
		
		if verbose:
			print "Fine iterations :"
			
		for it in range(fineit):
			if trace:
				pycs.gen.util.trace(lcs, [spline], tracedir = tracedir)
			pycs.spl.multiopt.opt_ts_indi(lcs, spline, optml=False, method="brute", brutestep=0.5, bruter=2, verbose = False, trace=False)
			pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 0, verbose=False, trace=False)
			pycs.spl.multiopt.opt_ts_indi(lcs, spline, optml=False, method="fmin", verbose = False, trace=False)
			pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 0, verbose=False, trace=False)
			pycs.spl.multiopt.opt_ml(lcs, spline, bokit=1, bokwindow=10.0, boktests=10, splflat=True, verbose=False)
			pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 1, bokwindow=5.0, boktests=5, bokeps=6.0, verbose=False, trace=False)
			
			if verbose:
				print "%s    (Iteration %2i, r2 = %8.1f)" % (pycs.gen.lc.getnicetimedelays(lcs, separator=" | "), it+1, spline.lastr2nostab)	
	else:
		if verbose:
			print "Skipping fine splines."
	
	# ====== Part 4 : stabilizing time shifts ========
	if verbose:
		print "Timeshift stabilization :"

	for it in range(5):
		
		pycs.spl.multiopt.opt_ts_indi(lcs, spline, optml=False, method="fmin", verbose = False, trace=False)
		pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 0, verbose=False, trace=False)
		if verbose:
			print "%s    (Iteration %2i, r2 = %8.1f)" % (pycs.gen.lc.getnicetimedelays(lcs, separator=" | "), it+1, spline.lastr2nostab)	
	"""
	if verbose:
		print "Setting border ML :"
		# Note that we don't do any timeshifts after this ...
	pycs.spl.multiopt.opt_ml(lcs, spline, bokit=0, splflat=False, verbose=False)
	pycs.spl.multiopt.opt_source(lcs, spline, method="extadj", bokit = 0, verbose=False, trace=False)
	if verbose:
			print "%s    (r2 = %8.1f)" % (pycs.gen.lc.getnicetimedelays(lcs, separator=" | "), spline.lastr2nostab)	
	"""
	
	if trace:
		pycs.gen.util.trace(lcs, [spline], tracedir = tracedir)
	return spline
	



