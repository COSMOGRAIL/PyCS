"""
Top level functions that use the src module to tweak microlensing and source splines.
These are the function to pass them to draw.draw or draw.multidraw.
"""

import pycs.sim.src
from numpy import arange, ones, array
#import pycs.gen.splml


def tweakml(lcs,spline, beta=-2.0, sigma=0.05, fmin=1/500.0, fmax=None, psplot=False, sampling=0.1):
	"""
	I tweak the SplineML of your curves by adding small scale structure.
	I DO modify your lcs inplace.
	
	
	"""
	for l in lcs:
		if l.ml == None:
			print "WARNING, curve %s has no ML to tweak !" % (str(l))
			continue
			# No ! Then we just add some flat spline
			#pycs.gen.splml.addtolc(l) # knotstep has no imortantce
			
		elif l.ml.mltype != "spline":
			print "WARNING, I can only tweak SplineML objects, curve %s has something else !" % (str(l))
			continue
		
		spline = l.ml.spline.copy()
		name = "ML(%s)" % (l.object)
		source = pycs.sim.src.Source(spline, name=name, sampling = sampling)
		
		if psplot==True:
			
			psspline = pycs.sim.src.PS(source, flux=False)
			psspline.plotcolour = "black"
			psspline.calcslope(fmin = 1/1000.0, fmax = 1/100.0)
			
		
		source.addplaw2(beta=beta, sigma=sigma, fmin=fmin, fmax=fmax, flux=False, seed=None)
		source.name += "_twk"
		newspline = source.spline()
		l.ml.replacespline(newspline)
		
		if psplot==True:
			
			psnewspline = pycs.sim.src.PS(source, flux=False)
			psnewspline.plotcolour = "red"
			psnewspline.calcslope(fmin = fmin, fmax = fmax)
			
			pycs.sim.src.psplot([psspline, psnewspline], nbins=50)



def tweakspl(spline, beta=-2.5, sigma=0.03, fmin=1/30.0, fmax=1/5.0, hann=False, psplot=False):
	"""
	Give me a spline, I return a tweaked version with added small scale structure.
	Note that the spline I return will have a LOT of knots.
	
	I DO NOT modify your spline, but return a new one.
	"""
	
	source = pycs.sim.src.Source(spline, sampling = 0.2) # No need to pass a copy, spline will not be modified, only evaluated.
		
	if psplot==True:
			
		psspline = pycs.sim.src.PS(source, flux=False)
		psspline.plotcolour = "black"
		psspline.calcslope(fmin = 1/1000.0, fmax = 1/30.0)
			
	source.addplaw2(beta=beta, sigma=sigma, fmin=fmin, fmax=fmax, hann=hann, flux=False, seed=None)
	source.name += "_twk"
	newspline = source.spline()
		
	if psplot==True:
			
		psnewspline = pycs.sim.src.PS(source, flux=False)
		psnewspline.plotcolour = "red"
		psnewspline.calcslope(fmin = fmin, fmax = fmax)

		pycs.sim.src.psplot([psspline, psnewspline], nbins=50)
	
	return newspline
	

def addspl(spl1, spl2, op='add'):
	"""
	Give me two splines, I return the addition/substraction of the two: spl1 +/- spl2

	Important: I assume that the two splines have the correct jds ! The resulting spline
	
	:param spl1: first spline 
	:param spl2: second spline
	:param op: 'add' for addition, 'sub' for subtraction
	:return: spline
	"""

	print "I'm NOT WORKING :("
	return

	src1 = pycs.sim.src.Source(spl1, sampling = 0.2)
	src2 = pycs.sim.src.Source(spl2, sampling = 0.2)

	newspline = spl1.copy() # deep copy, on which we will tweak the datapoints

	# defines the range on which the addition/subtraction happens
	minjd = max(src1.jdmin, src2.jdmin)
	maxjd = min(src1.jdmax, src2.jdmax)

	subjds = newspline.datapoints.jds
	subjdsbf = [jd for jd in subjds if jd<minjd]
	subjdsin = [jd for jd in subjds if jd>=minjd and jd <=maxjd]
	subjdsaf = [jd for jd in subjds if jd>maxjd]

	#subjds = arange(minjd, maxjd, step=0.2)
	# IMPORTANT ! I have to make sure that the two spline "overlap" in magnitude at their extrema points...and that the transition is smooth...?


	mags1 = src1.eval(subjdsin)
	mags2 = src2.eval(subjdsin)
	submagerrs = 0.001 * ones(subjds.shape) # default value, should not affect the following

	if op == 'add':
		submagsin = mags1 + mags2
	elif op == "sub":
		submagsin = mags1 - mags2

	submags = []
	for mags in [src1.eval(subjdsbf), submagsin, src1.eval(subjdsaf)]:
		for m in mags:
			submags.append(m)


	submags = array(submags)

	# Build a new spline from this new jds and mags
	newdatapoints = pycs.gen.spl.DataPoints(jds=subjds, mags=submags, magerrs=submagerrs)
	newspline = pycs.gen.spl.Spline(datapoints=newdatapoints, bokeps=0.2)
	#newspline.updatedp(newdatapoints=newdatapoints, dpmethod="stretch") # Nope

	newspline.uniknots(nint=0.2, n=False)
	#newspline.buildbounds()
	#newspline.optc()
	#newspline.bok()

	#spl1.display()
	newspline.display()
	pycs.gen.lc.display([], [spl1, spl2, newspline])

	print "allok"
	import sys
	sys.exit()
	

