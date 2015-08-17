"""
Top level functions that use the src module to tweak microlensing and source splines.
These are the function to pass them to draw.draw or draw.multidraw.
"""

import pycs.sim.src
#import pycs.gen.splml


def tweakml(lcs, beta=-2.0, sigma=0.05, fmin=1/500.0, fmax=None, psplot=False, sampling=0.1):
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
	

	
	
	

