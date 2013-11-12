Tips and Tricks, code snippets
==============================


.. contents::


Which PyCS am I using ?
-----------------------

If you have several copies, from SVN or installed...
::

	import pycs
	print pycs.__file__
	


"Faking" PyCS delay or error measurements
-----------------------------------------

This is useful if you want to include measurments from non-pycs techniques on the typical pycs plots
::

	# To fake a PyCS delay measurement :
	data = [{"label":"AB", "mean":-118.6, "med":0.0, "std":0.0}]

	# To fake a PyCS errorbar :
	"""
	data = [{
	"label":"AB",
	"sys":3.0,
	"ran":4.0,
	"tot":5.0,
	"bias":-3.0
	}]
	"""
	outname = "sims_copies_runresults_rathna_May07_delays.pkl"
	dc = pycs.sim.plot.delaycontainer(data = data, name = "Difference-smoothing technique", plotcolour = "darkorange", objects=["A", "B"])	
	pycs.gen.util.writepickle(dc, outname)



Full example to generate a nice big light curve plot
----------------------------------------------------

This was used for Fig. 4 of the RXJ1131 paper :

::

	# Loading the lists of lightcurves 
	EulerC2 = pycs.gen.util.readpickle("data/EulerC2_matched.pkl")
	ECAM = pycs.gen.util.readpickle("data/ECAM_matched.pkl")
	SMARTS = pycs.gen.util.readpickle("data/SMARTS_matched.pkl")
	Mercator = pycs.gen.util.readpickle("data/Mercator_matched.pkl")
	
	
	# Shifting the curves so to get reasonable mag scales :
	tels = [ECAM, SMARTS, Mercator, EulerC2]
	for tel in tels:
		for l in tel:
			l.shiftmag(14.5)
		tel[1].shiftmag(0.3)
		tel[3].shiftmag(0.3)
	
	
	# Tweaking the jdrange and magrange below so that it looks nice
	xw = 1700
	x1 = 52700
	x2 = 54400
	
	yw = 5
	y1= 0.85
	y2= -0.2
	
	# Text to be written on the plot :
	text = [
	(0.504, 0.89, r"$\mathrm{RX\,J1131-1231}$", {"fontsize":24, "horizontalalignment":"center"}),
	(0.03, 0.8, "A", {"fontsize":18}),
	(0.03, 0.67, "B", {"fontsize":18}),
	(0.03, 0.55, "C", {"fontsize":18}),
	(0.03, 0.20, "D", {"fontsize":18}),
	(0.052, 0.67, "+0.3 mag", {"fontsize":10}),
	(0.052, 0.20, "+0.3 mag", {"fontsize":10}),
	(0.052, 0.47, "Euler C2 : 265 epochs", {"fontsize":10, "color":"red"}),
	(0.052, 0.42, "Euler ECAM : 76 epochs", {"fontsize":10, "color":"#FF7700"}),
	(0.052, 0.37, "SMARTS 1.3-m : 288 epochs", {"fontsize":10, "color":"green"}),
	(0.052, 0.32, "Mercator : 78 epochs", {"fontsize":10, "color":"blue"})
	
	]
	
	pycs.gen.lc.display(ECAM + SMARTS + Mercator + EulerC2,
	nicefont=False, showdelays=False, showlegend=False, showdates=True, errorbarcolour="#777777",
	figsize=(10,4), plotsize=(0.05, 0.97, 0.05, 0.93), jdrange=(x1, x1+xw), magrange=(y1+yw, y1),
	filename="fig_lcpart1.pdf", markersize=3.0, capsize=0, text=text, jdmintickstep=50, magmintickstep=0.2, showgrid=True
	)
	
	pycs.gen.lc.display(ECAM + SMARTS + Mercator + EulerC2,
	nicefont=False, showdelays=False, showlegend=False, showdates=True, errorbarcolour="#777777",
	figsize=(10,4.25), plotsize=(0.05, 0.97, 0.11, 0.94), jdrange=(x2, x2+xw), magrange=(y2+yw, y2),
	filename="fig_lcpart2.pdf", markersize=3.0, capsize=0, jdmintickstep=50, magmintickstep=0.2, showgrid=True
	)
	




Building scrolling plots for long curves
----------------------------------------

::

	# Animated plot for talk :

	startjd = 52900.0
	width = 1000.0
	endjd = 55800.0
	n = 1000
	
	for i in range(n):
		
		a = startjd + i* (endjd - width - startjd)/(n-1)
		b = a + width
		
		filename = "mov/%i.png" % (i)
		pycs.gen.lc.display(lcs, nicefont=True, showdelays=False, showlegend=False, showdates=True, showgrid=True, magrange=(4.3, 0), jdrange=(a, b), filename=filename)
	

And then use ffmpeg (or any other similar tool) to turn this into a movie.

Tweaking magnitudes for individual seasons
------------------------------------------

For a lightcurve ``l``, ``l.mags`` is just a numpy array.
To *lower* the third season by 0.03 mags :
::
	
	seasons = pycs.gen.sea.autofactory(l)
	l.mags[seasons[2].indices] += 0.03
	



Playing with custom properties
------------------------------

You can perfectly create your own properties. It's just a list of dicts ...
::
	
	for i in range(len(l)):
		l.properties[i]["my_new_prop"] = "brocoli"
		
	# To see what properties a curve has :
	print l.longinfo()

"Common" properties are properties that all points of the curve have (this is usually the case). Only those "common" properties can be exported as columns in rdb files, for insance.


Splitting a curve by properties
-------------------------------

::
	
	def splitbyprop(l, prop = "telescope"):
		"""
		kills mask ...
		"""
		
		vals = sorted(list(set([l.properties[i][prop] for i in range(len(l))])))
		
		out = []
		for val in vals:
			lcp = l.copy()
			lcp.mask = np.array([l.properties[i][prop] == val for i in range(len(l))])
			lcp.cutmask()
			lcp.telescopename = val
			out.append(lcp)
			
		#pycs.gen.mrg.colourise(out)
		return out




Correcting for flux sharing
---------------------------

March 2012, only implemented for the spline method. Simple code works well, but quick tests on simulated data (HE2149) show degeneracies.
Need complete tests on simulated data with a little flux sharing, to see if it reduces systematic error.

::

	# draw fake curves :
	flcs = pycs.sim.draw.draw(lcs, spline, shotnoise="none", keepshifts=False)
	pycs.sim.draw.shareflux(flcs[0], flcs[1], frac=0.02)
	pycs.gen.lc.display(flcs)

	# then run pycs.spl.topopt.opt_fine, it has the option "redistribfluxes"
	
	
Plotting a structure function
-----------------------------


April 2012, implemented to see how it looks on residuals. Ugly. But might be intresting for curves ?

::

	(lcs, spline) = pycs.gen.util.readpickle("optspl.pkl")

	rls = pycs.gen.stat.subtract(lcs, spline)

	pycs.gen.lc.display(rls)
	
	pycs.gen.stat.sf(rls[0])



Playing with power spectra
--------------------------

It's what tweakml does :

::
	
	(lcs, spline) = pycs.gen.util.readpickle("optspl.pkl")
	#pycs.gen.lc.display(lcs, [spline])

	spline = l.ml.spline.copy()
	source = pycs.sim.src.Source(spline, sampling = 0.2)

	source.addplaw2(beta=-2.0, sigma=0.01, flux=False, fmin=1/1000.0, fmax=1/10.0, hann=False, seed=None)
	pycs.sim.src.sourceplot([source], filename=None, figsize=(12, 8), showlegend=True, marker=None)	

	newspline = source.spline()
	l.ml.replacespline(newspline)