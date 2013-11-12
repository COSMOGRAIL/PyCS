
.. _matchtels:

Matching telescopes
===================

Lightcurves from different telescopes usually need to be empirically "matched" in terms of (at least) magnitude offset. Even if you try to perform relative photometry with respect to the same reference stars, a small shift in magnitude is typically needed, as the detectors + filters of the telescopes are different. Often a shift in **flux** is also needed.

A note about flux shifts
------------------------


In pycs, lightcurves are encoded in terms of magnitudes, and thus also shown in this way on plots. Thus, multiplicative macrolensing amplification (i.e., the flux ratio between qso images), but also stationary microlensing amplifications correspond to a simple "vertical" shift of the lightcurves.

But aside of the attribute magshift, lightcurve objects also have an attribute fluxshift. **Note that a "shift" in flux does not correspond to a vertical shift of the curve plotted on a magnitude scale !** A shift in flux *deforms* the curves. That's why the flux shift represents a welcome parameter when mathing lightcurves from different telescopes. The physical origin of a shift in flux can be any additive contamination of the QSO image (like a constant foreground star or some structure of the lens galaxy). As these contaminations do not necessarily have the same colour as the QSO images, they can yield flux shifts not only between the QSO images, but also between curves of the same QSO image observed with different telescopes !

If your lightcurves comes from cosmouline/lcmanip, the "magnitudes" are simply related to the instrumental flux (in electrons) via :math:`m = -2.5 \log(f)`, i.e. :math:`f = 10^{(-0.4 \cdot m)}` : there is no zeropoint. Thus you can easily calculate yourself the fluxes of your curves.

Play a bit around with these fluxshifts (take a red curve, make a copy of it and set its plotcolour to blue, then set some fluxshifts, and plot both ...).


.. note:: Keep an eye on the values of the flux shifts, when optimizing them with the functions described below. There is probably no physical reason for a fluxshift to be larger than the actual measured flux of the source ...


Matching lightcurves from different telescopes
----------------------------------------------

Here we describe how to empirically adjust the magshift and fluxshift of different curves of the *same* QSO image, so that they match.

There's a module for this : :py:mod:`pycs.gen.mrg`.

Of course, this makes only sense if your curves show some overlap !

First of all, we need to define a *matching criterion*. A simple yet effective choice is a disperion measure.

::

	rawdispersionmethod = lambda lc1, lc2 : pycs.disp.disps.linintnp(lc1, lc2, interpdist = 30.0)
	dispersionmethod = lambda lc1, lc2 : pycs.disp.disps.symmetrize(lc1, lc2, rawdispersionmethod)

``dispersionmethod`` is now a function that takes two lightcurves as only arguments, and returns a value quantifying how well these two curves match.


The function :py:func:`pycs.gen.mrg.matchtels` will now optimize the magnitude shifts and flux shifts of some curves so that the match, given this criterion. It takes 3 mandatory arguments : a list of reference lightcurves (they will not be modified), a corresponding list of lightcurves to be adjusted so that they match to the reference curves, and the dispersionmethod defined above. Note that :py:func:`pycs.gen.mrg.matchtels` will attribute **one same magshift** and (if asked) **different individual fluxshifts** to all the curves ! If you want to adjust an individual maghift for each curve, simply call ``matchtels`` on dedicated lightcurve lists.

::

	# Let's assume we have a quad observed with telescopes Euler and Mercator.
	# This means we probably already have two corresponding lightcurves lists like e.g.
	# eulerlcs = [a_euler, b_euler, c_euler, d_euler]
	# mercatorlcs = [a_mercator, b_mercator, c_mercator, d_mercator]

	pycs.gen.mrg.matchtels(eulerlcs, mercatorlcs, dispersionmethod, fluxshifts=True)

	# Let's see if it worked :
	
	pycs.gen.lc.display(eulerlcs + mercatorlcs, showdates=True, showdelays=False)


Of course you can also set the fluxshifts and magshifts by hand ...
If you are happy with the match, merge the lightcurves. Instead of doing this one by one, there is a function to do this operating directly on lightcurve lists :py:func:`pycs.gen.mrg.merge`

::
	
	lcs = pycs.gen.mrg.merge([eulerlcs, mercatorlcs])
	pycs.gen.mrg.colourise(lcs)
	
	pycs.gen.lc.display(lcs, showdates=True, showdelays=False)
	
If you have another telescope, you can at this point call ``matchtels`` again to match this third telescope on th e merged curves ``lcs`` that you have just obtained.




	
	
