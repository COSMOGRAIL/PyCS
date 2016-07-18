Curve shifting, finally
=======================


PyCS makes it very easy to compare different point estimators. You can also add your own, without much code integration or inheritance. You would simply have to provide a python function that shifts instances of the lightcurve class.
High level functions are provided to run these point estimators on both real observations and the synthetic curves. This can be done in parallel on several CPUs, again using a very low tech approach.

For tuto :
If only the curve shifting methods where picklable, it would be trivial to use multiprocessing here.

This section describes how to do a first optimization of the time shifts (and of course the microlensing) so that the lightcurves of different QSO images "match".

Getting accurate time delay estimations (including the error analysis) is described in the next section.

We assume here that you have some curves (e.g. two, or four), nicely prepared according to the previous steps of this tutorial. Save these nicely merged curves into a pkl, as described earlier.


The dispersion optimizer
------------------------

*Optimizers* are functions that take a list of lightcurves as only argument, and shift these curves (potentially also adjusting the microlensing) so that they match.

Next to your scripts, make a new file called for instance ``myopt.py`` and fill it with the following content :

::

	import pycs.disp.disps
	import pycs.disp.topopt

	rawdispersionmethod = lambda lc1, lc2 : pycs.disp.disps.linintnp(lc1, lc2, interpdist = 30.0)
	dispersionmethod = lambda lc1, lc2 : pycs.disp.disps.symmetrize(lc1, lc2, rawdispersionmethod)

	def disp(lcs):
		return pycs.disp.topopt.opt_full(lcs, rawdispersionmethod, nit=5, verbose=True)


In this file you will set up all the "optimizers". Indeed some of them need to be adjusted for each individual lens.
From within your scripts, you can now ``import myopt`` and use these functions. This will make things easier later.

The function ``myopt.disp`` that you just defined is such an "optimizer". As you see it relies on existing optimizers, in this case :py:func:`pycs.disp.topopt.opt_full` (where *full* typically means that both time shifts and microlensing will get optimized). In fact your ``disp`` is just a one-line wrapper around this ``opt_full`` in this case, but it is useful to explicitly define it, as this will allow you to adjust some optimizer parameters so that it works best for your particular curve. OK for this simple first optimizer we don't have many parameters, just ``nit``, a number or iterations. But other optimizers might leave much more room for adjustments.

Let's try this optimizer.
First, we will have to equip the curves with some simple microlensing representations (start with simple polynoms across the full curves, to get a minimum amount of parameters).
Then, we shift them : 

::

	
	import pycs.gen.util
	import myopt
	
	lcs = pycs.gen.util.readpickle("merged.pkl") # Some cleaned and merged curves.
	# Let's assume you have a double lens. Thus, lcs is a list of two lightcurve objects.
	
	# We put some microlensing on one of the curves (simple polynoms in this case) :
	pycs.gen.polyml.addtolc(lcs[1], nparams=1, autoseasonsgap = 60.0)

	# And shift them !
	myopt.disp(lcs)
	
	pycs.gen.lc.display(lcs, showdates=True)


While they try to be robust and global, the optimizers are not made to find large delays by themselves, nor to discover delays that are not visible "by eye".
They might fail if the initial timeshifts are too far off. So if you see the delay, simply roughly set some time shifts before calling the optimizers. For instance, 

::
	
	lcs[3].shifttime(-90)


You can also perfectly do this already at the end of your importation / merging script, so that it gets saved into the pickle.

Now that some time shifts are set or optimized, you can also retry this :

::

	print pycs.gen.lc.getnicetimedelays(lcs, separator = " | ")
	


The spline optimizer
--------------------

The idea here is to fit one **single** spline (representing the "intrinsic" variation of the QSO) to all your curves, shifting the latter so that the chi2 between all the datapoints and this single spline is minimal.

This is a very parametric problem, and not a trivial task. The optimizer has to "simultaneously" adjust at least :

* the time/magnitude/flux -shifts of the curves
* the microlensing representation of the curves (polynom coefficients, or spline coeffs and knots)
* the intrinsic spline (both coefficients and knot positions)

In this first approach, we won't describe the internal details. In fact, the spline optimizer works in a very similar way than the dispersion optimizer described above : it shifts your curves, and adjusts their microlensing representations. But of course, it also involves one new object : the spline (that you will want to display on top of your lightcurves) !
Try it :

::

	lcs = pycs.gen.util.readpickle("merged.pkl")
	
	# We might need some microlensing representations for some of the curves :
	pycs.gen.splml.addtolc(lcs[1], knotstep=300) # Yep, let's try spline microlensing !
	
	# And now the optimizer. Note that it returns the spline object !
	spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=150)
	

To show this spline, we make use of :py:func:`pycs.gen.lc.display` that you've met many times before. Indeed, this function can also display an arbitrary number of splines !

As usual, it's first argument is simply a list of lightcurve objects (or an empty list, if you don't want to show any lightcurves at all). But you can also specify as optional second argument a **list of spline objects** (hence the ``[ ]``) :

::

	pycs.gen.lc.display(lcs, [spline])


Don't expect a perfect fit at this point -- this was just a first demo of the principles of the optimizer.
This particular spline is shown in black, with vertical ticks indicating its knot positions.

The spline optimizer seen above takes a few optional arguments. Some words about the arguments seen at this point :

* ``nit`` is a number of iterations, it's fine to leave it at 5 unless you dig into the details of these optimizers.
* ``knotstep`` sets the initial "spacing" (in days) of the knots. These knots will then move around, so the spacing will change... but the number of knots will not !

.. warning:: A *lower* knotstep corresponds to *more* knots !

When adding the spline microlensing to the curve, we specified ``knotstep=300`` : this is the same parameter, but for the knots of the microlensing. So choose a lower microlensing-``knotstep`` to get a more flexible microlensing.

The above example used the "rough" optimizer. This one is not made to get accurate time delays, but to roughly (and quickly) adjust the shifts and the microlensing so that the curves match somehow. Hence, for this *rough* part, leave a relatively high ``knotstep``.

Directly after this rough optimization, add a call to a finer optimizer :

::

	spline = pycs.spl.topopt.opt_fine(lcs, nit = 5, knotstep=100)
	

This optimizer will build a new spline from scratch (and return it), using a (usually finer) ``knotstep`` of your choice. Add this line just after the call to opt_rough, and play with the knotstep (e.g. 50) to see the effect. Also note that the knots are now effectively moving (the opt_rough didn't move them).


.. image:: ../_static/tutorial/spline.png
	:align: center
	:width: 800


It's now a good idea to add these optimizers to your ``myopt.py`` file, directly concatenating them ! This allows you to build a custom optimizer for your particular light curve. Here is an example (you should probably update the knotsteps, depending on the curves you want to process) :

::

	def spl(lcs):
		spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=150)
		for l in lcs:
			l.resetml()
		spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=40)
		spline = pycs.spl.topopt.opt_fine(lcs, nit=10, knotstep=40)

		return spline # Do not forget to return the spline !


You can now use ``myopt.spl(lcs)`` in the same way as ``myopt.disp(lcs)`` (except that myopt.spl(lcs) returns a spline, that you might want to "catch" by writing

::

	spline = myopt.spl(lcs)
	
As usual, after such an optimization, it might be convenient to save the shifted curves and in this case also the spline into a pickle file, so that you can work on them without rerunning the optimization. Tip : save both the curves and the spline into the same pickle file ! 

::
	
	pycs.gen.util.writepickle((lcs, spline), "optspl.pkl")
		
	# ...
		
	(lcs, spline) = pycs.gen.util.readpickle("optspl.pkl")
	pycs.gen.lc.display(lcs, [spline])
		



To learn more about the optional arguments of the spline optimizers, see the doc of :py:func:`pycs.spl.topopt.opt_rough` and :py:func:`pycs.spl.topopt.opt_fine`.

These spline optimizers also work with polynomial microlensing. You can mix the microlensing representations at will.

.. note:: Formally, the linear optimization of splines requires data points that are not only sorted, but also *strictly* increasing in jds : it cannot deal with lightcurves that have several data points taken at exactly the same epoch (which may happen as we shift the curves in time). This issue is automatically adressed by the class :py:class:`pycs.gen.spl.DataPoints`. As a user you don't have to worry about this in principle.


The regdiff optimizer
---------------------
This is de facto the easiest method to use, as it does not involve any explicit microlensing representation.

The idea is to shift the light curves so to minimize the variability of their "differences". To compute these difference curves, we need a regression, and in particular we use Gaussian Process Regression as provided by the ``pymc`` module.

.. note:: Therefore, to use the optimizer, you will have to **install** ``pymc`` first.
	
	Here is the website : `http://code.google.com/p/pymc/ <http://code.google.com/p/pymc/>`_



See the paper for a more detailed description of the idea.
In practice, as for the dispersion method and the splines, there is a simple top-level wrapper function, that you can add to your ``myopt.py`` :

::
	
	def regdiff(lcs):
		return pycs.regdiff.multiopt.opt_ts(lcs, pd=5, verbose=True)


But before blindly using the above optimizer, it is a good idea to test by yourself if the the Gaussian Process Regression (GPR) performs well on your lightcurve. The regressions are represented by "regularly sampled light curve" objects, implemented by the class :py:class:`pycs.regdiff.rslc.rslc`.
It is easy to perform a regression "manually", i.e. to obtain such a regularly sampled light curve starting from a usual light curve. The function that performs this GPR is :py:func:`pycs.regdiff.rslc.factory`, and you could for instance directly apply this directly to all your light curves :

::
	
	myrslcs = [pycs.regdiff.rslc.factory(l, pd=2) for l in lcs]
	# As this can take a minute, you might want to save the results :
	pycs.gen.util.writepickle(myrslcs, "myrslcs.pkl")


The parameter ``pd`` is a point density of the regression. Usually this is set to 5 (corresponding to one point every 0.2 days). Less points will give you a faster regression.

You can display these ``rslc`` objects with the usual display function, simply by putting them in the second argument list, as you would do for spline objects.

::
	
	myrslcs = pycs.gen.util.readpickle("myrslcs.pkl")
	pycs.gen.lc.display(lcs, myrslcs)


.. note:: These ``rslc`` have some attributes very similar to the usual ``lc`` objects, like ``jds``, ``mags``, ``magerrs``, ``plotcolour``. To shift an ``rslc`` in time, use ``myrslc.shifttime(12.3)``. To perform other operations, directly modify the attributes, for instance : ``myrslc.mags += 3.0``.


The reason why we want these finely sampled light curves is that we can easily subtract them from each other to get difference curves. This operation is implemented by :py:func:`pycs.regdiff.rslc.subtract`.  

::
	
	diffrslc = pycs.regdiff.rslc.subtract(myrslcs[0], myrslcs[1])
	# This diffrslc is the difference curve, and its again a rslc object.
	
	# Hence you can display the difference easily by putting it in the list, for instance :
	pycs.gen.lc.display([], myrslcs + [diffrslc])



Finally, the *WAV* of any ``rslc`` can be computed by calling the method :py:meth:`pycs.regdiff.rslc.rslc.wtv`.


