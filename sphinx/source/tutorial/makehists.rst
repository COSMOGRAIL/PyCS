Getting delay histograms and results
====================================


Generalities
------------

To obtain time-delay histograms and point and uncertainty estimates, several steps need to be done; they are presented in the following subsections.
We assume that you have generated some mock curves as described in the previous section.



Run curve shifting algorithms on the mock curves
------------------------------------------------


The wrapper function to run any method on any set of mock curves has the funny name :py:func:`pycs.sim.run.multirun`.
As input, it takes those pkl files containing simulated curves, made in the previous section. The ouput of this step are pkl files of runresult objects.

One important argument of ``multirun`` is ``tsrand``. It is the radius (in days) of a uniform randomization of the input time shifts.

.. note:: Something nice : you can launch several identical "calls" to this multirun function (for instance simply by launching your script on several CPUs), and this will indeed process the pkl files in parallel. For this to work, the ``multirun`` function stores a temporary file in its results directory as soon as it starts working on a pickle file, so that other scripts know that they should not run on this same pkl file as well. You can see those temporary files, of course. If something goes wrong and they don't get deleted automatically as the scirpt crashed, you might have to remove them by hand, otherwise ``multirun`` will just skip those pkl files.

Here is an example script running all 3 methods. We assume that you've defined your favorite optimizers in ``myopt.py``.


::
	
	lcs = pycs.gen.util.readpickle("merged.pkl")

	
	# Spline method
	"""
	
	# Initial conditions for the analysis :
	# (You might want to set some timeshifts, fluxshifts, as well)
	# In this case we'll just add some spline ML :
	pycs.gen.splml.addtolc(lcs[0], knotstep=100)
	pycs.gen.splml.addtolc(lcs[2], knotstep=300)
	pycs.gen.splml.addtolc(lcs[3], knotstep=300)
	
	#pycs.sim.run.multirun("copies", lcs, myopt.spl, optset="spl1", tsrand=10.0)
	#pycs.sim.run.multirun("sim1tsr5", lcs, myopt.spl, optset="spl1", tsrand=10.0)
	"""
	
	
	# Disp method
	"""
	# We set up some polynomial ML :
	pycs.gen.polyml.addtolc(lcs[0], nparams=3, autoseasonsgap = 60.0)
	pycs.gen.polyml.addtolc(lcs[2], nparams=4, autoseasonsgap = 1000.0)
	pycs.gen.polyml.addtolc(lcs[3], nparams=4, autoseasonsgap = 1000.0)
	
	#myopt.disp(lcs)
	#pycs.gen.lc.display(lcs)
	
	#pycs.sim.run.multirun("copies", lcs, myopt.disp, optset="disp1", tsrand=10.0)
	#pycs.sim.run.multirun("sim1tsr5", lcs, myopt.disp, optset="disp1", tsrand=10.0)
	"""
	
	
	# Regdiff method
	"""
	# No need for ml, but time delays should be set about right.
	
	#pycs.sim.run.multirun("copies", lcs, myopt.regdiff, optset="regdiff2", tsrand=10.0)
	#pycs.sim.run.multirun("sim1tsr5", lcs, myopt.regdiff, optset="regdiff2", tsrand=10.0)
	"""
	



Analysing the measurement results
---------------------------------


We read the "runresults" pickle files created at the previous step, and turn them into plots.
This is very flexible, as you might want to plot and analyse many things.

To start, we have the function :py:func:`pycs.sim.run.collect` that collects all the results from one directory::

	results = pycs.sim.run.collect(directory="./for/example/sims_copies_opt_spl")

The resulting object ``results`` is an instance of the class :py:class:`pycs.sim.run.runresults`. If you want to perform your own analysis of the results, you could directly access the following attributes::

	print results.labels # A list of the QSO image names (defines the order of QSO images with which the following results are given)
	print results.tsarray # A 2D array with the measured time shifts. Shape is (number of sets, number of QSO images)
	print results.truetsarray # Idem, for the TRUE time shifts, in case of simulated data
	print results.qs # A 1D array with the "chi2" or dispersion values. Shape is (number of sets).

Note that these "tsarrays" contain time shifts, not time delays. To get time delays between images "A" and "B" (i.e., ``results.labels[0]`` and ``results.labels[1]``), you would have to compute the differences yourself::

	measured_delays = results.tsarray[:,1] - results.tsarray[:,0]
	print measured_delays


If you want to go straight to some more or less automatic plots showing the results, here is a typical example:

::

		
	copiesres = [
		pycs.sim.run.collect("sims_copies_opt_disp1", "red", "Dispersion"),
		pycs.sim.run.collect("sims_copies_opt_spl1", "blue", "Spline"),
		pycs.sim.run.collect("sims_copies_opt_regdiff1", "green", "Regdiff")
	]
	
	pycs.sim.plot.hists(copiesres, r=30.0, nbins=100, dataout =True)
	
	
	simres = [
		pycs.sim.run.collect("sims_sim1tsr5_opt_disp1", "red", "Dispersion"),
		pycs.sim.run.collect("sims_sim1tsr5_opt_spl1", "blue", "Splines"),
		pycs.sim.run.collect("sims_sim1tsr5_opt_regdiff1", "green", "Regdiff")
	]
	
	
	pycs.sim.plot.hists(simres, r=30.0, nbins=100, dataout =True)
	
	pycs.sim.plot.measvstrue(simres, r=5.0, nbins = 10, plotpoints=True, ploterrorbars=True, sidebyside=True, errorrange=8, binclip=False, binclipr=20.0, dataout =True)
	
	delays_disp = pycs.gen.util.readpickle( "sims_copies_opt_disp1_delays.pkl")
	delays_spl = pycs.gen.util.readpickle( "sims_copies_opt_spl1_delays.pkl")
	delays_regdiff = pycs.gen.util.readpickle( "sims_copies_opt_regdiff1_delays.pkl")
	
	errorbar_disp = pycs.gen.util.readpickle( "sims_sim1tsr5_opt_disp1_errorbars.pkl")
	errorbar_spl = pycs.gen.util.readpickle( "sims_sim1tsr5_opt_spl1_errorbars.pkl")
	errorbar_regdiff = pycs.gen.util.readpickle( "sims_sim1tsr5_opt_regdiff_errorbars.pkl")
	
	totcontainer = [(delays_disp,errorbar_disp),(delays_spl,errorbar_spl),(delays_regdiff,errorbar_regdiff)]
	
	pycs.sim.plot.newdelayplot(totcontainer, rplot=6.0, displaytext=True)
	
	
	
	


	
