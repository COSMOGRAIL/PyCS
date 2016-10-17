import pycs

# We read the results obtained on the copies of the observed light curves :
dataresults = [
	pycs.sim.run.collect("sims_copies_opt_disp", "red", "Dispersion-like technique"),
	pycs.sim.run.collect("sims_copies_opt_regdiff", "green", "Regression difference technique"),
	pycs.sim.run.collect("sims_copies_opt_spl", "blue", "Free-knot spline technique")
]

# ... and turn this into simple histograms, that will give the intrinsic variance.
# The option dataout=True will save the delay point estimate, to be used below.
pycs.sim.plot.hists(dataresults, r=5.0, nbins=100, showqs=False,
	filename="fig_intrinsicvariance.pdf", dataout=True)
	
	
# We read the results obtained on the synthetic curves :
simresults = [
	pycs.sim.run.collect("sims_1Kset1_opt_disp", "red", "Dispersion-like technique"),
	pycs.sim.run.collect("sims_1Kset1_opt_regdiff", "green", "Regression difference technique"),
	pycs.sim.run.collect("sims_1Kset1_opt_spl", "blue", "Free-knot spline technique")
]

# ... and perform the error analysis :
# The option dataout=True will save the random and systematic error, to be used below
pycs.sim.plot.measvstrue(simresults, errorrange=3.5, r=5.0, nbins = 10, binclip=True, binclipr=20.0,
	plotpoints=False, filename="fig_measvstrue.pdf", dataout=True)

# With this same data we can also illustrate the correlations between measurements :
pycs.sim.plot.covplot(simresults, filename="fig_covplot.pdf")


# Finally we group the information saved by these steps to get the results in form of a summary plot :

dispres = (pycs.gen.util.readpickle("sims_copies_opt_disp_delays.pkl"),
	pycs.gen.util.readpickle("sims_1Kset1_opt_disp_errorbars.pkl"))

regdiffres = (pycs.gen.util.readpickle("sims_copies_opt_regdiff_delays.pkl"),
	pycs.gen.util.readpickle("sims_1Kset1_opt_regdiff_errorbars.pkl"))

splres = (pycs.gen.util.readpickle("sims_copies_opt_spl_delays.pkl"),
	pycs.gen.util.readpickle("sims_1Kset1_opt_spl_errorbars.pkl"))

pycs.sim.plot.newdelayplot([dispres, regdiffres, splres], rplot=6.0, displaytext=True,
	filename = "fig_delays.pdf", refshifts=[{"colour":"gray", "shifts":(0, -5, -20, -70)}])


