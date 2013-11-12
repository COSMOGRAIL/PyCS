import pycs
import myopt

# We will now run the curve-shifting methods on these data sets.

# To "define" the microlensing models and initial shifts,
# we make use of one single set of curves, e.g., the observed data.

lcs = pycs.gen.util.readpickle("data/trialcurves.pkl")

# We set some plausible guess delays (by eye, or maybe from the modelfit)
lcs[1].shifttime(-7.0)
lcs[2].shifttime(-22.0)
lcs[3].shifttime(-65.0)

# Now uncomment the calls to multirun below, and run them on the "copies" and on the synthetic data.
# These functions are designed so that several python instances can run them "in parallel"
# on the same data. They will not process pkl files which have already been measured.
# So when calling multirun by hand like we do here, I typically uncomment them one by one,
# and launch the script several times.


# Free-knot spline technique
"""
# We define the microlensing model. Here we use the same as in 3_modelfit :
pycs.gen.splml.addtolc(lcs[0], knotstep=150)
pycs.gen.splml.addtolc(lcs[1], knotstep=150)
pycs.gen.splml.addtolc(lcs[2], knotstep=150)
pycs.gen.splml.addtolc(lcs[3], knotstep=150)

pycs.sim.run.multirun("copies", lcs, myopt.spl, optset="spl", tsrand=10.0, keepopt=True)
pycs.sim.run.multirun("1Kset1", lcs, myopt.spl, optset="spl", tsrand=10.0, keepopt=True)

# When you run the spline technique using this multirun function, put keepopt=True as shown.
# This will allow to easily compare the residuals from the synthetic curves
# with the residuals from the observed data.

# You are free to change optset to whatever name you like. It should reflect the full method,
# including the settings of the microlensing. 

"""

# Dispersion-like technique
"""
pycs.gen.polyml.addtolc(lcs[0], nparams=2, autoseasonsgap = 60.0)
pycs.gen.polyml.addtolc(lcs[1], nparams=2, autoseasonsgap = 60.0)
pycs.gen.polyml.addtolc(lcs[2], nparams=2, autoseasonsgap = 60.0)
pycs.gen.polyml.addtolc(lcs[3], nparams=2, autoseasonsgap = 60.0)

pycs.sim.run.multirun("copies", lcs, myopt.disp, optset="disp", tsrand=10.0)
pycs.sim.run.multirun("1Kset1", lcs, myopt.disp, optset="disp", tsrand=10.0)

"""

# Regression difference technique
"""
pycs.sim.run.multirun("copies", lcs, myopt.regdiff, optset="regdiff", tsrand=10.0)
pycs.sim.run.multirun("1Kset1", lcs, myopt.regdiff, optset="regdiff", tsrand=10.0)
"""
