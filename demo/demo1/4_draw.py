import pycs
import myopt


# Two things to "draw" here : 

# 1) make say 200 plain copies of the data.
# These will be used to evaluate the intrinsic variance of the optimizers
# and compute the best point estimates for the delays.

lcs = pycs.gen.util.readpickle("data/trialcurves.pkl")
pycs.sim.draw.multidraw(lcs, onlycopy=True, n=10, npkl=20, simset="copies")


# 2) make say 1000 synthetic light curve sets with known true time delays,
# starting from the generative model.
# For this we don't need the "raw" data. Instead we use the "optimized" light curves,
# with their shifts and microlensing models, and the spline representing
# the intrinsic QSO varibility.

(modellcs, modelspline)  = pycs.gen.util.readpickle("optspline.pkl")

# At this place you could tweak the model, maybe adjust the ML or the delays by hand,
# if you know what you are doing. Otherwise, just go ahead :

pycs.sim.draw.saveresiduals(modellcs, modelspline)
pycs.sim.draw.multidraw(modellcs, modelspline, n=20, npkl=50, simset="1Kset1",
	truetsr=8.0, tweakml=[myopt.Atweakml, myopt.Btweakml, myopt.Ctweakml, myopt.Dtweakml])

# This "1Kset1" is a name that you can freely choose for your set of simulations.
# truetsr = 8.0 means that the synthetic curves will get random true time shifts
# in a range of 8.0 days around the time shifts of the modellcs.