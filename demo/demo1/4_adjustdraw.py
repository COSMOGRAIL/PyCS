import pycs
import myopt

# The aim of this script is to adjust the parameters of the randomized fast extrinsic variability
# (power law noise) that we will include in our synthetic curves.
# The functions that control this noise are called "tweakml" (in myopt.py). There's one for each curve.
# For a first overview you don't have to run this present script, as the current settings are fine.
# No output of this script is required for further steps -- just skip it and come back later.



# Scenario 1 : you have already drawn and analysed some synthetic curves (scripts 4 and 5)
# and want to check how similar this was to the observations.

"""
# The adjustement is based on residulas from a spline fit.
# We have to select a spline fit on the real data, that was obtained using the same parameters as when
# the spline technique is run on the synthetic curves. In this demo, this is the same spline as used by the 
# model fit :
(modellcs, modelspline)  = pycs.gen.util.readpickle("optspline.pkl")
# We compute the residuals of this fit, and save them as attributes to the light curves :
pycs.sim.draw.saveresiduals(modellcs, modelspline)

# You can now directly compare these residuals to the results from the simulations :
pycs.gen.stat.anaoptdrawn(modellcs, modelspline, simset="1Kset1", optset="spl")
# Which gives you the plots shown in the paper.
# (Change simset and optset to whatever name you have chosen)

"""
 

# Scenario 2 : you want to try a few tweakml settings on a small number of synthetic curves
# and see how it looks. This takes a few minutes to run for each settings.

"""
simset = "test1"

# Let's draw say 10 curves :
(modellcs, modelspline)  = pycs.gen.util.readpickle("optspline.pkl")
pycs.sim.draw.saveresiduals(modellcs, modelspline)

pycs.sim.draw.multidraw(modellcs, modelspline, n=10, npkl=1, simset=simset, truetsr=8.0, tweakml=[myopt.Atweakml, myopt.Btweakml, myopt.Ctweakml, myopt.Dtweakml])

# And fit them with the spline technique (see script 5...) :

lcs = pycs.gen.util.readpickle("data/trialcurves.pkl")
# We set some plausible guess delays :
lcs[1].shifttime(-7.0)
lcs[2].shifttime(-22.0)
lcs[3].shifttime(-65.0)
# We add the ML :
pycs.gen.splml.addtolc(lcs[0], knotstep=150)
pycs.gen.splml.addtolc(lcs[1], knotstep=150)
pycs.gen.splml.addtolc(lcs[2], knotstep=150)
pycs.gen.splml.addtolc(lcs[3], knotstep=150)
# And run with these settings on the simulated curves :
pycs.sim.run.multirun(simset, lcs, myopt.spl, optset="spl", tsrand=10.0, keepopt=True)

# Finally, we compare with residuals from the fit of the real observations, like for scenario 1 :
pycs.gen.stat.anaoptdrawn(modellcs, modelspline, simset=simset, optset="spl")

# Now change the tweakml, the simset name, and try again...

# Let's delete these test sims, so that if you rerun this script with different tweakml settings
# but omit to change the simset name, only the new sims are used :-)
import shutil
shutil.rmtree("sims_%s" % (simset))
shutil.rmtree("sims_%s_opt_spl" % (simset))

"""


# Scenario 3 : quick and dirty adjustment, without refitting the spline ...
# You can use this to quickly get some guess settings.
# Note that due to the fact that we do not refit the spline, you should aim here at
# getting simulated residuals that are slightly larger than the observed ones.
"""
(modellcs, modelspline)  = pycs.gen.util.readpickle("optspline.pkl")
pycs.sim.draw.saveresiduals(modellcs, modelspline)

real_residuals = pycs.gen.stat.subtract(modellcs, modelspline)

mocklcs = pycs.sim.draw.draw(modellcs, modelspline, tweakml=[myopt.Atweakml, myopt.Btweakml, myopt.Ctweakml, myopt.Dtweakml])
for l in mocklcs:
        l.plotcolour = "black"
sim_residuals = pycs.gen.stat.subtract(mocklcs, modelspline)

pycs.gen.stat.plotresiduals([real_residuals, sim_residuals], nicelabel=False, showsigmalines = True, errorbarcolour = "#999999")
"""