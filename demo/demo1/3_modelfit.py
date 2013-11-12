# We perform the free-knot spline fit that will be used as starting point
# to draw the synthetic curves.

import pycs
import myopt

lcs = pycs.gen.util.readpickle("data/trialcurves.pkl")

pycs.gen.splml.addtolc(lcs[0], knotstep=150)
pycs.gen.splml.addtolc(lcs[1], knotstep=150)
pycs.gen.splml.addtolc(lcs[2], knotstep=150)
pycs.gen.splml.addtolc(lcs[3], knotstep=150)

spline = myopt.spl(lcs)
pycs.gen.lc.display(lcs, [spline], knotsize=0.01, figsize=(20, 7), jdrange=(53900, 55500),filename="fig_modelfit.pdf")

# We save all this into a pkl,
# including the lcs and their optimized microlensing splines.
pycs.gen.util.writepickle((lcs, spline), "optspline.pkl")

