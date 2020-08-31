# We perform the free-knot spline fit that will be used as starting point
# to draw the synthetic curves.

import pycs
from pycs.regdiff.rslc import factory
from pycs.gen.lc import settimeshifts
import myopt

lcs = pycs.gen.util.readpickle("data/trialcurves.pkl")
lcs2 = [lc.copy() for lc in lcs]

#spline
# pycs.gen.splml.addtolc(lcs[0], knotstep=150)
# pycs.gen.splml.addtolc(lcs[1], knotstep=150)
# pycs.gen.splml.addtolc(lcs[2], knotstep=150)
# pycs.gen.splml.addtolc(lcs[3], knotstep=150)
#
# spline = myopt.spl(lcs)
# pycs.gen.lc.display(lcs, [spline], knotsize=0.01, figsize=(20, 7), jdrange=(53900, 55500),filename="fig_modelfit.pdf")
#
# # We save all this into a pkl,
# # including the lcs and their optimized microlensing splines.
# pycs.gen.util.writepickle((lcs, spline), "optspline.pkl")

#regdiff:
settimeshifts(lcs,shifts=[0, -5, -20, -60],includefirst=True)
myrslcs = [factory(l, pd=2, covkernel='matern',
                                     pow=1.7, amp=1., scale=200.0, errscale=2.) for l in lcs2]
pycs.regdiff.multiopt.opt_ts(lcs, pd=2, covkernel='matern',
                                     pow=1.7, amp=1., scale=200.0, errscale=2., verbose=True, method="weights")
for lc in lcs:
    print(lc.timeshift)
pycs.gen.lc.display(lcs2, myrslcs,filename="screen")
