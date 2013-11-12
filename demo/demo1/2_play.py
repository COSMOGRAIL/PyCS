import pycs

lcs = pycs.gen.util.readpickle("data/trialcurves.pkl")

# To see what we have :
for l in lcs:
	print l

# Let's try some curve shifting, without correcting for extrinsic variability...
# With the freek-knot spline technique :
import myopt
spline = myopt.spl(lcs)

# Show the result
pycs.gen.lc.display(lcs, [spline])
# (It doesn't look good, these curves don't overlap without a microlensing model)

# For humans, at any time we can print out the time delays between the curves,
# computed from the current time shifts of each curve :
print "Time delays:"
print pycs.gen.lc.getnicetimedelays(lcs, separator="\n", sorted=True)

# To get better results, we need to add microlensing models, e.g., polynomials
# to our light curves.
pycs.gen.polyml.addtolc(lcs[1], nparams=2, autoseasonsgap=60.0)
pycs.gen.polyml.addtolc(lcs[2], nparams=3, autoseasonsgap=600.0)
pycs.gen.polyml.addtolc(lcs[3], nparams=3, autoseasonsgap=600.0)
# (this choice is just an illustration)

# Let's try the free-knot spline optimization again :
spline = myopt.spl(lcs)
pycs.gen.lc.display(lcs, [spline])
# The result is obviously already better
