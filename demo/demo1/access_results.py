
import pycs

results = pycs.sim.run.collect(directory="./sims_copies_opt_spl")

# The object "results" is an instance of the class pycs.sim.run.runresults. You can directly access the following attributes:

print results.labels # A list of the QSO image names (defines the order of QSO images with which the following results are given)
print results.tsarray # A 2D array with the measured time shifts. Shape is (number of sets, number of QSO images)
print results.truetsarray # Idem, for the TRUE time shifts, in case of simulated data
print results.qs # A 1D array with the "chi2" or dispersion values. Shape is (number of sets).

# Note that these "tsarrays" contain time shifts, not time delays.
# To get time delays between images "A" and "B" (i.e., results.labels[0] and results.labels[1]),
# you will have to compute the differences yourself:

measured_delays = results.tsarray[:,1] - results.tsarray[:,0]
print measured_delays