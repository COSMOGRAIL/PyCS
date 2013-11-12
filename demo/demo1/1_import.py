# In this first script we "import" the data, in this case from a simple text file with
# headers (other formats are also supported, see doc)

import pycs

rdbfile = "data/trialcurves.txt"

lcs = [
	pycs.gen.lc.rdbimport(rdbfile, 'A', 'mag_A', 'magerr_A', "Trial"),
	pycs.gen.lc.rdbimport(rdbfile, 'B', 'mag_B', 'magerr_B', "Trial"),
	pycs.gen.lc.rdbimport(rdbfile, 'C', 'mag_C', 'magerr_C', "Trial"),
	pycs.gen.lc.rdbimport(rdbfile, 'D', 'mag_D', 'magerr_D', "Trial")
]

pycs.gen.mrg.colourise(lcs) # Gives each curve a different colour.

# Let's shift them by the "true" time shifts, for display purposes :
lcs[1].shifttime(-5.0)
lcs[2].shifttime(-20.0)
lcs[3].shifttime(-70.0)

# We show them :
pycs.gen.lc.display(lcs)

# Or if you prefer to save them into a file:
pycs.gen.lc.display(lcs, filename="fig_trialcurves.pdf")
# This function has many more options...

# We undo these shifts, as from now on we "forget" about these true delays.
for l in lcs:
	l.resetshifts()

# The main point of this script : we save the raw curves into a pkl file :
pycs.gen.util.writepickle(lcs, "data/trialcurves.pkl")

# Normally we would stop here.
# In any further scripts, you can now import the data by reading this pickle file :
lcs = pycs.gen.util.readpickle("data/trialcurves.pkl")

# ... and do something with it.
for l in lcs:
	print l.longinfo()
	
# For instance, we could export the data into a text file.
for l in lcs:
	l.resetshifts()
pycs.gen.util.multilcsexport(lcs, "out_trialcurves.txt", separator="\t", verbose=True, properties=None)
# Which gives in this case the same file as the one from which you read the data in first place.

