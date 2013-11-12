"""
Abandonned due to circularity.

Sample : related to the empirical model stuff.
These functions "select" runresults that are representative of the true delays.

"""

#import pycs.gen.lc
#import pycs.sim.draw
#import pycs.gen.util
import numpy as np
#import os
#import time




def keep(datarrs, simrrs, r=3.0, insigma=False):
	"""
	Selects those simrrs that did fall close to the datarrs.
	
	:param datarrs: runresults obtained on the actual data (i.e., draw = False).

	:param simrrs: runresults obtained by applying the black box on simulated curves, with a distribution of true shifts.
	
	:param r: radius around the mean of the delays of datarrs , in units of sigma if insigma is True
	
	
	:returns: a runresults object with only the selected results.
	
	"""
	
	print "Data : %s, %s, %i runs." % (datarrs.name, ", ".join(datarrs.labels), len(datarrs))
	print "Sims : %s, %s, %i runs." % (simrrs.name, ", ".join(simrrs.labels), len(simrrs))
	
	if datarrs.labels != simrrs.labels:
		raise RuntimeError("No the same curves !") # Important ! We rely on them being in the same order !
	
	# We build some rough ranges of delays to keep, using datarrs.
	
	n = len(datarrs.labels)
	coupleindices = [(i, j) for i in range(n) for j in range(n) if i < j]
	
	# We get the samples of delays, for each couple
	datamesdelays = [datarrs.tsarray[:,j] - datarrs.tsarray[:,i] for (i, j) in coupleindices]
	
	if insigma == True:
		coupleranges = [(np.mean(delays) - r*np.std(delays), np.mean(delays) + r*np.std(delays)) for delays in datamesdelays]
	elif insigma == False:
		coupleranges = [(np.mean(delays) - r, np.mean(delays) + r) for delays in datamesdelays]

	
	print "Ranges :"
	for i in range(len(coupleindices)):
		print "%s%s : (%8.2f, %8.2f)" % (datarrs.labels[coupleindices[i][0]], datarrs.labels[coupleindices[i][1]], coupleranges[i][0], coupleranges[i][1])
	
	
	# We go through the simrrs, and keep only those whose measured delays are within these ranges	
	simmesdelays = [simrrs.tsarray[:,j] - simrrs.tsarray[:,i] for (i, j) in coupleindices]
	
	# We process one couple at the time :
	keepmasklist = [np.logical_and(simmesdelays[i] > coupleranges[i][0], simmesdelays[i] < coupleranges[i][1]) for i in range(len(coupleindices))]
	# We collapse the conditions for each couple :
	keepmask = reduce(np.logical_and, keepmasklist)
	
	
	print "Keeping %i runs among %i." % (np.sum(keepmask), len(simrrs))
	#print keepmasks
	keeprrs = simrrs.copy()
	keeprrs.applymask(keepmask)
	keeprrs.name += "keep"
	keeprrs.plottrue = True
	keeprrs.plotgauss = True # Usually we want this, for the selected curves.
	keeprrs.plotcolour = "black"
	
	return keeprrs
	
	
	
	
	