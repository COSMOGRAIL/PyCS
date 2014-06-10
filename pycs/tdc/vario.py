"""
Variability analysis stuff
"""


import numpy as np

def vario(l, plot=False, filepath=None, nsamp=1000000, verbose=False):
	"""
	A simple stochatic variogram-like plot, and return the ratio
	delta_mag around 50-75 over delta_mag for the smallest step
	
	:param nsamp: the number of random samples to take

	"""
	
	if verbose:
		print "Starting vario analysis of %s" % (str(l))
	
	jds = l.getjds()
	mags = l.getmags()
	magerrs = l.magerrs
	
	stats = l.samplingstats(seasongap=60)
	sampling = float(stats["med"]) # The median sampling, in days
	seasonlength = 365 - stats["meansg"] # The mean season length, in days
	
	rsamp = int(np.floor(seasonlength / sampling)) # The typical number of datapoints we have per season
	#print "sampling", sampling, "seasonlength", seasonlength, "rsamp", rsamp
	
	# We now draw indices of points in lca and lcb to form pairs contributing to the variogram.
	# Instead of just drawing them randomly, we are smart and draw them relatively close by
	indas = np.random.random_integers(rsamp, len(jds)-1-rsamp, size=nsamp)
	#indbs = np.random.random_integers(0, len(jds)-1, size=nsamp)
	
	indbs = indas + np.random.random_integers(-rsamp, rsamp, size=nsamp)
	
	dts = np.fabs(jds[indas] - jds[indbs]) # time differences of the pairs
	dms = np.fabs(mags[indas] - mags[indbs]) # mag differences of the pairs
	#oneseasonindices = dts < seasonlength
	#dts = dts[oneseasonindices]
	#dms = dms[oneseasonindices]
	
	# We bin these
	binlims = np.arange(sampling/2.0, seasonlength, sampling) # from half the sampling to seasonlength, in days.
	indices = np.digitize(dts, binlims)
	#outofbins = np.logical_or(indices == 0, indices == len(binlims))
	binlows = binlims[:-1]
	binups = binlims[1:]
	bincenters = 0.5 * (binlows + binups)
	binmeans = np.array([np.mean(dms[indices == n]) for n in range(1, len(binlims))])
	binerrs = np.array([np.std(dms[indices == n]) / np.sqrt(len(dms[indices == n])) for n in range(1, len(binlims))])
	
	firstval = binmeans[0]
	others = np.arange(int(np.floor(50.0/sampling))-1, int(np.ceil(75.0/sampling)))
	# these are roughly the indices that correspond to the range from 50 to 75 days
	#print others
	#binmeans[others] = 0.0
	
	zoneval = np.mean(binmeans[others])
	vratio = zoneval / firstval
	#print vratio
	
	#plt.scatter(dts, dms, s=2, lw=0)
	#plt.plot(bincenters, binmeans, label=str(l))
	
	if plot:
		import matplotlib.pyplot as plt
		plt.clf()
		plt.errorbar(bincenters, binmeans, yerr=binerrs, label=str(l)+" %.3f"%(vratio))
		plt.title(l.object)
		plt.xlabel("Time separation")
		plt.ylabel("Average mag separation")
		if filepath != None:
			
			plt.savefig(filepath)
				
		else:
			plt.show()
		
	
	
	return {"vratio":vratio, "sampling":sampling, "seasonlength":seasonlength, "zoneval":zoneval, "firstval":firstval, "bincenters":bincenters, "binmeans":binmeans, "binerrs":binerrs}
		
	"""
	plt.xlim(-5, seasonlength)
	plt.xlabel("delta time [day]")
	plt.ylim(0, 0.3)
	plt.ylabel("delta mag")
	plt.show()
	"""



	
	
	
	

	

	