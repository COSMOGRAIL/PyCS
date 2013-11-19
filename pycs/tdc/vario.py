"""
Variability analysis stuff
"""


import numpy as np

def vario(l, plot=False, nsamp=10000):
	"""
	A simple stochatic variogram-like plot, and return the ratio
	delta_mag around 50-75 over delta_mag for the smallest step
	
	:param nsamp: the number of random samples to take

	"""
	
	jds = l.getjds()
	mags = l.getmags()
	magerrs = l.magerrs
	
	stats = l.samplingstats(seasongap=60)
	sampling = float(stats["med"])
	seasonlength = 365 - stats["meansg"]
	
	rsamp = int(np.floor(seasonlength / sampling))
	#print "sampling", sampling, "seasonlength", seasonlength, "rsamp", rsamp
	
	indas = np.random.random_integers(rsamp, len(jds)-1-rsamp, size=nsamp)
	#indbs = np.random.random_integers(0, len(jds)-1, size=nsamp)
	
	indbs = indas + np.random.random_integers(-rsamp, rsamp, size=nsamp)
	
	dts = np.fabs(jds[indas] - jds[indbs])
	dms = np.fabs(mags[indas] - mags[indbs])
	#oneseasonindices = dts < seasonlength
	#dts = dts[oneseasonindices]
	#dms = dms[oneseasonindices]
	
	# We bin these
	binlims = np.arange(sampling/2.0, seasonlength, sampling)
	indices = np.digitize(dts, binlims)
	#outofbins = np.logical_or(indices == 0, indices == len(binlims))
	binlows = binlims[:-1]
	binups = binlims[1:]
	bincenters = 0.5 * (binlows + binups)
	binmeans = np.array([np.mean(dms[indices == n]) for n in range(1, len(binlims))])
	binerrs = np.array([np.std(dms[indices == n]) / np.sqrt(len(dms[indices == n])) for n in range(1, len(binlims))])
	
	firstval = binmeans[0]
	others = np.arange(int(np.floor(50.0/sampling))-1, int(np.ceil(75.0/sampling)))
	#print others
	#binmeans[others] = 0.0
	
	zoneval = np.mean(binmeans[others])
	vratio = zoneval / firstval
	#print vratio
	
	#plt.scatter(dts, dms, s=2, lw=0)
	#plt.plot(bincenters, binmeans, label=str(l))
	
	if plot:
		import matplotlib.pyplot as plt
		plt.errorbar(bincenters, binmeans, yerr=binerrs, label=str(l)+" %.3f"%(vratio))
	
	
	return {"vratio":vratio, "sampling":sampling, "seasonlength":seasonlength}
		
	"""
	plt.xlim(-5, seasonlength)
	plt.xlabel("delta time [day]")
	plt.ylim(0, 0.3)
	plt.ylabel("delta mag")
	plt.show()
	"""



	
	
	
	

	

	