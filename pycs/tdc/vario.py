"""
Variability analysis stuff
"""


import numpy as np
import matplotlib.pyplot as plt

def vario(l):
	
	jds = l.getjds()
	mags = l.getmags()
	magerrs = l.magerrs
	
	stats = l.samplingstats(seasongap=60)
	sampling = stats["med"]
	seasonlength = 365 - stats["meansg"]
	
	
	
	nsamp = 1000
	rsamp = int(np.floor(seasonlength / sampling))
	
	print sampling, seasonlength, rsamp
	
	
	indas = np.random.random_integers(rsamp, len(jds)-1-rsamp, size=nsamp)
	#indbs = np.random.random_integers(0, len(jds)-1, size=nsamp)
	
	indbs = indas + np.random.random_integers(-rsamp, rsamp, size=nsamp)
	
	
	dts = np.fabs(jds[indas] - jds[indbs])
	dms = np.fabs(mags[indas] - mags[indbs])
	
	#oneseasonindices = dts < seasonlength
	#dts = dts[oneseasonindices]
	#dms = dms[oneseasonindices]
	
	plt.scatter(dts, dms, s=2, lw=0)
	plt.xlim(-5, 5000)
	plt.ylim(0, 0.3)
	plt.show()