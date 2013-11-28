"""
Stuff related to the TDC metrics
"""


import numpy as np
import matplotlib.pyplot as plt


def fN(estimates):
	return float(len(estimates))

def f(estimates, N):
	return float(len(estimates))/float(N)

def P(estimates):
	"""
	APPROXIMATION of the P metric...
	"""
	return (1.0/fN(estimates)) * np.sum(np.fabs(np.array([e.tderr/e.td for e in estimates])))
	#return (1.0/56.0)* np.sum(np.fabs(np.array([e.tderr/e.td for e in estimates])))
	
	
def sortbyP(estimates):
	"""
	I sort your estimates according to their claimed precision
	lowest precision first
	"""
	return sorted(estimates, key = lambda e: np.fabs(e.tderr/e.td))[::-1]


def maxPplot(estslist, N):
	"""
	Give me a list of estimate-lists, I plot P as a function of f
	N is the total number of curves (e.g. 56 for tdc0) (defines f=1)
	"""
	
	for ests in estslist:
		ests = sortbyP(ests)
		fs = []
		Ps = []
		for n in range(len(ests)):
			subests = ests[n:]
			fs.append(f(subests, N))
			Ps.append(P(subests))
		plt.plot(fs, Ps, ".-", label=ests[0].method)
	plt.xlabel("f")
	plt.ylabel("Approximation of P")
	plt.xlim(0.0, 0.6)
	plt.ylim(0.0, 0.5)
	plt.axvline(0.3, color="black")
	plt.axhline(0.15, color="black")
	plt.legend()
	plt.show()
		

	
	
	
	

	

	