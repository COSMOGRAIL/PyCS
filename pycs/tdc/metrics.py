"""
Stuff related to the TDC metrics
"""


import numpy as np
import matplotlib.pyplot as plt
import sys


############## Here, we work with estimates objects ##############

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
	lowest "precision" (== highest tderr/td) first -> select from end !
	"""
	return sorted(estimates, key = lambda e: np.fabs(e.tderr/e.td))[::-1]

def sortbyabstd(estimates):
	"""
	I sort your estimates according to the absolute value of their time delay.
	lowest "precision" (== lowest delays) first -> select from end !
	"""
	return sorted(estimates, key = lambda e: abs(e.td))



def maxPplot(estslist, N, filepath=None):
	"""
	Give me a list of estimate-lists, I plot P as a function of f
	N is the total number of curves (e.g. 56 for tdc0) (defines f=1)
	"""
	
	for ests in estslist:
		estscp = ests[:] # Important, we work on a copy, do not modify the original !
		for e in estscp:
			if e.td == 0.0: # We want to avoid that...
				#print e
				e.td = 0.1
		ests = sortbyP(estscp)
		fs = []
		Ps = []
		for n in range(len(ests)):
			subests = ests[n:]
			fs.append(f(subests, N))
			Ps.append(P(subests))
		plt.plot(fs, Ps, ".-", label=ests[0].method)
	plt.xlabel("f")
	plt.ylabel("Approximation of P")
	plt.xlim(0.0, 0.8)
	plt.ylim(0.0, 0.8)
	plt.axvline(0.3, color="black")
	plt.axhline(0.15, color="black")
	if len(estslist) > 1:
		plt.legend()
	plt.grid()
	if filepath:
		plt.savefig(filepath)
	else:
		plt.show()


		
############## Here, we work with the database from analyse_results ##############


def getP(db,method):
	"""
	Compute P for a given method stored in the database
	"""	
	subdb = [item for item in db if "%s_P" %(method) in item]
	#key = "%s_P" %method
	#print 'I RETURN: ',sum([item["%s_P" %(method)] for item in subdb]) / float(len(subdb))
	return sum([item["%s_P" %(method)] for item in subdb]) / float(len(subdb))
	
	
def getA(db,method):
	"""
	Compute A for a given method stored in the database
	"""	
	subdb = [item for item in db if "%s_A" %(method) in item]
	#print 'I RETURN: ',sum([item["%s_A" %(method)] for item in subdb]) / float(len(subdb))
	return sum([item["%s_A" %(method)] for item in subdb]) / float(len(subdb))
	
	
def getAmod(db,method):
	"""
	Compute A for a given method stored in the database
	"""	
	subdb = [item for item in db if "%s_Amod" %(method) in item]
	#print 'I RETURN: ',sum([item["%s_A" %(method)] for item in subdb]) / float(len(subdb))
	return sum([item["%s_Amod" %(method)] for item in subdb]) / float(len(subdb))		
	
def getchi2(db,method):
	"""
	Compute chi2 for a given method stored in the database
	"""	
	subdb = [item for item in db if "%s_chi2" %(method) in item]
	return sum([item["%s_chi2" %(method)] for item in subdb]) / float(len(subdb))
	
def getf(db,method,N):
	"""
	Compute f for a given method stored in the database
	N is the total number of curves	
	"""
	
	subdb = [item for item in db if "%s_td" %(method) in item]	
	return  float(len(subdb))/N


def Pplotall(db,methods,N):
	'''
	give me the db and a method, I plot P vs f for that method
	and chi2 vs f, A vs f for the same arrangement according to P
	Do it for one method only, then adapt to a list of methods...
	'''


	lfs =[]
	lPs=[]
	lAs=[]
	lchi2s=[]

	for method in methods:	
		db_P = [item for item in db if "%s_P" %(method) in item]	
		sorted_db_P = sorted(db_P, key = lambda item: item["%s_P" % method] )[::-1]

		fs = []
		Ps = []
		As = []
		chi2s = []
		for n in range(len(sorted_db_P)):
			subdb = sorted_db_P[n:]
			fs.append(getf(subdb,method=method,N=N))
			Ps.append(getP(subdb,method=method))
			As.append(getA(subdb,method=method))
			chi2s.append(getchi2(subdb,method=method))
		
		lfs.append(fs)
		lPs.append(Ps)
		lAs.append(As)
		lchi2s.append(chi2s)
			
	plt.subplot(3,1,1)
	for fs,Ps in zip(lfs,lPs):
		plt.plot(fs, Ps, ".-", label=method,color='indigo')
	plt.ylabel("P")
	plt.xlim(0.0, 0.8)
	plt.ylim(0.0, 0.1)
	plt.axvline(0.5, color="black")
	plt.axhline(0.03, color="black")

	plt.subplot(3,1,2)
	for fs,chi2s in zip(lfs,lchi2s):		
		plt.plot(fs, chi2s, ".-", label=method,color='chartreuse')				
	plt.ylabel("chi2")
	plt.xlim(0.0, 0.8)
	plt.ylim(0.0, 20)
	plt.axvline(0.5, color="black")
	plt.axhline(0.5, color="black")
	plt.axhline(1.5, color="black")

	plt.subplot(3,1,3)
	for fs,As in zip(lfs,lAs):
		plt.plot(fs, As, ".-", label=method,color='crimson')
	plt.xlabel("f")
	plt.ylabel("A")
	plt.xlim(0.0, 0.8)
	plt.ylim(-0.003, 0.003)
	plt.axvline(0.5, color="black")
	plt.axhline(0.03, color="black")
	plt.axhline(-0.03, color="black")
		
		
	plt.show()






