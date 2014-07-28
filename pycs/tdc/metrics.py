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
	and chi2 vs f, A vs f... for the same arrangement according to P

	'''


	lfs =[]
	lPs=[]
	lAs=[]
	lAmods=[]
	lchi2s=[]
	labels = []


	for method in methods:	
		db_P = [item for item in db if "%s_P" %(method) in item]	
		sorted_db_P = sorted(db_P, key = lambda item: item["%s_P" % method] )[::-1]
		
		fs = []
		Ps = []
		As = []
		Amods= []
		chi2s = []
		for n in range(len(sorted_db_P)):
			subdb = sorted_db_P[n:]
			fs.append(getf(subdb,method=method,N=N))
			Ps.append(getP(subdb,method=method))
			As.append(getA(subdb,method=method))
			Amods.append(getAmod(subdb,method=method))
			chi2s.append(getchi2(subdb,method=method))
		
		lfs.append(fs)
		lPs.append(Ps)
		lAs.append(As)
		lAmods.append(Amods)
		lchi2s.append(chi2s)
		labels.append(method)

	
	# And now, the plot
	
	
	colors = ["blue","green","red","chartreuse","crimson","black","magenta"]
	
	plt.figure('metrics vs f',(10,15))

	plt.subplot(4,1,1)
	for fs,Ps in zip(lfs,lPs):
		plt.plot(fs, Ps, ".-", label=labels[lfs.index(fs)], color=colors[lfs.index(fs)])
	plt.ylabel(r"$P$")
	plt.xlim(0.0, 0.5)
	plt.ylim(min([min(Ps) for Ps in lPs]),max([max(Ps) for Ps in lPs]))
	plt.axvline(0.5, color="black")
	plt.axhline(0.03, color="black")
	
	s = r"$ P = \frac{1}{fN} \sum_i \left( \frac{\sigma_i}{|\Delta t_i|} \right)$" 
	#s = 'P = 1/fN * sum(tderr_i/|td_i|)'
	plt.annotate(s=s, xy = (0.73,0.65), xycoords='axes fraction',
         		textcoords='axes fraction',size=18)	
	plt.legend(fontsize= 12,loc=2)
	
	
	plt.subplot(4,1,2)
	for fs,chi2s in zip(lfs,lchi2s):		
		plt.plot(fs, chi2s, ".-", label=labels[lfs.index(fs)],color=colors[lfs.index(fs)])				
	plt.ylabel(r"$\chi^2$")
	plt.xlim(0.0, 0.5)
	plt.ylim(min([min(chi2s) for chi2s in lchi2s]),max([max(chi2s) for chi2s in lchi2s]))
	plt.axvline(0.5, color="black")
	plt.axhline(0.5, color="black")
	plt.axhline(1.5, color="black")
	s = r"$ \chi^2 =\frac{1}{fN} \sum_i \left( \frac{\overline{\Delta} t_i - \Delta t_i}{\sigma_i} \right)$"

	plt.annotate(s=s, xy = (0.73,0.75), xycoords='axes fraction',
         		textcoords='axes fraction',size=18)

	plt.subplot(4,1,3)
	for fs,As in zip(lfs,lAs):
		plt.plot(fs, As, ".-", label=labels[lfs.index(fs)],color=colors[lfs.index(fs)])

	plt.ylabel(r"$A$")
	plt.xlim(0.0, 0.5)
	plt.ylim(min([min(As) for As in lAs]),max([max(As) for As in lAs]))
	plt.axvline(0.5, color="black")
	plt.axhline(0.03, color="black")
	plt.axhline(-0.03, color="black")
	
	s = r"$ A =\frac{1}{fN} \sum_i \left( \frac{\overline{\Delta} t_i - \Delta t_i}{|\Delta t_i|} \right)$"

	plt.annotate(s=s, xy = (0.73,0.75), xycoords='axes fraction',
         		textcoords='axes fraction',size=18)
	
	
	plt.subplot(4,1,4)
	for fs,Amods in zip(lfs,lAmods):
		plt.plot(fs, Amods, ".-", label=labels[lfs.index(fs)],color=colors[lfs.index(fs)])
	plt.xlabel(r"$f$")
	plt.ylabel(r"$A_{mod}$")
	plt.xlim(0.0, 0.5)
	plt.ylim(min([min(Amods) for Amods in lAmods]),max([max(Amods) for Amods in lAmods]))
	plt.axvline(0.5, color="black")
	plt.axhline(0.03, color="black")
	plt.axhline(-0.03, color="black")		
	s = r"$ A =\frac{1}{fN} \sum_i \left( \frac{\overline{\Delta} t_i - \Delta t_i}{\Delta t_i} \right)$"

	plt.annotate(s=s, xy = (0.73,0.75), xycoords='axes fraction',
         		textcoords='axes fraction',size=18)
			
	plt.show()



###### Towards a combined metric (experimental)




def combigauss(subtds, subtderrs, truetds, lensmodelsigma = 0.0):
	"""
	Give me submission delays and error bars, as well as the corresponding true delays, in form of numpy arrays.
	I compute the mean and sigma of the combined posterior on the fractional time delay distance error.
	"""

	from scipy.stats import norm

	subtdoffs = subtds - truetds
	centers = subtdoffs/truetds
	sigmas = subtderrs/np.fabs(truetds)
	
	# We convolve with the lensmodelsigma:
	sigmas = np.sqrt(sigmas**2 + lensmodelsigma**2)
	
	sigma_combi = 1.0 / np.sqrt(np.sum(1.0 / (sigmas**2)))
	center_combi = sigma_combi**2 * np.sum( centers/sigmas**2 )
	
	probazero = norm.pdf(0.0, center_combi, sigma_combi)
	
	return (center_combi, sigma_combi, probazero)
	
	# To plot this you could do:
	#plt.plot(x, norm.pdf(x, center_combi, sigma_combi), ls="--", color="black")
	



def Pplotcombi(db,methods,N,lensmodelsigma = 0.0):
	'''
	give me the db and a list of methods, I plot combigauss params (center, sigma, zerovalue) vs f
	f is selected according to bestP
	'''


	lfs =[]
	lcs=[]
	lss=[]
	lzs=[]
	labels = []


	for method in methods:	
		db_P = [item for item in db if "%s_P" %(method) in item]	
		sorted_db_P = sorted(db_P, key = lambda item: item["%s_P" % method] )[::-1]
		
		fs = []
		cs = []
		ss= []
		zs = []
		for n in range(len(sorted_db_P)):
			subdb = sorted_db_P[n:]
			fs.append(getf(subdb,method=method,N=N))

			subtds = []
			subtderrs = []
			truetds = []
			for entry in subdb:
				subtds.append(entry["%s_td" % method])
				subtderrs.append(entry["%s_tderr" % method])
				truetds.append(entry["truetd"])
			(c, s, z) = combigauss(np.array(subtds), np.array(subtderrs), np.array(truetds), lensmodelsigma = lensmodelsigma) 
			cs.append(c)
			ss.append(s)
			zs.append(z)


		lfs.append(fs)
		lcs.append(cs)
		lss.append(ss)
		lzs.append(zs)
		labels.append(method)

	
	# And now, the plot
	
	
	colors = ["blue","green","red","chartreuse","crimson","black","magenta"]
	
	plt.figure('combigauss vs f',(10,15))

	plt.subplot(3,1,1)
	for fs,cs in zip(lfs,lcs):
		plt.plot(fs, cs, ".-", label=labels[lfs.index(fs)], color=colors[lfs.index(fs)])
	plt.ylabel(r"$center_combi$")
	plt.xlim(0.0, 0.5)
	plt.ylim(min([min(cs) for cs in lcs]),max([max(cs) for cs in lcs]))
	#plt.axvline(0.5, color="black")
	#plt.axhline(0.03, color="black")
	
	'''
	s = r"$ P = \frac{1}{fN} \sum_i \left( \frac{\sigma_i}{|\Delta t_i|} \right)$" 
	#s = 'P = 1/fN * sum(tderr_i/|td_i|)'
	plt.annotate(s=s, xy = (0.73,0.65), xycoords='axes fraction',
         		textcoords='axes fraction',size=18)
	'''
				
	plt.legend(fontsize= 12,loc=2)
	
	
	plt.subplot(3,1,2)
	for fs,ss in zip(lfs,lss):		
		plt.plot(fs, ss, ".-", label=labels[lfs.index(fs)],color=colors[lfs.index(fs)])				
	plt.ylabel(r"$sigma_combi$")
	plt.xlim(0.0, 0.5)
	plt.ylim(min([min(ss) for ss in lss]),max([max(ss) for ss in lss]))
	#plt.axvline(0.5, color="black")
	#plt.axhline(0.5, color="black")
	#plt.axhline(1.5, color="black")
	'''
	s = r"$ \chi^2 =\frac{1}{fN} \sum_i \left( \frac{\overline{\Delta} t_i - \Delta t_i}{\sigma_i} \right)$"

	plt.annotate(s=s, xy = (0.73,0.75), xycoords='axes fraction',
         		textcoords='axes fraction',size=18)
	'''
	
	plt.subplot(3,1,3)
	for fs,zs in zip(lfs,lzs):
		plt.plot(fs, zs, ".-", label=labels[lfs.index(fs)],color=colors[lfs.index(fs)])

	plt.ylabel(r"$probazero_combi$")
	plt.xlim(0.0, 0.5)
	plt.ylim(min([min(zs) for zs in lzs]),max([max(zs) for zs in lzs]))
	#plt.axvline(0.5, color="black")
	#plt.axhline(0.03, color="black")
	#plt.axhline(-0.03, color="black")
	'''
	s = r"$ A =\frac{1}{fN} \sum_i \left( \frac{\overline{\Delta} t_i - \Delta t_i}{|\Delta t_i|} \right)$"

	plt.annotate(s=s, xy = (0.73,0.75), xycoords='axes fraction',
         		textcoords='axes fraction',size=18)
	
	'''
	
	'''
	plt.subplot(4,1,4)
	for fs,Amods in zip(lfs,lAmods):
		plt.plot(fs, Amods, ".-", label=labels[lfs.index(fs)],color=colors[lfs.index(fs)])
	plt.xlabel(r"$f$")
	plt.ylabel(r"$A_{mod}$")
	plt.xlim(0.0, 0.5)
	plt.ylim(min([min(Amods) for Amods in lAmods]),max([max(Amods) for Amods in lAmods]))
	plt.axvline(0.5, color="black")
	plt.axhline(0.03, color="black")
	plt.axhline(-0.03, color="black")		
	s = r"$ A =\frac{1}{fN} \sum_i \left( \frac{\overline{\Delta} t_i - \Delta t_i}{\Delta t_i} \right)$"

	plt.annotate(s=s, xy = (0.73,0.75), xycoords='axes fraction',
         		textcoords='axes fraction',size=18)
	'''		
	plt.show()




