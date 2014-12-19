"""
Stuff related to the TDC metrics
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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
	Compute P and Perr - std(P)/sqrt(len(P)) - for a given method stored in the database
	"""
	subdb = [item for item in db if "%s_P" %(method) in item]
	#key = "%s_P" %method
	#print 'I RETURN: ',sum([item["%s_P" %(method)] for item in subdb]) / float(len(subdb))
	Ps = [item["%s_P" %(method)] for item in subdb]
	nest = float(len(subdb))
	Perr = np.std(Ps)/np.sqrt(nest)
	P = sum(Ps)/nest
	return (P, Perr)

	
def getA(db,method):
	"""
	Compute A and Aerr - std(A)/sqrt(len(A)) - for a given method stored in the database
	"""	
	subdb = [item for item in db if "%s_A" %(method) in item]
	#print 'I RETURN: ',sum([item["%s_A" %(method)] for item in subdb]) / float(len(subdb))
	As = [item["%s_A" %(method)] for item in subdb]
	nest = float(len(subdb))
	Aerr = np.std(As)/np.sqrt(nest)
	A = sum(As)/nest
	return (A, Aerr)
	
	
def getAmod(db,method):
	"""
	Compute Amod and Amoderr - std(Amod)/sqrt(len(Amod)) - for a given method stored in the database
	"""
	subdb = [item for item in db if "%s_Amod" %(method) in item]
	Amods = [item["%s_Amod" %(method)] for item in subdb]
	nest = float(len(subdb))
	Amoderr = np.std(Amods)/np.sqrt(nest)
	Amod = sum(Amods)/nest
	return (Amod, Amoderr)

	
def getchi2(db,method):
	"""
	Compute chi2 and chi2err - std(chi2)/sqrt(len(chi2)) - for a given method stored in the database
	"""	
	subdb = [item for item in db if "%s_chi2" %(method) in item]
	chi2s = [item["%s_chi2" %(method)] for item in subdb]
	nest = float(len(subdb))
	chi2err = np.std(chi2s)/np.sqrt(nest)
	chi2 = sum(chi2s)/nest
	return (chi2, chi2err)


def getf(db,method,N):
	"""
	Compute f for a given method stored in the database
	N is the total number of curves	
	"""
	
	subdb = [item for item in db if "%s_td" %(method) in item]	
	return  float(len(subdb))/N


def Pplotall(db,methods,N,zoomchi2=False,zoomA=False):
	'''
	give me the db and a list of methods [method1,method2,...], I plot P vs f for each method
	and chi2 vs f, A vs f... for the same arrangement according to P

	If as a method_i you give me a tuple (method_i,chi2max), I will remove the values wich chi2 > chi2max from the db
	'''

	"""
	Plot fine tuning:

	# spl-vanilla-dou-XX
	fmaxplot = 0.32
	Pxy = (0.60,0.25)
	chi2xy = (0.60,0.65)
	Axy = (0.6,0.75)

	# spl-XX-dou-full
	fmaxplot = 0.35
	Pxy = (0.60,0.22)
	chi2xy = (0.60,0.68)
	Axy = (0.6,0.78)

	# XX-vanilla-dou-full
	fmaxplot = 0.32
	Pxy = (0.60,0.28)
	chi2xy = (0.60,0.73)
	Axy = (0.6,0.75)

	# d3cs-vanilla-XX-full
	fmaxplot = 0.8
	Pxy = (0.65,0.28)
	chi2xy = (0.65,0.73)
	Axy = (0.65,0.75)

	# spl-vanilla-XX-full
	fmaxplot = 0.67
	Pxy = (0.60,0.20)
	chi2xy = (0.60,0.65)
	Axy = (0.6,0.39)
	"""
	fmaxplot = 0.32
	Pxy = (0.60,0.20)
	chi2xy = (0.60,0.65)
	Axy = (0.6,0.39)

	colors = ["blue","green","red","chartreuse","crimson","black","magenta"] # color of the curves

	lfs =[]
	lPs=[]
	lPerrs=[]
	lAs=[]
	lAerrs=[]
	lAmods=[]
	lAmoderrs=[]
	lchi2s=[]
	lchi2errs=[]
	labels = []

	for method in methods:

		addcut = False
		# shoot the huge values in chi2
		if len(method) == 2:

			#db_P_reduced = [item for item in db if "%s_P" %method[0] in item and item["%s_chi2" %method[0]] < method[1]]
			db_P_reduced = [item for item in db if "%s_P" %method[0] in item and abs(item["%s_td" %method[0]]-item["truetd"]) < method[1]]
			sorted_db_P = sorted(db_P_reduced, key = lambda item: item["%s_P" % method[0]] )[::-1]
			method = method[0]
			print 'I shoot %i estimates...' %(len(db)-len(db_P_reduced))
			addcut = True

		else:
			db_P = [item for item in db if "%s_P" %(method) in item]
			sorted_db_P = sorted(db_P, key = lambda item: item["%s_P" % method] )[::-1]


		fs = []
		Ps = []
		Perrs = []
		As = []
		Aerrs = []
		Amods = []
		Amoderrs = []
		chi2s = []
		chi2errs = []
		for n in range(len(sorted_db_P)):
			subdb = sorted_db_P[n:]
			fs.append(getf(subdb,method=method,N=N))
			P = getP(subdb,method=method)
			Ps.append(P[0])
			Perrs.append(P[1])
			A = getA(subdb,method=method)
			As.append(A[0])
			Aerrs.append(A[1])
			Amod = getAmod(subdb,method=method)
			Amods.append(Amod[0])
			Amoderrs.append(Amod[1])
			chi2 = getchi2(subdb,method=method)
			chi2s.append(chi2[0])
			chi2errs.append(chi2[1])
		
		lfs.append(fs)
		lPs.append(Ps)
		lPerrs.append(Perrs)
		lAs.append(As)
		lAerrs.append(Aerrs)
		lAmods.append(Amods)
		lAmoderrs.append(Amoderrs)
		lchi2s.append(chi2s)
		lchi2errs.append(chi2errs)
		if addcut:
			method += '-cut'
		labels.append(method.split('_')[-1])



	# And now, the plot

	#ax.fill_between(x, y1, y2, where=y2>=y1, facecolor='green', interpolate=True) TO FILL BETWEEN CURVES


	fig = plt.figure('metrics vs f',(7.5,8.5))

	gs1 = gridspec.GridSpec(12, 1)
	gs1.update(left=0.141, right=0.99, top=0.98, bottom=0.09,hspace =0.1) # stupid hspace is not the same between gridspec and subplots adjust ! Crap !

	ax1 = plt.subplot(gs1[:3,0])

	if zoomchi2 == True:
		ax2 = plt.subplot(gs1[3:5,0])
		ax2s = plt.subplot(gs1[5,0])
	else:
		ax2 = plt.subplot(gs1[3:6,0])

	if zoomA == True:
		ax3 = plt.subplot(gs1[7:9,0])
		ax3s = plt.subplot(gs1[6,0])
	else:
		ax3 = plt.subplot(gs1[6:9,0])

	ax4 = plt.subplot(gs1[9:,0])



	ax1.xaxis.set_ticks([])
	ax1.axvline(0.5, color="black",linewidth=3,alpha=0.3)
	ax1.axhline(0.03, color="black",linewidth=3,alpha=0.3)
	for fs,Ps,Perrs in zip(lfs,lPs,lPerrs):
		errup = [P+Perr for P,Perr in zip(Ps,Perrs)]
		errdown = [P-Perr for P,Perr in zip(Ps,Perrs)]
		ax1.plot(fs, Ps, ".-", label=labels[lPs.index(Ps)], color=colors[lPs.index(Ps)],linewidth=3,ms=2)

		ax1.fill_between(fs, errup, errdown, facecolor=colors[lPs.index(Ps)], interpolate=True, alpha=0.3)

		print lPs.index(Ps)

	ax1.set_ylabel(r"$P$",fontsize=20)
	ax1.set_xlim(0.0, fmaxplot)
	ax1.set_ylim(min([min(Ps) for Ps in lPs]),max([max(Ps) for Ps in lPs]))

	ax1.legend(fontsize= 15,loc=2,framealpha=0)
	

	ax2.xaxis.set_ticks([])
	ax2.axvline(0.5, color="black",linewidth=3,alpha=0.3)
	ax2.axhline(0.5, color="black",linewidth=3,alpha=0.3)
	ax2.axhline(1.5, color="black",linewidth=3,alpha=0.3)
	ax2.axhline(1, color="black",linewidth=3,alpha=0.3,linestyle='--')
	for fs,chi2s,chi2errs in zip(lfs,lchi2s,lchi2errs):
		errup = [chi2+chi2err for chi2,chi2err in zip(chi2s,chi2errs)]
		errdown = [chi2-chi2err for chi2,chi2err in zip(chi2s,chi2errs)]
		ax2.plot(fs, chi2s, ".-", label=labels[lchi2s.index(chi2s)],color=colors[lchi2s.index(chi2s)],linewidth=3,ms=2)
		ax2.fill_between(fs, errup, errdown, facecolor=colors[lchi2s.index(chi2s)], interpolate=True, alpha=0.3)

	ax2.set_ylabel(r"$\chi^2$",fontsize=20) # to tweak label position : y=0.2,labelpad=+25
	ax2.set_xlim(0.0, fmaxplot)
	ax2.set_ylim(min([min(chi2s) for chi2s in lchi2s]),max([max(chi2s) for chi2s in lchi2s]))
	#plt.legend(fontsize= 15,loc=2)


	if zoomchi2 == True:
		ax2s.xaxis.set_ticks([])
		ax2s.axvline(0.5, color="black",linewidth=3,alpha=0.3)
		ax2s.axhline(0.5, color="black",linewidth=3,alpha=0.3)
		ax2s.axhline(1.5, color="black",linewidth=3,alpha=0.3)
		ax2s.axhline(1, color="black",linewidth=3,alpha=0.3,linestyle='--')
		for fs,chi2s,chi2errs in zip(lfs,lchi2s,lchi2errs):
			errup = [chi2+chi2err for chi2,chi2err in zip(chi2s,chi2errs)]
			errdown = [chi2-chi2err for chi2,chi2err in zip(chi2s,chi2errs)]
			ax2s.plot(fs, chi2s, ".-", label=labels[lchi2s.index(chi2s)],color=colors[lchi2s.index(chi2s)],linewidth=3,ms=2)
			ax2s.fill_between(fs, errup, errdown, facecolor=colors[lchi2s.index(chi2s)], interpolate=True, alpha=0.3)
		#ax2s.set_ylabel(r"$\chi^2$",fontsize=20)
		ax2s.set_xlim(0.0, fmaxplot)
		ax2s.set_ylim(0.5,1.0)
		ax2s.set_yticks([0.5, 0.75, 1])

		#plt.legend(fontsize= 15,loc=2)

		ax2.set_ylim(1.0,max([max(chi2s) for chi2s in lchi2s]))



	ax3.xaxis.set_ticks([])
	ax3.axvline(0.5, color="black",linewidth=3,alpha=0.3)
	ax3.axhline(0.03, color="black",linewidth=3,alpha=0.3)
	ax3.axhline(-0.03, color="black",linewidth=3,alpha=0.3)
	ax3.axhline(0, color="black",linewidth=3,alpha=0.3,linestyle='--')
	for fs,As,Aerrs in zip(lfs,lAs,lAerrs):
		errup = [A+Aerr for A,Aerr in zip(As,Aerrs)]
		errdown = [A-Aerr for A,Aerr in zip(As,Aerrs)]
		ax3.plot(fs, As, ".-", label=labels[lAs.index(As)],color=colors[lAs.index(As)],linewidth=3,ms=2)
		ax3.fill_between(fs, errup, errdown, facecolor=colors[lAs.index(As)], interpolate=False, alpha=0.3)

	ax3.set_ylabel(r"$A$",fontsize=20,y=0.8)
	ax3.yaxis.set_label_coords(-0.11, 0.5)
	#ax3.set_xlabel(r"$f$",fontsize=20)
	ax3.set_xlim(0.0, fmaxplot)
	ax3.set_ylim(min([min(As) for As in lAs]),max([max(As) for As in lAs]))

	if zoomA == True:
		ax3s.xaxis.set_ticks([])
		ax3s.axvline(0.5, color="black",linewidth=3,alpha=0.3)
		ax3s.axhline(0.03, color="black",linewidth=3,alpha=0.3)
		ax3s.axhline(-0.03, color="black",linewidth=3,alpha=0.3)
		ax3s.axhline(0, color="black",linewidth=3,alpha=0.3,linestyle='--')
		for fs,As,Aerrs in zip(lfs,lAs,lAerrs):
			errup = [A+Aerr for A,Aerr in zip(As,Aerrs)]
			errdown = [A-Aerr for A,Aerr in zip(As,Aerrs)]
			ax3s.plot(fs, As, ".-", label=labels[lAs.index(As)],color=colors[lAs.index(As)],linewidth=3,ms=2)
			ax3s.fill_between(fs, errup, errdown, facecolor=colors[lAs.index(As)], interpolate=False, alpha=0.3)

		#ax3s.set_ylabel(r"$A$",fontsize=20,y=0.2,labelpad=+25)
		ax3s.yaxis.set_label_coords(-0.11, 0.5)
		#ax3s.set_xlabel(r"$f$",fontsize=20)
		ax3s.set_xlim(0.0, fmaxplot)
		ax3s.set_ylim(-0.0015,0.004)
		ax3s.set_yticks([-0.001, 0.002])
		ax3.set_ylim(min([min(As) for As in lAs]),-0.001)



	ax4.axvline(0.5, color="black",linewidth=3,alpha=0.3)
	ax4.axhline(0.03, color="black",linewidth=3,alpha=0.3)
	ax4.axhline(-0.03, color="black",linewidth=3,alpha=0.3)
	ax4.axhline(0, color="black",linewidth=3,alpha=0.3,linestyle='--')
	for fs,Amods,Amoderrs in zip(lfs,lAmods,lAmoderrs):
		errup = [Amod+Amoderr for Amod,Amoderr in zip(Amods,Amoderrs)]
		errdown = [Amod-Amoderr for Amod,Amoderr in zip(Amods,Amoderrs)]
		ax4.plot(fs, Amods, ".-", label=labels[lAmods.index(Amods)],color=colors[lAmods.index(Amods)],linewidth=3,ms=2)
		ax4.fill_between(fs, errup, errdown, facecolor=colors[lAmods.index(Amods)], interpolate=False, alpha=0.3)

	ax4.set_ylabel(r"$Amod$",fontsize=20)
	ax4.yaxis.set_label_coords(-0.11, 0.5)
	ax4.set_xlabel(r"$f$",fontsize=20)
	ax4.set_xlim(0.0, fmaxplot)
	ax4.set_ylim(min([min(Amods) for Amods in lAmods]),max([max(Amods) for Amods in lAmods]))





	
	for ax in [ax1,ax2,ax3,ax4]:
		for tick in ax.xaxis.get_major_ticks():
  			tick.label.set_fontsize(15)
		for tick in ax.yaxis.get_major_ticks():
  			tick.label.set_fontsize(15)
		from matplotlib.ticker import MaxNLocator
		ax.yaxis.set_major_locator(MaxNLocator(5))

	if zoomchi2 == True:
		for tick in ax2s.xaxis.get_major_ticks():
  			tick.label.set_fontsize(15)
		for tick in ax2s.yaxis.get_major_ticks():
  			tick.label.set_fontsize(15)
	if zoomA == True:
		for tick in ax3s.xaxis.get_major_ticks():
  			tick.label.set_fontsize(15)
		for tick in ax3s.yaxis.get_major_ticks():
  			tick.label.set_fontsize(15)

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




