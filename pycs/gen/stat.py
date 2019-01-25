"""
Statistics related stuff.
"""

import sys
import os
import numpy as np
import math

import pycs.gen.util


def normal(x, mu, sigma):
	return (1.0/np.sqrt(2.0*np.pi*sigma*sigma)) * np.exp( - (x - mu)**2/(2*sigma*sigma))


def sf(l, binsize = 200, ssf=False):
	"""
	Structure function of a lightcurve
	
	ssf gives a 2D density plot, otherwise binned.
	
	For definition see for instance :
	De Vries, W.H. de, Becker, R., White, R., Loomis, C., 2005. Structure Function Analysis of Long-Term Quasar Variability. The Astronomical Journal 129, 615-615-629-629.
	
	"""
	import matplotlib.pyplot as plt
	
	mags = l.getmags()
	jds = l.getjds()
	n = len(l)
	
	#n = 1000
	#jds = np.arange(n)
	#mags = np.random.randn(n)*3.0
	
	
	ja = np.ones((n,n)) * jds
	jam = ja - ja.transpose()
	jam = jam.flatten()
	keep = jam > 0.0
	jam = jam[keep]
	
	ma = np.ones((n,n)) * mags
	mam = ma - ma.transpose()
	mam = mam.flatten()
	mam = mam[keep]

	if ssf: # stochastic structure function, we plot a 2d distribution
		
		jam = jam.flatten()
		mam = mam.flatten()
		
		plt.scatter(jam, mam, s=1.0)
		
		plt.xlabel("Delta t")
		plt.ylabel("Delta m")
		
		plt.show()
		
	
	else: # we do a normal structure function, "variance of increments versus delta t" :
	
		mam = np.square(mam)
	
		order = np.argsort(jam) # sorting according to the jd gaps
		jam = jam[order]
		mam = mam[order]
		
		m = len(jam)
		nbins = int(math.floor(m/binsize))
		jam = jam[0:nbins*binsize] # cutting to the nearest last bin
		mam = mam[0:nbins*binsize]
		
		jam = jam.reshape((nbins, binsize))
		mam = mam.reshape((nbins, binsize))
		
		cjs = np.mean(jam, axis=1)
		cms = np.sqrt(np.mean(mam, axis=1) / float(binsize))
		#cms2 = np.sqrt(np.median(mam, axis=1) / float(binsize))
		
	
		plt.scatter(cjs, cms)
		#plt.scatter(cjs, cms2, c="red")
		
		plt.xlabel("Delta t")
		plt.ylabel("SF")
		
		plt.show()
	
	
	

def erf(x):
	"""
	Error function. There is one in scipy, but this way we do it without scipy...
	
	scipy.special.erf(z)
	
	"""
	# save the sign of x
	sign = 1
	if x < 0: 
		sign = -1
	x = abs(x)

	# constants
	a1 =  0.254829592
	a2 = -0.284496736
	a3 =  1.421413741
	a4 = -1.453152027
	a5 =  1.061405429
	p  =  0.3275911

	# A&S formula 7.1.26
	t = 1.0/(1.0 + p*x)
	y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)
	return sign*y # erf(-x) = -erf(x)
	


def runstest(residuals, autolevel=False, verbose=True):
	"""
	One-sample runs test of randomness as presented in Practical Statistics for Astronomers by
	J. V. Wall and C. R. Jenkins, paragraph 5.3.3
	WARNING : ERROR IN THE BOOKS EXPECTATION FORMULA ! confirmed by author
	
	residuals is a numpy array of floats
	
	returns z (in units of sigmas)
	and p = assuming that the data are independent, the probability to get a result worse then what you got.
	
	"""

	medianlevel = np.median(residuals)
	if autolevel:
		if verbose:
			print "Leveling to median : %f" % (medianlevel)
		residuals -= medianlevel
	
	bools = residuals > 0 # So residuals = 0 would be set to False, but
	bools = bools[np.abs(residuals) > 0.000001] # we remove those very close to 0
	
	if verbose:
		print "Total : %i points / significant : %i points" % (len(residuals), len(bools))
	
	n = len(bools)
	if n <= 20:
		print "WARNING : to few points for a meaningful runs test (or too close to 0.0)"
	
	intbools = bools.astype(int)
	nplus = np.sum(bools)
	nminus = n - nplus
	
	# And the number of runs
	changes = np.abs(intbools[1:] - intbools[:-1])
	nruns = np.sum(changes) + 1
	
	if verbose:
		print "     + (m) : %i" % nplus
		print "     - (n) : %i" % nminus
		print "  Runs (r) : %i" % nruns
	
	# For large N, we can approximate the distribution of r by a gaussian :
	# Error from the book :
	#mur = (2.0 * nplus * nminus)/(nplus + nminus + 1.0)
	#sigmar = math.sqrt((2.0*nplus*nminus*(2.0*nplus*nminus - n))/(n*n*(n - 1.0)))
	# From corrected book and also wikipedia :
	mur = ((2.0 * nplus * nminus)/n) + 1.0
	sigmar = math.sqrt((mur-1.0)*(mur-2.0) / (n - 1.0))
	
	zruns = (nruns - mur) / sigmar
	
	# The probability to obtain a result worse than what you got :
	pruns = 1.0 - erf(abs(zruns)/math.sqrt(2.0))
	
	return {"zruns":zruns, "pruns":pruns, "nruns":nruns}
	

# def maprunstest(rl, seasons = None, autolevel=False, verbose=True):
# 	"""
# 	maps the function runstest above on the seasons of a residual lightcurve rl
# 	
# 	I DO NOT TAKE INTO ACCOUNT THE MASK, PLEASE CUT IT FIRST !
# 	AND THEN RECALCULATE YOUR SEASONS !!!
# 	
# 	>>> res = pycs.gen.stat.maprunstest(lca, seasons, autolevel=False, verbose = False)
# 	>>> out = [(season.name, res["z"]) for (season, res) in zip (seasons, res)]
# 
# 	
# 	"""
# 	
# 	if lc.hasmask():
# 		print "WARNING : I do not take into account the mask !"
# 	
# 	residualslist = []
# 	if seasons == None:
# 		residualslist.append(lc.getmags())
# 	else:
# 		mags = lc.getmags()
# 		for season in seasons:
# 			residuals = mags[season.indices]
# 			residualslist.append(residuals)
# 			
# 	return map(lambda x: runstest(x, autolevel=autolevel, verbose=verbose), residualslist)
			
			
			

def subtract(lcs, spline):
	"""
	I return a list of residual light curves ("lcs - spline").
	Technically, "residual" light curves are nothing but normal lightcurves objects.
	Of course, I take into account any shifts or microlensing of your lcs.
	I do not modify my input arguments.
	"""
	
	rls = []
	
	for l in lcs:
		
		if l.hasmask():
			print "WARNING : I do not take into account the mask !"
			
		lp = l.copy() # To avoid border effects
		
		lp.applyfluxshift()
		lp.applymagshift()
		if lp.ml != None:
			lp.applyml()

		#print lcs
		lp.mags -= spline.eval(lp.getjds())
			
		rls.append(lp)
		
	return rls



def resistats(rl):
	"""
	Give me a residual lightcurve, I return a dict with some descriptive stats about its magnitudes.
	"""
	
	meanmag = np.mean(rl.getmags())
	stdmag = np.std(rl.getmags())
	runs = runstest(rl.getmags(), autolevel=False, verbose=False)
	
	out = {"mean":meanmag, "std":stdmag}
	out.update(runs)
	
	return out
	
def mapresistats(rls):
	return [resistats(rl) for rl in rls]



def anaoptdrawn(optoriglcs, optorigspline, simset="simset", optset="optset", npkl=1000, plots=True, nplots=3, r=0.11, plotjdrange=None, plotcurveindexes=None, showplot=False, directory = "./", resihist_figsize = None):
	"""
	Not flexible but very high level function to analyse the spline-fit-residuals of drawn curves and comparing them to the
	real observations.
	This can be used to tune the parameters of the "drawing".
	.. warning:: The simset must have been optimized using spline fits, with option keepopt=True
	
	:param optoriglcs: optimized original curves
	
	:param optorigspline: spline that matches to these curves
	
	:param simset: Put this e.g. to plotcurveindex=(0,2) if you want to plot only the first 2 curves ...
	
	:param plotcurveindexes: allows you to plot only a subset of lcs (smaller plots). Give a tuple like eg (0, 2, 3)
	
	:param npkl: I read only the first npkl picke files.
	
	"""
	import matplotlib.pyplot as plt
	import glob
	
	print "Analysing the residuals of simset %s" % (simset)
	
	# For each light curve we make a dict that we will use to store stuff
	curves = [{"optoriglc":optoriglc} for optoriglc in optoriglcs]
	
	# We compute the residuals of the original curve
	optorigrlcs = subtract(optoriglcs, optorigspline)
	for (curve, optorigrlc) in zip(curves, optorigrlcs):
		curve["optorigrlc"] = optorigrlc
	
	# We read all the optimized mock curves :
	pkls = sorted(glob.glob(directory + "sims_%s_opt_%s/*_opt.pkl" % (simset, optset)))
	print directory + "sims_%s_opt_%s/*_opt.pkl" % (simset, optset)


	optmocksplinelist = []
	optmocklcslist = []

	for (i, pkl) in enumerate(pkls):
		if i >= npkl:
			break
		opttweak = pycs.gen.util.readpickle(pkl, verbose=False)
		optmocksplinelist.extend(opttweak["optfctoutlist"])
		optmocklcslist.extend(opttweak["optlcslist"])

	assert len(optmocksplinelist) == len(optmocklcslist)
	
	print "We have %i simulated curves" % (len(optmocksplinelist))

	# We compute all the residuals of the mock curves, and store them
	for curve in curves:
		curve["optmockrlclist"] = []
	
	for (optmocklcs, optmockspline) in zip(optmocklcslist, optmocksplinelist):
		assert len(optmocklcs) == len(optoriglcs)
		optmockrlcs = pycs.gen.stat.subtract(optmocklcs, optmockspline)
		
		for (curve, optmockrlc) in zip(curves, optmockrlcs):
			assert curve["optorigrlc"].object == optmockrlc.object
			curve["optmockrlclist"].append(optmockrlc)
		
	
	# We want to return the displayed statistics
	stats = []
	for curve in curves:
		curve["origresistats"] = resistats(curve["optorigrlc"])
		
		curve["mockresistats"] = map(resistats, curve["optmockrlclist"])
		curve["meanmockresistats"] = dict([[key, np.mean(np.array([el[key] for el in curve["mockresistats"]]))] for key in curve["origresistats"].keys()])
		curve["medmockresistats"] = dict([[key, np.median(np.array([el[key] for el in curve["mockresistats"]]))] for key in curve["origresistats"].keys()])
		curve["stdmockresistats"] = dict([[key, np.std(np.array([el[key] for el in curve["mockresistats"]]))] for key in curve["origresistats"].keys()])
	
	
		print "++++++ %s ++++++" % (curve["optorigrlc"].object)
		curve["zrunstxt"] = "zruns : %.2f (obs) vs %.2f +/- %.2f (sim)" % (curve["origresistats"]["zruns"], curve["meanmockresistats"]["zruns"], curve["stdmockresistats"]["zruns"])
		curve["sigmatxt"] = "sigma : %.4f (obs) vs %.4f +/- %.4f (sim)" % (curve["origresistats"]["std"], curve["meanmockresistats"]["std"], curve["stdmockresistats"]["std"])
		print curve["zrunstxt"]
		print curve["sigmatxt"]

		# return the original, mean and std of mocks zruns, then original, mean and std of mocks of sigma
		stats.append([curve["origresistats"]["zruns"], curve["meanmockresistats"]["zruns"], curve["stdmockresistats"]["zruns"], curve["origresistats"]["std"], curve["meanmockresistats"]["std"], curve["stdmockresistats"]["std"]])
	
		#for item in curve["mockresistats"]:
		#	print item["zruns"]



	# Now we proceed with making plots.
	
	# The plots of the residuals statistics, one for each curve :
	"""
	if plots:
		for curve in curves:
			
			plt.figure(figsize=(12, 6))
			
			plt.figtext(0.03, 0.5, curve["optorigrlc"].object, fontsize=30)
			
			# Histo of zruns
			plt.subplot(121)
			plt.hist(np.array([el["zruns"] for el in curve["mockresistats"]]), 30, color="gray")
			plt.axvline(curve["origresistats"]["zruns"], color="green", linewidth=3.0)
			plt.title(curve["zrunstxt"])
			plt.xlabel("zruns")
			plt.xlim(-10.0, 1.0)
			
			# Histo of residuals
			r = 0.1
			plt.subplot(122)
			plt.hist(np.concatenate([rlc.mags for rlc in curve["optmockrlclist"]]), 50, range=(-r, r), facecolor='gray', normed=True)
			# Gaussian for the mock hist :
			plt.plot(np.linspace(-r, r, 100), normal(np.linspace(-r, r, 100), curve["meanmockresistats"]["mean"], curve["meanmockresistats"]["std"]), color="gray")
			plt.hist(curve["optorigrlc"].mags, 50, facecolor='green', alpha=0.5, range=(-r, r), normed=True)
			plt.xlabel("Residuals [mag]")
			plt.title(curve["sigmatxt"])
			plt.xlim(-r, r)
			
			plt.savefig("anaoptdrawn_%s_%s_%s.pdf" % (simset, optset, curve["optorigrlc"].object))
	"""	
	
	# Resi histos combined into one nicer figure :
	"""
	if plots:
		r = 0.11
		plt.figure(figsize=(3*len(curves), 3))
		plt.subplots_adjust(left=0.03, bottom=0.19, right=0.97, top=None, wspace=None, hspace=None)
		for (i,curve) in enumerate(curves):
			#print (1, len(curves), i+1)
			plt.subplot(1, len(curves), i+1)
			plt.hist(np.concatenate([rlc.mags for rlc in curve["optmockrlclist"]]), 50, range=(-r, r), facecolor='gray', normed=True, histtype="stepfilled")
			# Gaussian for the mock hist :
			#plt.plot(np.linspace(-r, r, 100), normal(np.linspace(-r, r, 100), curve["origresistats"]["mean"], curve["origresistats"]["std"]), color="green")
			plt.hist(curve["optorigrlc"].mags, 50, facecolor='green', alpha=0.4, range=(-r, r), normed=True, histtype="stepfilled")
			plt.xlabel("Spline fit residuals [mag]")
			#plt.title(sigmatxt)
			#print plt.gca().get_ylim()
			plt.text(-r+0.1*r, 0.85*plt.gca().get_ylim()[1], curve["optorigrlc"].object, fontsize=20)
			plt.xlim(-r, r)
			plt.gca().get_yaxis().set_ticks([])

		#plt.show()
		plt.savefig("anaoptdrawn_%s_%s_resihists.pdf" % (simset, optset))	
	"""
	# Resi and zruns histos combined into one nicer figure :
	
	if plots:
		if resihist_figsize == None :
			plt.figure(figsize=(3*len(curves), 4))
		else :
			plt.figure(figsize=resihist_figsize)
		plt.subplots_adjust(left=0.02, bottom=0.12, right=0.98, top=0.98, wspace=0.08, hspace=0.37)
		
		# Resi histos :
		for (i,curve) in enumerate(curves):
			#print (1, len(curves), i+1)
			plt.subplot(2, len(curves), i+1)
			plt.hist(np.concatenate([rlc.mags for rlc in curve["optmockrlclist"]]), 50, range=(-r, r), facecolor='black', alpha=0.4, normed=True, histtype="stepfilled")
			# Gaussian for the mock hist :
			#plt.plot(np.linspace(-r, r, 100), normal(np.linspace(-r, r, 100), curve["origresistats"]["mean"], curve["origresistats"]["std"]), color="green")
			plt.hist(curve["optorigrlc"].mags, 50, facecolor='green', alpha=0.4, range=(-r, r), normed=True, histtype="stepfilled")
			plt.xlabel("Spline fit residuals [mag]")
			
			#print plt.gca().get_ylim()
			plt.text(-r+0.1*r, 0.8*plt.gca().get_ylim()[1], curve["optorigrlc"].object, fontsize=18)
			plt.xlim(-r, r)
			plt.gca().get_yaxis().set_ticks([])

		# zruns histos :
		for (i,curve) in enumerate(curves):
			#print (1, len(curves), i+1)
			plt.subplot(2, len(curves), len(curves)+i+1)
			
			plt.hist(np.array([el["zruns"] for el in curve["mockresistats"]]), 20, facecolor="black", alpha=0.4, normed=True, histtype="stepfilled")
			plt.axvline(curve["origresistats"]["zruns"], color="green", linewidth=2.0, alpha=0.7)
			
			plt.xlabel(r"$z_{\mathrm{r}}$", fontsize=18)
			# plt.xlim(-5.0, 5.0)

			#plt.text(-9.0, 0.85*plt.gca().get_ylim()[1], curve["optorigrlc"].object, fontsize=20)
			plt.gca().get_yaxis().set_ticks([])

		if showplot:
			plt.show()
		plt.savefig("fig_anaoptdrawn_%s_%s_resihists.png" % (simset, optset))
	
	
	# A detailed plot of some residuals, just for a few drawn curves
	
	if plots:		
		for i in range(nplots):
			
			optmockrlcs = [curve["optmockrlclist"][i] for curve in curves]
			for l in optmockrlcs:
				l.plotcolour = "black"
			
			optorigrlcs = [curve["optorigrlc"] for curve in curves]
			
			if plotcurveindexes != None:
				optorigrlcs = [optorigrlcs[index] for index in plotcurveindexes]
				optmockrlcs = [optmockrlcs[index] for index in plotcurveindexes]
			plotresiduals([optorigrlcs, optmockrlcs], jdrange=plotjdrange, nicelabel=False, showlegend=False, showsigmalines = False, errorbarcolour = "#999999", filename="fig_anaoptdrawn_%s_%s_resi_%i.png" % (simset, optset, i+1))

	
	return stats
	
	

def plotresiduals(rlslist, jdrange=None, magrad=0.1, errorbarcolour = "#BBBBBB", 
				  showerrorbars=True, showlegend=True, nicelabel=True, 
				  showsigmalines=True, filename=None, ax = None):
	"""
	We plot the residual lightcurves in separate frames.
	
	The arguement rlslist is a *list* of *lists* of lightcurve objects.
	Ths sublists should have the same length, I'll choose my number of panels accordingly.
	The structure is : [[lca, lcb], [lca_sim1, lcb_sim1], ...]
	If you have only one lightcurve object, you can of course pass [[l]] ...
	
	:param rlslist: 
	
	I disregard the timeshift of the curves !
	"""
	import matplotlib.pyplot as plt
	import matplotlib.ticker
	
	minorLocator = matplotlib.ticker.MultipleLocator(50)
	majorLocator = matplotlib.ticker.MultipleLocator(200)
	
	
	resminorLocator = matplotlib.ticker.MultipleLocator(0.01)
	resmajorLocator = matplotlib.ticker.MultipleLocator(0.05)
	
	eps = 0.001
	
	npanels = len(rlslist[0])
	if ax == None:
		fig = plt.figure(figsize=(12,1.6*npanels))	# sets figure size
		fig.subplots_adjust(left=0.07, right=0.99, top=0.95, bottom=0.14, hspace=0.05)
		ax = plt.gca()
		ihaveax = False
	else :
		ihaveax = True

	
	# fig = plt.figure(figsize=(12,1.6*npanels))
	# fig.subplots_adjust(left = 0.07, right=0.99, top=0.95, bottom=0.14, hspace=0.05)


	#plt.rc('font', family = 'serif', serif = 'STIXGeneral')

	
	for i in range(npanels): # i is the panel index
		
		rls = [rlslist[j][i] for j in range(len(rlslist))] # j is the curve index.
		
		if ihaveax : 
			ax0 = ax
		else : 
			if i > 0:
				ax = plt.subplot(npanels, 1, i+1, sharex=ax0, sharey=ax0)
			else:
				ax = plt.subplot(npanels, 1, i+1)
				ax0 = ax
			
		for (j, rl) in enumerate(rls):
		
			stats = resistats(rl)
			
			label = "[%s/%s] (std: %.4f, zruns : %.3f)" % (rl.telescopename, rl.object, stats["std"], stats["zruns"])
			#print label
			if nicelabel:
				label = "%s" % (rl.object)
				#label = "%s (std: %.4f, zruns : %.3f)" % (rl.object, stats["std"], stats["zruns"])
			
			if showerrorbars:
				ax.errorbar(rl.jds, rl.getmags(), rl.magerrs, fmt=".", color=rl.plotcolour, markeredgecolor=rl.plotcolour, ecolor=errorbarcolour, label=label, elinewidth=0.5)
			else:
				ax.plot(rl.jds, rl.getmags(), marker=".", markersize=3.0, linestyle="None", markeredgecolor=rl.plotcolour, color=rl.plotcolour, label=label)
			
			
			if showsigmalines:
				ax.axhline(y = stats["std"], lw=0.5, color=rl.plotcolour)
				ax.axhline(y = -stats["std"], lw=0.5, color=rl.plotcolour)
			
			if nicelabel:
				ax.text(0.04 + (0.087 * j), 0.82 , label, transform=ax.transAxes, color = rl.plotcolour)
			else:
				if not showlegend:
					if j == 0:
						ax.text(0.01 , 0.81 , rl.object, transform=ax.transAxes, color = rl.plotcolour, fontsize=17)
			
			
		ax.axhline(0, color="gray", dashes=(3,3))
		
		ax.xaxis.set_minor_locator(minorLocator)
		ax.xaxis.set_major_locator(majorLocator)
		
		ax.yaxis.set_minor_locator(resminorLocator)
		ax.yaxis.set_major_locator(resmajorLocator)
		
		ax.set_ylim(-magrad+eps, magrad-eps)
		ax.set_ylim(ax.get_ylim()[::-1])
		
		if showlegend:
			ax.legend(numpoints = 1, prop={'size':10})
	
		
		#ax.set_ylabel("Residual [mag]")
		
		
		ax.set_xlabel("HJD - 2400000.5 [day]")
		#a.set_xlim(52750, 55400)
		
		if i != npanels-1:
			plt.setp(ax.get_xticklabels(), visible=False)
			ax.set_xlabel("")
	
	ax.text(0.01, 0.5, 'Residuals [mag]', rotation=90, verticalalignment="center", horizontalalignment="center")
		
	if jdrange != None:
		plt.xlim(jdrange[0], jdrange[1])
	else:
		plt.xlim(np.min(rlslist[0][0].jds)-50, np.max(rlslist[0][0].jds)+50)
	
	if filename:
		plt.savefig(filename)	
	else:
		plt.show()
	

	
	
