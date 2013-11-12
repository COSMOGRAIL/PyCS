"""
Subpackage with functions to plot all kind of results from runs.
"""


import numpy as np
import math

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import scipy.ndimage

#import pycs.sim.run



def normal(x, mu, sigma):
	"""
	Plain normal distribution.
	You can directly apply me on numpy arrays x, mu, sigma.
	"""
	
	return (1.0/np.sqrt(2.0*np.pi*sigma*sigma)) * np.exp( - (x - mu)**2/(2*sigma*sigma))


def hists(rrlist, r=10.0, nbins=100, showqs=True, showallqs=False, qsrange=None, title=None, niceplot=False, displaytext=True, figsize=(16, 9), left = 0.06, right=0.95, bottom=0.065, top=0.95, wspace=0.2, hspace=0.2, filename=None):
	"""
	Comparing the delay distributions from different run result objects.
	
	:param rrlist: a list of runresults object.		
	:param r: a range radius for the hists
	
	:param showqs: If True, I overplot the qs as scatter points.
	
	"""
		
	n = rrlist[0].nimages()
	
	labels = rrlist[0].labels
	
	# To get some fixed ranges for the histograms, we will use the center of the histos :
	#reftrueshifts = (1.0/len(rrlist)) * np.sum(np.array([rr.getts()["center"] for rr in rrlist]), axis=0)
	reftrueshifts = 0.5 * (np.max(np.array([rr.getts()["center"] for rr in rrlist]), axis=0) + np.min(np.array([rr.getts()["center"] for rr in rrlist]), axis=0))
	#print reftrueshifts
	
	for rr in rrlist:
		if rr.labels != labels:
			raise RuntimeError("Don't ask me to overplot runresults of different curves !")
		#if not np.all(rr.gettruets()["center"] == reftrueshifts):
		#	print "Warning : I use the trueshift of the first rr to set the ranges."
		
		rr.trues = rr.gettruets() # To store this, avoids calculating it again and again.
	
	
	fig = plt.figure(figsize=figsize)
	#fig.subplots_adjust(left = 0.03, right=0.95, bottom=0.05, top=0.95, wspace=0.2, hspace=0.2)
	#Looks good :
	#fig.subplots_adjust(left = 0.06, right=0.95, bottom=0.065, top=0.95, wspace=0.2, hspace=0.2)
	fig.subplots_adjust(left = left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)

	axisNum = 0
	for i in range(n): # [A, B, C, D]
		for j in range(n):
			
			if (i == 0) or (j == n-1) :
				continue # No plot
				
			
			axisNum += 1
			if j >= i:
				continue
			
			ax = plt.subplot(n-1, n-1, axisNum)
			#ax = plt.gca()
			if i == n-1:
				plt.xlabel("Delay [day]")
			if showqs:
				axscatter = ax.twinx()
				
						
			# Hide the y ticks :
			#ax.get_yaxis().set_ticks([])
			
			# Ranges to plot
			reftruedelay = reftrueshifts[i] - reftrueshifts[j]
			histrange = (reftruedelay - r, reftruedelay + r)
			
			# Delay label
			delaylabel="%s%s" % (labels[j], labels[i])
			
			for irr, rr in enumerate(rrlist):
				# We will express the delays "i - j"
				if rr.plottrue == True:
					delays = rr.truetsarray[:,i] - rr.truetsarray[:,j]
				else:
					delays = rr.tsarray[:,i] - rr.tsarray[:,j]
					
				#(counts, bins, patches) = ax.hist(delays, bins=nbins, range=histrange, histtype="step", color=colours[irr % len(colours)], normed=True)
				(counts, bins, patches) = ax.hist(delays, bins=nbins, range=histrange, histtype="bar", color=rr.plotcolour, alpha = 0.4, lw=0, normed=True)
				
				if niceplot:
					majorLocator = MultipleLocator(1)
					minorLocator = MultipleLocator(0.5)
					ax.xaxis.set_major_locator(majorLocator)
					ax.xaxis.set_minor_locator(minorLocator)
					ax.yaxis.set_ticks([])
				
				if showqs and not rr.plottrue :
					if showallqs:
						axscatter.scatter(delays, rr.qs, s=1, facecolor=rr.plotcolour, lw = 0)
					else:
						axscatter.scatter(delays[::5], rr.qs[::5], s=1, facecolor=rr.plotcolour, lw = 0)
						#cmap = colors.LinearSegmentedColormap.from_list('custom',['white', rr.plotcolour],gamma=1.0)
						#axscatter.hexbin(delays, rr.qs, gridsize=(5, 2), mincnt=1, cmap=cmap, edgecolor="none")
						# extent=(histrange[0], histrange[1], -r, r)
					
					if qsrange:
						axscatter.set_ylim(qsrange)
					if niceplot:
						majorLocator = MultipleLocator(500)
						axscatter.yaxis.set_major_locator(majorLocator)
						if axisNum == 1:
							axscatter.set_ylabel(r"$\chi^2$", fontsize=18)
					
				# We plot the true shifts (if available) as a straight line individually for each rr :
				if rr.trues["type"] == "same":
					truedelay = rr.trues["center"][i] - rr.trues["center"][j]
					plt.axvline(x=truedelay, linewidth=1, linestyle="--", color=rr.plotcolour)
				
				# We compute and display the mean and std of the hist :
				#if getattr(rr, "plotgauss", False) == True:
				if displaytext == True:
					meddelay = np.median(delays)
					meandelay = np.mean(delays)
					stddelay = np.std(delays)
					
					if getattr(rr, "plotgauss", False) == True:
						x = np.linspace(histrange[0], histrange[1], 100)
						y = normal(x, meandelay, stddelay)
						ax.plot(x, y, linestyle="-", color = rr.plotcolour)
					
					delaytext = r"%+.1f $\pm$ %.1f" % (meddelay, stddelay)
					#print rr.name
					#print delaylabel + "   " + delaytext
					#ax.text(meddelay, np.max(y)/2.0, "%.1f +/- %.1f" % (meddelay, stddelay), horizontalalignment = "center", color = rr.plotcolour)
					ax.annotate(delaytext, xy=(0.04, 0.75 - 0.08*irr),  xycoords='axes fraction', color = rr.plotcolour)
					
					
			plt.xlim(histrange)
			# Looked ok on big plots :
			#plt.annotate(delaylabel, xy=(0.03, 0.88),  xycoords='axes fraction', fontsize=12, color="black")
			plt.annotate(delaylabel, xy=(0.05, 0.84),  xycoords='axes fraction', fontsize=12, color="black")
			
		
		
	for irr, rr in enumerate(rrlist):
	
		if niceplot:
			labeltxt = "%s" % (getattr(rr, 'name', 'NoName'))
			plt.figtext(x = 0.95, y = 0.93 - 0.05*irr, s = labeltxt, horizontalalignment="right", color=rr.plotcolour, fontsize=16)	
		else:
			labeltxt = "%s (%s, %i) " % (getattr(rr, 'name', 'NoName'), "Truth" if rr.plottrue else "Measured", rr.tsarray.shape[0])
			plt.figtext(x = 0.95, y = 0.93 - 0.03*irr, s = labeltxt, horizontalalignment="right", color=rr.plotcolour)	

		print 'Plotting "%s"' % labeltxt
		print "     Labels : %s" % (", ".join(rr.labels))
		print "     Median shifts : %s" % (", ".join(["%.2f" % (np.median(rr.tsarray[:,i])) for i in range(len(rr.labels))]))
		print "     Std shifts : %s" % (", ".join(["%.2f" % (np.std(rr.tsarray[:,i])) for i in range(len(rr.labels))]))
	
	
	if title != None:
		plt.figtext(x = 0.5, y = 0.95, s = title, horizontalalignment="center", color="black", fontsize=18)	

	if filename == None:
		plt.show()
	else:
		plt.savefig(filename)



def measvstrue(rrlist, r=10.0, nbins = 20, plotpoints=True, ploterrorbars=False, sidebyside=False, errorrange=None, binclip=False, binclipr=10.0, title=None, figsize=(10, 6), left = 0.06, right=0.97, top=0.99, bottom=0.08, wspace=0.15, hspace=0.3, txtstep=0.04, majorticksstep=2, displayn=True, filename=None):
	"""
	
	Plots measured delays versus true delays
	
	:param r: radius of simulation input delays to plot (x axis range)
	:param nbins: number of bins for the bar plot within this range.
	:param plotpoints: should I plot the points (scatter plot) ?
	:param ploterrorbars: should I add errorbars upon the bar plot ?
	:param sidebyside: should I plot bars side by side, or overplot them ?
	:param errorrange: radius of measurement errors to plot (y axis range)
	:param binclip: should I clip errors larger than binclipr days (catastrophic failures of methods) ?
	:param binclipr: see binclip ...
	
	"""
		
	n = rrlist[0].nimages()
	
	labels = rrlist[0].labels
	
	# To get some fixed ranges for the histograms, we will use the first element of rrlist.
	reftrueshifts = np.round(rrlist[0].gettruets()["center"])

	for rr in rrlist:
		if rr.labels != labels:
			raise RuntimeError("Don't ask me to overplot runresults of different curves !")
		#if not np.all(rr.gettruets()["center"] == reftrueshifts):
		#	print "Warning : I use the trueshift of the first rr to set the ranges."
		
		rr.trues = rr.gettruets() # To store this, avoids calculating it again and again.
	
	
	fig = plt.figure(figsize=figsize)
	fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)
	
	axisNum = 0
	for i in range(n): # [A, B, C, D]
		for j in range(n):
			#print i, j
			if (i == 0) or (j == n-1) :
				continue # No plot
			axisNum += 1
			if j >= i:
				continue
			
			ax = plt.subplot(n-1, n-1, axisNum)
			
			minorLocator = MultipleLocator(1.0)
			majorLocator = MultipleLocator(majorticksstep)
			ax.xaxis.set_minor_locator(minorLocator)
			ax.xaxis.set_major_locator(majorLocator)
		
			
			reftruedelay = reftrueshifts[i] - reftrueshifts[j]
			plotrange = (reftruedelay - r, reftruedelay + r)
			
			# Identity line :
			line = np.linspace(plotrange[0], plotrange[1], 100)
			zeros = np.zeros(100)
			plt.plot(line, zeros, color="black", lw=0.5)
			
			# Delay label
			delaylabel="%s%s" % (labels[j], labels[i])
			
			for irr, rr in enumerate(rrlist): # We go through the different runresult objects
				# We will express the delays "i - j"
				truedelays = rr.truetsarray[:,i] - rr.truetsarray[:,j]
				measdelays = rr.tsarray[:,i] - rr.tsarray[:,j]
				
				resis = measdelays-truedelays
				
				# A simple scatter plot of the residues :
				if plotpoints:
					ax.scatter(truedelays, resis, s=2, facecolor=rr.plotcolour, lw = 0)
				
				# We bin those :
				binlims = np.linspace(plotrange[0], plotrange[1], nbins + 1)
				digitized = np.digitize(truedelays, binlims)
				
				binvals = [resis[digitized == bini] for bini in range(1, len(binlims))]
				binstds = map(np.std, binvals)
				binmedians = map(np.median, binvals)
				binmeans = map(np.mean, binvals)
				
				#print binstds
				#print binmedians
				
				if binclip:
					for (bini, binvalarray) in enumerate(binvals):
						
						#keep = np.logical_and(binvalarray < (binmedians[bini] + 1*binstds[bini]), binvalarray > (binmedians[bini] - 1*binstds[bini]))
						#keep = np.logical_and(binvalarray < np.max(binvalarray), binvalarray > np.min(binvalarray))
						keep = np.logical_and(binvalarray < binclipr, binvalarray > -binclipr)
						if np.sum(keep == False) != 0:
							print "Kicking %i points." % (np.sum(keep == False))
						binvals[bini] = binvalarray[keep]
					binstds = map(np.std, binvals)
					binmedians = map(np.median, binvals)
					binmeans = map(np.mean, binvals)
				
				#binmeans = [np.median(resis[digitized == bini]) for bini in range(1, len(binlims))]
				#binstds = [np.std(resis[digitized == bini]) for bini in range(1, len(binlims))]
				
				width = binlims[1] - binlims[0]
				
				if not sidebyside:
					if ploterrorbars:
						ax.bar(binlims[:-1], binmeans, yerr=binstds, width=width, color=rr.plotcolour, ecolor=rr.plotcolour, edgecolor=rr.plotcolour, alpha = 0.2)
					else:
						ax.bar(binlims[:-1], binmeans, width=width, color=rr.plotcolour, edgecolor=rr.plotcolour, alpha = 0.2)
				else:
					width = width/len(rrlist)
					if ploterrorbars:
						ax.bar(binlims[:-1] + irr*width, binmeans, yerr=binstds, width=width, color=rr.plotcolour, ecolor=rr.plotcolour, edgecolor=rr.plotcolour, alpha = 0.3)
					else:
						ax.bar(binlims[:-1] + irr*width, binmeans, width=width, color=rr.plotcolour, edgecolor=rr.plotcolour, alpha = 0.3)
					
				
				# That's it for the different runresult objects, back to the common stuff for this particular panel :
			
			# on all border plots :	
			#if i == n-1:
			#	plt.xlabel("Synthetic input delay [day]")
			#if j == 0:
			#	plt.ylabel("Delay measurement error [day]")
			
			# Just on 2 plots :
			if i == n-1:
				plt.xlabel("True delay [day]")
			if j == 0 and i == int(math.floor(n/2.0)):
				plt.ylabel("Delay measurement error [day]")
			
			
			plt.xlim(plotrange)
			#plt.ylim(plotrange)
			if errorrange != None:
				plt.ylim((-errorrange, errorrange))
			
			plt.annotate(delaylabel, xy=(0.03, 0.88-txtstep),  xycoords='axes fraction', fontsize=12, color="black")
			
			# That's it for this panel, back to the figure :
			
	for irr, rr in enumerate(rrlist):
	
		if displayn:
			labeltxt = "%s (%i) " % (getattr(rr, 'name', 'NoName'), rr.tsarray.shape[0])
		else:
			labeltxt = "%s" % (getattr(rr, 'name', 'NoName'))
		plt.figtext(x = right, y = top - txtstep*irr, s = labeltxt, verticalalignment="top", horizontalalignment="right", color=rr.plotcolour)	

	if title != None:
		plt.figtext(x = 0.5, y = 0.95, s = title, horizontalalignment="center", color="black", fontsize=18)	

	if filename==None:
		plt.show()
	else:
		plt.savefig(filename)


def covplot(rrlist, showpoints=False, showcontour=False, showdensity=True, bins=50, smoothing=0.0, figsize=(12, 12), left=0.02, right=0.98, bottom=0.02, top=0.98, wspace=0.05, hspace=0.05, r=5.0, title=None, txtstep=0.04, filename=None):
	"""
	Covariance scatter of all measurement errors.
	Give me a single runresults object (from a sim, with known true delays).
	"""
	import scipy.stats
	import matplotlib.colors as colors
	
	nimages = rrlist[0].nimages()
	imginds = np.arange(nimages)
	#nruns = len(rr[0])
	labels = rrlist[0].labels
	
	couplelist = [(i, j) for j in imginds for i in imginds if i > j]
	ncouples = len(couplelist)
	
	fig = plt.figure(figsize=figsize)
	fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)
	
	axisNum = 0
	for ii, i in enumerate(couplelist): # (0, 1), (0, 2) ...
		for jj, j in enumerate(couplelist):
			
			if (ii == 0) or (jj == ncouples-1) :
				continue # No plot
			axisNum += 1
			if jj >= ii:
				continue
			
			#print i, j, axisNum
			ax = plt.subplot(ncouples-1, ncouples-1, axisNum, aspect='equal')
			ax.axhline(0, color="black")
			ax.axvline(0, color="black")
				
			
			for rr in rrlist:
			
				
				#print idelaylabel, " vs ", jdelaylabel
				itruedelays = rr.truetsarray[:,i[0]] - rr.truetsarray[:,i[1]]
				imeasdelays = rr.tsarray[:,i[0]] - rr.tsarray[:,i[1]]
				iresis = imeasdelays - itruedelays
				
				jtruedelays = rr.truetsarray[:,j[0]] - rr.truetsarray[:,j[1]]
				jmeasdelays = rr.tsarray[:,j[0]] - rr.tsarray[:,j[1]]
				jresis = jmeasdelays - jtruedelays
				
				
				if showdensity or "diff" in rr.name:
					cmap = colors.LinearSegmentedColormap.from_list('custom',['white', rr.plotcolour],gamma=1.0)
					#cmap = colors.LinearSegmentedColormap.from_list('custom',[rr.plotcolour, rr.plotcolour],gamma=1.0)
					#cmap._init()
					#alphas = np.abs(np.linspace(0.0, 0.5, cmap.N))
					#cmap._lut[:-3,-1] = alphas
					ax.hexbin(iresis, jresis, gridsize=bins, extent=(-r, r, -r, r), mincnt=1, cmap=cmap, edgecolor="none")
				
				
				if showpoints:
					ax.scatter(iresis, jresis, s=2, facecolor=rr.plotcolour, lw = 0)
					#ax.hexbin(iresis, jresis, gridsize=20, extent=(-r, r, -r, r))
				
				if showcontour:
					"""
					H, xedges, yedges = np.histogram2d(iresis, jresis, range=[[-r,r], [-r,r]], bins=(bins, bins))
					H = H.transpose()
					if smoothing > 0.01:
						H = scipy.ndimage.filters.gaussian_filter(H, smoothing, mode='constant', cval=0.0)		
					extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
					#levels = [np.mean(H), np.max(H)/2.0]
					#levels = [2.0*np.mean(H), 6.0*np.mean(H)]
					#levels = (1.0e4, 1.0e3, 1.0e2, 2.0e1)
					levels = [scipy.stats.scoreatpercentile(H.flatten(), 95.45), scipy.stats.scoreatpercentile(H.flatten(), 68.27)]
					#levels = [scipy.stats.scoreatpercentile(H.flatten(), 68.27)]
					
					cset = ax.contour(H, levels=levels, origin="lower", colors=rr.plotcolour, extent=extent, linewidth=0.5)
					"""
					
					H, xedges, yedges = np.histogram2d(iresis, jresis, range=[[-r,r], [-r,r]], bins=(bins, bins))
					extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
					
					data = np.vstack((iresis, jresis))
					#print data.shape
					kde = scipy.stats.kde.gaussian_kde(data)
					
					grid = np.mgrid[-r:r:1j*bins, -r:r:1j*bins]
					grid_coords = np.append(grid[0].reshape(-1,1),grid[1].reshape(-1,1),axis=1) 
					
					z = kde(grid_coords.T)
					z = z.reshape(bins,bins)
					
					#levels = [scipy.stats.scoreatpercentile(z.flatten(), 95.45)]
					levels = [np.max(z)*0.45]
				
					cset = ax.contour(grid[0], grid[1], z, levels=levels, origin="lower", colors=rr.plotcolour, extent=extent, linewidth=0.5)
					
					
					
					
				
		
			
			idelaylabel="%s%s" % (labels[i[1]], labels[i[0]])
			jdelaylabel="%s%s" % (labels[j[1]], labels[j[0]])
			#ax.set_xlabel(idelaylabel)
			#ax.set_ylabel(jdelaylabel)
			
			if figsize[0] > 8:
				ax.annotate(idelaylabel, xy=(0.9, 0.05),  xycoords='axes fraction', ha="center") # x axis
				ax.annotate(jdelaylabel, xy=(0.06, 0.85),  xycoords='axes fraction', ha="left", rotation=90.0) # y axis
			else:
				ax.annotate(idelaylabel, xy=(0.78, 0.08),  xycoords='axes fraction', ha="center") # x axis
				ax.annotate(jdelaylabel, xy=(0.08, 0.76),  xycoords='axes fraction', ha="left", rotation=90.0) # y axis

			ax.set_xlim(-r, r)
			ax.set_ylim(-r, r)
			majorLocator = MultipleLocator(1.0)
			ax.xaxis.set_major_locator(majorLocator)
			majorLocator = MultipleLocator(1.0)
			ax.yaxis.set_major_locator(majorLocator)
			ax.set_xticklabels([])
			ax.set_yticklabels([])
			
			#ax.annotate(delaytext, xy=(0.03, 0.78 - 3*txtstep*(irr+0.5)),  xycoords='axes fraction', color = datarr.plotcolour)
	if title != None:
		plt.figtext(x = 0.5, y = 0.97, s = title, horizontalalignment="center", color="black", fontsize=18)	
	
	#for (irr, rr) in enumerate(rrlist):
	#	plt.figtext(x = left + 0.25*irr, y = 0.96, s = getattr(rr, 'name', 'NoName'), horizontalalignment="left", color=rr.plotcolour)	
	
	for irr, rr in enumerate(rrlist):
			labeltxt = "%s" % (getattr(rr, 'name', 'NoName'))
			plt.figtext(x = right, y = top - txtstep*irr, s = labeltxt, verticalalignment="top", horizontalalignment="right", color=rr.plotcolour)	

		
		
	if filename==None:
		plt.show()
	else:
		plt.savefig(filename)


def delayplot(datarrlist, simrrlist, type="gaussians", rplot=7.0, rbins=7, nbins = 20, binclipr=10.0, displaytext=True, total=True, title=None, figsize=(10, 6), left = 0.06, right=0.97, top=0.99, bottom=0.08, wspace=0.15, hspace=0.3, txtstep=0.04, majorticksstep=2, filename=None, trueshifts=None):
	"""
	Plots final delay measurements for several techniques
	delays are the mean delays as measured on the datarrlist
	errors are statistical + systematic ones, as measured from the simrrlist
	
	
	:param type: either "errorbars" or "gaussians"
	
	:param trueshifts: I will draw vertical lines for the delays corresponding to these shifts.
		Give me a tuple, like (0.0, -5.0, -20.0, -70.0) (for a quad in this case)
	
	"""
	
	n = datarrlist[0].nimages()
	
	labels = datarrlist[0].labels
	
	# To get some fixed centers for the histograms AND the bins, we will use the first element of simrrlist.
	# If you change this, keep using the sims for the bins ...
	reftrueshifts = np.round(simrrlist[0].gettruets()["center"])
	
	# Checks that we compare the same lenses and curves ...
	assert len(datarrlist) == len(simrrlist)
	
	for rr in datarrlist:
		if rr.labels != labels:
			raise RuntimeError("Don't ask me to overplot runresults of different curves !")
	for rr in simrrlist:
		if rr.labels != labels:
			raise RuntimeError("Don't ask me to overplot runresults of different curves !")
	for (datarr, simrr) in zip(datarrlist, simrrlist):
		if datarr.plotcolour != simrr.plotcolour:
			raise RuntimeError("Hmm, plotcolours of data and sim don't correspond !")
		print "Data : %s <-> Simulations : %s" % (datarr.name, simrr.name)
		
	
	fig = plt.figure(figsize=figsize)
	fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)

	axisNum = 0
	for i in range(n): # A, B, C, D and so on
		for j in range(n):
			
			#print i, j
			if (i == 0) or (j == n-1) :
				continue # No plot
				
			axisNum += 1
			if j >= i:
				continue
		
			ax = plt.subplot(n-1, n-1, axisNum)
		
			# Hide the y ticks :
			ax.get_yaxis().set_ticks([])
			minorLocator = MultipleLocator(1.0)
			majorLocator = MultipleLocator(majorticksstep)
			ax.xaxis.set_minor_locator(minorLocator)
			ax.xaxis.set_major_locator(majorLocator)
		
		
			# We determine the plot range :
			reftruedelay = reftrueshifts[i] - reftrueshifts[j]
			plotrange = (reftruedelay - rplot, reftruedelay + rplot)
			
			# Delay label
			# We will express the delays "i - j"
			delaylabel="%s%s" % (labels[j], labels[i])
			
			# We go through all the methods, simultaneously for the data and the corresponding sims :
				
			pdfs = []
			plotx = np.linspace(plotrange[0], plotrange[1], 200)
				
			for (irr, (datarr, simrr)) in enumerate(zip(datarrlist, simrrlist)):
				
				#print getattr(datarr, 'name', 'NoName'), getattr(simrr, 'name', 'NoName')
				# Getting the actual mean delays as measured on the data :
				
				datameasdelays = datarr.tsarray[:,i] - datarr.tsarray[:,j]
				meandatameasdelay = np.mean(datameasdelays)
				
				# Check plot :
				#(counts, bins, patches) = ax.hist(datameasdelays, bins=200, range=plotrange, histtype="bar", color=datarr.plotcolour, alpha = 0.4, lw=0, normed=True)
				#plt.axvline(x=meandatameasdelay, linewidth=1, linestyle="-", color=datarr.plotcolour)
			
				# Getting the corresponding errors from the sims
				
				
				simtruedelays = simrr.truetsarray[:,i] - simrr.truetsarray[:,j]
				simmeasdelays = simrr.tsarray[:,i] - simrr.tsarray[:,j]
					
				simresis = simmeasdelays-simtruedelays
				
				# Check plot :
				#ax.scatter(simtruedelays, simresis, s=2, facecolor=simrr.plotcolour, lw = 0)
				
				# We bin those :		
				binlims = np.linspace(reftruedelay-rbins, reftruedelay+rbins, nbins + 1)
				digitized = np.digitize(simtruedelays, binlims)
				binvals = [simresis[digitized == bini] for bini in range(1, len(binlims))]
				for (bini, binvalarray) in enumerate(binvals):	
					keep = np.logical_and(binvalarray < binclipr, binvalarray > -binclipr)
					if np.sum(keep == False) != 0:
						print "Kicking %i points (%s, %s)" % (np.sum(keep == False), delaylabel, simrr.name)
					binvals[bini] = binvalarray[keep]
				
				binstds = map(np.std, binvals)
				binmedians = map(np.median, binvals)
				binmeans = map(np.mean, binvals)
				
				syserror = np.max(np.fabs(binmeans))
				randerror = np.max(binstds)
				toterror = np.sqrt(syserror*syserror + randerror*randerror)
				
				#randpdf = normal(plotx, mu = meandatameasdelay, sigma = syserror)
				#randpdf /= np.sum(randpdf)
				totpdf = normal(plotx, mu = meandatameasdelay, sigma = toterror)
				totpdf /= np.sum(totpdf)
				
				pdfs.append(totpdf)
				
				#delaytext = r"$%+.1f \pm %.1f (\mathrm{ran}) \pm %.1f (\mathrm{sys})$" % (meandatameasdelay, randerror, syserror)
				delaytext = r"$%+.1f \pm (%.1f,\, %.1f)$" % (meandatameasdelay, randerror, syserror)
				
				
				print "%s : %.2f +/- %.2f (ran) +/- %.2f (sys) using %s" % (delaylabel, meandatameasdelay, randerror, syserror, getattr(datarr, 'name', 'NoName'))
			
				if type == "gaussians":			
					plt.plot(plotx, totpdf, linewidth=1, linestyle="-", color=datarr.plotcolour)
					if displaytext:
						ax.annotate(delaytext, xy=(0.03, 0.78 - 3*txtstep*(irr+0.5)),  xycoords='axes fraction', color = datarr.plotcolour)
			
					
				if type == "errorbars":	
					plt.errorbar([meandatameasdelay], [len(datarrlist) - irr], yerr=None, xerr=toterror, fmt='-', ecolor=datarr.plotcolour, elinewidth=1.5, capsize=3, barsabove=False)
					plt.plot([meandatameasdelay], [len(datarrlist) - irr], marker='o', markeredgecolor=datarr.plotcolour, color=datarr.plotcolour)
					
					if displaytext:
						ax.annotate(delaytext, xy=(meandatameasdelay, len(datarrlist) - irr + 0.3), color = datarr.plotcolour, horizontalalignment="center")
			
					
				# Just to see, we try to get the acutal bias :
				#sysbias = np.mean(binmeans)
				#debiaspdf = normal(plotx, mu = meandatameasdelay-sysbias, sigma = randerror)
				#debiaspdf /= np.sum(debiaspdf)
				#plt.plot(plotx, debiaspdf, linewidth=1, linestyle="--", color=datarr.plotcolour)
			
			
				#if displaytext and type == "gaussians":
					#delaytext = r"$%+.1f \pm %.1f$" % (meandatameasdelay, toterror)
					#print rr.name
					#print delaylabel + "   " + delaytext
					#ax.text(meddelay, np.max(y)/2.0, "%.1f +/- %.1f" % (meddelay, stddelay), horizontalalignment = "center", color = rr.plotcolour)
					
		
			if total and type == "gaussians":
				multpdf = np.ones(len(pdfs[0]))
				for pdf in pdfs:
					multpdf *= pdf
				multpdf /= np.sum(multpdf)
				plt.plot(plotx, multpdf, linewidth=1, linestyle="-", color="black")
				
				if displaytext:
					
					meanmultpdf = np.sum(multpdf * plotx)
					stdmultpdf = np.sqrt(np.sum(multpdf * (plotx - meanmultpdf)**2))
					
					#testpdf = normal(plotx, mu = meanmultpdf, sigma = stdmultpdf)
					#testpdf /= np.sum(testpdf)
					#plt.plot(plotx, testpdf, linewidth=1, linestyle="-", color="yellow")
					#plt.axvline(x=meanmultpdf, linewidth=1, linestyle="-", color="black")
				
					delaytext = r"$%+.1f \pm %.1f$" % (meanmultpdf, stdmultpdf)
					ax.annotate(delaytext, xy=(0.03, 0.78 - 0.1*(irr+1.5)),  xycoords='axes fraction', color = "black")
				
			if i == n-1:
				plt.xlabel("Delay [day]")
			if j == 0 and type == "gaussians":
				plt.ylabel("PDF")
			
			plt.xlim(plotrange)
			if type == "errorbars":	
				plt.ylim((0.5, len(datarrlist)+1.5))
			
			#if errorrange != None:
			#	plt.ylim((-errorrange, errorrange))
			
			plt.annotate(delaylabel, xy=(0.03, 0.88-txtstep),  xycoords='axes fraction', fontsize=12, color="black")
			
			if trueshifts != None:
				truedelay = trueshifts[i] - trueshifts[j]
				plt.axvline(truedelay, color="gray", linestyle="--", dashes=(3, 3), zorder=-20)
			
			
	# The "legend" :
	for irr, datarr in enumerate(datarrlist):
	
		labeltxt = "%s" % (getattr(datarr, 'name', 'NoName'))
		plt.figtext(x = right, y = top - txtstep*irr, s = labeltxt, verticalalignment="top", horizontalalignment="right", color=datarr.plotcolour)	

	if title != None:
		plt.figtext(x = 0.5, y = 0.95, s = title, horizontalalignment="center", color="black", fontsize=18)	

	if filename==None:
		plt.show()
	else:
		plt.savefig(filename)



# 
# 
# def qsscatter(rr, r=10.0, nbins=100):
# 	"""
# 	Delay distribution vs qs values
# 	"""
# 	
# 	n = rr.tsarray.shape[1]
# 	labels = rr.labels
# 	trueshifts = rr.trueshifts
# 	qs = rr.qs
# 	
# 	fig = plt.figure(figsize=(16, 9))
# 	
# 	fig.subplots_adjust(left = 0.03, right=0.95, bottom=0.05, top=0.9, wspace=0.2, hspace=0.2)
# 
# 	axisNum = 0
# 	for i in range(n): # [A, B, C, D]
# 	    for j in range(n):
# 	        if (i == 0) or (j == n-1) :
# 	        	continue # No plot
# 		axisNum += 1
# 		if j >= i:
# 			continue
# 		
# 		ax = plt.subplot(n-1, n-1, axisNum)
# 	
# 		truedelay = trueshifts[i] - trueshifts[j]
# 		plotrange = (truedelay - r, truedelay + r)
# 		delaylabel="%s%s" % (labels[j], labels[i])
# 		print delaylabel
# 		
# 	    	delays = rr.tsarray[:,i] - rr.tsarray[:,j]
# 	    	assert len(delays) == len(qs)
# 	    	
# 		(counts, bins, patches) = ax.hist(delays, bins=nbins, histtype="bar", range=plotrange, normed=True, facecolor='0.5', lw = 0)
# 		# Hide the y ticks :
# 		ax.get_yaxis().set_ticks([])
# 		
# 		axscatter = ax.twinx()
# 	
# 		axscatter.scatter(delays, qs, s=1, facecolor='black', lw = 0)
# 		#ax2.set_ylabel('sin', color='r')
# 		#for tl in ax2.get_yticklabels():
#    		#tl.set_color('r')
# 		#plt.scatter(delays, qs, s=1, facecolor='black', lw = 0)
# 		
# 		# We plot the true shift :
# 		plt.axvline(x=truedelay, linewidth=1, color="black")
# 	
# 		plt.xlim(plotrange)
# 		#plt.ylim(0, np.max(counts))
# 		plt.annotate(delaylabel, xy=(0.03, 0.88),  xycoords='axes fraction', fontsize=12, color="black")
# 	
# 	# A title for the figure
# 	labeltxt = "%s (%i)" % (getattr(rr, 'simset', 'NoName'), rr.tsarray.shape[0])
# 	plt.figtext(x = 0.4, y = 0.93, s = labeltxt, color="black")	
# 
# 	
# 	print 'Plotting simset "%s"' % labeltxt
# 	print "     Labels : %s" % (", ".join(rr.labels))
# 	print "     Median shifts : %s" % (", ".join(["%.2f" % (np.median(rr.tsarray[:,i])) for i in range(len(rr.labels))]))
# 	print "     Std shifts : %s" % (", ".join(["%.2f" % (np.std(rr.tsarray[:,i])) for i in range(len(rr.labels))]))
# 	
# 	plt.show()
# 		
# 
# 
# 
# 
# def qshist(rr, r=10.0, nbins=50, stat="mean"):
# 	"""
# 	Delay distribution histograms, coloured according to median qs values.
# 	
# 	stat = "mean", "med", "min"
# 	
# 	"""
# 	
# 	if rr.qs == None:
# 		#rr.qs = rr.tsarray[:, 2]
# 		raise RuntimeError("Your rr has no qs !")
# 	
# 	n = rr.tsarray.shape[1]
# 	labels = rr.labels
# 	trueshifts = rr.trueshifts
# 	qs = rr.qs
# 	statnorm = colors.normalize(np.min(qs), np.max(qs)) # Important, this way we get the same colours for all histograms
# 	
# 	print "Colorbar min = %.2f, max = %.2f" % (np.min(qs), np.max(qs))
# 	
# 	plt.figure(figsize=(16, 9))
# 
# 	axisNum = 0
# 	for i in range(n): # [A, B, C, D]
# 	    for j in range(n):
# 	        if (i == 0) or (j == n-1) :
# 	        	continue # No plot
# 		axisNum += 1
# 		if j >= i:
# 			continue
# 		
# 		ax = plt.subplot(n-1, n-1, axisNum)
# 	
# 		truedelay = trueshifts[i] - trueshifts[j]
# 		histrange = (truedelay - r, truedelay + r)
# 		delaylabel="%s%s" % (labels[j], labels[i])
# 		print delaylabel
# 		
# 	    	delays = rr.tsarray[:,i] - rr.tsarray[:,j]
# 	    	assert len(delays) == len(qs)
# 	    	
# 		(counts, bins, patches) = plt.hist(delays, bins=nbins, histtype="bar", normed=True, range=histrange)
# 		
# 		assert nbins+1 == len(bins)
# 		# bins are nbins+1 floats
# 		
# 		# We construct the mean q for each bin
# 		statbinqs = []
# 		spanbinqs = []
# 		
# 		for bini in range(nbins):
# 			a = bins[bini]
# 			b = bins[bini+1]
# 			qsofbin = []
# 			for qsi in range(len(qs)):
# 				if delays[qsi] < b and delays[qsi] > a:
# 					qsofbin.append(qs[qsi])
# 			if len(qsofbin) == 0:
# 				qsofbin.append(0)
# 			
# 			qsofbin = np.array(qsofbin)
# 			#print "Bin %i, %.1f -> %.1f, relat = %f, min(qs) = %.1f, max(qs) = %.1f" % (bini, a, b, len(qsofbin)/counts[bini], np.min(qsofbin), np.max(qsofbin))
# 			if stat == "mean":
# 				statbinqs.append(np.mean(qsofbin))
# 			elif stat == "min":
# 				statbinqs.append(np.min(qsofbin))
# 			else:
# 				raise RuntimeError("Unknown stat")
# 			
# 			spanbinqs.append(np.max(qsofbin) - np.min(qsofbin))
# 			
# 		assert len(statbinqs) == nbins
# 			
# 		
# 		statbinqs = np.array(statbinqs)	
# 		spanbinqs = np.array(spanbinqs)
# 		
# 		spannorm = colors.normalize(np.min(spanbinqs[spanbinqs>0.0]), np.max(spanbinqs))
# 		#statnorm = colors.normalize(np.min(statbinqs[spanbinqs>0.0]), np.max(statbinqs))
# 		#print "Warning I redefine statnorm for every plot"
# 		
# 		for statbinq, spanbinq, patch in zip(statbinqs, spanbinqs, patches):
# 			patch.set_facecolor(cm.jet(statnorm(statbinq)))
# 			patch.set_edgecolor("None")
# 			#patch.set_edgecolor(cm.jet(spannorm(spanbinq)))
# 			
# 		# We plot the true shift :
# 		plt.axvline(x=truedelay, linewidth=1, color="black")
# 	
# 		plt.xlim(histrange)
# 		#plt.ylim(0, np.max(counts))
# 		plt.annotate(delaylabel, xy=(0.03, 0.88),  xycoords='axes fraction', fontsize=12, color="black")
# 	
# 	# A title for the figure
# 	labeltxt = "%s (%i)" % (getattr(rr, 'simset', 'NoName'), rr.tsarray.shape[0])
# 	plt.figtext(x = 0.4, y = 0.93, s = labeltxt, color="black")	
# 
# 	
# 	print 'Plotting simset "%s"' % labeltxt
# 	print "     Labels : %s" % (", ".join(rr.labels))
# 	print "     Median shifts : %s" % (", ".join(["%.2f" % (np.median(rr.tsarray[:,i])) for i in range(len(rr.labels))]))
# 	print "     Std shifts : %s" % (", ".join(["%.2f" % (np.std(rr.tsarray[:,i])) for i in range(len(rr.labels))]))
# 	
# 	plt.show()
# 		
# 
# 
# 
# 
# 
# def couplehists(rrlist, r=10.0, nbins=50):
# 	"""
# 	Comparing the delay distributions from different run result objects.
# 	
# 	rrlist : a list of runresults object, where each object is usually provided by collect().
# 	I will plot each rr of this list in a different colour.
# 		
# 	r = a radius (we make plots around the delays from the truth ?)
# 	"""
# 		
# 	n = rrlist[0].tsarray.shape[1]
# 	
# 	labels = rrlist[0].labels
# 	trueshifts = rrlist[0].trueshifts
# 	
# 	for rr in rrlist:
# 		if rr.labels != labels:
# 			raise RuntimeError("Don't ask me to overplot runresults of different curves !")
# 		if rr.trueshifts != trueshifts:
# 			print "Warning : I use the trueshift of the first rr to set the ranges."
# 			#raise RuntimeError("Don't ask me to overplot runresults with different true shifts !")
# 	
# 	colours = ["blue", "red", "#008800", "cyan", "gray"]
# 	
# 	fig = plt.figure(figsize=(16, 9))
# 	fig.subplots_adjust(left = 0.03, right=0.95, bottom=0.05, top=0.9, wspace=0.2, hspace=0.2)
# 
# 	axisNum = 0
# 	for i in range(n): # [A, B, C, D]
# 	    for j in range(n):
# 	        #print i, j
# 	        if (i == 0) or (j == n-1) :
# 	        	continue # No plot
# 		axisNum += 1
# 		if j >= i:
# 			continue
# 		
# 		ax = plt.subplot(n-1, n-1, axisNum)
# 		axscatter = ax.twinx()
# 		
# 		# Labels only for some plots :
# 		#if j == 0:
# 		#	plt.ylabel(labels[i])
# 		#if i == n-1:
# 		#	plt.xlabel(labels[j])
# 			
# 		# Hide the y ticks :
# 		ax.get_yaxis().set_ticks([])
# 		
# 		# Ranges to plot
# 		truedelay = trueshifts[i] - trueshifts[j]
# 		histrange = (truedelay - r, truedelay + r)
# 		
# 		# Delay label
# 		delaylabel="%s%s" % (labels[j], labels[i])
# 		
# 	    	for irr, rr in enumerate(rrlist):
# 	    		# We will express the delays "i - j"
# 	    		delays = rr.tsarray[:,i] - rr.tsarray[:,j]
# 	    		
# 			#(counts, bins, patches) = ax.hist(delays, bins=nbins, range=histrange, histtype="step", color=colours[irr % len(colours)], normed=True)
# 			(counts, bins, patches) = ax.hist(delays, bins=nbins, range=histrange, histtype="bar", color=colours[irr % len(colours)], alpha = 0.4, lw=0, normed=True)
# 			
# 			axscatter.scatter(delays, rr.qs, s=1, facecolor=colours[irr % len(colours)], lw = 0)
# 
# 			# We plot the true shifts individually for each rr :
# 			#truedelay = rr.trueshifts[i] - rr.trueshifts[j]
# 			#plt.axvline(x=truedelay, linewidth=1, color=colours[irr % len(colours)])
# 			
# 		plt.axvline(x=truedelay, linewidth=1, color='black')
# 		
# 		plt.xlim(histrange)
# 		plt.annotate(delaylabel, xy=(0.03, 0.88),  xycoords='axes fraction', fontsize=12, color="black")
# 		
# 	for irr, rr in enumerate(rrlist):
# 	
# 		labeltxt = "%s (%i)" % (getattr(rr, 'simset', 'NoName'), rr.tsarray.shape[0])
# 		plt.figtext(x = 0.1 + 0.2*irr, y = 0.93, s = labeltxt, color=colours[irr % len(colours)])	
# 
# 		
# 		print 'Plotting simset "%s"' % labeltxt
# 		print "     Labels : %s" % (", ".join(rr.labels))
# 		print "     Median shifts : %s" % (", ".join(["%.2f" % (np.median(rr.tsarray[:,i])) for i in range(len(rr.labels))]))
# 		print "     Std shifts : %s" % (", ".join(["%.2f" % (np.std(rr.tsarray[:,i])) for i in range(len(rr.labels))]))
# 		
# 	plt.show()
# 		