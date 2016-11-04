"""
Subpackage with functions to plot all kind of results from runs.
"""


import numpy as np
import math, sys

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator

import scipy.ndimage

import pycs.gen.util

class delaycontainer:
	"""
	Stores the delay or error bar measurement(s) (one for each curve pair).
	This object is usually produced by the plot-functions ``hists`` or ``meanvstrue`` below.
	
	markers : [ 7 | 4 | 5 | 6 | 'o' | 'D' | 'h' | 'H' | '_' | '' | 'None' | ' ' | None | '8' | 'p' | ',' | '+' | '.' | 's' | '*' | 'd' | 3 | 0 | 1 | 2 | '1' | '3' | '4' | '2' | 'v' | '<' | '>' | '^' | '|' | 'x' | '$...$' | tuple | Nx2 array ]
	
	"""
	
	def __init__(self, data, name="No name", objects=None, plotcolour = "black", marker=None):
		"""
		self.data is a list of dicts, the fields depend on if it's delays or errorbars
		
			* delays : label, mean, med, std
			* errorbars : label, tot, sys, ran, bias
		
		
		"""	
		self.data = data
		self.name = name
		self.objects = objects
	
	
		self.plotcolour = plotcolour
		self.marker = marker
		self.markersize = 5.0
		self.yshift = 0.0 # allows to "group" measurements


def newdelayplot(plotlist, rplot=7.0, displaytext=True, hidedetails=False, showbias=True, showran=True, showlegend=True, text=None, figsize=(10, 6), left = 0.06, right=0.97, top=0.99, bottom=0.08, wspace=0.15, hspace=0.3, txtstep=0.04, majorticksstep=2, filename=None, refshifts=None, refdelays=None, legendfromrefdelays=False, hatches=None, centershifts=None, ymin=0.2, hlines=None):
	"""
	Plots delay measurements from different methods, telescopes, sub-curves, etc in one single plot.
	For this I use only ``delaycontainer`` objects, i.e. I don't do any "computation" myself.
	
	:param plotlist: Give me a list of tuples (delays, errorbars), where delays and errorbars are delaycontainer objects as written into pkl files by ``hists`` and ``meanvstrue``.
	:type plotlist: list
	
	:param rplot: radius of delay axis, in days.
	
	:param displaytext: Show labels with technique names and values of delays
	:type displaytext: boolean
	
	:param hidedetails: Do not show (ran, sys) in labels
	:type hidedetails: boolean
	
	
	:param refshifts: This is a list of dicts like {"colour":"gray", "shifts":(0, 0, 0, 90)}. Will be plotted as dashed vertical lines.
	:type refshifts: list
	
	:param refdelays: a list of tuples (delays, errorbars) to be plotted as shaded vertical zones.
	:type refdelays: list

	:param legendfromrefdelays: if you want to display the refdelays name in the legend panel
	:type refdelays: boolean

	:param hatches: list of hatch keyword for the refdelays plotting
	:type refdelays: list
	
	:param showbias: draws a little cross at the position of the delay "corrected" for the bias.
	:type showbias: boolean
	
	:param showran: draws "minor" error bar ticks using the random error only.
	:type showran: boolean
	
	:param text:
		Text that you want to display, in the form : [line1, line2, line3 ...]
		where line_i is (x, y, text, kwargs) where kwargs is e.g. {"fontsize":18} and x and y are relative positions (from 0 to 1). 
	:type text: list

	
	"""
	
	# Some checks :
	objects = plotlist[0][0].objects
	for (delays, errors) in plotlist:
		if delays.objects != objects or errors.objects != objects:
			raise RuntimeError("Don't ask me to overplot stuff from different objects !")
	
	n = len(objects)
	nmeas = len(plotlist)
	print "Objects : %s" % (", ".join(objects))
	
	for (delays, errors) in plotlist:
		if delays.plotcolour != errors.plotcolour:
			raise RuntimeError("Hmm, plotcolours of delays and errors don't correspond !")
		print "Delays : %s <-> Errors : %s" % (delays.name, errors.name)
		
	
	fig = plt.figure(figsize=figsize)
	fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)

	axisNum = 0
	print "#"*80
	for i in range(n): # A, B, C, D and so on
		for j in range(n):
			
			#print i, j
			if (i == 0) or (j == n-1) :
				continue # No plot
				
			axisNum += 1
			if j >= i:
				continue
		
			ax = plt.subplot(n-1, n-1, axisNum)
		
			# We will express the delays "i - j"
			delaylabel="%s%s" % (objects[j], objects[i])
			print "           Delay %s" % (delaylabel)
			
			# General esthetics :
			ax.get_yaxis().set_ticks([])
			minorLocator = MultipleLocator(1.0)
			majorLocator = MultipleLocator(majorticksstep)
			ax.xaxis.set_minor_locator(minorLocator)
			ax.xaxis.set_major_locator(majorLocator)
		
			# To determine the plot range :
			paneldelays = []
			
			# Going throuh plotlist :
			for (ipl,(delays, errors)) in enumerate(plotlist):
				
				# Getting the delay for this particular panel
				delay = [meas for meas in delays.data if meas["label"] == delaylabel][0]
				error = [meas for meas in errors.data if meas["label"] == delaylabel][0]
				
				paneldelays.append(delay["mean"])
				
				ypos = nmeas - ipl + delays.yshift
				
				plt.errorbar([delay["mean"]], [ypos], yerr=None, xerr=error["tot"], fmt='-', ecolor=delays.plotcolour, elinewidth=1.5, capsize=3, barsabove=False)
				if showran:
					plt.errorbar([delay["mean"]], [ypos], yerr=None, xerr=error["ran"], fmt='-', ecolor=delays.plotcolour, elinewidth=0.5, capsize=2, barsabove=False)

				
				if delays.marker == None or delays.marker == ".":
					plt.plot([delay["mean"]], [ypos], marker='o', markersize=delays.markersize, markeredgecolor=delays.plotcolour, color=delays.plotcolour)
				else:
					plt.plot([delay["mean"]], [ypos], marker=delays.marker, markersize=delays.markersize, markeredgecolor=delays.plotcolour, color="white")
				
				if showbias:
					plt.plot([delay["mean"] - error["bias"]], [ypos], marker="x", markersize=delays.markersize, markeredgecolor=delays.plotcolour, color=delays.plotcolour)

				if hidedetails or (error["ran"] < 0.001 and error["sys"] < 0.001): # Then we ommit to write them.
					delaytext = r"$%+.1f \pm %.1f$" % (delay["mean"], error["tot"])
				else:
					delaytext = r"$%+.1f \pm %.1f\,(%.1f, %.1f)$" % (delay["mean"], error["tot"], error["ran"], error["sys"])
				
				if n==2: # For doubles, we include the technique name into the txt :
					delaytext = r"%s : " % (delays.name) + delaytext	
			
				if displaytext:
					ax.annotate(delaytext, xy=(delay["mean"], ypos + 0.3), color = delays.plotcolour, horizontalalignment="center", fontsize=14)
				
				print "%45s : %+6.2f +/- %.2f (%.2f, %.2f)" % (delays.name, delay["mean"], error["tot"], error["ran"], error["sys"])
		
			print "#"*80
			# Now this panel is done. Some general settings :
			
			if centershifts != None:
				centerdelay = centershifts[i] - centershifts[j]
			else:
				centerdelay = np.median(paneldelays)
		
			plt.xlim((centerdelay - rplot, centerdelay + rplot))
			plt.ylim((ymin, nmeas+1.5))
			
			if i == n-1:
				plt.xlabel("Delay [day]", fontsize=14)
			if n != 2: # otherwise only one panel, no need
				plt.annotate(delaylabel, xy=(0.03, 0.88-txtstep),  xycoords='axes fraction', fontsize=14, color="black")
	
			if refshifts != None:
				for item in refshifts:
					refdelay = item["shifts"][i] - item["shifts"][j]
					plt.axvline(refdelay, color=item["colour"], linestyle="--", dashes=(3, 3), zorder=-20)
			
			if refdelays != None:
				for (ipl,(delays, errors)) in enumerate(refdelays):
				
					# Getting the delay for this particular panel
					delay = [meas for meas in delays.data if meas["label"] == delaylabel][0]
					error = [meas for meas in errors.data if meas["label"] == delaylabel][0]

					if hatches!=None:
						plt.axvspan(delay["mean"]-error["tot"], delay["mean"]+error["tot"], facecolor=delays.plotcolour, alpha=0.25, zorder=-20, edgecolor="none", linewidth=0, hatch=hatches[ipl])
					else:
						plt.axvspan(delay["mean"]-error["tot"], delay["mean"]+error["tot"], facecolor=delays.plotcolour, alpha=0.25, zorder=-20, edgecolor="none", linewidth=0)
					plt.axvline(delay["mean"], color=delays.plotcolour, linestyle="--", dashes=(5, 5), lw=1.0, zorder=-20)
					#plt.axvline(delay["mean"], color=item.plotcolour, linestyle="-", lw=2, alpha=0.5, zorder=-20)
			
			if hlines != None:
				for hline in hlines:
					plt.axhline(hline, lw=0.5, color="gray", zorder=-30)
					
	
	
	# The "legend" :
	if showlegend:
		for (ipl,(delays, errors)) in enumerate(plotlist):
			line = "%s" % (delays.name)
			plt.figtext(x = right, y = top - txtstep*ipl, s = line, verticalalignment="top", horizontalalignment="right", color=delays.plotcolour, fontsize=14)
		if legendfromrefdelays:
			for (ipl,(delays, errors)) in enumerate(refdelays):
				line = "%s" % (delays.name)
				plt.figtext(x = right, y = top - txtstep*(ipl+len(plotlist)), s = line, verticalalignment="top", horizontalalignment="right", color=delays.plotcolour, fontsize=14)

	
	# Generic text :
	if text != None:
		for line in text:
			plt.figtext(x=line[0], y=line[1], s=line[2], **line[3])


	if filename==None:
		plt.show()
	else:
		plt.savefig(filename)

	
	


def newdelayplot2(plotlist, rplot=7.0, displaytext=True, hidedetails=False, showbias=True, showran=True, showlegend=True, text=None, figsize=(10, 6), left = 0.06, right=0.97, top=0.99, bottom=0.08, wspace=0.15, hspace=0.3, txtstep=0.04, majorticksstep=2, filename=None, refshifts=None, refdelays=None, legendfromrefdelays=False, hatches=None, centershifts=None, ymin=0.2, hlines=None):
	"""
	Plots delay measurements from different methods, telescopes, sub-curves, etc in one single plot.
	For this I use only ``delaycontainer`` objects, i.e. I don't do any "computation" myself.

	Difference from newdelayplot is that the previously hatched/shaded regions are plotted as smaller points, without infos on the time-delay

	:param plotlist: Give me a list of tuples (delays, errorbars), where delays and errorbars are delaycontainer objects as written into pkl files by ``hists`` and ``meanvstrue``.
	:type plotlist: list

	:param rplot: radius of delay axis, in days.

	:param displaytext: Show labels with technique names and values of delays
	:type displaytext: boolean

	:param hidedetails: Do not show (ran, sys) in labels
	:type hidedetails: boolean


	:param refshifts: This is a list of dicts like {"colour":"gray", "shifts":(0, 0, 0, 90)}. Will be plotted as dashed vertical lines.
	:type refshifts: list

	:param refdelays: a list of tuples (delays, errorbars) to be plotted as shaded vertical zones.
	:type refdelays: list

	:param legendfromrefdelays: if you want to display the refdelays name in the legend panel
	:type refdelays: boolean

	:param hatches: list of hatch keyword for the refdelays plotting
	:type refdelays: list

	:param showbias: draws a little cross at the position of the delay "corrected" for the bias.
	:type showbias: boolean

	:param showran: draws "minor" error bar ticks using the random error only.
	:type showran: boolean

	:param text:
		Text that you want to display, in the form : [line1, line2, line3 ...]
		where line_i is (x, y, text, kwargs) where kwargs is e.g. {"fontsize":18} and x and y are relative positions (from 0 to 1).
	:type text: list


	"""

	# Some checks :
	objects = plotlist[0][0].objects
	for (delays, errors) in plotlist:
		if delays.objects != objects or errors.objects != objects:
			raise RuntimeError("Don't ask me to overplot stuff from different objects !")

	n = len(objects)
	nmeas = len(plotlist)+len(refdelays)/2 +1
	print "Objects : %s" % (", ".join(objects))

	for (delays, errors) in plotlist:
		if delays.plotcolour != errors.plotcolour:
			raise RuntimeError("Hmm, plotcolours of delays and errors don't correspond !")
		print "Delays : %s <-> Errors : %s" % (delays.name, errors.name)


	fig = plt.figure(figsize=figsize)
	fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)

	axisNum = 0
	print "#"*80
	for i in range(n): # A, B, C, D and so on
		for j in range(n):

			#print i, j
			if (i == 0) or (j == n-1) :
				continue # No plot

			axisNum += 1
			if j >= i:
				continue

			ax = plt.subplot(n-1, n-1, axisNum)

			# We will express the delays "i - j"
			delaylabel="%s%s" % (objects[j], objects[i])
			print "           Delay %s" % (delaylabel)

			# General esthetics :
			ax.get_yaxis().set_ticks([])
			minorLocator = MultipleLocator(1.0)
			majorLocator = MultipleLocator(majorticksstep)
			ax.xaxis.set_minor_locator(minorLocator)
			ax.xaxis.set_major_locator(majorLocator)

			# To determine the plot range :
			paneldelays = []

			# Going throuh plotlist :
			for (ipl,(delays, errors)) in enumerate(plotlist):

				# Getting the delay for this particular panel
				delay = [meas for meas in delays.data if meas["label"] == delaylabel][0]
				error = [meas for meas in errors.data if meas["label"] == delaylabel][0]

				paneldelays.append(delay["mean"])

				ypos = nmeas - ipl*1.3 + delays.yshift


				plt.errorbar([delay["mean"]], [ypos], yerr=None, xerr=error["tot"], fmt='-', ecolor=delays.plotcolour, elinewidth=delays.markersize/5.0*1.5, capsize=3, barsabove=False)
				if showran:
					plt.errorbar([delay["mean"]], [ypos], yerr=None, xerr=error["ran"], fmt='-', ecolor=delays.plotcolour, elinewidth=0.5, capsize=2, barsabove=False)


				if delays.marker == None or delays.marker == ".":
					plt.plot([delay["mean"]], [ypos], marker='o', markersize=delays.markersize, markeredgecolor=delays.plotcolour, color=delays.plotcolour)
				else:
					plt.plot([delay["mean"]], [ypos], marker=delays.marker, markersize=delays.markersize, markeredgecolor=delays.plotcolour, color="white")

				if showbias:
					plt.plot([delay["mean"] - error["bias"]], [ypos], marker="x", markersize=delays.markersize, markeredgecolor=delays.plotcolour, color=delays.plotcolour)

				if hidedetails or (error["ran"] < 0.001 and error["sys"] < 0.001): # Then we ommit to write them.
					delaytext = r"$%+.1f \pm %.1f$" % (delay["mean"], error["tot"])
				else:
					delaytext = r"$%+.1f \pm %.1f\,(%.1f, %.1f)$" % (delay["mean"], error["tot"], error["ran"], error["sys"])

				if n==2: # For doubles, we include the technique name into the txt :
					delaytext = r"%s : " % (delays.name) + delaytext

				if displaytext:
					if delays.markersize>5:
						ax.annotate(delaytext, xy=(delay["mean"], ypos + 0.3), color = delays.plotcolour, horizontalalignment="center", fontsize=16)
					else:
						ax.annotate(delaytext, xy=(delay["mean"], ypos + 0.3), color = delays.plotcolour, horizontalalignment="center", fontsize=14)

				print "%45s : %+6.2f +/- %.2f (%.2f, %.2f)" % (delays.name, delay["mean"], error["tot"], error["ran"], error["sys"])


			# Going throuh plotlist :
			for (ipl,(delays, errors)) in enumerate(refdelays):

				# Getting the delay for this particular panel
				delay = [meas for meas in delays.data if meas["label"] == delaylabel][0]
				error = [meas for meas in errors.data if meas["label"] == delaylabel][0]

				paneldelays.append(delay["mean"])

				if ipl in [0, 1]:
					ypos = nmeas - (ipl/2.5+4.2) + 0.6
				elif ipl in [2, 3]:
					ypos = nmeas - (ipl/2.5+4.2) + 0.6 -0.4
				elif ipl in [4, 5]:
					ypos = nmeas - (ipl/2.5+4.2) + 0.6 -0.8


				plt.errorbar([delay["mean"]], [ypos], yerr=None, xerr=error["tot"], fmt='-', ecolor=delays.plotcolour, elinewidth=1.0, capsize=3, barsabove=False)
				if showran:
					plt.errorbar([delay["mean"]], [ypos], yerr=None, xerr=error["ran"], fmt='-', ecolor=delays.plotcolour, elinewidth=0.33, capsize=2, barsabove=False)


				if delays.marker == None or delays.marker == ".":
					plt.plot([delay["mean"]], [ypos], marker='o', markersize=delays.markersize/1.5, markeredgecolor=delays.plotcolour, color=delays.plotcolour)
				else:
					plt.plot([delay["mean"]], [ypos], marker=delays.marker, markersize=delays.markersize/1.5, markeredgecolor=delays.plotcolour, color=delays.plotcolour)

				if showbias:
					plt.plot([delay["mean"] - error["bias"]], [ypos], marker="x", markersize=delays.markersize, markeredgecolor=delays.plotcolour, color=delays.plotcolour)

				if hidedetails or (error["ran"] < 0.001 and error["sys"] < 0.001): # Then we ommit to write them.
					delaytext = r"$%+.1f \pm %.1f$" % (delay["mean"], error["tot"])
				else:
					delaytext = r"$%+.1f \pm %.1f\,(%.1f, %.1f)$" % (delay["mean"], error["tot"], error["ran"], error["sys"])

				if n==2: # For doubles, we include the technique name into the txt :
					delaytext = r"%s : " % (delays.name) + delaytext

				if displaytext:
					pass
					#ax.annotate(delaytext, xy=(delay["mean"], ypos + 0.3), color = delays.plotcolour, horizontalalignment="center", fontsize=14)

				print "%45s : %+6.2f +/- %.2f (%.2f, %.2f)" % (delays.name, delay["mean"], error["tot"], error["ran"], error["sys"])

			if axisNum ==1:
				"""
				ax.annotate(r"0", xy=(-13.3, 4.1), color = "crimson", horizontalalignment="center", fontsize=16)
				ax.annotate(r"1", xy=(-13.3, 2.9), color = "crimson", horizontalalignment="center", fontsize=16)
				ax.annotate(r"2", xy=(-13.33, 1.7), color = "crimson", horizontalalignment="center", fontsize=16)
				ax.annotate(r"3", xy=(-13.37, 0.5), color = "crimson", horizontalalignment="center", fontsize=16)
				"""

				"""
				ax.annotate(r"$\diamond$", xy=(-13.3, 3.0), color = "crimson", horizontalalignment="center", fontsize=18)
				ax.annotate(r"$\dag$", xy=(-13.33, 1.8), color = "crimson", horizontalalignment="center", fontsize=18)
				ax.annotate(r"$\bowtie$", xy=(-13.37, 0.6), color = "crimson", horizontalalignment="center", fontsize=18)
				"""


			print "#"*80
			# Now this panel is done. Some general settings :

			if centershifts != None:
				centerdelay = centershifts[i] - centershifts[j]
			else:
				centerdelay = np.median(paneldelays)

			plt.xlim((centerdelay - rplot, centerdelay + rplot))
			plt.ylim((ymin, nmeas+1.5))

			if i == n-1:
				plt.xlabel("Delay [day]", fontsize=14)
			if n != 2: # otherwise only one panel, no need
				plt.annotate(delaylabel, xy=(0.03, 0.88-txtstep),  xycoords='axes fraction', fontsize=14, color="black")

			if refshifts != None:
				for item in refshifts:
					refdelay = item["shifts"][i] - item["shifts"][j]
					plt.axvline(refdelay, color=item["colour"], linestyle="--", dashes=(3, 3), zorder=-20)


			if hlines != None:
				for hline in hlines:
					plt.axhline(hline, lw=0.5, color="gray", zorder=-30)



	# The "legend" :
	if showlegend:
		for (ipl,(delays, errors)) in enumerate(plotlist):
			line = "%s" % (delays.name)
			plt.figtext(x = right, y = top - txtstep*ipl, s = line, verticalalignment="top", horizontalalignment="right", color=delays.plotcolour, fontsize=16)
		"""
		if legendfromrefdelays:
			for (ipl,(delays, errors)) in enumerate(refdelays):
				line = "%s" % (delays.name)
				plt.figtext(x = right, y = top - txtstep*(ipl+len(plotlist)), s = line, verticalalignment="top", horizontalalignment="right", color=delays.plotcolour, fontsize=12)
		"""
		"""
		plt.figtext(x = right-0.123, y = top - txtstep*len(plotlist) - 0.025, s = r"$\diamond$", verticalalignment="top", horizontalalignment="right", color="crimson", fontsize=18)
		plt.figtext(x = right-0.125, y = top - txtstep*(len(plotlist)+1) - 0.023 , s = r"$\dag$", verticalalignment="top", horizontalalignment="right", color="crimson", fontsize=18)
		plt.figtext(x = right-0.12, y = top - txtstep*(len(plotlist)+2) - 0.025, s = r"$\bowtie$", verticalalignment="top", horizontalalignment="right", color="crimson", fontsize=18)

		plt.figtext(x = right, y = top - txtstep*len(plotlist) - 0.03, s = "- 2003-2007", verticalalignment="top", horizontalalignment="right", color="crimson", fontsize=13)
		plt.figtext(x = right, y = top - txtstep*(len(plotlist)+1) - 0.03 , s = "- 2008-2012", verticalalignment="top", horizontalalignment="right", color="crimson", fontsize=13)
		plt.figtext(x = right, y = top - txtstep*(len(plotlist)+2) - 0.03, s = "- 2013-2016", verticalalignment="top", horizontalalignment="right", color="crimson", fontsize=13)
		"""

	# Generic text :
	if text != None:
		for line in text:
			plt.figtext(x=line[0], y=line[1], s=line[2], **line[3])


	if filename==None:
		plt.show()

	else:
		plt.savefig(filename)


def normal(x, mu, sigma):
	"""
	Plain normal distribution.
	You can directly apply me on numpy arrays x, mu, sigma.
	"""
	
	return (1.0/np.sqrt(2.0*np.pi*sigma*sigma)) * np.exp( - (x - mu)**2/(2*sigma*sigma))


def hists(rrlist, r=10.0, nbins=100, showqs=True, showallqs=False, qsrange=None, title=None, niceplot=False, displaytext=True, figsize=(16, 9), left = 0.06, right=0.95, bottom=0.065, top=0.95, wspace=0.2, hspace=0.2, txtstep=0.04, majorticksstep=2, trueshifts=None, filename=None, dataout=False):
	"""
	Comparing the delay distributions from different run result objects.
	
	:param rrlist: a list of runresults object.		
	:param r: a range radius for the hists
	
	:param showqs: If True, I overplot the qs as scatter points.
	
	:param dataout: True means that I'll write the pkl file needed to make the delayplot.
	
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
	
		rr.tmpdata = []
	
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
			
			# Delay label, used not only for display purposes, but also for the output pkl.
			delaylabel="%s%s" % (labels[j], labels[i])
			
			if i == n-1:
				if n == 2: # Only one panel -> we write the object names into the xlabel
					plt.xlabel("Delay %s%s [day]" % (labels[j], labels[i]))			
				else:
					plt.xlabel("Delay [day]")
			if showqs:
				axscatter = ax.twinx()
				
						
			# Hide the y ticks :
			#ax.get_yaxis().set_ticks([])
			
			# Ranges to plot
			reftruedelay = reftrueshifts[i] - reftrueshifts[j]
			histrange = (reftruedelay - r, reftruedelay + r)
			
			
			for irr, rr in enumerate(rrlist):
				# We will express the delays "i - j"
				if rr.plottrue == True:
					delays = rr.truetsarray[:,i] - rr.truetsarray[:,j]
				else:
					delays = rr.tsarray[:,i] - rr.tsarray[:,j]
				
				meddelay = np.median(delays)
				meandelay = np.mean(delays)
				stddelay = np.std(delays)
				
				# We save these :
				rr.tmpdata.append({"label":delaylabel, "mean":meandelay, "med":meddelay, "std":stddelay})
				
				#(counts, bins, patches) = ax.hist(delays, bins=nbins, range=histrange, histtype="step", color=colours[irr % len(colours)], normed=True)
				(counts, bins, patches) = ax.hist(delays, bins=nbins, range=histrange, histtype="bar", color=rr.plotcolour, alpha = 0.4, lw=0, normed=True)
				
				if niceplot:
					majorLocator = MultipleLocator(majorticksstep)
					minorLocator = MultipleLocator(1.0)
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
					
					if getattr(rr, "plotgauss", False) == True:
						x = np.linspace(histrange[0], histrange[1], 100)
						y = normal(x, meandelay, stddelay)
						ax.plot(x, y, linestyle="-", color = rr.plotcolour)
					
					delaytext = r"%+.1f $\pm$ %.1f" % (meandelay, stddelay)
					#print rr.name
					#print delaylabel + "   " + delaytext
					#ax.text(meddelay, np.max(y)/2.0, "%.1f +/- %.1f" % (meddelay, stddelay), horizontalalignment = "center", color = rr.plotcolour)
					ax.annotate(delaytext, xy=(0.04, 0.75 - 0.1*irr),  xycoords='axes fraction', color = rr.plotcolour, fontsize=8)
					
					
			plt.xlim(histrange)
			
			# We increase ylim by 30% if it
			ylims = list(ax.get_ylim())
			if n == 2: # single panel
				ylims[1] *= 1.4
			else:
				ylims[1] *= 1.1
			ax.set_ylim(ylims)
			
			
			# Looked ok on big plots :
			#plt.annotate(delaylabel, xy=(0.03, 0.88),  xycoords='axes fraction', fontsize=12, color="black")
			if n != 2: # otherwise we have only one single panel
				plt.annotate(delaylabel, xy=(0.05, 0.84),  xycoords='axes fraction', fontsize=12, color="black")
			
			if trueshifts != None:
				truedelay = trueshifts[i] - trueshifts[j]
				plt.axvline(truedelay, color="gray", linestyle="--", dashes=(3, 3), zorder=-20)
	
	if dataout:
		for rr in rrlist:
			
			dc = delaycontainer(data = rr.tmpdata, name = rr.name, plotcolour = rr.plotcolour, objects=labels[:])	
			pycs.gen.util.writepickle(dc, "%s_delays.pkl" % (rr.autoname))
			rr.tmpdata = None
	
	labelspacetop = 0.0
	labelspaceright = 0.0
	if n == 2:
		labelspacetop = 0.04
		labelspaceright = 0.04
	
	for irr, rr in enumerate(rrlist):
	
		if niceplot:
			labeltxt = "%s" % (getattr(rr, 'name', 'NoName'))
			plt.figtext(x = right - labelspaceright, y = top - labelspacetop - txtstep*irr, s = labeltxt, verticalalignment="top", horizontalalignment="right", color=rr.plotcolour)	
		else:
			labeltxt = "%s (%s, %i) " % (getattr(rr, 'name', 'NoName'), "Truth" if rr.plottrue else "Measured", rr.tsarray.shape[0])
			plt.figtext(x = right - labelspaceright, y = top - labelspacetop - txtstep*irr, s = labeltxt, verticalalignment="top", horizontalalignment="right", color=rr.plotcolour)	

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


def newcovplot(rrlist, r=5, rerr=3, nbins = 10, nbins2d=3, binclip=False, binclipr=10.0, figsize=(20, 12), left=0.06, right=0.97, top=0.99, bottom=0.08, wspace=0.3, hspace=0.3, method='indepbin', detailplots=True, verbose=True):

	assert (method in ['depbin', 'indepbin'])

	nimages = rrlist[0].nimages()
	imginds = np.arange(nimages)
	labels = rrlist[0].labels

	couplelist = [(i, j) for j in imginds for i in imginds if i > j]
	ncouples = len(couplelist)
	# print couplelist

	tderrsdicts = []
	# rrlist is just a list of rr, we treat them one after the other
	for rr in rrlist:
		# for each rr, we compute the error from the true delay
		truetsslist = rr.truetsarray
		tsslist = rr.tsarray-truetsslist
		for ind, tss in enumerate(tsslist):
			tderrs = []
			truetds = []
			for (i, j) in couplelist:
				tderrs.append(tss[i]-tss[j])
				truetds.append(truetsslist[ind][i]-truetsslist[ind][j])
			tderrsdicts.append({"tderrs": tderrs, "truetds": truetds})

	#tderrsdict contains the errors on the true delays, as well as the true delays for each simulation

	#todo: we want to make a control plot out of this procedure"
	# figure 1: general covariance plot for each pair of delays. Diagonal elements are the same than newdelayplot, off-diagonal elements are covariance for all the runresults
	allcovplot = plt.figure(figsize=figsize)
	allcovplot.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)
	# figure 2: covariance computed in each bin, for each pair of delays. Diagonal elements are the same than newdelayplot, off diagonal elements are colored tiles of covariance per true delays bins, with points overplotted.
	bincovplot = plt.figure(figsize=figsize)
	bincovplot.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)
	# figure 3: for each pair, plot the covariance in each bin. One figure per pair
	# todo: define me below !!

	axisNum = 0

	# create the empty covariance matrix
	covmat = []
	for ind in range(len(couplelist)):
		covmat.append(np.zeros(len(couplelist)))
	indepbins = np.zeros(len(couplelist))
	depbins = np.zeros(len(couplelist))
	rranges = np.zeros(len(couplelist))
	for ii, i in enumerate(couplelist): # (0, 1), (0, 2) ...

		xtderrs = [tderrsdict["tderrs"][ii] for tderrsdict in tderrsdicts]
		xtruetds = [tderrsdict["truetds"][ii] for tderrsdict in tderrsdicts]
		maxx = np.max(xtruetds)
		minx = np.min(xtruetds)


		### fill the diagonal element

		ax1 = allcovplot.add_subplot(ncouples, ncouples, 6*ii + (ii+1))
		ax2 = bincovplot.add_subplot(ncouples, ncouples, 6*ii + (ii+1))
		majorLocator = MultipleLocator(1.0)

		for ax in [ax1, ax2]:
			ax.yaxis.set_major_locator(majorLocator)
			ax.xaxis.set_major_locator(MaxNLocator(5))

		# way 1 - binning independent of xtruedelays distribution. User choose the plot range. Similar to newdelayplot()
		reftrueshifts = np.round(rrlist[0].gettruets()["center"])
		reftruedelay = reftrueshifts[i[0]] - reftrueshifts[i[1]]
		plotrange = (reftruedelay - r, reftruedelay + r)
		binlims = np.linspace(plotrange[0], plotrange[1], nbins + 1)

		# If we want to compare to newdelayplot():
		# xtruetds = truedelays
		# xtderrs = resis

		# needed for binvals:
		xtderrs = np.array(xtderrs)

		digitized = np.digitize(xtruetds, binlims)
		binvals = [xtderrs[digitized == bini] for bini in range(1, len(binlims))]
		binstds = map(np.std, binvals)
		binmeans = map(np.mean, binvals)

		if binclip:
				for (bini, binvalarray) in enumerate(binvals):

					keep = np.logical_and(binvalarray < binclipr, binvalarray > -binclipr)
					if np.sum(keep == False) != 0:
						print "Kicking %i points." % (np.sum(keep == False))
					binvals[bini] = binvalarray[keep]
				binstds = map(np.std, binvals)
				binmeans = map(np.mean, binvals)

		syserror = np.max(np.fabs(binmeans))
		randerror = np.max(binstds)
		toterror = np.sqrt(syserror*syserror + randerror*randerror)
		indepbins[ii] = toterror

		# Plot the result !
		line = np.linspace(plotrange[0], plotrange[1], 100)
		zeros = np.zeros(100)
		delaylabel="%s%s" % (labels[i[1]], labels[i[0]])
		width = binlims[1] - binlims[0]
		errorrange=1.5
		for ax in [ax1, ax2]:
			ax.plot(line, zeros, color="black", lw=0.5)
			ax.bar(binlims[:-1], binmeans, yerr=binstds, width=width, color=rr.plotcolour, ecolor=rr.plotcolour, error_kw={"capsize":2.5, "capthick":0.5, "markeredgewidth":0.5}, edgecolor=rr.plotcolour, alpha = 0.2)
			ax.set_ylim((-errorrange, errorrange))
			if figsize[0] > 8:
				ax.annotate(delaylabel, xy=(0.9, 0.05),  xycoords='axes fraction', ha="center") # x axis
			else:
				ax.annotate(delaylabel, xy=(0.78, 0.08),  xycoords='axes fraction', ha="center")
			ax.set_xlim(plotrange)


		# way 2 - binning dependent on the xtruedelays samples: min and max vals corresponds to the extremas of xtruedelays distribution

		xbinvals = np.linspace(minx, maxx, num=nbins+1, endpoint=True)
		rranges[ii] = maxx-minx

		binmeans = []
		binstds = []
		for indx, xbinval in enumerate(xbinvals[:nbins]):
			subsamples = []
			for (ind, xtruetd) in enumerate(xtruetds):
				if xtruetd > xbinval and xtruetd < xbinvals[indx+1]:
					subsamples.append(xtderrs[ind])
			binmeans.append(np.mean(subsamples))
			binstds.append(np.std(subsamples))
		syserror = np.max(np.fabs(binmeans))
		randerror = np.max(binstds)
		toterror = np.sqrt(syserror*syserror + randerror*randerror)
		depbins[ii] = toterror

		# We let the user choose which method he prefers
		# Dear user, be CAREFUL with your choice !

		if method == 'depbin':
			if ii == 0 and verbose : print "You chose a binning depending on the sample values"
			covmat[ii][ii] = depbins[ii]
		elif method == 'indepbin':
			if ii == 0 and verbose : print "You chose a binning independent of the sample values"
			covmat[ii][ii] = indepbins[ii]


		### fill the off-diagonal elements
		for jj, j in enumerate(couplelist):

			axisNum += 1
			if (ii == 0) or (jj == ncouples-1) :
				continue # No plot
			if jj >= ii:
				continue

			if detailplots:
				# figure 3: for each pair, plot the covariance in each bin. One figure per pair
				bincovplot2 = plt.figure(figsize=figsize)
				bincovplot2.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)

			ytderrs = [tderrsdict["tderrs"][jj] for tderrsdict in tderrsdicts]
			ytruetds = [tderrsdict["truetds"][jj] for tderrsdict in tderrsdicts]

			ax1 = allcovplot.add_subplot(ncouples, ncouples, axisNum)
			ax2 = bincovplot.add_subplot(ncouples, ncouples, axisNum)
			majorLocator = MultipleLocator(2.0)
			for ax in [ax1]:
				ax.set_xlim(-rerr, rerr)
				ax.set_ylim(-rerr, rerr)
				ax.xaxis.set_major_locator(majorLocator)
				ax.yaxis.set_major_locator(majorLocator)
				ax.axhline(0, color="black")
				ax.axvline(0, color="black")


			## binning independent of xtrudelays and ytruedelays distribution.
			xbinlims2d = np.linspace(plotrange[0], plotrange[1], nbins2d + 1)

			yreftruedelay = reftrueshifts[j[0]] - reftrueshifts[j[1]]
			yplotrange = (yreftruedelay - r, yreftruedelay + r)
			ybinlims2d = np.linspace(yplotrange[0], yplotrange[1], nbins2d + 1)

			covsindep=[]
			for indx, xbinlim in enumerate(xbinlims2d[:nbins2d]):
				for indy, ybinlim in enumerate(ybinlims2d[:nbins2d]):
					subsamples = []
					for (ind, xtruetd), ytruetd in zip(enumerate(xtruetds), ytruetds):
						if xtruetd > xbinlim and xtruetd < xbinlims2d[indx+1] and ytruetd > ybinlim and ytruetd < ybinlims2d[indy+1]:
							subsamples.append((xtderrs[ind], ytderrs[ind]))

					#print len(subsamples), len(subsamples[0]), subsamples[0]
					covval = np.cov(subsamples, rowvar=False)[0][1]
					xcoord = xbinlim + (xbinlims2d[indx+1]-xbinlim)/2
					ycoord = ybinlim + (ybinlims2d[indy+1]-ybinlim)/2
					ax2.annotate("%.2f" % covval, xy=(xcoord, ycoord), ha="center")
					maxval=0.5
					alpha = min(np.abs(covval/maxval), 1.0)
					from matplotlib.patches import Rectangle
					rect = Rectangle((xbinlim, ybinlim), xbinlims2d[indx+1]-xbinlim, ybinlims2d[indy+1]-ybinlim, color=rrlist[0].plotcolour, alpha=alpha)
					ax2.add_patch(rect)

					if detailplots:
						# add a axes on the new figure for each bin, and plot the errors
						# mapping the maptlotlib indice is a bit tricky:
						# if we use nbins2dx and nbins2dy: nbins2dx*nbins2dy - (nbins2dx-1-indx) - (nbins2dy*indy)
						spind = nbins2d*nbins2d - (nbins2d-1-indx) - (nbins2d*indy)
						ax3 = bincovplot2.add_subplot(nbins2d, nbins2d, spind)
						ax3.set_xlim(-rerr, rerr)
						ax3.set_ylim(-rerr, rerr)
						showdensity = True
						bins = 10
						if showdensity:
							cmap = colors.LinearSegmentedColormap.from_list('custom', ['white', rrlist[0].plotcolour],gamma=1.0)
							ax3.hexbin([s[0] for s in subsamples], [s[1] for s in subsamples], gridsize=bins, extent=(-rerr, rerr, -rerr, rerr), mincnt=1, cmap=cmap, edgecolor="none")
						showpoints=True
						if showpoints:
							ax3.scatter([s[0] for s in subsamples], [s[1] for s in subsamples], s=5, facecolor=rrlist[0].plotcolour, lw=0)
						showcontour=True
						if showcontour:
							H, xedges, yedges = np.histogram2d([s[0] for s in subsamples], [s[1] for s in subsamples], range=[[-r, r], [-r, r]], bins=(bins, bins))
							extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
							data = np.vstack((xtderrs, ytderrs))
							kde = scipy.stats.kde.gaussian_kde(data)
							grid = np.mgrid[-r:r:1j*bins, -r:r:1j*bins]
							grid_coords = np.append(grid[0].reshape(-1,1),grid[1].reshape(-1,1),axis=1)
							z = kde(grid_coords.T)
							z = z.reshape(bins,bins)
							levels = [np.max(z)*0.45]
							cset = ax3.contour(grid[0], grid[1], z, levels=levels, origin="lower", colors=rrlist[0].plotcolour, extent=extent, linewidth=0.5)

						xdelaylabel="%s%s [%.1f , %.1f]" % (labels[i[1]], labels[i[0]], xbinlim, xbinlims2d[indx+1])
						ydelaylabel="%s%s [%.1f , %.1f]" % (labels[j[1]], labels[j[0]], ybinlim, ybinlims2d[indy+1])

						if figsize[0] > 8:
							ax3.annotate(xdelaylabel, xy=(0.85, 0.05),  xycoords='axes fraction', ha="center")
							ax3.annotate(ydelaylabel, xy=(0.04, 0.85),  xycoords='axes fraction', ha="left", rotation=90.0)

					covsindep.append(covval)
			detailplots = False
			mincovindep = np.min(covsindep)
			maxcovindep = np.max(covsindep)

			if abs(mincovindep) > maxcovindep:
				covindep = mincovindep
			else:
				covindep = maxcovindep

			#plotting ax2 uses the 2d binning
			for ind, xbinlim in enumerate(xbinlims2d):
				ax2.axvline(xbinlim, linestyle='--', color='black', alpha=0.5)
				ax2.axhline(ybinlims2d[ind], linestyle='--', color='black', alpha=0.5)
			showpoints=True
			if showpoints:
				ax2.scatter(xtruetds, ytruetds, s=2, facecolor=rrlist[0].plotcolour, lw=0)

			ax2.set_xlim(plotrange)
			ax2.set_ylim(yplotrange)

			# plotting ax1 is pretty basic, that's only the points
			showdensity = True
			bins = 10
			if showdensity:
				cmap = colors.LinearSegmentedColormap.from_list('custom', ['white', rrlist[0].plotcolour],gamma=1.0)
				ax1.hexbin(xtderrs, ytderrs, gridsize=bins, extent=(-rerr, rerr, -rerr, rerr), mincnt=1, cmap=cmap, edgecolor="none")

			showpoints=False
			if showpoints:
				ax1.scatter(xtderrs, ytderrs, s=2, facecolor=rrlist[0].plotcolour, lw=0)

			showcontour=True
			if showcontour:
				H, xedges, yedges = np.histogram2d(xtderrs, ytderrs, range=[[-r, r], [-r, r]], bins=(bins, bins))
				extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]

				data = np.vstack((xtderrs, ytderrs))
				kde = scipy.stats.kde.gaussian_kde(data)

				grid = np.mgrid[-r:r:1j*bins, -r:r:1j*bins]
				grid_coords = np.append(grid[0].reshape(-1,1),grid[1].reshape(-1,1),axis=1)

				z = kde(grid_coords.T)
				z = z.reshape(bins,bins)

				levels = [np.max(z)*0.45]
				cset = ax1.contour(grid[0], grid[1], z, levels=levels, origin="lower", colors=rrlist[0].plotcolour, extent=extent, linewidth=0.5)

			xdelaylabel="%s%s" % (labels[i[1]], labels[i[0]])
			ydelaylabel="%s%s" % (labels[j[1]], labels[j[0]])

			if figsize[0] > 8:
				ax1.annotate(xdelaylabel, xy=(0.9, 0.05),  xycoords='axes fraction', ha="center") # x axis
				ax1.annotate(ydelaylabel, xy=(0.06, 0.85),  xycoords='axes fraction', ha="left", rotation=90.0) # y axis
			else:
				ax1.annotate(xdelaylabel, xy=(0.78, 0.08),  xycoords='axes fraction', ha="center") # x axis
				ax1.annotate(ydelaylabel, xy=(0.08, 0.76),  xycoords='axes fraction', ha="left", rotation=90.0) # y axis

			## binning dependent of true delays, for comparision
			xbinvals = np.linspace(minx, maxx, num=nbins2d+1, endpoint=True)
			maxy = np.max(ytruetds)
			miny = np.min(ytruetds)
			ybinvals = np.linspace(miny, maxy, num=nbins2d+1, endpoint=True)

			covsdep=[]
			for indx, xbinval in enumerate(xbinvals[:nbins2d]):
				for indy, ybinval in enumerate(ybinvals[:nbins2d]):
					subsamples = []
					for (ind, xtruetd), ytruetd in zip(enumerate(xtruetds), ytruetds):
						if xtruetd > xbinval and xtruetd < xbinvals[indx+1] and ytruetd > ybinval and ytruetd < ybinvals[indy+1]:
							subsamples.append((xtderrs[ind], ytderrs[ind]))

					#TODO: due to the non-uniform sampling of the simulated true tds, some regions of the truetd_x vs truetd_y are rather empty (less than 10 samples). Should we i) increase the number of simulated samples, ii) discard these regions from the analysis, iii) transfer these samples to the nearest bin  ?
					#print len(subsamples), len(subsamples[0]), subsamples[0]

					covsdep.append(np.cov(subsamples, rowvar=False)[0][1])
			mincovdep = np.min(covsdep)
			maxcovdep = np.max(covsdep)

			if abs(mincovdep) > maxcovdep:
				covdep = mincovdep
			else:
				covdep = maxcovdep

			if method == "depbin":
				covmat[ii][jj] = covdep
				covmat[jj][ii] = covdep
			elif method == "indepbin":
				covmat[ii][jj] = covindep
				covmat[jj][ii] = covindep

			print "-"*15
			print i, j
			print covdep, covindep

	plt.show()
	sys.exit()
	# now let's compare indepbins and depbins
	if verbose:
		print "-"*35
		print "nbins = %i" % nbins
		print "indepbins - r = %.1f" % r
		print "depbins - r(max-min) =", np.mean(rranges)
		print "-"*35
		print "pair - indepbins - depbins - diff"
		print "-"*35
		print "AB - %.2f - %.2f - %.1f%%" % (indepbins[0], depbins[0], (max(indepbins[0], depbins[0])-min(indepbins[0], depbins[0])) / max(indepbins[0], depbins[0])*100)
		print "AC - %.2f - %.2f - %.1f%%" % (indepbins[1], depbins[1], (max(indepbins[1], depbins[1])-min(indepbins[1], depbins[1])) / max(indepbins[1], depbins[1])*100)
		print "BC - %.2f - %.2f - %.1f%%" % (indepbins[3], depbins[3], (max(indepbins[3], depbins[3])-min(indepbins[3], depbins[3])) / max(indepbins[3], depbins[3])*100)
		print "AD - %.2f - %.2f - %.1f%%" % (indepbins[2], depbins[2], (max(indepbins[2], depbins[2])-min(indepbins[2], depbins[2])) / max(indepbins[2], depbins[2])*100)
		print "BD - %.2f - %.2f - %.1f%%" % (indepbins[4], depbins[4], (max(indepbins[4], depbins[4])-min(indepbins[4], depbins[4])) / max(indepbins[4], depbins[4])*100)
		print "CD - %.2f - %.2f - %.1f%%" % (indepbins[5], depbins[5], (max(indepbins[5], depbins[5])-min(indepbins[5], depbins[5])) / max(indepbins[5], depbins[5])*100)
		print "-"*35

	return covmat




def measvstrue(rrlist, r=10.0, nbins = 10, plotpoints=True, plotrods=True, ploterrorbars=True, sidebyside=True, errorrange=None, binclip=False, binclipr=10.0, title=None, figsize=(10, 6), left = 0.06, right=0.97, top=0.99, bottom=0.08, wspace=0.15, hspace=0.3, txtstep=0.04, majorticksstep=2, displayn=True, filename=None, dataout=False):
	"""
	
	Plots measured delays versus true delays
	
	:param r: radius of simulation input delays to plot (x axis range)
	:param nbins: number of bins for the bar plot within this range.
	:param plotpoints: should I plot the points (scatter plot) ?
	:param plotrods: should I plot the avg within each bin ?
	:param ploterrorbars: should I add errorbars upon the bar plot ?
	:param sidebyside: should I plot bars side by side, or overplot them ?
	:param errorrange: radius of measurement errors to plot (y axis range). You can also give a tuple (low, high), to make asymetric plots.
	:param binclip: should I clip errors larger than binclipr days (catastrophic failures of methods) ?
	:param binclipr: see binclip ...
	
	"""
		
	n = rrlist[0].nimages()
	
	labels = rrlist[0].labels
	
	# To get some fixed ranges for the histograms, we will use the first element of rrlist.

	reftrueshifts = np.round(rrlist[0].gettruets()["center"])
	#@todo: WAAARNING ! Depending on the shape your rrlist (is it a 1x1000 runresults or 50x20 runresults), reftrueshift will have different values, impacting the final determination of the systematic and random error you compute. This can lead to a variation >10% on the final error !!!! DO SOMETHING !!!
	#print len(rrlist), rrlist[0].gettruets()["center"]
	#sys.exit()

	for rr in rrlist:
		if rr.labels != labels:
			raise RuntimeError("Don't ask me to overplot runresults of different curves !")
		#if not np.all(rr.gettruets()["center"] == reftrueshifts):
		#	print "Warning : I use the trueshift of the first rr to set the ranges."
		
		rr.trues = rr.gettruets() # To store this, avoids calculating it again and again.
		
		rr.tmpdata = []
	
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
			ax.yaxis.set_minor_locator(MultipleLocator(1.0))
			
			reftruedelay = reftrueshifts[i] - reftrueshifts[j]
			plotrange = (reftruedelay - r, reftruedelay + r)
			
			# Identity line :
			line = np.linspace(plotrange[0], plotrange[1], 100)
			zeros = np.zeros(100)
			plt.plot(line, zeros, color="black", lw=0.5)
			
			# Delay label
			delaylabel="%s%s" % (labels[j], labels[i])
			
			# Preparing the bins :
			binlims = np.linspace(plotrange[0], plotrange[1], nbins + 1)


			for irr, rr in enumerate(rrlist): # We go through the different runresult objects
				# We will express the delays "i - j"
				truedelays = rr.truetsarray[:,i] - rr.truetsarray[:,j]
				measdelays = rr.tsarray[:,i] - rr.tsarray[:,j]

				resis = measdelays-truedelays

				# A simple scatter plot of the residues :
				if plotpoints:
					ax.scatter(truedelays, resis, s=2, facecolor=rr.plotcolour, lw = 0)

				# We bin those :
				digitized = np.digitize(truedelays, binlims)
				
				binvals = [resis[digitized == bini] for bini in range(1, len(binlims))]
				binstds = map(np.std, binvals)
				binmedians = map(np.median, binvals)
				binmeans = map(np.mean, binvals)

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

				# We save the maximum sys and ran error :

				syserror = np.max(np.fabs(binmeans))
				randerror = np.max(binstds)
				toterror = np.sqrt(syserror*syserror + randerror*randerror)

				bias = np.mean(binmeans) # The signed bias
				rr.tmpdata.append({
					"label":delaylabel,
					"sys":syserror,
					"ran":randerror,
					"tot":toterror,
					"bias":bias
					})
				
				#binmeans = [np.median(resis[digitized == bini]) for bini in range(1, len(binlims))]
				#binstds = [np.std(resis[digitized == bini]) for bini in range(1, len(binlims))]
				
				width = binlims[1] - binlims[0]
				
				if plotrods:
					if not sidebyside:
						if ploterrorbars:
							ax.bar(binlims[:-1], binmeans, yerr=binstds, width=width, color=rr.plotcolour, ecolor=rr.plotcolour, error_kw={"capsize":2.5, "capthick":0.5, "markeredgewidth":0.5}, edgecolor=rr.plotcolour, alpha = 0.2)
						else:
							ax.bar(binlims[:-1], binmeans, width=width, color=rr.plotcolour, edgecolor=rr.plotcolour, alpha = 0.2)
					else:
						width = width/len(rrlist)
						squeezefactor = 1.0
						plotwidth = squeezefactor * width
						offset = width * (1.0-squeezefactor)/2.0
						
						if ploterrorbars:
							ax.bar(binlims[:-1] + offset + irr*plotwidth, binmeans, yerr=binstds, width=plotwidth, color=rr.plotcolour, ecolor=rr.plotcolour, error_kw={"capsize":2.5, "capthick":0.5, "markeredgewidth":0.5}, edgecolor=rr.plotcolour, alpha = 0.35, linewidth=0)
						else:
							ax.bar(binlims[:-1] + offset + irr*plotwidth, binmeans, width=plotwidth, color=rr.plotcolour, edgecolor=rr.plotcolour, alpha = 0.35)
					
				# That's it for the different runresult objects, back to the common stuff for this particular panel :
			
			if sidebyside:
				for binlim in binlims:
					plt.axvline(binlim, lw=0.5, color="#AAAAAA", zorder=-30)
				
			
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
				if hasattr(errorrange, '__iter__'): # then its a tuple or list
					plt.ylim((errorrange[0], errorrange[1]))
				else:
					plt.ylim((-errorrange, errorrange))
			
			if n != 2: # otherwise we have only 1 delay and panel
				plt.annotate(delaylabel, xy=(0.03, 0.88-txtstep),  xycoords='axes fraction', fontsize=12, color="black")
			
			# That's it for this panel, back to the total figure :
	
	if dataout:
		for rr in rrlist:
			
			dc = delaycontainer(data = rr.tmpdata, name = rr.name, plotcolour = rr.plotcolour, objects=labels[:])	
			pycs.gen.util.writepickle(dc, "%s_errorbars.pkl" % (rr.autoname))
			rr.tmpdata = None
	
	labelspacetop = 0.0
	labelspaceright = 0.0
	if n == 2:
		labelspacetop = 0.04
		labelspaceright = 0.04

	for irr, rr in enumerate(rrlist):
	
		if displayn:
			labeltxt = "%s (%i) " % (getattr(rr, 'name', 'NoName'), rr.tsarray.shape[0])
		else:
			labeltxt = "%s" % (getattr(rr, 'name', 'NoName'))
		plt.figtext(x = right - labelspaceright, y = top - labelspacetop - txtstep*irr, s = labeltxt, verticalalignment="top", horizontalalignment="right", color=rr.plotcolour)	

	if title != None:
		plt.figtext(x = 0.5, y = 0.95, s = title, horizontalalignment="center", color="black", fontsize=18)	

	if filename==None:
		plt.show()
	else:
		plt.savefig(filename)


def covplot(rrlist, showpoints=False, showcontour=True, showdensity=False, fractionalresiduals=False, bins=50, smoothing=0.0, figsize=(12, 12), left=0.02, right=0.98, bottom=0.02, top=0.98, wspace=0.05, hspace=0.05, r=5.0, title=None, txtstep=0.04, filename=None):
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
				if fractionalresiduals:
					iresis = (imeasdelays - itruedelays)/itruedelays
				else:
					iresis = imeasdelays - itruedelays
				
				jtruedelays = rr.truetsarray[:,j[0]] - rr.truetsarray[:,j[1]]
				jmeasdelays = rr.tsarray[:,j[0]] - rr.tsarray[:,j[1]]
				if fractionalresiduals:
					jresis = (jmeasdelays - jtruedelays)/jtruedelays
				else:
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
	print "I'M OBSOLETE, PLEASE USE NEWDELAYPLOT !"

	n = datarrlist[0].nimages()
	
	labels = datarrlist[0].labels
	
	# To get some fixed centers for the bins, we will use the first element of simrrlist.
	# If you change this, check also meanvstrue !!!
	reftrueshifts = np.round(simrrlist[0].gettruets()["center"])
	
	# Center for the delay plot :
	#plotshifts = reftrueshifts
	plotshifts = 0.5 * (np.max(np.array([rr.getts()["center"] for rr in datarrlist]), axis=0) + np.min(np.array([rr.getts()["center"] for rr in datarrlist]), axis=0))
	
	
	
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
		
		
			# We determine the bin range :
			reftruedelay = reftrueshifts[i] - reftrueshifts[j]
			
			
			# and for the plot :
			plotdelay = plotshifts[i] - plotshifts[j]
			plotrange = (plotdelay - rplot, plotdelay + rplot)
			
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

	print "I'M OBSOLETE, PLEASE USE NEWDELAYPLOT !"


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