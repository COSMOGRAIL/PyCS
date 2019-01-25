"""
Plot functions. Partly rewritting some of the  sim.plot functions, as this will end up replacing it
"""

import numpy as np
import math, sys, os

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator

import scipy.ndimage

import pycs.gen.util




def delayplot(plotlist, rplot=7.0, autoobj= None, displaytext=True, hidedetails=False, showbias=True, showran=True, showerr=True, showlegend=True, text=None, figsize=(10, 6), left=0.06, right=0.97, top=0.99, bottom=0.08, wspace=0.10, hspace=0.15, txtstep=0.03, majorticksstep=2, filename=None, refgroup=None, legendfromrefgroup=False, centerdelays=None, ymin=0.2, hlines=None, blindness=False, horizontaldisplay=False, showxlabelhd=True, update_group_style = True, auto_radius =False,
			  tick_step_auto = True, legendx = 0.85, legendy_offset = 0.12 ):

	"""
	Plots delay measurements from different methods, telescopes, sub-curves, etc in one single plot.

	For this I use only ``Group`` objects, i.e. I don't do any "computation" myself. I can of course handle asymmetric errors.


	@param plotlist: list of Groups
	@param rplot: radius of delay axis, in days.

	@param displaytext: Show labels with technique names and values of delays

	@param hidedetails: Do not show (ran, sys) in labels

	@param refgroup: a group to be plotted as shaded vertical zones.

	@param legendfromrefdelays: if you want to display the refdelays name in the legend panel


	@param showbias: draws a little cross at the position of the delay "corrected" for the bias. Only if your group has a sys_errors attribute.


	@param showran: draws "minor" error bar ticks using the random error only. Only if your group has a ran_errors attribute.


	@param text:
		Text that you want to display, in the form : [line1, line2, line3 ...]
		where line_i is (x, y, text, kwargs) where kwargs is e.g. {"fontsize":18} and x and y are relative positions (from 0 to 1).

	@param blindness: Shift the measurements by their mean, so the displayed value are centered around 0

	@param horizontaldisplay: display the delay panels on a single line. Works only for three-delay containers.

	@param showxlabelhd: display or not the x label when horizontal display is True

	@param auto_radius : automatically adjust the xlim, if true, radius won't be used.

	@param tick_step_auto : automatically adjust the tickstep, if true, radius won't be used.



	"""

	"""
	# Some checks :
	objects = plotlist[0][0].objects
	for (delays, errors) in plotlist:
		if delays.objects != objects or errors.objects != objects:
			raise RuntimeError("Don't ask me to overplot stuff from different objects !")
	"""
	pairs = plotlist[0].labels
	if autoobj == None :
		objects = sorted(list(set("".join(pairs))))
	else :
		objects = autoobj

	n = len(objects)
	nmeas = len(plotlist)
	print "Objects : %s" % (", ".join(objects))

	if horizontaldisplay and n != 3:
		print "Horizontal display works only for three delays, you have %i" % n
		print "Switching back to regular display"
		horizontaldisplay = False

	"""
	for (delays, errors) in plotlist:
		if delays.plotcolour != errors.plotcolour:
			raise RuntimeError("Hmm, plotcolours of delays and errors don't correspond !")
		print "Delays : %s <-> Errors : %s" % (delays.name, errors.name)
	"""


	fig = plt.figure(figsize=figsize)
	fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)

	axisNum = 0
	print "#" * 80
	for i in range(n):  # A, B, C, D and so on
		for j in range(n):

			#print i, j

			if (i == 0) or (j == n - 1):
				continue  # No plot

			if not horizontaldisplay:
				axisNum += 1

			if j >= i:
				continue


			if horizontaldisplay:
				axisNum += 1
				ax = plt.subplot(1, n, axisNum)
			else:
				ax = plt.subplot(n - 1, n - 1, axisNum)

			# We will express the delays "i - j"
			delaylabel = "%s%s" % (objects[j], objects[i])
			print "           Delay %s" % (delaylabel)

			# To determine the plot range :
			paneldelays = []
			errors_up_list = []
			errors_down_list = []

			# Going throuh plotlist :

			if blindness:
				blinddelays = []
				for ipl, group in enumerate(plotlist):
					blindmedians.append([med for med, label in zip(group.medians, group.labels) if label == delaylabel][0])
					blindshift = np.mean(blinddelays)

			for ipl, group in enumerate(plotlist):

				# Getting the delay for this particular panel
				labelindex = group.labels.index(delaylabel)
				median = group.medians[labelindex]
				if blindness:
					median -= blindshift

				#error_up = [e for e, l in zip(group.errors_up, group.labels) if l == delaylabel][0]
				#error_down = [e for e, l in zip(group.errors_down, group.labels) if l == delaylabel][0]

				error_up = group.errors_up[labelindex]
				error_down = group.errors_down[labelindex]

				paneldelays.append(median)
				errors_up_list.append(error_up)
				errors_down_list.append(error_down)
				ypos = nmeas - ipl

				xerr = np.array([[error_down, error_up]]).T


				# define style attributes:

				# error line width
				if not hasattr(group, 'elinewidth'):
					group.elinewidth = 1.5

				# color
				if not hasattr(group, 'plotcolor'):
					group.plotcolor = "royalblue"

				# marker
				if not hasattr(group, 'marker'):
					group.marker = None

				# size of the marker
				if not hasattr(group, 'markersize') or update_group_style:
					if n > 2 :
						group.markersize = math.floor(-(4./12.) * (nmeas - 4) + 8)
					else :
						group.markersize = 8


				# size of the delay annoted on top of the measurement
				if not hasattr(group, 'labelfontsize')or update_group_style:
					if n > 2:
						group.labelfontsize = max(math.floor(-(8./12.)*(nmeas-4) + 18), 8)
					else :
						group.labelfontsize = 18

				# size of the legend
				if not hasattr(group, 'legendfontsize')or update_group_style:
					if n > 2  :
						group.legendfontsize = math.floor(-(4/12.)*(nmeas-4) + 16)
						txtstep = -(0.01/12.)*(nmeas-4) + 0.03
					else :
						group.legendfontsize = 16

				# extra properties: elinewidth, plotcolor, marker, markersize, labelfontsize, legendfontsize

				plt.errorbar([median], [ypos], yerr=None, xerr=xerr, fmt='-', ecolor=group.plotcolor, elinewidth=group.elinewidth, capsize=3, barsabove=False)

				if showran and hasattr(group, "ran_errors"):
					plt.errorbar([median], [ypos], yerr=None, xerr=group.ran_errors[labelindex], fmt='-', ecolor=group.plotcolor, elinewidth=0.5, capsize=2, barsabove=False)

				if showbias:
					plt.plot([median - group.sys_errors[labelindex]], [ypos], marker="x", markersize=group.markersize, markeredgecolor=group.plotcolor, color=group.plotcolor)

				if group.marker == None or group.marker == ".":
					plt.plot([median], [ypos], marker='o', markersize=group.markersize, markeredgecolor=group.plotcolor, color=group.plotcolor)
				else:
					plt.plot([median], [ypos], marker=group.marker, markersize=group.markersize,markeredgecolor=group.plotcolor, color=group.plotcolor)



				if hidedetails or (group.ran_errors[labelindex] < 0.001 and group.sys_errors[labelindex] < 0.001):  # Then we ommit to write them.
					delaytext = r"$%+.1f^{+%.1f}_{-%.1f}$" % (median, error_up, error_down)
				else:
					delaytext = r"$%+.1f^{+%.1f}_{-%.1f}\,(%.1f, %.1f)$" % (median, error_up, error_down, group.ran_errors[labelindex], group.sys_errors[labelindex])

				# if you want to hide the error...
				if not showerr:
					delaytext = r"$%+.1f$" % median

				if n == 2:  # For doubles, we include the technique name into the txt :
					delaytext = r"%s : " % (group.name) + delaytext

				if displaytext:
					ax.annotate(delaytext, xy=(median, ypos + 0.3), color=group.plotcolor,
					            horizontalalignment="center", fontsize=group.labelfontsize)

					print "%45s : %+6.2f + %.2f - %.2f" % (group.name, median, error_up, error_down)

			print "#" * 80

			# Now this panel is done. Some general settings :
			if centerdelays != None:
				centerdelay = centerdelays[delaylabel]
			else:
				centerdelay = np.median(paneldelays)

			if auto_radius :
				rplot = (np.max(errors_up_list) + np.max(errors_down_list)) / 2.0 * 2.5
			plt.xlim((centerdelay - rplot, centerdelay + rplot))
			plt.ylim((ymin, nmeas + 1.5))

			# General esthetics :
			if tick_step_auto :
				majorticksstep = max(2.0, int(rplot/5.0) )
			ax.get_yaxis().set_ticks([])
			minorLocator = MultipleLocator(1.0)
			majorLocator = MultipleLocator(majorticksstep)
			ax.xaxis.set_minor_locator(minorLocator)
			ax.xaxis.set_major_locator(majorLocator)

			# Blindness display options
			if blindness:
				xlabel = "Blind delay [day]"
			else:
				xlabel = "Delay [day]"

			plt.xticks(fontsize=15 - 2*(n-2))
			xlabelfontsize = 18

			if i == n - 1 and not horizontaldisplay:
				plt.xlabel(xlabel, fontsize=xlabelfontsize)
			elif horizontaldisplay:
				if showxlabelhd:
					plt.xlabel(xlabel, fontsize=xlabelfontsize)
				else:
					ax.get_xaxis().set_ticks([])

			if n != 2:  # otherwise only one panel, no need
				plt.annotate(delaylabel, xy=(0.03, 0.88 - txtstep), xycoords='axes fraction', fontsize=14, color="black")

			if refgroup != None:

				reflabelindex = refgroup.labels.index(delaylabel)

				refmedian = refgroup.medians[reflabelindex]
				referror_up = refgroup.errors_up[reflabelindex]
				referror_down = refgroup.errors_down[reflabelindex]

				plt.axvspan(refmedian - referror_down, refmedian + referror_up,facecolor="grey", alpha=0.25, zorder=-20, edgecolor="none", linewidth=0)

				plt.axvline(refmedian, color="grey", linestyle="--", dashes=(5, 5), lw=1.0, zorder=-20)

			if hlines != None:
				for hline in hlines:
					plt.axhline(hline, lw=0.5, color="gray", zorder=-30)

	# The "legend" :
	if showlegend:
		for ipl, group in enumerate(plotlist):
			line = "%s" % (group.name)

			plt.figtext(x=legendx, y=top - txtstep * ipl - legendy_offset, s=line, verticalalignment="top",
			            horizontalalignment="center", color=group.plotcolor,
			            fontsize=group.legendfontsize)  # for 3-delay plots
		if legendfromrefgroup and refgroup is not None:
			if not hasattr(refgroup, 'legendfontsize'):
				refgroup.legendfontsize = 16

			line = "%s" % (refgroup.name)
			plt.figtext(x=legendx, y=top - txtstep * len(plotlist) - legendy_offset, s=line, verticalalignment="top", horizontalalignment="center", color="grey", fontsize=refgroup.legendfontsize)

	# Generic text :
	if text != None:
		for line in text:
			plt.figtext(x=line[0], y=line[1], s=line[2], **line[3])

	if filename == None:
		plt.show()
	else:
		plt.savefig(filename)


def write_delays(group, write_dir=None, mode="GLEE"):
	"""
	Write the group linarized distributions into a txt file, to be used seamlessly by various modeling code. So far, only GLEE is implemented.

    Parameters
    ----------

	write_dir : string or None (default=None)
		Directory path of the output files. If None, use the current working directory (cwd)

	mode: string or None (default="GLEE")
		Defines how the output has to be written.
		GLEE: no header, positively increasing values, constant step
	"""

	if write_dir is None:
		_wd = "."
	else:
		_wd = write_dir

	for ind, (label, bins, vals) in enumerate(zip(group.labels, group.binslist, group.lins)):

		# label: AB, AC, BC, etc...
		# bins: edges of the histogram bins
		# vals: value of the measured pdf at the middle of each bin

		# rebase bins into a new vector whose values are the middle of the bins
		xs = [(bins[i] + bins[i+1])/2. for i, _ in enumerate(bins[:-1])]

		if mode is "GLEE":
			# force a positively increasing delays values
			if xs[1] < xs[0]:
				xs = xs[::-1]
				vals = vals[::-1]


			# assert constant step
			try:
				assert len(list(set(["%.3f" % (xs[i + 1] - xs[i]) for i in range(len(xs) - 1)]))) == 1
			except:
				print list(set(["%.3f" % (xs[i + 1] - xs[i]) for i in range(len(xs) - 1)]))
				raise AssertionError("The delay values step is not constant. This might be due to a rounding error: either change the digit precision of your delays, or use a different binning when linearizing.")


		# save the data in a txt file for easier use
		f = open(os.path.join(_wd, "%s_%s.txt" % (group.name, label)), "w")
		if mode is not "GLEE":
			f.write("Dt\tprob\n")
			f.write("==\t====\n")
		for x, val in zip(xs, vals):
			f.write("%.3f\t%.8f\n" % (x, val))
		f.close()

		# also save the revert delays (i.e. BA instead of AB)
		_reverted_xs = np.array(xs) * -1.0
		if len(label) == 2:
			newlabel = label[::-1]
		elif len(label) == 3:
			newlabel = label[2:]+label[:2]
		elif len(label) == 4:
			newlabel = label[2:]+label[:2]

		if mode is "GLEE":
			# force a positively increasing delays values
			if _reverted_xs[1] < _reverted_xs[0]:
				_reverted_xs = _reverted_xs[::-1]
				vals = vals[::-1]


		f = open(os.path.join(_wd, "%s_%s.txt" % (group.name, newlabel)), "w")
		if mode is not "GLEE":
			f.write("Dt\tprob\n")
			f.write("==\t====\n")
		for x, val in zip(_reverted_xs, vals):
			f.write("%.3f\t%.8f\n" % (x, val))
		f.close()
