"""
Stuff to manipulate time delay estimations from different techniques, specific to the TDC.
"""

import pycs.gen.lc
import numpy as np
import sys
import csv
import glob
import os
import datetime


class Estimate:
	"""
	Class to hold delay estimates for TDC curves, as obtained by various delay estimators.
	"""
	
	def __init__(self, set="tdc0", rung=0, pair=1, method="None", methodpar="None", td=0.0, tderr=0.0, ms=0.0, confidence=0, timetaken=0.0):
		"""

		:param set: What set (e.g., "tdc0") the curves are from
		:param rung: rung
		:param pair: pair
		
		:param method: Name of the method
		:param methodpar: String containing parameters of this method, or username...

		:param td: time delay point estimate
		:param tderr: 1sigma error estimate on the time delay
		:param ms: magnitude shift estimate
		
		:param confidence: confidence level. 0 = not estimated, 1 = doubtless, 2 = plausible, ...
		:param timetaken: seconds it took to get this estimate

		"""
	
		self.set = set
		self.rung = rung
		self.pair = pair
		
		self.method = method
		self.methodpar = methodpar

		self.td = td
		self.tderr = tderr
		self.ms = ms
		
		self.confidence = confidence
		self.timetaken = timetaken

	def __str__(self):
		return "%s %s (%s rung %i pair %i): %.2f +/- %.2f, conf. %i" % (self.method, self.methodpar, self.set, self.rung, self.pair, self.td, self.tderr, self.confidence)
	
	def aslist(self):
		return [self.set, self.rung, self.pair, self.method, self.methodpar, self.td, self.tderr, self.ms, self.confidence, self.timetaken]
	
	def check(self):
		
		#assert self.tderr >= 0.0
		if self.tderr < 0.0:
			raise RuntimeError("Negative error...")
	

	
def readcsv(filepath):
	"""
	Read a CSV file of estimates.
	You can specify a directory instead of a file, and I will read in all .csv files form there.
	"""
	
	if os.path.isdir(filepath):
		files = sorted(glob.glob(os.path.join(filepath, "*.csv")))
		estimates = []
		for file in files:
			estimates.extend(readcsv(file))
		return estimates
			
	else:
		f = open(filepath, 'rb')
		reader = csv.reader(f, delimiter=',', quotechar='"')
		estimates = []
		for row in reader:
			for i in [1, 2, 8]:
				row[i] = int(row[i])
			for i in [5, 6, 7, 9]:
				row[i] = float(row[i])
			estimates.append(Estimate(*row))
		f.close()
		print "Read %i estimates from %s" % (len(estimates), filepath)
		return estimates
        
	
def writecsv(estimates, filepath):
	"""
	Write a CSV file of estimates
	"""
	f  = open(filepath, "wb")
	writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
  	for est in estimates:
 		writer.writerow(est.aslist()) 
	f.close()
	print "Wrote %i estimates into %s" % (len(estimates), filepath)


def importfromd3cs(filepath, set="tdc0"):
	"""
	Reads a d3cs log file and returns the list of estimates
	"""
	
	logfile = open(filepath, 'r')
	lines = logfile.readlines()
	logfile.close()
	
	lines = [line.split(',') for line in lines]
	lines = [[element.strip() for element in line] for line in lines]
	
	estimates = [
		Estimate(set=set,
			rung = int(line[3]),
			pair = int(line[4]),
			method = "D3CS", 
			methodpar = line[2], 
			td = float(line[5]), 
			tderr = float(line[6]), 
			ms = float(line[7]),  
			confidence = int(line[9]), 
			timetaken = float(line[10]))
		for line in lines]

	return estimates
	

def group(estimates):
	"""
	Groups estimates by "quasar"
	In other words : takes a messy list of mixed estimates and returns a list of lists of estimates for a pair.
	"""
	for est in estimates:
		est.id = "%s_%i_%i" % (est.set, est.rung, est.pair)
	estids = sorted(list(set([est.id for est in estimates])))
	
	return [[est for est in estimates if est.id == estid] for estid in estids]
		
	
def checkunique(estimates):
	"""
	Checks that there is only one estimate per pair
	"""
	if len(estimates) != len(group(estimates)):
		raise RuntimeError("Your estimates are not unique !")
	
	
def writesubmission(estimates, filepath):
	"""
	Write a submissible TDC file from a list of estimates
	It takes td and tderr of estimates
		
	The delay td is in positive units (?)
	"""
	
	#if os.path.exists(filepath):
	#	print "WARNING, THAT FILE EXISTS !"
	#	return
	
	checkunique(estimates)
	tdcfile = open(filepath, "w")
	
	tdcfile.write("# TDC submission written by PyCS\n")
	tdcfile.write("# %s\n" % (datetime.datetime.now()))
	tdcfile.write("# \n")
	tdcfile.write("# (Some room for your comment...)\n")
	tdcfile.write("# \n")
	tdcfile.write("# \n")
	tdcfile.write("# datafile              dt      dterr\n")
		
	for est in estimates:
		est.check()
		name = "%s_rung%i_pair%i.txt" % (est.set, est.rung, est.pair)
		tdcfile.write("%s\t%.2f\t%.2f\n" % (name, est.td, est.tderr))
		
	tdcfile.close()	
	



def bigplot(estimates, plotpath = None, wholeset=False):
	import matplotlib.pyplot as plt
	for est in estimates:
		est.tmpid = "(%i, %i)" % (est.rung, est.pair)
	estids = sorted(list(set([est.tmpid for est in estimates])))
	
	#estids = estids[:40]
	
	def colour(est):
		if est.confidence==0: return "black"
		if est.confidence==1: return "blue"
		if est.confidence==2: return "green"
		if est.confidence==3: return "orange"
		if est.confidence==4: return "red"
		
	fig, axes = plt.subplots(nrows=len(estids), figsize=(10, 1.0*(len(estids))))
	fig.subplots_adjust(top=0.99, bottom=0.05, hspace=0.32)
   
	for ax, estid in zip(axes, estids):
		thisidests = [est for est in estimates if est.tmpid == estid]
		n = len(thisidests)
		
		tds = np.array([est.td for est in thisidests])
		tderrs = np.array([est.tderr for est in thisidests])
 		ys = np.arange(n)
 		colours = map(colour, thisidests)
		#ax.scatter(tds, ys)
		ax.errorbar(tds, ys, yerr=None, xerr=tderrs, fmt='.', ecolor="gray", capsize=3)
		ax.scatter(tds, ys, s = 50, c=colours, linewidth=0, zorder=20)#, cmap=plt.cm.get_cmap('jet'), vmin=0, vmax=4)
		
		for (est, y) in zip(thisidests, ys):
			ax.text(est.td, y+0.3, "  %s (%s)" % (est.method, est.methodpar), va='center', ha='left', fontsize=6)
		
		
		ax.set_ylim(-1, n)
		
		meantd = np.mean(tds)
		maxdist = np.max(np.fabs(tds - meantd))
		if maxdist < 200:
			tdr = 200
		else:
			tdr = maxdist*1.2
		ax.set_xlim(meantd - tdr, meantd + tdr)
		
  		pos = list(ax.get_position().bounds)
   		x_text = pos[0] - 0.01
		y_text = pos[1] + pos[3]/2.
		fig.text(x_text, y_text, estid, va='center', ha='right', fontsize=14)
		
	for ax in axes:
		ax.set_yticks([])

	#cbar = plt.colorbar(sc, cax = axes[0], orientation="horizontal")
	#cbar.set_label('Confidence')

	if plotpath:
		plt.savefig(plotpath)
	else:
		plt.show()



def interactivebigplot(estimates, plotpath = None):
	
	import matplotlib.pyplot as plt
	import matplotlib.axes as maxes
	from matplotlib.widgets import Button
	
	for est in estimates:
		est.tmpid = "(%i, %i)" % (est.rung, est.pair)
	estids = sorted(list(set([est.tmpid for est in estimates])))
	
	
	# init lists for interactive plotting...
				
	buttonstod3cs=[]
	buttonstoshow=[]
	
				
	def colour(est):
		if est.confidence==0: return "black"
		if est.confidence==1: return "blue"
		if est.confidence==2: return "green"
		if est.confidence==3: return "orange"
		if est.confidence==4: return "red"
		
	fig, axes = plt.subplots(nrows=len(estids), figsize=(10, 1.0*(len(estids))))	
	fig.subplots_adjust(top=0.99, bottom=0.05, hspace=0.32)

	
	for ax, estid in zip(axes, estids):
	
		# resize ax and add a new box we will fill with other informations (see below)
		# ugly stuff (but fast), sorry Malte
		
		bbox = ax.get_position()
		points = bbox.get_points()
		
		xright = points[1][0] 
		points[1][0] = 0.73
		bbox.set_points(points)
		ax.set_position(bbox)
		
		hspace = 0.02
		width  = xright-hspace-points[1][0]
		height = points[1][1]-points[0][1]
		
		rect = points[1][0]+hspace, points[0][1], width, height
		ax2 = fig.add_axes(rect)
		
		
				
		### arange and compute values to fill the left box (ax)	
		
		thisidests = [est for est in estimates if est.tmpid == estid]
		n = len(thisidests)
		
		tds = np.array([est.td for est in thisidests])
		tderrs = np.array([est.tderr for est in thisidests])
		meantd = np.mean(tds)

 		ys = np.arange(n)
		
		colours = map(colour, thisidests)
		
		
		#ax.scatter(tds, ys)
		ax.errorbar(tds, ys, yerr=None, xerr=tderrs, fmt='.', ecolor="grey", capsize=3)
		ax.scatter(tds, ys, s = 50, c=colours, linewidth=0, zorder=20)#, cmap=plt.cm.get_cmap('jet'), vmin=0, vmax=4)
		
		for (est, y) in zip(thisidests, ys):
			ax.text(est.td, y+0.3, "  %s (%s)" % (est.method, est.methodpar), va='center', ha='left', fontsize=6)
		
		
		ax.set_ylim(-1, n)
		
		maxdist = np.max(np.fabs(tds - meantd))
		if maxdist < 200:
			tdr = 200
		else:
			tdr = maxdist*1.2
		ax.set_xlim(meantd - tdr, meantd + tdr)
		
  		pos = list(ax.get_position().bounds)
   		x_text = pos[0] - 0.01
		y_text = pos[1] + pos[3]/2.
		fig.text(x_text, y_text, estid, va='center', ha='right', fontsize=14)
		ax.set_yticks([])
		
		
		
		### OK, now fill the right box (ax2)
		
		# some inits...
		conflevelmin = 3
		maxtolerr = 8
		
		thisidests_disc = [est for est in estimates if est.tmpid == estid and est.confidence <= conflevelmin]
		tds_disc = np.array([est.td for est in thisidests_disc])
	
		meantd_disc = np.mean(tds_disc)
		meantderr1 = np.std(tds_disc)
		meantderr2 = np.std(tds_disc)/np.sqrt(len(tds_disc))
				
				
		if meantderr1 < maxtolerr:
			mcolor1 = 'black'
		else:	
			mcolor1 = 'red'
			
		if meantderr2 < maxtolerr:
			mcolor2 = 'black'
		else:	
			mcolor2 = 'red'			
		
		ax2.errorbar(meantd_disc,1,yerr=None, xerr=[meantderr1], fmt='.', ecolor=mcolor1, capsize=3)
		ax2.errorbar(meantd_disc,1,yerr=None, xerr=[meantderr2], fmt='.', ecolor=mcolor2, capsize=3)					
		ax2.scatter(meantd_disc,1, s=50, c=mcolor2, linewidth=0, zorder=20)
		
		
		from matplotlib.ticker import MaxNLocator
		ax2.xaxis.set_major_locator(MaxNLocator(4))
		
		ax2.set_xlim(meantd_disc - maxtolerr, meantd_disc + maxtolerr)		
		ax2.set_yticks([])
		

		### Interactive plotting options
		# We create buttons
		
		sizescale = 3.3
		
		axtod3cs = plt.axes([points[0][0], points[0][1], width/sizescale, height/sizescale])
		axtoshow = plt.axes([points[0][0], points[0][1]+height-height/sizescale, width/sizescale, height/sizescale])
				
		buttonstod3cs.append(Button(axtod3cs, 'D3CS'))
		buttonstoshow.append(Button(axtoshow, 'Show'))

		

	for buttontoshow, buttontod3cs, estid in zip(buttonstoshow, buttonstod3cs, estids):
	
		
		# weird way of doing things... but didn't find a cleaner way to have everything working... 
						
		rung = int(estid[1]) # Done this way to avoid mix up if I select my estimates in a weird order...
		pair = int(estid[4])
		
		class Goto:
	    		myrung = rung
			mypair = pair
	    		def show(self, event):
				show(estimates,self.myrung,self.mypair)
				
			def d3cs(self, event):
				d3cs(self.myrung,self.mypair)	
				
		goto = Goto()				
		buttontoshow.on_clicked(goto.show)
		buttontod3cs.on_clicked(goto.d3cs)
	
		
	
	#cbar = plt.colorbar(sc, cax = axes[0], orientation="horizontal")
	#cbar.set_label('Confidence')
	
	# add next/previous rung buttons
	
	''' Work in progress !
	
	axnext = plt.axes()
	axprev = plt.axes()
	
	buttonnext = Button(axnext,'Next Rung')
	buttonprev = Button(axprev,'Previous Rung')
	
		
	class ChangeRound:
    		ind = 0
    		def next(self, event):
        		self.ind += 1
        		i = self.ind % len(freqs)
        		ydata = np.sin(2*np.pi*freqs[i]*t)
        		l.set_ydata(ydata)
        		plt.draw()
	'''
	
	

	if plotpath:
		plt.savefig(plotpath)
	else:
		plt.show()

def show(estimates, rung, pair):

	'''
	display the rung / pair curve for each corresponding estimate 
	'''
	
	
	# I select only the estimates with the rung and pair desired
	estimates = [est for est in estimates if est.rung == rung and est.pair == pair]

	setlist=[]
		
	# import the curve from data
	filepath = 'data/rung%0i/tdc0_rung%0i_pair%0i.txt' % (rung,rung,pair)
		
	for est in estimates:	
		lcs = pycs.tdc.util.read(filepath)
		lcs[1].shifttime(est.td)
		lcs[1].shiftmag(est.ms)
		setlist.append([lcs,est.methodpar])
	
	pycs.gen.lc.multidisplay(setlist, showlegend = False, showdelays = True)		
		
	
def d3cs(rung,pair):
	'''
	Open d3cs with the rung/pair curve in your default browser
		
	'''
	import webbrowser
	cmd='http://www.astro.uni-bonn.de/~mtewes/d3cs/index.php?user=display&loadrung=%i&loadpair=%i' %(rung,pair)
	webbrowser.open(cmd)

	
	


