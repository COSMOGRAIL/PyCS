"""
Stuff to manipulate time delay estimations from different techniques, specific to the TDC.
"""

import pycs.gen.lc
import numpy as np


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
	
	

def bigplot(estimates, plotpath = None):
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



def displaymeplease(estimates, rung, pair):

	'''
	display the rung / pair curve for each estimate corresponding
	'''
	
	
	# I select only the estimates with the rung and pair desired
	estimates = [est for est in estimates if est.rung == rung and est.pair == pair]

	setlist=[]
		
	# import the curve from data
	datapath = '/archive/vbonvin/disk1/LENSES/TDC0/data/'
	filepath = datapath + 'rung%0i/tdc0_rung%0i_pair%0i.txt' % (rung,rung,pair)
		
	for est in estimates:	
		lcs = pycs.tdc.util.read(filepath)
		lcs[1].shifttime(est.td)
		lcs[1].shiftmag(est.ms)
		setlist.append([lcs,est.methodpar])
	
	pycs.gen.lc.multidisplay(setlist, showlegend = False, showdelays = True)		
		
	


	
	


