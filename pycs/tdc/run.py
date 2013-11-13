"""
Wrapper stuff to run PyCS on TDC data
"""

import os
import pycs
import datetime
import numpy as np
import matplotlib.pyplot as plt

class Run:
	"""
	Class to contain a run of PyCS methods on a given lens
	"""
	
	def __init__(self, iniest, lca, lcb, plots=True, outdir=None):
		"""

		"""
			
		self.iniest = iniest
		self.lca = lca
		self.lcb = lcb
		
		self.plots = plots
		
		if outdir == None:
			outdir = os.getcwd()
		self.outdir = outdir
		
		if not os.path.isdir(self.outdir):
			os.mkdir(self.outdir)
			
		self.plotdir = os.path.join(self.outdir, "diagnostics")
		if not os.path.isdir(self.plotdir):
			os.mkdir(self.plotdir)
		self.datadir = os.path.join(self.outdir, "data")
		if not os.path.isdir(self.datadir):
			os.mkdir(self.datadir)
			
		self.logpath = os.path.join(self.plotdir, "log.txt")
		
		self.gogogo = True
		
		
	def check(self):
		if self.gogogo == False:
			return
			
	def rmlog(self):
		if os.path.exists(self.logpath):
			os.remove(self.logpath)

	def log(self, message):
		logfile = open(self.logpath, "a")
		logfile.write(message + "\n")
		logfile.close()
		print "### LOG ###: %s" % (message)

	def setup(self):
		
		self.rmlog()
		self.log("Setup at %s" % (str(datetime.datetime.now())))
		self.log("Data: %s" % ", ".join([str(self.lca), str(self.lcb)]))
		self.log("Initial estimate: %s" % str(self.iniest))
		
		# Set the initial shift
		self.lca.timeshift = 0.0
		self.lcb.timeshift = self.iniest.td
		pycs.tdc.util.setnicemagshift([self.lca, self.lcb])
		
				
	def fitsourcespline(self, sploptfct):
		self.check()
		# Adding some simple ML to B (required even just for the magshift !)
		pycs.gen.polyml.addtolc(self.lcb, nparams=1, autoseasonsgap = 100000.0)

		# And run a first spline
		self.sourcespline = sploptfct([self.lca, self.lcb])
		self.log("Fitted sourcespline, delay %.2f" % self.lcb.timeshift)
		
		relshift = np.fabs(self.lcb.timeshift - self.iniest.td) / self.iniest.tderr
		self.log("Relative time shift: %.2f" % relshift)
		
		if relshift > 5.0:
			self.log("That's too much, I stop")
			self.gogogo = False
		
		#pycs.gen.util.writepickle(self.sourcespline, os.path.join(self.datadir, "sourcespline.pkl"))
		
		# Visu
		pycs.gen.lc.display([self.lca, self.lcb], [self.sourcespline], figsize=(18, 6), filename=os.path.join(self.plotdir, "sourcespline.png"), verbose=False)
		
		

	def runobs(self, optfct, n=5):
		"""
		Run the optimizer n times on the real data.
		We draw the copies, run on them
		"""
		self.check()
		
		#simdir = os.path.join(self.datadir, "copies")
		#if os.path.isdir(simdir):
		#	os.remove(simdir)
		#pycs.sim.draw.multidraw([self.lca, self.lcb], onlycopy=True, n=n, npkl=1, simset="copies", simdir=simdir)
		
		self.log("Running %i times on the observations..." % n)
		lcsclist = []
		for i in range(n):
			
			# We make some copies to work with
			lcac = self.lca.copy()
			lcbc = self.lcb.copy()
			lcsc = [lcac, lcbc]
			
			# Random initial shifts
			lcac.shifttime(float(np.random.uniform(low=-self.iniest.tderr, high=self.iniest.tderr, size=1)))
			lcbc.shifttime(float(np.random.uniform(low=-self.iniest.tderr, high=self.iniest.tderr, size=1)))
			
			# Randomize order
			pycs.gen.lc.shuffle(lcsc)
			
			# And add it
			lcsclist.append(lcsc)
		
		
		# We run the optimizer on these copies
		#optfctouts = pycs.sim.run.applyopt(optfct, lcsclist)
		
		# We run the optimizer and unshuffle
		for lcsc in lcsclist:
			optfct(lcsc)
			pycs.gen.lc.objsort(lcsc, verbose=False)
		
		# We collect the results
		self.obsres = np.array([([lcsc[1].timeshift - lcsc[0].timeshift]) for lcsc in lcsclist])
		self.log("Done! I measured %.2f +/- %.2f (intrinsic)" % (np.median(self.obsres), np.std(self.obsres)))
	
	
	def runobsshow():
		pass
				
		
		
			
			
			
			
			
			
	
	
