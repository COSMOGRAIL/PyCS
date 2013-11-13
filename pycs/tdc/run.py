"""
Wrapper stuff to run PyCS on TDC data
"""

import os
import pycs
import est
import copy
import datetime
import numpy as np
import matplotlib.pyplot as plt

class Run:
	"""
	Class to contain a run of PyCS methods on a given lens
	"""
	
	def __init__(self, iniest, lca, lcb, method="", methodpar="", plots=True, outdir=None):
		"""

		"""
			
		self.iniest = iniest
		
		# The output estimate of this run :
		self.outest = copy.deepcopy(iniest)
		self.outest.method = method
		self.outest.methodpar = methodpar
		self.outest.td = 0.0
		self.outest.tderr = 0.0
		self.outest.ms = 0.0
		self.outest.confidence = 0
		self.outest.timetaken = 0.0
	
		self.lca = lca
		self.lcb = lcb
		(self.jdrange, self.magrange) = pycs.gen.lc.displayrange([self.lca, self.lcb]) # Just to be sure that they are defined.
		
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
		# Adding some simple ML to B (required for the magshift ...)
		#pycs.gen.polyml.addtolc(self.lca, nparams=1, autoseasonsgap = 100000.0)
		pycs.gen.polyml.addtolc(self.lcb, nparams=1, autoseasonsgap = 100000.0)

		# And run a first spline
		self.sourcespline = sploptfct([self.lca, self.lcb])
		self.log("Fitted sourcespline, delay %.2f" % self.lcb.timeshift)
		
		relshift = np.fabs(self.lcb.timeshift - self.iniest.td) / self.iniest.tderr
		self.log("Time shift: %.2f sigma" % relshift)
		
		if relshift > 5.0:
			self.log("That's too much, I stop")
			self.gogogo = False
		
		#pycs.gen.util.writepickle(self.sourcespline, os.path.join(self.datadir, "sourcespline.pkl"))
		
		# Visu
		
		(self.jdrange, self.magrange) = pycs.gen.lc.displayrange([self.lca, self.lcb])
		
		pycs.gen.lc.display([self.lca, self.lcb], [self.sourcespline], jdrange=self.jdrange, magrange = self.magrange, figsize=(18, 6), filename=os.path.join(self.plotdir, "sourcespline.png"), verbose=False)
		
		# And we save the residuals, to be used later when drawign simulated curves.
		pycs.sim.draw.saveresiduals([self.lca, self.lcb], self.sourcespline)

		

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
		lcscoplist = []
		for i in range(n):
			
			# We make some copies to work with
			lcacop = self.lca.copy()
			lcbcop = self.lcb.copy()
			lcscop = [lcacop, lcbcop]
			
			# Random initial shifts
			lcacop.shifttime(float(np.random.uniform(low=-self.iniest.tderr, high=self.iniest.tderr, size=1)))
			lcbcop.shifttime(float(np.random.uniform(low=-self.iniest.tderr, high=self.iniest.tderr, size=1)))
			
			# Randomize order
			pycs.gen.lc.shuffle(lcscop)
			lcscoplist.append(lcscop)
		
		# We run the optimizer and unshuffle
		for lcscop in lcscoplist:
			optfct(lcscop)
			pycs.gen.lc.objsort(lcscop, verbose=False)
		
		# We collect the results
		self.obsmesdelays = np.array([([lcscop[1].timeshift - lcscop[0].timeshift]) for lcscop in lcscoplist])
		self.log("Done! I measured %.2f +/- %.2f (intrinsic)" % (np.median(self.obsmesdelays), np.std(self.obsmesdelays)))
		
		self.outest.td = np.median(self.obsmesdelays)
		self.intrinsicratio = np.std(self.obsmesdelays) / self.iniest.tderr
		self.log("Intrinsic/initial error ratio: %.2f" % self.intrinsicratio)
		
	
	def runobsplot(self):
		"""
		A quick histogram to see that intrinsic variance compared to the initial estimate
		"""
		tdmin = self.iniest.td - 3.0*self.iniest.tderr
		tdmax = self.iniest.td + 3.0*self.iniest.tderr
		
		fig = plt.figure(figsize=(10, 3))
		fig.subplots_adjust(top=0.95, bottom=0.2)
		plt.hist(self.obsmesdelays, range=(tdmin, tdmax), bins=200, color="green", lw=0)
		plt.xlim(tdmin, tdmax)
		plt.xlabel("Delay [day]")
		plt.ylabel("Counts")
		#ax = plt.gca()
		plt.figtext(0.15, 0.8, "Intrinsic/initial error ratio: %.2f" % self.intrinsicratio)
		plt.axvline(self.iniest.td - self.iniest.tderr, color="gray", linestyle="-", zorder=20)
		plt.axvline(self.iniest.td + self.iniest.tderr, color="gray", linestyle="-", zorder=20)
		plt.axvline(self.outest.td, color="red", linestyle="-", zorder=20)
		plt.savefig(os.path.join(self.plotdir, "intrinsic_variance.png"))
		
	
	
	def runsim(self, optfct, n=5, plots=True):
		"""
		Drawing the simulated curves and runing on them
		"""
		self.check()
		
		self.log("Drawing %i simulations..." % n)
		
		timeshiftoriga = self.lca.timeshift # So this is what the single spline fit gave us.
		timeshiftorigb = self.lcb.timeshift
		
		lcssimlist = []
		for i in range(n):
			
			# Set some random "true delays"
			truetsr = self.iniest.tderr
			self.lca.timeshift = timeshiftoriga + float(np.random.uniform(low = -truetsr, high = truetsr, size=1))
			self.lcb.timeshift = timeshiftorigb + float(np.random.uniform(low = -truetsr, high = truetsr, size=1))
			
			# Draw them
			lcssim = pycs.sim.draw.draw([self.lca, self.lcb], self.sourcespline, shotnoise="sigma")
	
			if plots:
				pycs.gen.lc.display(lcssim, [self.sourcespline], jdrange=self.jdrange, magrange = self.magrange, figsize=(18,6), filename = os.path.join(self.plotdir, "sim_%i.png" % (i+1)))
			
			# Set some wrong "initial delays" for the analysys, around these "true delays".
			lcssim[0].shifttime(float(np.random.uniform(low=-self.iniest.tderr, high=self.iniest.tderr, size=1)))
			lcssim[1].shifttime(float(np.random.uniform(low=-self.iniest.tderr, high=self.iniest.tderr, size=1)))
			
			# Randomize order
			pycs.gen.lc.shuffle(lcssim)
			lcssimlist.append(lcssim)
		
		# We run the optimizer and unshuffle
		self.log("Running on %i simulations..." % n)
		for lcssim in lcssimlist:
			optfct(lcssim)
			pycs.gen.lc.objsort(lcssim, verbose=False)
		
		# We collect the results
		self.simtruedelays = np.array([([lcssim[1].truetimeshift - lcssim[0].truetimeshift]) for lcssim in lcssimlist])
		self.simmesdelays = np.array([([lcssim[1].timeshift - lcssim[0].timeshift]) for lcssim in lcssimlist])
		
		syserr = np.fabs(np.mean(self.simmesdelays - self.simtruedelays))
		ranerr = np.std(self.simmesdelays - self.simtruedelays)
		toterr = np.hypot(syserr, ranerr)
		
		self.log("Done! Systematic %.2f day, random %.2f day, total error %.2f day." % (syserr, ranerr, toterr))
		
		self.outest.tderr = toterr
		self.totalratio = toterr / self.iniest.tderr
		self.log("Total/initial error ratio: %.2f" % self.totalratio)
			

	def runsimplot(self):
		"""
		Viz of the MC results
		"""
		
		tdmin = self.iniest.td - 3.0*self.iniest.tderr
		tdmax = self.iniest.td + 3.0*self.iniest.tderr
		
		fig = plt.figure(figsize=(10, 3))
		fig.subplots_adjust(top=0.95, bottom=0.2)
		plt.scatter(self.simtruedelays, self.simmesdelays - self.simtruedelays)
		plt.xlim(tdmin, tdmax)
		plt.xlabel("True delay [day]")
		plt.ylabel("Measurement error [day]")
		plt.axvline(self.iniest.td - self.iniest.tderr, color="gray", linestyle="-", zorder=20)
		plt.axvline(self.iniest.td + self.iniest.tderr, color="gray", linestyle="-", zorder=20)
		plt.axhline(0, color="black", linestyle="-", zorder=20)
		plt.axhline(- self.iniest.tderr, color="gray", linestyle="-", zorder=20)
		plt.axhline(+ self.iniest.tderr, color="gray", linestyle="-", zorder=20)
		plt.axhspan(-self.outest.tderr, self.outest.tderr, xmin=0, xmax=1, lw=0, color="red", alpha=0.2)
		#plt.ylim(- 2.0*self.iniest.tderr, + 2.0*self.iniest.tderr)
		
		
		plt.figtext(0.15, 0.8, "Total/initial error ratio: %.2f" % self.totalratio)
		
		#ax = plt.gca()
		#plt.figtext(0.15, 0.8, "Intrinsic/initial error ratio: %.2f" % self.intrinsicratio)
		#plt.axvline(self.iniest.td - self.iniest.tderr, color="gray", linestyle="-", zorder=20)
		#plt.axvline(self.iniest.td + self.iniest.tderr, color="gray", linestyle="-", zorder=20)
		#plt.axvline(self.outest.td, color="red", linestyle="-", zorder=20)
		
		plt.savefig(os.path.join(self.plotdir, "measvstrue.png"))
			
			
	
	
