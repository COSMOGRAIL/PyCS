"""
Wrapper stuff to run PyCS on TDC data
"""

import os
import pycs
import est
import datetime
import numpy as np
import matplotlib.pyplot as plt

class Run:
	"""
	Class to contain a run of PyCS methods on a given lens
	"""
	
	def __init__(self, iniest, lca, lcb, method="", methodpar="", outdir=None):
		"""

		"""
			
		self.iniest = iniest
		
		# The output estimate of this run :
		self.outest = iniest.copy()
		self.outest.method = method
		self.outest.methodpar = methodpar
		self.outest.td = 0.0
		self.outest.tderr = -1.0
		self.outest.ms = 0.0
		self.outest.confidence = 0
		self.outest.timetaken = 0.0
	
		self.lca = lca
		self.lcb = lcb
		(self.jdrange, self.magrange) = pycs.gen.lc.displayrange([self.lca, self.lcb]) # Just to be sure that they are defined.
		
		if outdir == None:
			outdir = os.getcwd()
		self.outdir = outdir
		
		if not os.path.isdir(self.outdir):
			os.mkdir(self.outdir)
		
		self.plotdir = self.outdir
		#self.plotdir = os.path.join(self.outdir, "diagnostics")
		#if not os.path.isdir(self.plotdir):
		#	os.mkdir(self.plotdir)
		#self.datadir = os.path.join(self.outdir, "data")
		#if not os.path.isdir(self.datadir):
		#	os.mkdir(self.datadir)
			
		self.logpath = os.path.join(self.plotdir, "log.txt")
	
	
		
	def check(self):
		pass
			
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
		
	
	def varioanalysis(self):
		"""
		Computes the vario analysis, for instance to guess best knot steps
		"""
		self.log("Starting vario analysis...")
		self.lca.vario = pycs.tdc.vario.vario(self.lca)
		self.lcb.vario = pycs.tdc.vario.vario(self.lcb)
		self.log("vario results for A : %s" % str(self.lca.vario))
		self.log("vario results for B : %s" % str(self.lcb.vario))
		
		# To keep the coding simpler, we also write this into the lc objects :
		#self.lca.knotstep = self.knotstep
		#self.lcb.knotstep = self.knotstep
		self.knotstep = pycs.tdc.splopt.calcknotstep([self.lca.vario, self.lcb.vario])
		self.log("Computed knotstep: %.2f" % self.knotstep)
	
		
	
	def show(self, interactive=False):
		"""
		A simple display, with some text
		If I have a source spline, I show it !
		"""
		# Updating the ranges
		(self.jdrange, self.magrange) = pycs.gen.lc.displayrange([self.lca, self.lcb])
		
		stats = self.lca.samplingstats()
		statstxt = "%s: %i points, %i seasons, gap length : %.1f +/- %.1f days, sampling : median %.1f / mean %.1f days" % (
			self.iniest.niceid, len(self.lca), stats["nseas"], stats["meansg"], stats["stdsg"], stats["med"], stats["mean"])
		text = [(0.01, 0.93, statstxt, {"fontsize":12, "color":"black"})]
		
		if interactive:
			filename = "screen"
		else:
			filename = os.path.join(self.plotdir, "iniest.png")
		
		if hasattr(self, "sourcespline"):
			splines = [self.sourcespline]
		else:
			splines = []
		
		
		pycs.gen.lc.display([self.lca, self.lcb], splines, legendloc=1,
			filename = filename,# capsize=0, markersize=2,
			knotsize=0.005,
			showdelays=True,
			figsize=(18, 6), plotsize=(0.05, 0.98, 0.09, 0.98),
			jdrange=self.jdrange, magrange = self.magrange,
			showgrid=True, text=text)	
		
			
			
	def fitsourcespline(self, sploptfct, addmlfct=None, saveplot=True):
		# Adding some simple ML to B (required for the magshift ...)
		#pycs.gen.polyml.addtolc(self.lca, nparams=1, autoseasonsgap = 100000.0)
		#pycs.gen.polyml.addtolc(self.lcb, nparams=1, autoseasonsgap = 100000.0)
		
		# Adding the ML
		self.lca.magshift = 0.0
		self.lcb.magshift = 0.0
		if addmlfct != None:
			addmlfct([self.lca, self.lcb])

		# And run a first spline
		self.sourcespline = sploptfct([self.lca, self.lcb])
		self.log("Fitted sourcespline, delay %.2f" % self.lcb.timeshift)
		
		relshift = np.fabs(self.lcb.timeshift - self.iniest.td) / self.iniest.tderr
		self.log("Time shift: %.2f sigma" % relshift)
		
		if relshift > 5.0:
			self.log("That's too much, I stop")
			raise RuntimeError("Source spline fit failed.")
		
		#pycs.gen.util.writepickle(self.sourcespline, os.path.join(self.datadir, "sourcespline.pkl"))
		
		# We update the plot ranges
		(self.jdrange, self.magrange) = pycs.gen.lc.displayrange([self.lca, self.lcb])
		
		if saveplot:
			pycs.gen.lc.display([self.lca, self.lcb], [self.sourcespline], showdelays=True, jdrange=self.jdrange, magrange = self.magrange, figsize=(18, 6), plotsize=(0.05, 0.98, 0.09, 0.98), filename=os.path.join(self.plotdir, "sourcespline.png"), verbose=False, legendloc=1)
			pycs.gen.lc.display([self.lca, self.lcb], [self.sourcespline], showdelays=True, jdrange=(0, 1300), magrange = self.magrange, figsize=(18, 6), plotsize=(0.05, 0.98, 0.09, 0.98), filename=os.path.join(self.plotdir, "zoom_sourcespline.png"), verbose=False, legendloc=1)
		
		# And we save the residuals, to be used later when drawign simulated curves.
		pycs.sim.draw.saveresiduals([self.lca, self.lcb], self.sourcespline)

		

	def runobs(self, optfct, addmlfct=None, n=5, saveplots=True):
		"""
		Run the optimizer n times on the real data.
		We draw the copies, run on them
		"""
		
		#simdir = os.path.join(self.datadir, "copies")
		#if os.path.isdir(simdir):
		#	os.remove(simdir)
		#pycs.sim.draw.multidraw([self.lca, self.lcb], onlycopy=True, n=n, npkl=1, simset="copies", simdir=simdir)
		
		self.log("Running %i times on the observations..." % n)
		if n == 0:
			self.intrinsicratio = 0.0
			self.obsmesdelays = np.array([])
			return
		
		lcscoplist = []
		tsr = np.max([3.0, self.iniest.tderr])
		
		for i in range(n):
			
			# We make some copies to work with
			lcacop = self.lca.copy()
			lcbcop = self.lcb.copy()
			lcacop.rmml()
			lcbcop.rmml()
			lcscop = [lcacop, lcbcop]
			
			# Resetting magshifts
			lcacop.magshift = 0.0
			lcbcop.magshift = 0.0
			
			# Adding ML
			if addmlfct != None:
				addmlfct(lcscop)
			
			# Random initial shifts
			lcacop.shifttime(float(np.random.uniform(low=-tsr, high=tsr, size=1)))
			lcbcop.shifttime(float(np.random.uniform(low=-tsr, high=tsr, size=1)))
			
			# Randomize order
			pycs.gen.lc.shuffle(lcscop)
			lcscoplist.append(lcscop)
		
		# We run the optimizer and unshuffle
		for (i, lcscop) in enumerate(lcscoplist):
			out = optfct(lcscop)
			pycs.gen.lc.objsort(lcscop, verbose=False)
			
			if hasattr(out, "bokeps"): # then its a spline :)
				displaysplines = [out]
			else:
				displaysplines = []
			if saveplots:
				(jdrange, magrange) = pycs.gen.lc.displayrange(lcscop)
				pycs.gen.lc.display(lcscop, displaysplines, showdelays=True, jdrange=jdrange, magrange=magrange, figsize=(18,6), plotsize=(0.05, 0.98, 0.09, 0.98), legendloc=1, filename = os.path.join(self.plotdir, "obsopt_%i.png" % (i+1)))
				
		
		# We collect the results
		self.obsmesdelays = np.array([([lcscop[1].timeshift - lcscop[0].timeshift]) for lcscop in lcscoplist])
		self.log("Done! I measured %.2f +/- %.2f (intrinsic)" % (np.median(self.obsmesdelays), np.std(self.obsmesdelays)))
		
		self.outest.td = np.median(self.obsmesdelays)
		self.intrinsicratio = np.std(self.obsmesdelays) / self.iniest.tderr
		self.log("Intrinsic/initial error ratio: %.2f" % self.intrinsicratio)
		
		# We also estimate the magnitude shift, assuming that the method optimized it.
		self.obsmesms = np.array([np.median(lcscop[1].getmags() - lcscop[1].mags) - np.median(lcscop[0].getmags() - lcscop[0].mags) for lcscop in lcscoplist])
		self.outest.ms = np.median(self.obsmesms)
		
		
	
	def runobsplot(self):
		"""
		A quick histogram to see that intrinsic variance compared to the initial estimate
		"""
		
		tdmin = self.iniest.td - 3.0*self.iniest.tderr
		tdmax = self.iniest.td + 3.0*self.iniest.tderr
		
		fig = plt.figure(figsize=(6, 3))
		fig.subplots_adjust(top=0.95, bottom=0.2)
		if len(self.obsmesdelays) != 0:
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
		plt.close()
		
	
	
	def runsim(self, optfct, addmlfct=None, n=5, saveplots=True):
		"""
		Drawing the simulated curves and runing on them
		"""
		
		self.log("Drawing %i simulations..." % n)
		if n == 0:
			self.totalratio = 0.0
			self.simtruedelays = np.array([])
			self.simmesdelays = np.array([])
			return

		timeshiftoriga = self.lca.timeshift # So this is what the single spline fit gave us.
		timeshiftorigb = self.lcb.timeshift
		truetsr = np.max([3.0, self.iniest.tderr])
		tsr = np.max([3.0, self.iniest.tderr])
		
		lcssimlist = []
		for i in range(n):
			
			# Set some random "true delays"
			self.lca.timeshift = timeshiftoriga + float(np.random.uniform(low = -truetsr, high = truetsr, size=1))
			self.lcb.timeshift = timeshiftorigb + float(np.random.uniform(low = -truetsr, high = truetsr, size=1))
			
			# Draw them
			lcssim = pycs.sim.draw.draw([self.lca, self.lcb], self.sourcespline, shotnoise="sigma")
	
			if saveplots:
				pycs.gen.lc.display(lcssim, [self.sourcespline], showdelays=True, jdrange=self.jdrange, magrange = self.magrange, figsize=(18,6), plotsize=(0.05, 0.98, 0.09, 0.98), legendloc=1, filename = os.path.join(self.plotdir, "sim_%i.png" % (i+1)))
				pycs.gen.lc.display(lcssim, [self.sourcespline], showdelays=True, jdrange=(0, 1300), magrange = self.magrange, figsize=(18,6), plotsize=(0.05, 0.98, 0.09, 0.98), legendloc=1, filename = os.path.join(self.plotdir, "zoom_sim_%i.png" % (i+1)))
			
			# Remove ML and add a new one :
			lcssim[0].magshift = 0.0
			lcssim[1].magshift = 0.0
			lcssim[0].rmml()
			lcssim[1].rmml()
			if addmlfct != None:
				addmlfct(lcssim)
			
			# Set some wrong "initial delays" for the analysys, around these "true delays".
			lcssim[0].shifttime(float(np.random.uniform(low=-tsr, high=tsr, size=1)))
			lcssim[1].shifttime(float(np.random.uniform(low=-tsr, high=tsr, size=1)))
			
			# Propagate the "pirate" attributes:
			lcssim[0].vario = self.lca.vario
			lcssim[1].vario = self.lcb.vario
			
			# Randomize order
			pycs.gen.lc.shuffle(lcssim)
			lcssimlist.append(lcssim)
			
			
		# We run the optimizer and unshuffle
		self.log("Running on %i simulations..." % n)
		for (i, lcssim) in enumerate(lcssimlist):
			out = optfct(lcssim)
			pycs.gen.lc.objsort(lcssim, verbose=False)
			
			if hasattr(out, "bokeps"): # then its a spline :)
				displaysplines = [out]
			else:
				displaysplines = []
			if saveplots:
				(jdrange, magrange) = pycs.gen.lc.displayrange(lcssim)
				pycs.gen.lc.display(lcssim, displaysplines, jdrange=jdrange, magrange=magrange, figsize=(18,6), plotsize=(0.05, 0.98, 0.09, 0.98), legendloc=1, filename = os.path.join(self.plotdir, "simopt_%i.png" % (i+1)))

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
		
		fig = plt.figure(figsize=(6, 3))
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
		plt.close()
	
	def getoutest(self):
		"""
		Called at the end...
		"""
		return self.outest
	
	
	def save(self):
		"""
		I write myself into a pickle
		"""
		pycs.gen.util.writepickle(self, os.path.join(self.outdir, "run.pkl"))
	
	
def multirun(iniests, 
	sploptfct, optfct, addmlfct,
	nobs = 5, nsim = 5,
	ncpu=0,
	diagnostics = True,
	saveplots = False,
	tdcpath = "./tdc0", outdir="./multirun", method="", methodpar=""):

	"""
	Running a technique on several TDC curves.
	
	:param iniests: a list of (unique) initial estimates. I will run on those quasars.
	
	:param sploptfct: a spline optimizer function to draw simulations
	:param optfct: the optimizer function that you want to use
	:param addmlfct: a function that adds some ML to the curves passed as argument.
	
	
	:param nobs: how often should I run on the observations (intrinsic variance)
	:param nsim: how many Monte Carlo runs should I do ?
	
	:param ncpu: how many cpus should I use (0 = automatic)
	
	:param diagnostics: do you want me to output some summary diagnostic plots / data ?
	:param saveplots: do you want me to output even more detaild plots (for debugging)
	
	:param tdcpath: path to the data directory.
	
	:param outdir: where I'll write some output
	
	:param method: free string with a name for this run, e.g. the optimizer that you selected
	:param methodpar: free string with some method parameters
	
	"""
	
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
		
	outcsv = os.path.join(outdir, "estimates.csv")
	crashcsv = os.path.join(outdir, "crash.csv")
	
	for iniest in iniests:
		
		try:
			lcspath = pycs.tdc.util.tdcfilepath(set=iniest.set, rung=iniest.rung, pair=iniest.pair)
			(lca, lcb) = pycs.tdc.util.read(lcspath, shortlabel=False)
		
			starttime = datetime.datetime.now()
			
			r = Run(iniest, lca, lcb, method=method, methodpar=methodpar, outdir=os.path.join(outdir, iniest.id))
			r.setup()
			
			# Computing some vario stats
			r.varioanalysis()
			
			# Display of the input light curves
			if diagnostics:
				r.show()
			
			# First fit
			r.fitsourcespline(sploptfct = sploptfct, addmlfct = addmlfct, saveplot = diagnostics)
			
			# Getting the delay
			r.runobs(optfct = optfct, addmlfct = addmlfct, n = nobs, saveplots=saveplots)
			if diagnostics:
				r.runobsplot()
	
			# Getting the error bar
			r.runsim(optfct = optfct, addmlfct = addmlfct, n = nsim, saveplots=saveplots)
			if diagnostics:
				r.runsimplot()
			
			endtime = datetime.datetime.now()
			r.outest.timetaken = (endtime - starttime).total_seconds()
			r.log("This took me %.0f seconds" % (r.outest.timetaken))

			# Writing the output
			est.writecsv([r.outest], outcsv, append=True)
			
			# And saving myself, so that you can inspect me afterwards !
			r.save()
			
			
		except RuntimeError as error:
			r.log("Damn, a RuntimeError !") # Be polite !
			r.log(str(error))
			est.writecsv([iniest], crashcsv, append=True)
			print error
			raise
			
		
	

	
	
	
	
	
	
