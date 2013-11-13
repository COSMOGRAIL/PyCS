"""
Wrapper stuff to run PyCS on TDC data
"""

import os
import pycs

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
			
		self.plotdir = os.path.join(self.outdir, "plots")
		if not os.path.isdir(self.plotdir):
			os.mkdir(self.plotdir)
		self.datadir = os.path.join(self.outdir, "data")
		if not os.path.isdir(self.datadir):
			os.mkdir(self.datadir)
	
		
		# Gather some stats about the lcs :
		stats = lca.samplingstats(seasongap=100)
		self.sampling = stats["med"]


	def setup(self):
		
		# Set the initial shift
		self.lca.timeshift = 0.0
		self.lcb.timeshift = self.iniest.td
		pycs.tdc.util.setnicemagshift([self.lca, self.lcb])
		

	# We define a spline optimiser
	def spl(self):
		spline = pycs.spl.topopt.opt_rough([self.lca, self.lcb], nit = 5, shifttime=False, crit="r2", 
			knotstep=2.0*self.sampling, stabext=300.0, stabgap=100.0, stabstep=10.0, stabmagerr=-2.0,
			verbose=True)
		spline = pycs.spl.topopt.opt_fine([self.lca, self.lcb], nit=2, shifttime=True,
			knotstep=10.0*self.sampling, stabext=300.0, stabgap=100.0, stabstep=10.0, stabmagerr=-2.0,
			bokeps=2.0*self.sampling, boktests=5, bokwindow=None,
			splflat=True, verbose=True)
		return spline
		
		
	def fitsourcespline(self):
		# Adding some simple ML to B (required even just for the magshift !)
		pycs.gen.polyml.addtolc(self.lcb, nparams=1, autoseasonsgap = 100000.0)

		# And run a first spline
		self.sourcespline = self.spl()
		
		pycs.gen.util.writepickle(self.sourcespline, os.path.join(self.datadir, "sourcespline.pkl"))
		
		# Visu
		pycs.gen.lc.display([self.lca, self.lcb], [self.sourcespline], figsize=(18, 6), filename=os.path.join(self.plotdir, "sourcespline.png"))


	def draw(self):
		"""
		We draw the copies
		"""

	
	
