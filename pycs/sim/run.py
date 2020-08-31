"""
Functions to run curve shifting techniques on lightcurves produced by sim.multidraw.
We define the class runresults that "holds" results obtained by curve shifting techniques (saved into pickles).
"""

#import pycs.gen.lc
import pycs.sim.draw
import pycs.gen.util
import pycs.gen.lc
import numpy as np
import os
import time
#import pycs.sim.frk
import copy as pythoncopy

from glob import glob
import random
import os, sys


def applyopt(optfct, lcslist, **kwargs):
	"""
	Applies optfct (an optimizing function that takes a list of lightcurves as single argument)
	to all the elements (list of lightcurves)
	of lcset (a list of lists of lightcurves).
	
	Optimizes the lightcuves themselves, in place, and returns a list of the ouputs of the optimizers, corresponding to the lcslist.
	For instance, if the optfct output is a spline, it also contains the final r2s, that I will later save into the pkls !
	
	About multi cpu :
	First try using the multiprocessing module -> failed, as the optfct cannot be pickled. Perhaps need to rethink the 
	strategy.
	Second try using good old forkmap... it works !
	
	ncpu : None = I will use all CPUs, -1 = I will use all - 1 CPUs, and otherwise I will use ncpu CPUs.
	
	"""
	
	# Ouch, we cannot pickle optfct, hence multiprocessing module does not work...
	# Then I tried forkmap
	# Big bug ... border effect ? It works, but messes up the delays.
	# Could not solve this ... back to 1 cpu for now.
	
	
	
#	def retfct(lcs):
# 		"""
# 		Used only for forkmap stuff, as we cannot modifiy the lcslist without
# 		making crap, it seems.
# 		"""
# 		mylcs = [l.copy() for l in lcs] # Is this needed ? Not even sure, but better logic.
# 		optfctout = optfct(mylcs)
# 		return (optfctout, mylcs)
	
	ncpu = 1
	verbose = True
	if ncpu == 1:
		if verbose:
			print "Starting the curve shifting on a single CPU, no multiprocessing..."
		start = time.time()
		kwargs_vec = [kwargs for k in lcslist]
		# optfctouts = [optfct(lcs, **kwargs_vec[i]) for i,lcs in enumerate(lcslist)] # Ok to use optfct directly
		optfctouts = []
		sucess_dic = {'success':True, 'failed_id':[], 'error_list':[]}
		for i, lcs in enumerate(lcslist):
			try :
				optout = optfct(lcs, **kwargs_vec[i])
			except Exception as e:
				print "WARNING : I have a probleme with the curve number %i."%(i)
				sucess_dic['failed_id'].append(i)
				sucess_dic['success'] = False
				sucess_dic['error_list'].append(e)
			else :
				optfctouts.append(optout)

		print "Shifted %i simulations, using 1 CPU, time : %s" % (len(lcslist), pycs.gen.util.strtd(time.time() - start))

# 	else:
# 		ncpuava = pycs.sim.frk.nprocessors()
# 		if ncpu == None:
# 			ncpu = ncpuava
# 		if ncpu == -1:
# 			ncpu = ncpuava - 1
# 		if verbose:
# 			print "Starting the curve shifting on %i/%i CPUs." % (ncpu, ncpuava)
# 		start = time.time()
# 		# Bug :
# 		#optfctouts = pycs.sim.frk.map(optfct, lcslist, n=ncpu)
# 		
# 		retfctouts = pycs.sim.frk.map(retfct, lcslist, n=ncpu)
# 		optfctouts = [retfctout[0] for retfctout in retfctouts]
# 		optlcslist = [retfctout[1] for retfctout in retfctouts]
# 		# And now we use the fact that lcslist is mutable ...
# 		assert len(lcslist) == len(optlcslist)
# 		for lcs, optlcs in zip(lcslist, optlcslist):
# 			lcs = optlcs
# 		
# 		print "Shifted %i simulations on %i/%i CPUs, time : %s" % (len(lcslist), ncpu, ncpuava, pycs.gen.util.strtd(time.time() - start))
	if len(optfctouts) ==0 :
		print("### WARNING : it seems that your optfct does not return anything ! ###")
	# if optfctouts[0] == None:
	# 	print("### WARNING : it seems that your optfct does not return anything ! ###")
	
	return optfctouts, sucess_dic
# 		"""
# 		sys.exit("Sorry, mp does not work yet...")
# 		import multiprocessing
# 		# And using multiprocessing :
# 		ncpuava = multiprocessing.cpu_count()
# 		print "We have %i CPUs on this machine." % (ncpuava)
# 		if ncpu == None:
# 			ncpu = ncpuava
# 		print "I will use %i CPUs." % (ncpu)
# 		start = time.time()
# 		pool = multiprocessing.Pool(processes=ncpu)
# 		optlcslist = pool.map(optfct, lcslist) # Here we use this retfct
# 		lcslist = optlcslist
# 		print "Done on %i CPUs, time : %s" % (ncpu, pycs.gen.util.strtd(time.time() - start))
# 		"""

class runresults:
	"""
	Summarizes the huge list of list of lightcurves as a numpy array of timeshifts and some further info,
	to serve as input for plots, and actual time delay determinations.
	This replaces the old "boot.pkl" files ...
	The TRUE shifts are also saved (if available)
	
	All this is not related to a particular optimization technique.
	Please also provide the success_dic to remove the curves where the optimiser failed.
	
	Note the mask funcitonallity.
	"""
	
	def __init__(self, lcslist, qs=None, name="None", plotcolour = "#008800", success_dic = None):
		"""
		lcslist may or may not have "truetimeshifts". If not, I will put 0.0 as true shifts.
		
		qs should be a numpy array (as long as lcslist) that contains some chi2, r2, or d2 stuff to quantify how good the fit was.
		
		All the lcs in lcslist should be in the same order (I will check this, even if slow).
		I will not sort these lcs (as you might want them in an unsorted order).
		
		"""

		if qs is not None:
			self.qs = qs
			if qs.shape[0] != len(lcslist):
				raise RuntimeError("These qs don't have the right length !")
		else:
			# We put zeros...
			self.qs = np.zeros(len(lcslist))
				
		if len(lcslist) == 0:
			raise RuntimeError("Should this happen ?")
		
		self.tsarray = np.vstack(np.array([l.timeshift for l in lcs]) for lcs in lcslist)
		# First index selects the simulation, second index selects the timeshifts of the curves in each lcs.
		# We build a similar array for the true shifts (value = 0.0 if the curve was not drawn)
		self.truetsarray = np.vstack(np.array([getattr(l, "truetimeshift", 0.0) for l in lcs]) for lcs in lcslist)
			
		
		# We check the ordering of the lcs in lcslist
		objectstringsasset = set(["/".join([l.object for l in lcs]) for lcs in lcslist])
		if len(objectstringsasset) != 1:
			raise RuntimeError("Ouch, your lcs in lcslist are not identical/ordered !")
		
		self.labels = [l.object for l in lcslist[0]]
		
		self.name = name
		self.autoname = name
		
		self.plottrue = False # By default we plot the measured delays.
		self.plotgauss = False 
		self.plotcolour = plotcolour
		self.success_dic = success_dic
		
		self.check()

	def __len__(self):
		"""
		The number of runs
		"""
		return self.tsarray.shape[0]
	
	def nimages(self):
		"""
		The number of images (4 for a quad, 2 for a double) ...
		"""
		return self.tsarray.shape[1]
	
	def __str__(self):
		return "Runresults '%s' (%i)" % (getattr(self, "name", "untitled"), len(self))
	
	def copy(self):
		return pythoncopy.deepcopy(self)
	
	def check(self):
		
		if self.qs.shape[0] != len(self):
			raise RuntimeError("qs length error")
		#print self.truetsarray
		#print self.tsarray
		if self.tsarray.shape != self.truetsarray.shape:	
			raise RuntimeError("tsarray shape error")
	
	def applymask(self, mask):
		"""
		Removes some of the runresults according to your mask.
		"""
		
		self.tsarray = self.tsarray[mask]
		self.truetsarray = self.truetsarray[mask]
		self.qs = self.qs[mask]
		self.check()
		
	
	def gettruets(self):
		"""
		Returns some summary stats about the true delays.
		Used to find histogram ranges for plots, etc.
		"""
		#print self.truetsarray
		ret = {}
		ret["center"] = np.median(self.truetsarray, axis=0)
		ret["max"] = np.max(self.truetsarray, axis=0)
		ret["min"] = np.min(self.truetsarray, axis=0)
		spans = ret["max"] - ret["min"]
		ret["type"] = "distribution"
		if np.all(spans < 0.00001): # Then all true delays are identical
			ret["type"] = "same"
			if np.all(np.absolute(ret["center"]) < 0.00001): # Then we have no true delays (i.e., they are all 0).
				ret["type"] = "none"
		
		return ret
		
	def getts(self):
		"""
		A bit similar to gettruets, we return the median of the measured ts...
		Used for plots etc, not for calculations.
		"""
		ret = {}
		ret["center"] = np.median(self.tsarray, axis=0)
		ret["max"] = np.max(self.tsarray, axis=0)
		ret["min"] = np.min(self.tsarray, axis=0)
		ret["type"] = "distribution"
		return ret

	def get_delays_from_ts(self):
		"""
        Return the time delays, from the timeshifts. I do not account for the true timeshift.
        :return: dictionary containing the median, max, and min delays + delay labels
        """
		n = len(self.labels)
		couples = [(self.tsarray[:, i], self.tsarray[:, j]) for i in range(n) for j in range(n) if i < j]
		label_couple = [self.labels[i] + self.labels[j] for i in range(n) for j in range(n) if i < j]
		ret = {"center": [np.median(lcs2 - lcs1) for (lcs1, lcs2) in couples]}
		ret["max"] = [np.max(lcs2 - lcs1) for (lcs1, lcs2) in couples]
		ret["min"] = [np.min(lcs2 - lcs1) for (lcs1, lcs2) in couples]
		ret["delay_label"] = label_couple
		ret["type"] = "delay distribution"
		return ret
		
		
		
def joinresults(rrlist):
	"""
	Give me a list of runresults objects, I join those into a single one an return the latter.
	"""
	
	if len(rrlist) == 0:
		raise RuntimeError("Your rrlist is empty !")
	
	joined = rrlist[0].copy() # Just to get an object, with labels from the first rr.
	# Perform lots of test if it is ok to join these results ...
	for rr in rrlist:
		if rr.labels != joined.labels:
			raise RuntimeError("Don't ask me to join runresults of different objects !")
		
	#joined.name = "+".join(list(set([getattr(rr, 'simset', 'NoName') for rr in rrlist])))
	joined.name = "+".join(list(set([getattr(rr, 'name', 'NoName') for rr in rrlist])))
	joined.autoname = "%s" % (joined.name)
	joined.tsarray = np.vstack([rr.tsarray for rr in rrlist])
	joined.truetsarray = np.vstack([rr.truetsarray for rr in rrlist])
	joined.qs = np.concatenate([rr.qs for rr in rrlist])

	joined.check()
	return joined

def collect(directory = "./test", plotcolour="#008800", name=None):
	"""
	Collects the runresult objects from a directory (typically from a multirun run),
	and returns the joined runresults.
	
	"""
	#directory = "simset_%s" % (simset)
	if not os.path.isdir(directory):
		raise RuntimeError("I cannot find the directory %s" % directory)
	pklfiles = sorted(glob(os.path.join(directory, "*_runresults.pkl")))
	if len(pklfiles) == 0:
		raise RuntimeError("I couldn't find pkl files in directory %s" % directory)
	print "Reading %i runresult pickles..." % (len(pklfiles))
	rrlist = [pycs.gen.util.readpickle(pklfile, verbose=False) for pklfile in pklfiles]
	#for rr in rrlist:
	#	if not hasattr(rr, "qs"):
	#		rr.qs = None
	jrr = pycs.sim.run.joinresults(rrlist)
	jrr.plotcolour = plotcolour
	if name is not None:
		jrr.name = name
	print "OK, I have collected %i runs from %s" % (len(jrr), jrr.name)
	return jrr
	



def multirun(simset, lcs, optfct, kwargs_optim=None, optset="multirun", tsrand=10.0, analyse = True, shuffle=True, keepopt=False, trace=False, verbose=True, destpath = "./"):
	"""
	Top level wrapper to get delay "histograms" : I will apply the optfct to optimize the shifts
	between curves that you got from :py:func:`pycs.sim.draw.multidraw`, and save the results in
	form of runresult pickles.
	
	.. note: Remove my ".workingon" file and I will finish the current pkl and skip the remaining ones !
		This is useful to stop we cleanly.
	
	It is perfectly ok to launch several instances of myself on the same simset, to go faster.
	I will process every pkl of the simset only once, and prevent other instances from processing the same files.
	
	You can use me for a lot of different tasks. (note from VB : not to make coffee apparently)
	
	:param simset: The name of the simulations to run on. Those are in a directory called ``sims_name``.
	
	:param lcs: Lightcurves that define the initial shifts and microlensings you want to use.
		I will take the lightcurves from the simset, and put these shifts and ML on them.

	:param kwargs_optim: kwargs to be passed to your optfct
	
	:param optset: A new name for the optimisation.

	:param optfct: The optimizing function that takes lcs as single argument, fully optimizes the curves,
		and returns a spline, or a d2 value.
		Can be None if argument analyse is False (used for tests).
	:type optfct: function
	
	:param tsrand: I will randomly shift the simulated curves before running the optfct
		This randomizes the initial conditions.
		(uniform distrib from -tsrand to tsrand)
	
	:param shuffle: if True, I will shuffle the curves before running optc on them, and then sort them immediatly afterwards.
		
	:param keepopt: a bit similar to Trace, but simpler : we write the optimized lightcurves as well as the output of the optimizers into one pickle file per input pickle file.
		{"optfctoutlist":optfctouts, "optlcslist":simlcslist}
	
	"""
	if kwargs_optim is None :
		kwargs_optim = {}
	
	# We look for the sims directory OH GOD THIS IS SO UGLY !
	simdir = destpath + "sims_%s" % (simset)
	if not os.path.isdir(simdir):
		raise RuntimeError("Sorry, I cannot find the directory %s" % simset)
		
	simpkls = sorted(glob(os.path.join(simdir, "*.pkl")))
	if verbose:
		print "I have found %i simulation pickles in %s." % (len(simpkls), simdir)
	
	
	# We prepare the destination directory	
	destdir = destpath+"sims_%s_opt_%s" % (simset, optset)
	if verbose:
		print "I'll write my results into the directory %s." % (destdir)
	
	if not os.path.isdir(destdir):
		os.mkdir(destdir)
	else:
		if verbose:
			print "(The latter already exists.)"
	
	# The initial conditions that I will set to the sims
	if verbose:
		print "Initial conditions : "
		for l in lcs:
			print l

	success_dic = {'success': True, 'failed_id': [], 'error_list': []}
	for simpkl in simpkls:
		
		# First we test if this simpkl is already processed (or if another multirun is working on it).
		simpklfilebase = os.path.splitext(os.path.basename(simpkl))[0]
		
		workingonfilepath = os.path.join(destdir, simpklfilebase+".workingon")
		resultsfilepath = os.path.join(destdir, simpklfilebase+"_runresults.pkl")
		optfilepath = os.path.join(destdir, simpklfilebase+"_opt.pkl")
		
		if os.path.exists(workingonfilepath) or os.path.exists(resultsfilepath):
			continue
		
		# Ok, we start, hence we want to avoid other instances to work on the same pkl ...
		os.system("date > %s" % workingonfilepath)
		
		print "--- Casino running on simset %s, optset %s ---" % (simset, optset)
		simlcslist = pycs.gen.util.readpickle(simpkl)
		print "Working for %s, %i simulations." % (resultsfilepath, len(simlcslist))
		
		
		# We set the initial conditions for the curves to analyse, based on the lcs argument as reference.
		
		for simlcs in simlcslist:
			pycs.sim.draw.transfershifts(simlcs, lcs)
		
		# Now we add uniform noise to the initial time shifts 
		if tsrand != 0.0:
			for simlcs in simlcslist:
				for simlc in simlcs:
					simlc.shifttime(float(np.random.uniform(low=-tsrand, high=tsrand, size=1)))
		else:
			if verbose:
				print "I do NOT randomize initial contidions for the time shifts !"
		
		# And to the actual shifting, that will take most of the time
		if analyse:
			if shuffle:
				for simlcs in simlcslist:
					pycs.gen.lc.shuffle(simlcs)
			optfctouts, success_dic = applyopt(optfct, simlcslist, **kwargs_optim)
			if shuffle: # We sort them, as they will be passed the constructor of runresuts.
				for simlcs in simlcslist:
					pycs.gen.lc.objsort(simlcs, verbose=False)
		else:
			# Else, we just skip this step, and save the results anyway.
			if verbose:
				print "I do NOT analyse the curves !"
			optfctouts = [None]*len(simlcslist)
		
		# And now we want to save the results.
		
		# If the optfct was a spline optmization, this optfctouts is a list of splines.
		# Else it might be something different, we deal with this now.
		if hasattr(optfctouts[0], "lastr2nostab"): # then it's a spline, and we will collect these lastr2nostab values.
			tracesplinelists = [[optfctout] for optfctout in optfctouts] # just for the trace
			qs = np.array([s.lastr2nostab for s in optfctouts])
			if np.all(qs < 1.0):
				print "### WARNING : qs values are very small, did you fit that spline ? ###"
		
		else:
			try:
				qs = np.array(map(float, optfctouts)) # Then it's some kind of chi2, or a d2 : easy !
				tracesplinelists = [[]]*len(simlcslist) # just for the trace
			except Exception as e :
				print e
				tracesplinelists = [[]]*len(simlcslist) # just for the trace
				qs = None
				if analyse == True:
					print type(optfctouts[0]), optfctouts
					print "Oh no, I don't know what to do with the optfctouts !"
		
		# Trace after shifting
		if trace:
			print "Saving trace of optimized curves ..."
			tracedir = "trace_sims_%s_opt_%s" % (simset, optset)
			for (simlcs, tracesplinelist) in zip(simlcslist, tracesplinelists):
				pycs.gen.util.trace(lclist=simlcs, splist=tracesplinelist, tracedir = tracedir)

		clean_simlcslist = clean_simlist(simlcslist, success_dic)
		if keepopt:
			# A bit similar to trace, we save the optimized lcs in a pickle file.
			outopt = {"optfctoutlist":optfctouts, "optlcslist":clean_simlcslist}
			pycs.gen.util.writepickle(outopt, optfilepath)


		# Saving the results
		rr = runresults(clean_simlcslist, qs = qs, name="sims_%s_opt_%s" % (simset, optset), success_dic = success_dic)
		pycs.gen.util.writepickle(rr, resultsfilepath)

		# We remove the lock for this pkl file.
		# If the files does not exist we stop !
		if not os.path.exists(workingonfilepath):
			print "WORKINGON FILE REMOVED -> I STOP HERE"
			break
		
		else:
			os.remove(workingonfilepath)

	return success_dic
			
def clean_simlist(simlcslist, success_dic):
	for i in reversed(success_dic['failed_id']):
		print "remove simlcs ",i
		del simlcslist[i]

	return simlcslist


