"""
Wrapper stuff to run PyCS on TDC data

We change the philosophy here: we make a clear separation between functions that draw lcs,
and functions that analyse them.

The goal is to create "once and for all" a huge amount of simulated curves, and then
run the optimizer on a limited number of simulated curves, randomly chosen.

The copied and simulated lcs are stored each into an individual pkl (one per lcs)
"""

import os, sys
import pycs
import est
import datetime
import numpy as np
from copy import copy
import glob


def createdir(estimate, path):
	
	dirpath = os.path.join(os.getcwd(),path)

	# 1) Create main directory
	
	if not os.path.isdir(dirpath):
		os.mkdir(dirpath)
		
	# 2) Create sub-directories

	subdirpath = os.path.join(dirpath,estimate.id)
	if not os.path.isdir(subdirpath):
		os.mkdir(subdirpath)	
	

def drawcopy(estimate, path, n=1, maxrandomshift = None, datadir=''):
	"""
	Draw n times one single copy of lcs (the ones your estimate is about) into the path directory.
	NEW: set can be the name of a COSMOGRAIL lens.
	
	:estimates	- list of estimate objects, with id, td and tderr
	
	:path		- where your copy will be written
	
	:addmlfct	- If you want to add different microlensing to your obslcs
	
	:n		- Number of time you will run drawobs (each run produce new copycurves)
	
	:maxrandomshift	- Maximum range for the random time shift added to the copycurves
	"""

	copydir = os.path.join(os.getcwd(),path)

	'''
	if os.path.isfile(copydir)==False:
		try:
			os.system('mkdir %s' % copydir)
		except:
			print "going on..."	
	'''

	if estimate.set in ['tdc0', 'tdc1']:
		lcspath = os.path.join(datadir,pycs.tdc.util.tdcfilepath(set=estimate.set, rung=estimate.rung, pair=estimate.pair))
		(lca, lcb) = pycs.tdc.util.read(lcspath, shortlabel=False)

	else:
		lcspath = os.path.join(datadir,pycs.tdc.util.tdcfilepath(set=estimate.set, rung=estimate.rung, pair=estimate.pair))
		(lca, lcb) = pycs.tdc.util.read(lcspath, pair = estimate.pair)


	# compute the vario analysis of these lcs
	if not (hasattr(lca, "vario") and hasattr(lcb, "vario")):
		lca.vario = pycs.tdc.vario.vario(lca, verbose=True)
		lcb.vario = pycs.tdc.vario.vario(lcb, verbose=True)
		print '---Vario Analysis Done---'


	ind=0

	while ind < n:


		# Reset ML and shifts (get them "as read" from scratch...)
		lca.magshift = 0.0
		lcb.magshift = 0.0	
		lca.timeshift = 0.0
		lcb.timeshift = 0.0

		lcs = [lca,lcb] 

		# Time shift around the initial value		
		tsr = np.max([3.0, estimate.tderr]) ## NEED TO GIVE A MAXIMUM VALUE ON THIS (10?), otherwise optimizer fails

		if not maxrandomshift == None:
			if tsr >= maxrandomshift:
				tsr = maxrandomshift


		lcb.shifttime(estimate.td +float(np.random.uniform(low=-tsr, high=tsr, size=1)))


		# Now, we write that copycurve

		if os.path.exists(os.path.join(copydir,estimate.id)) == False:
			try:
				os.system('mkdir %s' % os.path.join(copydir,estimate.id))
				print ' \n ----- New copy/sim directory created: %s ' % os.path.join(copydir,estimate.id)

			except:
				print "going on..."

		find_max_ind = 'ls %s | wc -l' % os.path.join(copydir,estimate.id,'c*')
		next_copy_ind = 'c'+str(int(os.popen(find_max_ind).readline())).rjust(3,'0')+'.pkl'
		copypath = os.path.join(copydir,estimate.id,next_copy_ind)
		pycs.gen.util.writepickle(lcs,copypath)	
		ind += 1
	

	
def drawsim(estimate, path, sploptfct, n=1, maxrandomshift = None, datadir='') :
	"""
	Draw n times one single sim curves of lcs (the ones your esimate is about) into the path directory.
	
	:estimates	- list of estimate objects, with id, td and tderr
	
	:path		- where your copy will be written
	
	:addmlfct	- If you want to add different microlensing to your obslcs
	
	:n		- Number of time you will run drawsim (each run produce new simcurves)	

	:maxrandomshift	- Maximum range for the random "true" and "guessed" time shift added to the simcurves
	"""

	simdir = os.path.join(os.getcwd(),path)
	
	'''
	if os.path.isfile(simdir)==False:
		try:
			os.system('mkdir %s' % simdir)
		except:
			print "going on..."		
	'''	
				
	
	# get the model lightcurves
	if estimate.set in ['tdc0', 'tdc1']:
		lcspath = os.path.join(datadir,pycs.tdc.util.tdcfilepath(set=estimate.set, rung=estimate.rung, pair=estimate.pair))
		(lca, lcb) = pycs.tdc.util.read(lcspath, shortlabel=False)

	else:
		lcspath = os.path.join(datadir,pycs.tdc.util.tdcfilepath(set=estimate.set, rung=estimate.rung, pair=estimate.pair))
		(lca, lcb) = pycs.tdc.util.read(lcspath, pair = estimate.pair)
	#pycs.gen.lc.display([lca,lcb])
	#sys.exit()
	ind=0

	# shift lcs as estimated by d3cs
	lcb.shifttime(estimate.td)
	#pycs.gen.lc.display([lca,lcb])
	#sys.exit()

	# fit a spline on the data and save residuals
	sourcespline = sploptfct([lca, lcb])
	pycs.sim.draw.saveresiduals([lca, lcb], sourcespline)


	# define 'initial timeshifts', given by fitsourcespline
	timeshifta = lca.timeshift
	timeshiftb = lcb.timeshift		
	#pycs.gen.lc.display([lca,lcb],[sourcespline])
	#sys.exit()

	#TODO : propagate knotstep attribute from varioanalysis to every copy, to avoid running varioanalysis for every simcurve...

	# now, draw n copycurves
	while ind < n:

		# set  a random "true" delay 
		truetsr = np.max([3.0, estimate.tderr])
		if not maxrandomshift == None:
			if truetsr >= maxrandomshift:
				truetsr = maxrandomshift

		lca.timeshift = timeshifta
		lcb.timeshift = timeshiftb + (float(np.random.uniform(low = -truetsr, high = truetsr, size=1))) 


		'''
		# Reset ML and shifts
		lca.rmml()
		lcb.rmml()	
		lca.magshift = 0.0
		lcb.magshift = 0.0			
		'''


		# Use the magic draw function to create simulated lcs

		lcssim = pycs.sim.draw.draw([lca, lcb], sourcespline, shotnoise="sigma")

		# add the vario attribues from original lcs (instead of recomputing it everytime...)
		lcssim[0].vario = lca.vario
		lcssim[1].vario = lcb.vario

		#pycs.gen.lc.display(lcssim)			
		#sys.exit()

		# Remove magshift, remove ML and add a new one :
		lcssim[0].magshift = 0.0
		lcssim[1].magshift = 0.0
		lcssim[0].rmml()
		lcssim[1].rmml()
		
		tsr = np.max([3.0, estimate.tderr])
		if not maxrandomshift == None:
			if tsr>=maxrandomshift:
				tsr=maxrandomshift

		# Set some wrong "initial delays" for the analysis, around the "true delays".
		lcssim[1].shifttime(float(np.random.uniform(low=-tsr, high=tsr, size=1)))

		#print 'TRUE DELAY  : ',lcssim[1].truetimeshift-lcssim[0].truetimeshift
		#print 'INITIAL DELAY: ',lcssim[1].timeshift-lcssim[0].timeshift						
		#pycs.gen.lc.display(lcssim)
		#sys.exit()

		# And write that simcurve

		if os.path.exists(os.path.join(simdir,estimate.id)) == False:
			try:
				os.system('mkdir %s' % os.path.join(simdir,estimate.id))
				print ' \n ----- New copy/sim directory created: %s ' % os.path.join(simdir,estimate.id)
			except:
				print "going on..."	

		find_max_ind = 'ls %s | wc -l' % os.path.join(simdir,estimate.id,'s*')
		next_sim_ind = 's'+str(int(os.popen(find_max_ind).readline())).rjust(3,'0')+'.pkl'
		simpath = os.path.join(simdir,estimate.id,next_sim_ind)

		pycs.gen.util.writepickle(lcssim,simpath)

		ind += 1	
			

def runcopy(estimate, path, optfct, n=1, clist=None):
	"""
	Run the optimizer (optfct) on n different copycurves
	Return the optimised timeshift and magnitude for each copycurves
	
	The n copycurves are chosen randomly in the copydir, unless you give a clist of numbers
	i.e. if clist = [1,3] runobs will run on c001.pkl and c003.pkl
	"""
	
	index_list = copy(clist)

	copydir = os.path.join(os.getcwd(),path)	

	# Check that there is enough copycurves in the copydir:

	find_max_ind = 'ls %s | wc -l' % os.path.join(copydir,estimate.id,'c*')
	number_of_files = int(os.popen(find_max_ind).readline())

	if number_of_files < n:
		print 'Not enough copyfiles in %s ' % os.path.join(copydir,estimate.id)		
		raise RuntimeError('not enough copycurves !!')
		sys.exit()


	if clist == None:

		mylist = np.arange(number_of_files)
		np.random.shuffle(mylist)
		index_list = mylist[:n] 	

	copytds  = []
	copymags = []
	optlcs = []
	splines = []


	for ind in index_list:
		try:
			copyfile = 'c%s' % str(ind).rjust(3,'0')+'.pkl'
			copypath = os.path.join(copydir,estimate.id,copyfile)

			lcs = pycs.gen.util.readpickle(copypath)

			#pycs.gen.lc.display(lcs)
			#sys.exit()
			spline = optfct(lcs)

			copytds.append(lcs[1].timeshift-lcs[0].timeshift)
			copymags.append(np.median(lcs[1].getmags() - lcs[1].mags) - np.median(lcs[0].getmags() - lcs[0].mags))

			optlcs.append(lcs)
			splines.append(spline)
		except:
			print 'Running on c%s.pkl' % str(ind).rjust(3,'0')+' failed'


	return (copytds, copymags, optlcs, splines)
		
		
def runsim(estimate, path, optfct, n=1, slist=None):		
	"""
	Run the optimizer (optfct) on n different simcurves
	Return the optimised timeshift and magnitude for each copycurves
	Also return the true delays of each simcurve
	
	The n simcurves are chosen randomly in the simdir, unless you give a slist of numbers
	i.e. if slist = [1,3] runobs will run on s001.pkl and s003.pkl
	"""		

	index_list = copy(slist)

	simdir = os.path.join(os.getcwd(),path)		

	# Check that there is enough simcurves in the copydir:

	find_max_ind = 'ls %s | wc -l' % os.path.join(simdir,estimate.id,'s*')
	number_of_files = int(os.popen(find_max_ind).readline())

	if number_of_files < n:
		print 'Not enough simfiles in %s ' % os.path.join(simdir,estimate.id)
		raise RuntimeError('not enough simcurves !!')
		sys.exit()


	if slist == None:

		mylist = np.arange(number_of_files)
		np.random.shuffle(mylist)
		index_list = mylist[:n] 	

	simtds  = []
	simttds = []
	simmags = []
	optlcs = []
	inlcs = []
	splines = []

	# run the optimizer on each simcurve

	for ind in index_list:

		try:
			simfile = 's%s' % str(ind).rjust(3,'0')+'.pkl'
			simpath = os.path.join(simdir,estimate.id,simfile)

			lcs = pycs.gen.util.readpickle(simpath)
			inlcs.append([lcs[0].copy(), lcs[1].copy()])

			out = optfct(lcs)
			print 'TRUE DELAY: ',lcs[1].truetimeshift-lcs[0].truetimeshift
			print 'MEASURED DELAY: ',lcs[1].timeshift-lcs[0].timeshift
			#pycs.gen.lc.display(lcs)
			#sys.exit()


			simttds.append(lcs[1].truetimeshift-lcs[0].truetimeshift)
			simtds.append(lcs[1].timeshift-lcs[0].timeshift)
			simmags.append(np.median(lcs[1].getmags() - lcs[1].mags) - np.median(lcs[0].getmags() - lcs[0].mags))
			#print "simttds",lcs[1].truetimeshift-lcs[0].truetimeshift
			#print "simtds",lcs[1].timeshift-lcs[0].timeshift
			#sys.exit()

			optlcs.append(lcs)
			splines.append(out)
		except:
			print 'Running on s%s.pkl' % str(ind).rjust(3,'0')+' failed'

		
	return (simtds, simmags, simttds, inlcs, optlcs, splines)
		

def multirun(estimate, path, optfct, ncopy, nsim, clist=None, slist=None):

	"""
	Wrapper around runsim and runobs
	Run the optimizer ncopy and nsim times on the copy and sim
	Return the results in a list of estimates objects
	"""

	copyout = runcopy(estimate, path, optfct = optfct, n=ncopy, clist=clist)
	pycs.gen.util.writepickle(copyout, os.path.join(path, "copyout_%s.pkl" % estimate.id))
	
	simout = runsim(estimate, path, optfct = optfct, n=nsim, slist=slist)
	pycs.gen.util.writepickle(simout, os.path.join(path, "simout_%s.pkl" % estimate.id))
	


def viz(estimate, path, datadir):
	"""
	Look at some light curves. change in place...
	"""
	
	copytds,copymags,copyoptlcs,copyoptsplines = pycs.gen.util.readpickle(os.path.join(path, "copyout_%s.pkl" % estimate.id))
	simtds,simmags,simttds,siminlcs,simoptlcs,simoptsplines = pycs.gen.util.readpickle(os.path.join(path, "simout_%s.pkl" % estimate.id))
	
	lcspath = os.path.join(datadir,pycs.tdc.util.tdcfilepath(set=estimate.set, rung=estimate.rung, pair=estimate.pair))

	if estimate.set in ['tdc0', 'tdc1']:
		origlcs = pycs.tdc.util.read(lcspath, shortlabel=False)
	else:
		origlcs = pycs.tdc.util.read(lcspath, pair=estimate.pair)

	for l in origlcs:  
		l.plotcolour = "black"
		
	# See the fit on the copies, and compare copies with originial curve:	
	"""
	for (spline, lcs) in zip(copyoptsplines, copyoptlcs):
		pycs.gen.lc.display(origlcs + lcs, [spline])
	"""
	for (spline, lcs) in zip(copyoptsplines, copyoptlcs):
		pycs.gen.lc.display(lcs, [spline])
	
	
	
	# Compare the sims to the data

	for lcs in siminlcs:
		pycs.gen.lc.display(lcs + origlcs, figsize=(22, 10))

	
	# See the fit on the optimized sims:

	for (spline, lcs) in zip(simoptsplines, simoptlcs):
		pycs.gen.lc.display(lcs, [spline], figsize=(22, 10))

	
	
	
	
def summarize(estimate, path, makefig=False, skipdone=True):
	"""
	I will silently skip non-yet-ready pairs.
	If skipdone, I will silently skip those pairs that have already been summarized.
	
	"""

	resultpklpath = os.path.join(path,'%s.pkl' % estimate.id)
	if skipdone and os.path.exists(resultpklpath):
		return
	
	copyoutpklpath = os.path.join(path, "copyout_%s.pkl" % estimate.id)
	simoutpklpath = os.path.join(path, "simout_%s.pkl" % estimate.id)
	if not os.path.exists(copyoutpklpath):
		return
	if not os.path.exists(simoutpklpath):
		return

	# We start by creating Estimate objects that will contain our final result:

	method = 'PyCS'
	methodpar = 'PyCS'
	
	outest = pycs.tdc.est.Estimate(set=estimate.set, rung=estimate.rung, pair=estimate.pair,
		 				method = method, methodpar = methodpar)

	
	copytds,copymags,copyoptlcs,copyoptsplines = pycs.gen.util.readpickle(copyoutpklpath)
	simtds,simmags,simttds,siminlcs,simoptlcs,simoptsplines = pycs.gen.util.readpickle(simoutpklpath)
	
	
	# And we put the results in the output estimates, and save each individual output estimate
	
	print '='*30
	print outest.niceid
	outest.td = np.median(copytds)
	outest.tdvarint = np.std(copytds)
	print 'time delay:  ',outest.td
	syserr = np.fabs(np.mean([simtds[i]-simttds[i] for i in np.arange(len(simtds))]))
	ranerr = np.std([simtds[i]-simttds[i] for i in np.arange(len(simtds))])
	outest.tderr  = np.hypot(syserr, ranerr)
	outest.tdranerr = ranerr
	outest.tdsyserr = syserr
	print '   systematic error:  ',syserr
	print '       random error:  ',ranerr
	print '        total error:  ',outest.tderr		

	
	pycs.gen.util.writepickle(outest,resultpklpath)
	
	# We are done.
	if not makefig:	
		return
		
	# If we can, we make a fig
	try:
		import matplotlib.pyplot as plt
	except:
		print "can't import matplotlib"
		return

	mind3cs = estimate.td - estimate.tderr
	maxd3cs = estimate.td + estimate.tderr
	minrange = estimate.td - 2.0*estimate.tderr
	maxrange = estimate.td + 2.0*estimate.tderr
	
	plt.figure(figsize=(14, 5))
	
	plt.subplot(121)
	plt.hist(np.array(copytds), bins=30) # Do not specify range, so that we see outliers !
	plt.axvline(x=outest.td, linewidth=4, color='r')
	plt.axvline(x=mind3cs, linewidth=4, color='g')
	plt.axvline(x=maxd3cs, linewidth=4, color='g')
	plt.xlabel("copytds")
	
	plt.title(estimate.id)
	
	plt.subplot(122)
	plt.scatter(np.array(simttds), np.array(simtds)-np.array(simttds))
	plt.xlabel("simttds")
	plt.ylabel("simtds - simttds")
	#plt.xlim(minrange, maxrange) # No, do not use xlim here -- we want to see all points, if something goes wrong.
	plt.axhline(y=0, color='black')
	plt.axhline(y=outest.tderr, linewidth=4, color='red')  # The output error
	plt.axhline(y=-outest.tderr, linewidth=4, color='red')
	plt.axhline(y=estimate.tderr, linewidth=4, color='g')  # The input error
	plt.axhline(y=-estimate.tderr, linewidth=4, color='g')
	plt.axvline(x=mind3cs, linewidth=4, color='g')  # The input interval
	plt.axvline(x=maxd3cs, linewidth=4, color='g')
	plt.title("Uncertainty [d]: %.1f (%.1f ran, %.1f sys), varint: %.1f (%.0f %%)" % (outest.tderr, outest.tdranerr, outest.tdsyserr, outest.tdvarint, 100.0*outest.tdvarint/outest.tderr)  )
	
	plt.savefig(os.path.join(path, "copytds_%s.png" % estimate.id))



def collect(estimates, path):
	"""
	I gather the estimates.
	"""
	
	print "You ask me to search for %i estimates" % (len(estimates))
	print "in %s" % (path)
	
	foundestimates = glob.glob(os.path.join(path, "*.pkl"))
	
	print "I have found %i files that look like estimates that are ready." % (len(foundestimates))
	
	outests = []
	for estimate in estimates:
		resultpklpath = os.path.join(path,'%s.pkl' % estimate.id)
		if resultpklpath in foundestimates:
			outest = pycs.gen.util.readpickle(resultpklpath, verbose=False)
			outests.append(outest)
	
	print "Collected %i estimates" % (len(outests))
	if len(outests) != len(estimates):
		print "WARNING, some estimates are still missing, I'm returning only %i estimates." % len(outests)
	else:
		print "OK! I found every estimate."
	
	return outests
	
	
def displayresults(estimates, paths):
	"""
	I gather all the estimates of a single set/rung combination and plot the delays + errors for all the corresponding pairs.

	I use the delaycontainer class and the newdelayplot function from pycs.sim.plot

	WARNING : This is not designed to be used for tdc0 and tdc1 estimates, but for cosmograil data analysed in a PyCS@TDC way. If set == 'tdc1' or 'tdc0', I do not execute and return an error.
	"""

	'''
	# Check that every pack of estimates have a path
	try:
		len(lestimates) != len(paths)
	except:
		print "list of estimates and number of related paths do not match!"
		sys.exit()
	'''

	ests = []
	colours = ['black', 'blue', 'seagreen', 'crimson', 'orchid', 'dodgerblue', 'gold']


	# First step, collect all the estimates and sort them by set/rung:
	for path in paths:
		# This is utterly ugly, but it works.
		# Here we replace the rung attribute by the absolute path of the estimates. Thus, we take into account the simulation parameters, like nsim, ncopy or maxshift. As this is encoded in the path, we simply replace est.rung by current path basename.
		myests = collect(estimates, path)
		for myest in myests:
			myest.rung = os.path.basename(path)
		ests += myests
		#print ests[-1].td,ests[-1].pair, ests[-1].rung


	estsets = sorted(list(set([est.set for est in ests])))
	groupbysets = [[est for est in ests if est.set == estset] for estset in estsets]
	print "Grouped %i estimates of %i different lenses" % (len(ests), len(groupbysets))

	for group in groupbysets:
		estrungs = sorted(list(set([est.rung for est in group])))
		groupbyrungs = [[est for est in group if est.rung == estrung] for estrung in estrungs]
		print "Grouped %i estimates of %i different rungs" % (len(group), len(groupbyrungs))
		print estrungs

		toplot = []
		for ind, group in enumerate(groupbyrungs):
			# Then, create a delaycontainer object with all pycs plots.
			name = group[0].rung
			colour = colours[ind]
			nests = len(group)
			cops = []
			sims = []
			print 'WARNING !! estimates delay are negative if positive in PyCS convention ! I do not take that into account here.'
			for est in group:
				cops.append({'mean': est.td, 'label': est.pair, 'std': 0, 'med': 0 })
				sims.append({'tot': est.tderr,  'label': est.pair, 'sys': est.tdsyserr, 'ran': est.tdranerr, 'bias': 0 })

			results=(pycs.sim.plot.delaycontainer(data=cops,objects=['A','B','C','D'],name=name, plotcolour = colour),pycs.sim.plot.delaycontainer(data=sims,objects=['A','B','C','D'],name=name, plotcolour = colour))

			toplot.append(results)

		pycs.sim.plot.newdelayplot(toplot, rplot=5.0, displaytext=True)