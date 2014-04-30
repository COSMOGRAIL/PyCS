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
import matplotlib.pyplot as plt
from copy import copy

def fitsourcespline(lcs, sploptfct, addmlfct=None):
	"""
	Return a spline fitting the given lightcurves
	You can add microlensing to your lcs before fitting the spline
	"""

	lca = lcs[0]
	lcb = lcs[1]
	
	# Adding the ML
	lca.magshift = 0.0
	lcb.magshift = 0.0
	if addmlfct != None:
		addmlfct([lca, lcb])

	# And run a first spline
	sourcespline = sploptfct([lca, lcb])
	
	
	return sourcespline




def createdir(estimates, path):
	
	dirpath = os.path.join(os.getcwd(),path)


	# 1) Create main directory
	
	if not os.path.isdir(dirpath):
		os.mkdir(dirpath)
		
	# 2) Create sub-directories for each estimate
	
	for est in estimates:
		subdirpath = os.path.join(dirpath,est.id)
		if not os.path.isdir(subdirpath):
			os.mkdir(subdirpath)	
	

def drawcopy(estimates, path, addmlfct=None, n=1):
	"""
	Draw n times one single copy of lcs (the ones your esimate is about) into the path directory.
	
	:estimates	- list of estimate objects, with id, td and tderr
	
	:path		- where your copy will be written
	
	:addmlfct	- If you want to add different microlensing to your obslcs
	
	:n		- Number of time you will run drawobs (each run produce new copycurves)
	"""
	

	
	copydir = os.path.join(os.getcwd(),path)


	'''
	if os.path.isfile(copydir)==False:
		try:
			os.system('mkdir %s' % copydir)
		except:
			print "going on..."	
	'''
	
	for estimate in estimates:
	
		lcspath = pycs.tdc.util.tdcfilepath(set=estimate.set, rung=estimate.rung, pair=estimate.pair)
		(lca, lcb) = pycs.tdc.util.read(lcspath, shortlabel=False)

		ind=0

		while ind < n:

			
			# Reset ML and shifts (get them "as read" from scratch...)
			lca.rmml()
			lcb.rmml()	
			lca.magshift = 0.0
			lcb.magshift = 0.0	
			lca.timeshift = 0.0
			lcb.timeshift = 0.0
			
			
			lcs = [lca,lcb] 

			# Add new ML
			if addmlfct != None:
				addmlfct(lcs)

			# Time shift around the initial value		
			tsr = np.max([3.0, estimate.tderr])

			lca.shifttime(float(np.random.uniform(low=-tsr, high=tsr, size=1)))
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
	
	
def drawsim(estimates, path, sploptfct, addmlfct=None, n=1) :
	"""
	Draw n times one single sim curves of lcs (the ones your esimate is about) into the path directory.
	
	:estimates	- list of estimate objects, with id, td and tderr
	
	:path		- where your copy will be written
	
	:addmlfct	- If you want to add different microlensing to your obslcs
	
	:n		- Number of time you will run drawsim (each run produce new simcurves)	

	"""

	simdir = os.path.join(os.getcwd(),path)
	
	'''
	if os.path.isfile(simdir)==False:
		try:
			os.system('mkdir %s' % simdir)
		except:
			print "going on..."		
	'''	
	for estimate in estimates:				
	
		# get the model lightcurves	
		lcspath = pycs.tdc.util.tdcfilepath(set=estimate.set, rung=estimate.rung, pair=estimate.pair)
		(lca, lcb) = pycs.tdc.util.read(lcspath, shortlabel=False)
		#pycs.gen.lc.display([lca,lcb])
		#sys.exit()
		ind=0
		
		# shift lcs as estimated by d3cs
		lcb.shifttime(estimate.td)
		#pycs.gen.lc.display([lca,lcb])
		#sys.exit()
		
		# fit a spline on the data and save residuals
		sourcespline = fitsourcespline([lca, lcb], sploptfct, addmlfct)
		pycs.sim.draw.saveresiduals([lca, lcb], sourcespline)
		
		# define 'initial timeshifts', given by fitsourcespline
		timeshifta = lca.timeshift
		timeshiftb = lcb.timeshift		
		#pycs.gen.lc.display([lca,lcb],[sourcespline])
		#sys.exit()

		#TODO : propagate knotstep attribute from varioanalysis to every copy, to avoid running varioanalysis for every simcurve...

		# now, draw n copycurves
		while ind < n:
	
			print 'IND---->',ind
						
			# set  a random "true" delay 
			truetsr = np.max([3.0, estimate.tderr])
			lca.timeshift = timeshifta + (float(np.random.uniform(low = -truetsr, high = truetsr, size=1)))
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

			print 'TRUE DELAY: ',lcssim[1].truetimeshift-lcssim[0].truetimeshift
			#pycs.gen.lc.display(lcssim)			
			#sys.exit()

			# Remove magshift, remove ML and add a new one :
			lcssim[0].magshift = 0.0
			lcssim[1].magshift = 0.0
			lcssim[0].rmml()
			lcssim[1].rmml()
			if addmlfct != None:
				addmlfct(lcssim)	

			tsr = np.max([3.0, estimate.tderr])

			# Set some wrong "initial delays" for the analysis, around the "true delays".
			lcssim[0].shifttime(float(np.random.uniform(low=-tsr, high=tsr, size=1)))
			lcssim[1].shifttime(float(np.random.uniform(low=-tsr, high=tsr, size=1)))
			
			print 'TRUE DELAY  : ',lcssim[1].truetimeshift-lcssim[0].truetimeshift
			print 'WRONG INITIAL DELAY: ',lcssim[1].timeshift-lcssim[0].timeshift						
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
			

def runcopy(estimates, path, optfct, n=1, clist=None):
	"""
	Run the optimizer (optfct) on n different copycurves
	Return the optimised timeshift and magnitude for each copycurves
	
	The n copycurves are chosen randomly in the copydir, unless you give a clist of numbers
	i.e. if clist = [1,3] runobs will run on c001.pkl and c003.pkl
	"""
	copytdslist  = []
	copymagslist = []
	
	
	for estimate in estimates:
	
		index_list = copy(clist)
	
		copydir = os.path.join(os.getcwd(),path)	

		# Check that there is enough copycurves in the copydir:

		find_max_ind = 'ls %s | wc -l' % os.path.join(copydir,estimate.id,'c*')
		number_of_files = int(os.popen(find_max_ind).readline())

		if number_of_files < n:
			print 'Not enough copyfiles in %s ' % os.path.join(copydir,estimate.id)
			print 'I will run on %i files only' % number_of_files
			n = number_of_files


		if clist == None:

			mylist = np.arange(number_of_files)
			np.random.shuffle(mylist)
			index_list = mylist[:n] 	

		copytds  = []
		copymags = []

		for ind in index_list:

			copyfile = 'c%s' % str(ind).rjust(3,'0')+'.pkl'
			copypath = os.path.join(copydir,estimate.id,copyfile)

			lcs = pycs.gen.util.readpickle(copypath)

			#pycs.gen.lc.display(lcs)
			#sys.exit()
			out = optfct(lcs)

			if hasattr(out, "bokeps"): # then its a spline !!
				displaysplines = [out]
			else:
				displaysplines = []


			copytds.append(lcs[1].timeshift-lcs[0].timeshift)
			copymags.append(np.median(lcs[1].getmags() - lcs[1].mags) - np.median(lcs[0].getmags() - lcs[0].mags))

		copytdslist.append(copytds)
		copymagslist.append(copymags)
	
	return copytdslist, copymagslist		
		
		
def runsim(estimates, path, optfct, n=1, slist=None):		
	"""
	Run the optimizer (optfct) on n different simcurves
	Return the optimised timeshift and magnitude for each copycurves
	Also return the true delays of each simcurve
	
	The n simcurves are chosen randomly in the simdir, unless you give a slist of numbers
	i.e. if slist = [1,3] runobs will run on s001.pkl and s003.pkl
	"""		
	
	simtdslist = []
	simmagslist = []
	simttdslist = []
	
	for estimate in estimates:
	
		index_list = copy(slist)
		
		simdir = os.path.join(os.getcwd(),path)		

		# Check that there is enough simcurves in the copydir:

		find_max_ind = 'ls %s | wc -l' % os.path.join(simdir,estimate.id,'s*')
		number_of_files = int(os.popen(find_max_ind).readline())

		if number_of_files < n:
			print 'Not enough simfiles in %s ' % os.path.join(simdir,estimate.id)
			print 'I will run on %i files only' % number_of_files
			n = number_of_files


		if slist == None:

			mylist = np.arange(number_of_files)
			np.random.shuffle(mylist)
			index_list = mylist[:n] 	

		simtds  = []
		simttds = []
		simmags = []

		# run the optimizer on each simcurve
		
		for ind in index_list:

			simfile = 's%s' % str(ind).rjust(3,'0')+'.pkl'
			simpath = os.path.join(simdir,estimate.id,simfile)

			lcs = pycs.gen.util.readpickle(simpath)

			
			out = optfct(lcs)
			print 'TRUE DELAY: ',lcs[1].truetimeshift-lcs[0].truetimeshift
			print 'WRONG GUESSED DELAY: ',lcs[1].timeshift-lcs[0].timeshift			
			#pycs.gen.lc.display(lcs)
			#sys.exit()
			if hasattr(out, "bokeps"): # then its a spline !!
				displaysplines = [out]
			else:
				displaysplines = []


			simttds.append(lcs[1].truetimeshift-lcs[0].truetimeshift)
			simtds.append(lcs[1].timeshift-lcs[0].timeshift)
			simmags.append(np.median(lcs[1].getmags() - lcs[1].mags) - np.median(lcs[0].getmags() - lcs[0].mags))
			#print "simttds",lcs[1].truetimeshift-lcs[0].truetimeshift
			#print "simtds",lcs[1].timeshift-lcs[0].timeshift
			#sys.exit()
		 			
		simtdslist.append(simtds)
		simmagslist.append(simmags)
		simttdslist.append(simttds)
		
	return simtdslist,simmagslist,simttdslist
		

def multirun(estimates, path, optfct, ncopy, nsim, clist=None, slist=None):
	"""
	Wrapper around runsim and runobs
	"""
	
	# We start by creating Estimate objects that will contain our final result:
	
	outests = []
	

	method = 'PyCS-Run'
	methodpar = str(optfct)
	
	
	for estimate in estimates:	
		outests.append(pycs.tdc.est.Estimate(set=estimate.set, rung=estimate.rung, pair=estimate.pair,
		 				method = method, methodpar = methodpar))

	
	
	# Then, we run the optimizers on the copy and the sims
	
	copytdslist,copymagslist = runcopy(estimates, path, optfct = optfct, n=ncopy, clist=clist)
	simtdslist,simmagslist,simttdslist = runsim(estimates, path, optfct = optfct, n=nsim, slist=slist)
	

	
	# And we put the results in the output estimates
	
	for outest,copytds,simtds,simttds in zip(outests,copytdslist,simtdslist,simttdslist):
		print '='*30
		print outest.niceid
		outest.td = np.median(copytds)
		print 'time delay:  ',outest.td
		syserr = np.fabs(np.mean([simtds[i]-simttds[i] for i in np.arange(len(simtds))]))
		ranerr = np.std([simtds[i]-simttds[i] for i in np.arange(len(simtds))])
		outest.tderr  = np.hypot(syserr, ranerr)
		print '   systematic error:  ',syserr
		print '       random error:  ',ranerr
		print '        total error:  ',outest.tderr		

	
	return outests
	
					
