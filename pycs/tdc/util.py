"""
General purpose functions related to the TDC.
"""

import os
import numpy as np
import pycs.gen.lc



def read(filepath, verbose=True, shortlabel=True):
	"""
	Imports TDC light curves.
	So far we expect exactly 2 curves in each file, as TDC simulates only doubles.
	Can be generalized later ...
	"""
	
	tdcid = os.path.splitext(os.path.basename(filepath))[0]
	
	if shortlabel == False:
		lc1 = pycs.gen.lc.flexibleimport(filepath, magcol=2, errcol=3, telescopename="TDC", object=tdcid+"_A", plotcolour="red", verbose = verbose)
		lc2 = pycs.gen.lc.flexibleimport(filepath, magcol=4, errcol=5, telescopename="TDC", object=tdcid+"_B", plotcolour="blue", verbose = verbose)
	else:
		lc1 = pycs.gen.lc.flexibleimport(filepath, magcol=2, errcol=3, telescopename="TDC", object="A", plotcolour="red", verbose = verbose)
		lc2 = pycs.gen.lc.flexibleimport(filepath, magcol=4, errcol=5, telescopename="TDC", object="B", plotcolour="blue", verbose = verbose)			

	lcs = [lc1, lc2]

	for l in lcs:
		#l.jds += 56586.0 # Starts on 21 October 2013 :)
		
		# First we compute the error bars :
		l.magerrs = 2.5 * np.log10(1.0 + (l.magerrs / l.mags))
		# Then the magnitudes, assumging nanomaggies :
		l.mags = 22.5 - 2.5 * np.log10(l.mags)
			
	return lcs



def setnicemagshift(lcs):
	"""
	Sets a "nice" magshift to the n-1 last curves of lcs, so that they appear
	nicely below each other when displayed.
	Also works if you curves contain ML models or are alredy shifted.
	"""
	
	bottom = np.max(lcs[0].getmags())
	for l in lcs[1:]:
		top = np.min(l.getmags())
		l.shiftmag(bottom-top)
		bottom = np.max(l.getmags())

def cutlcs(lcs, nseasons=3.0,overlapfrac=0.2):
	"""
	From a given set of lcs, return cutted lcs, only part of the original lcs
	The goal is then to run the optimisation on these cutted lcs to estimate the error
	by taking the scatter between these estimations (not very smart, but fast)
	Generalised to n lightcurves in lcs, but only 3 cutted lcs returned (should be generalised later)
	
	This function return 3 cutted lcs for each lc 
	
	"""
	import sys
	
	
	# some inits...
	
	n = len(lcs)
	
	jds=[]
	for j in lcs[0].jds:
		jds.append(j)
		
	lseason = len(jds)/nseasons
	noverlap = int(overlapfrac*lseason)

	# we define the jds range...
			
	sinit  = [int(0) , int(lseason+noverlap)]
	smid   = [int(lseason-noverlap/2) , int(2*lseason+noverlap/2)]
	sfinal = [int(len(jds)-lseason-noverlap) , int(len(jds))] 
	
	# now, let's create the truncated lightcurves

	cuts=[]
	
	for l in lcs:
	
		jd1 = jds[sinit[0]:sinit[1]]
		jd2 = jds[smid[0]:smid[1]]
		jd3 = jds[sfinal[0]:sfinal[1]]
		
		mags = l.getmags()
		mag1 = mags[sinit[0]: sinit[1]]
		mag2 = mags[smid[0]: smid[1]]
		mag3 = mags[sfinal[0]:sfinal[1]]		
		
		magerrs = l.getmagerrs()
		magerr1 = magerrs[sinit[0]: sinit[1]]
		magerr2 = magerrs[smid[0]: smid[1]]
		magerr3 = magerrs[sfinal[0]:sfinal[1]]		
		
		lc1 = pycs.gen.lc.factory(jd1,mag1,magerr1)
		lc1.plotcolour='blue'
		lc2 = pycs.gen.lc.factory(jd2,mag2,magerr2)
		lc2.plotcolour='black'
		lc3 = pycs.gen.lc.factory(jd3,mag3,magerr3)
		lc3.plotcolour='red'		
		
		cuts.append([lc1,lc2,lc3])

	# Ok, now the cutted lcs are created. We want to return something like
	# [lcA_1, lcB_1] (first season),[lcA_2, lcB_2] (second season)... 
		
	return map(None,*cuts)


		
def importfromd3cs(filepath='d3cslog.txt', user='vivien-08/11', maxconflevel=4, writeto=None):	

	'''
	This function read the d3cs database located at filepath, and put the values
	added by user in a database. If you give a writeto path, the function will write
	a .pkl file there containing the database. 
	
	If a pkl of the db already exists and you give a writeto path, the function
	will update the existing database
	
	WARNING : for the same username+rung+pair, only the latest entry in the db is considered
	
	Return a database: db[rung][pair] gives you the corresponding delay
	'''
		
	print 'Looking for a database in: %s ...' %filepath

	eyedbfile = open(filepath, 'r')

	lines = eyedbfile.readlines()
	lines = [line.split(', ') for line in lines]
	

	### structure of the d3cslogfile

	date	 = 0    	# Date and time of submission
	ip	 = 1		# IP address 
	set	 = 2		# Username
	rung	 = 3		# Rung
	pair	 = 4		# Pair
	delay	 = 5		# Time shift [day]
	detdelay = 6		# Time shift uncertainty [day] 
	mag	 = 7		# Mag shift 
	flux     = 8		# Flux shift [fraction of median flux] 
	cl	 = 9		# Confidence category
	elapt 	 = 10		# Time between page load and submission [s] 
	
	
	# create database
	
	if writeto != None:
		if os.path.isfile(writeto):
			print 'database located in %s already exists' % writeto
		 	print 'going on...'

		else:
			print 'The database does not exist yet'
	  		print 'I will create a database on %s' % writeto

	  		'''
	  		structure of the database:
	  		[rung][pair][val]
				[val=0] = delay
				[val=1]	= deltadelay
	  		'''				

	  		# init size, hardcoded (must be changed ! booh! )

	  		nrungs   = 7
	  		npairs   = 8
	  		nvals 	 = 2
		
		
			print '%i rungs, with %s pairs of lightcurves each' %(nrungs,npairs) 
		
			eyedb=[]
			for nrung in np.arange(nrungs):
  				eyedb.append([])
				for npair in np.arange(npairs):
					eyedb[nrung].append([])
					for nval in np.arange(nvals):
						eyedb[nrung][npair].append([])		
			
		
		
			pycs.gen.util.writepickle(eyedb,writeto)
		

		eyedb = pycs.gen.util.readpickle(writeto)
	
	else:
	
	  	# init size, hardcoded (must be changed ! booh! )

	  	nrungs   = 7
	  	npairs   = 8
	  	nvals 	 = 2	
	

		print '%i rungs, with %s pairs of lightcurves each' %(nrungs,npairs) 
		
		eyedb=[]
		for nrung in np.arange(nrungs):
  			eyedb.append([])
			for npair in np.arange(npairs):
				eyedb[nrung].append([])
				for nval in np.arange(nvals):
					eyedb[nrung][npair].append([])	
											
	
	# Ok, now we fill the database (eyedb)  
	
	
	for line in lines:

		if (str(line[set]) == str(user) and int(line[cl]) < maxconflevel):
			print 'rung %i pair %i : eyedelay is %.2f' %(int(line[rung]), int(line[pair]), float(line[delay]) )
			eyedb[int(line[rung])][int(line[pair])-1][0].append(float(line[delay]))
			eyedb[int(line[rung])][int(line[pair])-1][1].append(float(line[detdelay]))					

	if writeto != None:
		pycs.gen.util.writepickle(eyedb,writeto)
		
	return eyedb	
