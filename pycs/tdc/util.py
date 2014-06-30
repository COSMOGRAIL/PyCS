"""
General purpose functions related to the TDC.
"""

import os
import numpy as np
import math
import pycs.gen.lc
import pycs.tdc.est
import datetime


def listtdc1v2pairs():
	"""
	A list of all pairs that are available in TDC1 v2 (useful, as some pairs have been removed...)
	"""
	tdc1v2doubles = np.arange(1, 721)
	tdc1v2quads = np.array(
		[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 
		34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 
		65, 66, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 95, 
		96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 
		122, 123, 124, 125, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 139, 140, 141, 142, 143, 144, 145, 146, 
		147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158]) # Removed are: 11,15,67,94,126,138
	return list(np.sort(np.concatenate([tdc1v2doubles, 720+(2*tdc1v2quads-1), 720+(2*tdc1v2quads)])))


def tdcfilepath(set, rung, pair, skipset=False):
	"""
	
	For tdc1 :
	rung is in [0, 4]
	pair is in [1, 1036]... but some are missing.
	
	"""
	
	if set=='tdc0':
		return "%s/rung%i/%s_rung%i_pair%i.txt" % (set, rung, set, rung, pair)


	if set=='tdc1':
	
		if pair<=720:
			if skipset == False:	
				return "%s/rung%i/%s_rung%i_double_pair%i.txt" % (set, rung, set, rung, pair)
			else:
				return "rung%i/%s_rung%i_double_pair%i.txt" % (rung, set, rung, pair)
		else:
			modpair  = int(pair-720) 		# 1, 2, 3, 4
			quadpair = (modpair+1) // 2 	# 1, 1, 2, 2, ...
			quadcode = "A" if modpair % 2 == 1 else "B"
			 
			"""
			intdiv   = modpair/2    # gives the integer portion
			floatdiv = modpair/2.0  # gives the float value
			
			if abs(intdiv-floatdiv)> 1e-5:
				quadpair = str(intdiv+1)+str('A')
			else:
				quadpair = str(intdiv)+str('B')	
			"""
			
			if skipset == False:
				return "%s/rung%i/%s_rung%i_quad_pair%s%s.txt" % (set, rung, set, rung, quadpair, quadcode)	
			else:
				return "rung%i/%s_rung%i_quad_pair%s%s.txt" % (rung, set, rung, quadpair, quadcode)



def pogmag(flux, fluxerr, m0 = 22.5):
	"""
	Computes a "normal" Pogson magnitude from TDC info.
	Used for TDC0

	flux and fluxerr are two numpy arrays
	
	"""
	# This is the code used for TDC0 :
# 	# First we compute the error bars :
# 	l.magerrs = 2.5 * np.log10(1.0 + (l.magerrs / l.mags))
# 	# Then the magnitudes, assumging nanomaggies :
# 	l.mags = 22.5 - 2.5 * np.log10(l.mags)
	
	mag = m0 - 2.5 * np.log10(flux)
	magerr = 2.5 * np.log10(1.0 + (fluxerr / flux))

	return (mag, magerr)
	
	
	
def asinhmag(flux, fluxerr,  m0 = 22.5, f0=1.0, b=0.01):
	"""
	Implements
	http://ssg.astro.washington.edu/elsst/opsim.shtml?lightcurve_mags
	"""

	mag = m0 -(2.5/np.log(10.)) * ( np.arcsinh( flux / (f0 * 2.0 * b)) + np.log(b) )
	
	magplu = m0 -(2.5/np.log(10.)) * ( np.arcsinh( (flux+fluxerr) / (f0 * 2.0 * b)) + np.log(b) )
	magmin = m0 -(2.5/np.log(10.)) * ( np.arcsinh( (flux-fluxerr) / (f0 * 2.0 * b)) + np.log(b) )
	magerr = 0.5*(magmin - magplu)
	
	return (mag, magerr)


def read(filepath, mag="asinh", verbose=True, shortlabel=True):
	"""
	Imports TDC light curves, WITH ASINH MAGNITUDES 
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
		
		if mag == "pog":
			(l.mags, l.magerrs) = pogmag(l.mags, l.magerrs, m0 = 22.5)
		if mag == "asinh":
			(l.mags, l.magerrs) = asinhmag(l.mags, l.magerrs, m0 = 22.5, b=0.01)
			
	return lcs



def godtweak(estimates):
	"""
	I return a list of estimates in which I have modified some by divin inspiration.
	This should be systematically done before writing a submission.
	"""
	
	outestimates = estimates[:] # make a copy
	pycs.tdc.est.checkunique(outestimates)
	
	print "#"*80
	print "This is god speaking, I am fixing your estimates..."
	print "#"*80
	
	for est in outestimates:
		
		# The doubtless ones with Thibault going mad:
		if (est.rung, est.pair) == (3, 505):
			est.td = 91.6
			est.tderr = 2.9
			print "Fixed", est		
		if (est.rung, est.pair) == (3, 593):
			est.tderr = 4.0
			print "Fixed", est	
		if (est.rung, est.pair) == (0, 557):
			est.tderr = 3.5
			print "Fixed", est
		if (est.rung, est.pair) == (3, 612):
			est.tderr = 3.0
			print "Fixed", est
	
	
	
	return outestimates
	


def writesubmission(estimates, filepath, commentlist=None, theseonly=False):
	"""
	Writes a list of estimates into a TDC1 submission
	
	Accordign to the doc, filepath should be something like pycs_tdc1_algorithm.dt
	
	if theseonly is set to True, I will not write -99 flags for those estimates that are missing in your list !

	"""
	
	pycs.tdc.est.checkunique(estimates)
	#est.sort(estimates) # We do not rely on this sorting anymore. We just go through all curves anyway...

	
	tdcfile = open(filepath, "w")
	
	tdcfile.write("# TDC submission written by PyCS\n")
	tdcfile.write("# %s\n" % (datetime.datetime.now()))
	tdcfile.write("# \n")
	tdcfile.write("# Malte Tewes <mtewes@astro.uni-bonn.de>, Uni Bonn\n")
	tdcfile.write("# Vivien Bonvin <vivien.bonvin@epfl.ch>, EPFL\n")
	tdcfile.write("# Frederic Courbin <frederic.courbin@epfl.ch>, EPFL\n")
	tdcfile.write("# \n")
	tdcfile.write("# Orignial filename: %s\n" % (os.path.basename(filepath)))
	tdcfile.write("# %i estimates have been selected\n" % (len(estimates)))
	tdcfile.write("# Comments:\n")
	for comment in commentlist:
		tdcfile.write("#   %s\n" % comment)
	tdcfile.write("# \n")
	tdcfile.write("# datafile                         dt          dterr\n")
	
	
	estids = [est.id for est in estimates]
	notfound = 0
	for rung in np.arange(5):
		for pair in listtdc1v2pairs():	
		
			# We reconstruct the filename of this pair. Carefully.
			name = tdcfilepath(set="tdc1", rung=rung, pair=pair)
			name = os.path.basename(name)
						
			searchid = 'tdc1_%i_%i' % (rung, pair)
					
			try:				
				ind = estids.index(searchid)
				tdcfile.write("%s\t%8.2f\t%8.2f\n" % (name, -1.0 * estimates[ind].td, estimates[ind].tderr))
				# --- WARNING --- The TDC convention for the delays is the inverse of PyCS, thus the "-" sign above
				
				# Just a little check to make Vivien happy: **Woohoo** !
				'''
				if rung==3 and pair==505:
					assert estimates[ind].td == 91.6 and estimates[ind].tderr == 2.9
				'''
			except ValueError:	
				notfound += 1
				if not theseonly:
					tdcfile.write("%s\t%8.2f\t%8.2f\n" % (name, -99.0, -99.0))
	
	
	#tdcfile.write("\n")
	tdcfile.close()	

	print "Wrote %s" % (filepath)
	print "Number of estimates: %i" % (len(estimates))
	print "Number of light curve pairs that I could not find among your estimates: %i" % (notfound)
	if theseonly:
		print "WARNING, I did not write -99 flags for the missing estimates !"
	


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

'''
def cutlcs(lcs, nseasons=3.0,overlapfrac=0.2):

	# APPARENTLY UNUSED, TO BE DELETED...? (10.01.2014)

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

'''
'''
def cutdb(filepath,newpath,season):

	import sys	
	"""
	small hand-tuned function to cut the d3cs database in two
	"""

	db=open(filepath,'r')
	newdb=open(newpath,'w')
	
	
	lines = db.readlines()

	#for ind,line in enumerate(lines):
		#print ind,line
		
	if season == 2:
		cut = 7518 
	
	for line in lines[cut:]:
		newdb.write(line)	
	 	 
	db.close()
	newdb.close()
'''
		
def goingon():
	"""
	Asks the user if he wants to proceed. If not, exits python.
	"""
	import sys
	answer = raw_input("Tell me, do you want to go on ? (yes/no) ")
	if answer[:3] != "yes":
		sys.exit("Ok, bye.")
	print ""	# to skip one line after the question.

		
	
