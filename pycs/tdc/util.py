"""
General purpose functions related to the TDC.
"""

import os, sys
import numpy as np
import pycs.gen.lc
import pycs.tdc.est
import datetime


def listtdc1v2pairs():
	"""
	A list of all pairs that are available in TDC1 v2 (useful, as some pairs have been removed...)

	@return: a list with all the existing pairs
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
	Return the path to a TDC light curve in the form of "set/rung/set_rung_pair.txt"

	@param set: "tdc0" or "tdc1"
	@param rung: in [0, 4] for tdc1
	@param pair: in [1, 1036] for tdc1 (some are missing)
	@param skipset: boolean. If true, skip the first "set" in the path

	@return: string. Path "set/rung/set_rung_pair.txt"

	"""
	
	if set=='tdc0':
		return "%s/rung%i/%s_rung%i_pair%i.txt" % (set, rung, set, rung, pair)


	elif set=='tdc1':
	
		if pair<=720:
			if skipset == False:	
				return "%s/rung%i/%s_rung%i_double_pair%i.txt" % (set, rung, set, rung, pair)
			else:
				return "rung%i/%s_rung%i_double_pair%i.txt" % (rung, set, rung, pair)
		else:
			modpair  = int(pair-720) 		# 1, 2, 3, 4
			quadpair = (modpair+1) // 2 	# 1, 1, 2, 2, ...
			quadcode = "A" if modpair % 2 == 1 else "B"
			
			if skipset == False:
				return "%s/rung%i/%s_rung%i_quad_pair%s%s.txt" % (set, rung, set, rung, quadpair, quadcode)	
			else:
				return "rung%i/%s_rung%i_quad_pair%s%s.txt" % (rung, set, rung, quadpair, quadcode)

	else: # if this is a cosmograil target
		return "%s/%s_%s.rdb" % (set, set, rung)



def pogmag(flux, fluxerr, m0 = 22.5):
	"""
	Computes a "normal" Pogson magnitude from TDC info.
	Used for TDC0

	@param flux: numpy array of floats
	@param fluxerr: numpy array of floats
	@param m0: float. Zero-point of the magnitude system.
	@return: tuple containing two lists of fluxes converted into magnitudes.
	
	"""
	
	mag = m0 - 2.5 * np.log10(flux)
	magerr = 2.5 * np.log10(1.0 + (fluxerr / flux))

	return (mag, magerr)
	
	
	
def asinhmag(flux, fluxerr,  m0 = 22.5, f0=1.0, b=0.01):
	"""
	Implements asinh magnitudes, following
	http://ssg.astro.washington.edu/elsst/opsim.shtml?lightcurve_mags

	@param flux: numpy array of floats
	@param fluxerr: numpy array of floats
	@param m0: float. Zero-point of the magnitude system.
	@param f0: float. asinh mag normalisation constant.
	@param b: float. Asinh mag normalisation constant.
	@return: tuple containing two lists of fluxes converted into magnitudes

	"""

	mag = m0 -(2.5/np.log(10.)) * ( np.arcsinh( flux / (f0 * 2.0 * b)) + np.log(b) )
	
	magplu = m0 -(2.5/np.log(10.)) * ( np.arcsinh( (flux+fluxerr) / (f0 * 2.0 * b)) + np.log(b) )
	magmin = m0 -(2.5/np.log(10.)) * ( np.arcsinh( (flux-fluxerr) / (f0 * 2.0 * b)) + np.log(b) )
	magerr = 0.5*(magmin - magplu)
	
	return (mag, magerr)


def read(filepath, mag="asinh", verbose=True, shortlabel=True):
	"""
	Imports TDC light curves as a PyCS Lightcuve object.
	So far we expect exactly 2 curves in each file, as TDC simulates only doubles.

	@param filepath: string. Path to the tdc lightcurves, typically the output of tdcfilepath
	@param mag: can be "asinh" of "pog". Defines the type of magnitude considered
	@param verbose: boolean. Defines verbosity
	@param shortlabel: boolean. If True, gives a shorter object name to each light curve, omitting its TDC id.

	@return: Return two PyCS Lightcurve objects, one per tdc light curve
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
		#l.jds += 56586.0 # Starts on 21 October 2013

		if mag == "pog":
			(l.mags, l.magerrs) = pogmag(l.mags, l.magerrs, m0 = 22.5)
		if mag == "asinh":
			(l.mags, l.magerrs) = asinhmag(l.mags, l.magerrs, m0 = 22.5, b=0.01)

	return lcs





def writesubmission(estimates, filepath, commentlist=None, theseonly=False):
	"""
	Writes a list of estimates into a TDC1 submission


	@param estimates: list of Estimate objects, for which you want to write the time delay and error estimations.
	@param filepath: path of the file to be written. Accordign to the doc, filepath should be something like pycs_tdc1_algorithm.dt
	@param commentlist: list of srings. Additionnal comments to be written in the header of the submission file.
	@param theseonly: boolean. If set to True, I will not write -99 flags for those estimates that are missing in your list !

	"""
	
	pycs.tdc.est.checkunique(estimates)

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

			except ValueError:	
				notfound += 1
				if not theseonly:
					tdcfile.write("%s\t%8.2f\t%8.2f\n" % (name, -99.0, -99.0))

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

	@param lcs: list of Lightcurve objects
	"""
	
	bottom = np.max(lcs[0].getmags())
	for l in lcs[1:]:
		top = np.min(l.getmags())
		l.shiftmag(bottom-top)
		bottom = np.max(l.getmags())

		
def goingon():
	"""
	Asks the user if he wants to proceed. If not, exits python.
	"""
	import sys
	answer = raw_input("Tell me, do you want to go on ? (yes/no) ")
	if answer[:3] != "yes":
		sys.exit("Ok, bye.")
	print ""	# to skip one line after the question.

		
	
def readsubmission(filepath):
	"""
	Read a TDC1 submission file, and return its values in PyCS standards.

	@param filepath: string. Path to the TDC1 submission

	@return: list of tuples. Each tuple are (estimate id, estimated delay, estimated error)
	"""
	
	subinfos = []

	# we read the submission file and remove the comments
	lines = open(filepath,'r').xreadlines()
	cleanlines = [line for line in lines if line[0] != '#']
	
	# now, we extract the values we are interested in (id, td, tderr)	
	# Typical lines of a .dt submission file :
	#  	tdc1_rung0_double_pair1.txt	  -70.45	   11.70
	#	tdc1_rung4_quad_pair142A.txt	   17.55	    5.70

	for line in cleanlines:
	
		infos = [info for info in line.split(' ') if info != '']
		
		# get the id
		filename = infos[0].split('\t')[0]
		
		rung	 = filename.split('rung')[1].split('_')[0]
		pair_tdc = filename.split('pair')[1].split('.txt')[0]
		
		
		if 'A' in pair_tdc or 'B' in pair_tdc: # then it's a quad
		
			quadcode  = pair_tdc[-1] # get the A or B
			quadadd	  = -1 if quadcode == 'A' else 0
			pair_tdc  = pair_tdc.split('A')[0].split('B')[0]
			pair_pycs = 2*int(pair_tdc) + 720 + quadadd
		
		else: # then it's a double
			pair_pycs = pair_tdc

		estid   = 'tdc1_%s_%s' % (rung,pair_pycs)

		
		# get the delay td and the error tderr
		try:
			td    = infos[1].split('\t')[0]
			tderr = infos[2].split('\n')[0]
		except:
			print infos
			sys.exit()	

		if td != '-99.00' and tderr != '-99.00':
			subinfos.append((estid,-float(td),float(tderr))) # The - sign is as conventions change between tdc and pycs 

	print "I read submission %s, with %i entries" %(os.path.basename(filepath), len(subinfos))
	return subinfos		



def godtweak(estimates):

	"""
	I return a list of estimates in which I have modified some by divine inspiration.
	This should be systematically done before writing a TDC1-D3CS submission.

	@param estimates: list of Estimate objects
	@return: tweaked list of estimates
	"""

	print "WARNING - I was designed for a very specific purpose, and I will probably screw up your list of estimates. Don't proceed unless your are REALLY sure you know what you are doing"
	goingon()

	outestimates = estimates[:] # make a copy
	pycs.tdc.est.checkunique(outestimates)

	print "#"*80
	print "This is god speaking, I am fixing your estimates..."
	print "#"*80

	for est in outestimates:

		# The doubtless ones where we disagreed:
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
