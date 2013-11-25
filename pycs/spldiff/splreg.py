"""
spline regression module...

"""
import os,sys
os.environ['OMP_NUM_THREADS'] = "1"

import pycs.gen.lc as lc
import pycs.gen.spl as spl
import pycs.gen.sea as sea


import numpy as np

#from scipy.interpolate import UnivariateSpline
#import matplotlib.pyplot as plt




def splreg(jds, mags, magerrs,knotstep=20.0, n=None, stab=True,stabext=300.0, stabgap=20.0, stabstep=5.0,
		stabmagerr=-2.0, stabrampsize=0, stabrampfact=1.0, bokit=1, bokeps=2.0, boktests=5,
		bokwindow=None, k=3, verbose=True):
		
	"""
	Give me datapoints x y, with errorbars on y
	I return a function : you pass an array of new x, the func returns (newy, newyerr)
	
	I return a spline as well, so I can use it for magshifts later 
	performed with a BOK spline fitting on the datas and the errors.
	The errors on y take into account the seasons in y 
	
	"""
	
	l     = lc.factory(jds,mags,np.zeros(len(jds))+1e-5)
	l_err = lc.factory(jds,magerrs,np.zeros(len(jds))+1e-5)
	
	spl_mag = spl.fit([l],knotstep=knotstep, n=n, stab=stab, stabext=stabext, stabgap=stabgap, stabstep=stabstep,
				stabmagerr=stabmagerr, stabrampsize=stabrampsize, stabrampfact=stabrampfact, bokit=bokit,
				bokeps=bokeps,boktests=boktests,bokwindow=bokwindow, k=k, verbose=verbose) 
				
	spl_magerr = spl.fit([l_err],knotstep=knotstep, n=n, stab=stab, stabext=stabext, stabgap=stabgap, stabstep=stabstep,
				stabmagerr=stabmagerr, stabrampsize=stabrampsize, stabrampfact=stabrampfact, bokit=bokit,
				bokeps=bokeps,boktests=boktests,bokwindow=bokwindow, k=k, verbose=verbose)
	
	interseasons = sea.autofactory(l,seasongap=30, tpe='interseasons')
	
	
	def normaldistrib(x, mu, sig):
    		return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
		
	def errorinterseason(minintsea,maxintsea,rsjd):
		
		# we modelise the error enveloppe as a gaussian curve, with sigma adapted to the length of interseason
	
		normrsjd = (rsjd-minintsea+1e-5)/(maxintsea-minintsea)		
		return normaldistrib(normrsjd,0.5,0.2)**1.5		      
	        
	def outfct(rsjds):
				
		newy = spl_mag.eval(rsjds) 
		newyerr = spl_magerr.eval(rsjds)
		
		#let's modify the error enveloppe for the rsjds between seasons:
		
		for interseason in interseasons:
			for ind,rsjd in enumerate(rsjds):
				
				#the height of the enveloppe will be related to the length of the interseason
				ndays = interseason[1]-interseason[0]
				multfact = ndays/8.0  #arbitrary factor, may be changed...
				
				if rsjd < interseason[1] and rsjd > interseason[0]:
					newyerr[ind] += newyerr[ind]*errorinterseason(interseason[0],interseason[1],rsjd)*multfact
		
		return (newy, newyerr)	
	
	return outfct, spl_mag
