"""
Wrapper around pymc's GP module

"""
import os
os.environ['OMP_NUM_THREADS'] = "1"

import pymc.gp
import pymc.gp.GPutils
from pymc.gp.cov_funs import matern
from pymc.gp.cov_funs import gaussian
from pymc.gp.cov_funs import pow_exp

import numpy as np

#from scipy.interpolate import UnivariateSpline
#import matplotlib.pyplot as plt



def regression(x, y, yerr, priormeanfct):
	"""
	Give me data points
	
	yerr is the 1sigma error of each y
	
	I return a function : you pass an array of new x, the func returns (newy, newyerr)
	"""
	obs_mesh = x
	obs_vals = y
	obs_V = yerr*yerr # Converting std to variance
	

	# pymc GP objects :
	M = pymc.gp.Mean(priormeanfct)
	C = pymc.gp.Covariance(eval_fun = matern.euclidean, diff_degree=1.5, amp=2.0, scale=200.)
	
	#C = pymc.gp.Covariance(eval_fun = matern.euclidean, diff_degree=1.5, amp = 1.0, scale =120.)
	#C = pymc.gp.Covariance(eval_fun = pow_exp.euclidean, pow=1.7, amp = 0.4, scale = 150.)
	#C = pymc.gp.Covariance(eval_fun = gaussian.euclidean, amp = 1.0, scale = 20.)

	# Impose observations on the GP
	pymc.gp.GPutils.observe(M, C, obs_mesh=obs_mesh, obs_V = obs_V, obs_vals = obs_vals)

	def outfct(jds):
		(m_out, v_out) = pymc.gp.GPutils.point_eval(M, C, jds)
		
		newy = m_out
		newyerr = np.sqrt(v_out)
		
		return (newy, newyerr)
	
	#print "Done"
	return outfct


