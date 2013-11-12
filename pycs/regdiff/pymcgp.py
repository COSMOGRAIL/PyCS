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



def regression(x, y, yerr, priormeanfct, pow=1.5, amp=2.0, scale=200.0, errscale=5.0, verbose=True):
	"""
	Give me data points
	
	yerr is the 1sigma error of each y
	
	I return a function : you pass an array of new x, the func returns (newy, newyerr)
	
	pow, amp and scale are params for the covariance function.
	
	"""
	obs_mesh = x
	obs_vals = y
	obs_V = yerr*yerr # Converting std to variance
	
	if verbose:
		print "Computing GPR with params pow=%.1f, amp=%.1f, scale=%.1f, errscale=%.1f" % (pow, amp, scale, errscale)
	# pymc GP objects :
	M = pymc.gp.Mean(priormeanfct)
	
	# v1, best restults for RXJ1131 -> default values of these options
	#C = pymc.gp.Covariance(eval_fun = matern.euclidean, diff_degree=1.5, amp=2.0, scale=200.)
	#obs_V *= 5.0
	
	C = pymc.gp.Covariance(eval_fun = matern.euclidean, diff_degree=pow, amp=amp, scale=scale)
	obs_V *= errscale
	
	#C = pymc.gp.Covariance(eval_fun = matern.euclidean, diff_degree=1.5, amp=1.0, scale=100.)

	# v2 :
	#C = pymc.gp.Covariance(eval_fun = matern.euclidean, diff_degree=1.5, amp=0.5, scale=100.)
	#obs_V *= 3.0
	
	# v3 (?) :
	#C = pymc.gp.Covariance(eval_fun = pow_exp.euclidean, pow=1.7, amp = 2.0, scale = 300.)
	
	# sandbox :
	#C = pymc.gp.Covariance(eval_fun = matern.euclidean, diff_degree=1.5, amp = 1.0, scale =120.)
	#C = pymc.gp.Covariance(eval_fun = pow_exp.euclidean, pow=1.7, amp = 2.0, scale = 200.)
	#C = pymc.gp.Covariance(eval_fun = gaussian.euclidean, amp = 0.6, scale = 100.)

	# Impose observations on the GP
	pymc.gp.GPutils.observe(M, C, obs_mesh=obs_mesh, obs_V = obs_V, obs_vals = obs_vals)

	def outfct(jds):
		(m_out, v_out) = pymc.gp.GPutils.point_eval(M, C, jds)
		
		newy = m_out
		newyerr = np.sqrt(v_out)
		
		return (newy, newyerr)
	
	#print "Done"
	return outfct


