"""
The top of the food chain :-)
Global optimizers that make use of the building blocks from multiopt.py
The functions take raw lightcurves (any microlensing will be deleted),
eventually ML, and optimize all this.
They return a float, the dispersion value.
They can be used with :py:func:`pycs.sim.run.multirun`

"""
import sys, os
import pycs.gen.lc
import pycs.disp.multiopt

	


def opt_full(lcs, rawdispersionmethod, iter=5, verbose=True):
	"""
	A first optimization of ml and ts, alternatively doing one after the other.
	Works great, seems to be a keeper.
	Note that I do **keep** the ML.
	
	You can put a rawdispersionmethod (i.e. non symmetric) into these guys, as they will anyway apply it on AB and BA for every pair AB.
	"""
	
	if verbose:
		
		print "Starting dispersion optimization of :\n%s" % ("\n".join([str(lc) for lc in lcs]))
		print "Initial delays :"
		print "%s" % (pycs.gen.lc.getnicetimedelays(lcs, separator=" | "))	
		
	# We start with some magnitude shifts :
	pycs.disp.multiopt.opt_magshift(lcs, rawdispersionmethod, verbose=False)
	
	for i in range(iter):
	
		d2 = pycs.disp.multiopt.opt_ml(lcs, rawdispersionmethod, maxit = 1, verbose=False)
		#if verbose:
		#	print "Iteration %i ml done, d2 = %8.3f" % (i+1, d2)
		d2 = pycs.disp.multiopt.opt_ts_mix(lcs, rawdispersionmethod, movefirst=False, verbose=False)
		#if verbose:
		#	print "Iteration %i ts done, d2 = %8.3f" % (i+1, d2)
		
		if verbose:
			print "Iteration %i done, d2 = %8.3f" % (i+1, d2)
			print "%s" % (pycs.gen.lc.getnicetimedelays(lcs, separator=" | "))	
		
	if verbose:
		print "Done with optimization of :\n%s" % ("\n".join([str(lc) for lc in lcs]))
	return d2



def opt_full_addML(lcs, rawdispersionmethod, **kwargs):
	"""
	We include the microlensing addition into the optimizer, to be similar to the spline optimizers.
	This will by principle allow to randomize also the ML.
	For instance if shuffle of multirun is True
	"""
	for l in lcs:
		l.rmml()
	for l in lcs[1:]:
		pycs.gen.splml.addtolc(l, knotstep=150)
	return pycs.disp.opt1.opt_full(lcs, rawdispersionmethod, **kwargs)
	
	
	
	
	