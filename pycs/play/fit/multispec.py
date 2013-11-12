# """
# 
# Stuff to make multi-D dispersion spectra, for instance for Mayavi ...
# For now this works only for 4 curves, to give 3 dimensions.
# 
# """
# 
# import sys
# import numpy as np
# #import matplotlib.pyplot as plt
# 
# import pycs.gen.util as util
# import pycs.gen.lc as lc
# import pycs.gen.ml as ml
# 
# 
# 
# 
# 
# def cube(lcs, fitmethod, verbose=True, timewidth=30, timestep=1.0, filename="chi2cube.pkl"):
# 	"""
# 	3D specplot, calculates the chi2 over a cube of time-delays. And writes the result in a pickle.
# 	
# 	This pickle can then be looked at with Mayavi (see example below)
# 	"""
# 	
# 	if len(lcs)!=4:
# 		raise RuntimeError, "I want 4 lightcurves."
# 	
# 	lcsc = [l.copy() for l in lcs]
# 	
# 	# We apply microlensing
# 	for l in lcsc:
# 		if l.ml != None:
# 			l.applyml()
# 			
# 	def chi2(delays):
# 		lc.multisettimedelays(lcsc, delays)
# 		chi2 = fitmethod(lcsc)["chi2n"]
# 		return chi2
# 	
# 	initparams = np.concatenate([lc.multigettimedelays(lcsc)])
# 	print "Initial shifts : ", initparams
# 	
# 	timeshifts = np.arange(-(timewidth)*timestep/2.0, (timewidth+1)*timestep/2.0, timestep)
# 	cubeindexes = np.arange(timewidth + 1)
# 	
# 	print "Points to calculate :", len(timeshifts)**3
# 	chi2cube = np.zeros((timewidth+1, timewidth+1, timewidth+1))
# 	
# 	xshifts = timeshifts + initparams[0]
# 	yshifts = timeshifts + initparams[1]
# 	zshifts = timeshifts + initparams[2]
# 	
# 	for ix in cubeindexes:
# 		print "Slice %i of %i" % (ix + 1, timewidth+1)
# 		for iy in cubeindexes:
# 			for iz in cubeindexes:
# 				chi2cube[ix, iy, iz] = chi2([xshifts[ix], yshifts[iy], zshifts[iz]])
# 	
# 	
# 	beg = -(timewidth)*timestep/2.0
# 	end = (timewidth+1)*timestep/2.0
# 	step = timestep
# 	
# 	x, y, z = np.mgrid[beg:end:step, beg:end:step, beg:end:step]
# 	#print x, y, z
# 	x += initparams[0]
# 	y += initparams[1]
# 	z += initparams[2]
# 	#print x, y, z
# 	
# 	util.writepickle({"lcs":lcs, "x":x, "y":y, "z":z, "chi2":chi2cube}, filename)
# 	
# 	
# 	
# #   To give an idea how to plot such a data cube with Mayavi2/mlab :
# #		import sys
# #		sys.path.append("../")
# #		from pycs.gen import *
# #		import numpy as np
# #		from enthought.mayavi import mlab
# #		
# #		pkldict = util.readpickle("dispcube50.pkl")
# #		
# #		maxval = 1.5
# #		minval = 1.43
# #		
# #		x = pkldict["x"]
# #		y = pkldict["y"]
# #		z = pkldict["z"]
# #		d2 = pkldict["d2"]
# #		
# #		lcs = pkldict["lcs"]
# #		
# #		minpos = np.argmin(d2)
# #		minpos = np.unravel_index(minpos, d2.shape)
# #		min_x = x[minpos]
# #		min_y = y[minpos]
# #		min_z = z[minpos]
# #		
# #		
# #		mlab.clf()
# #		
# #		src = mlab.pipeline.scalar_field(x, y, z, d2)
# #		
# #		# in green, the minimum
# #		mlab.points3d([min_x], [min_y], [min_z], color=(0,1,0), mode="cube", scale_mode="none", resolution=14, scale_factor=0.15)
# #		
# #		mlab.pipeline.scalar_cut_plane(src, vmin=minval, vmax=maxval)
# #		
# #		mlab.colorbar(title='Dispersion', orientation='vertical')
# #		
# #		mlab.xlabel("%s%s"% (lcs[0].object, lcs[1].object))
# #		mlab.ylabel("%s%s"% (lcs[0].object, lcs[2].object))
# #		mlab.zlabel("%s%s"% (lcs[0].object, lcs[3].object))
# #		
# #		
# #		mlab.show()
# 
# 
# 
# 
# 
# 
# 
# 
# 
