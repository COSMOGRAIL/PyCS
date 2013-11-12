"""
Functions to study r2 landscapes (equivalent of dispersion spectra) and similar stuff.

"""

import sys, os
import pycs.gen.lc
import pycs.gen.spl
import pycs.gen.util

import pycs.spl.multiopt

import numpy as np
import scipy.optimize as spopt




def explore(lcs, sourcespline, tss):
	"""
	We explore a volume of time shifts, and calculate the r2 at each point, by
	optimizing the source spline coeffs only.
	
	tss (timeshifts) is a list of arrays of timeshifts to try.
	>>> tss = np.mgrid[0:0:1j, -10:10:5j, -10:10:5j, -10:10:5j]
	
	todo : see if optimizing the ML coeffs changes something
	"""
		
	def calc(x, y, z):
		
		
		# We make a local copy of the lcs and the spline :
		mys = sourcespline.copy()
		mylcs = [l.copy() for l in lcs]
		# Set the shifts :
		for (l, ts) in zip(mylcs[1:], [x, y, z]):
			l.shifttime(ts) # So this is relative, starting from the beginning.
		# Optimize what we want :
		r2 = pycs.spl.multiopt.opt_source(mylcs, mys, verbose=False)
	
		#print shifts, r2
		#return np.array([l.timeshift for l in mylcs] + [r2])
		
		return r2
	
	veccalc = np.vectorize(calc, otypes=[np.ndarray])
	return veccalc(*tss).astype(float)
	
	#r2s = np.zeros(tss[-1].shape)
	#print tss.shape
	
	"""
	vectors = [shifts.flatten() for shifts in tss]
	results = []
	for shifts in zip(*vectors):
		results.append(calc(shifts))
	
	return np.vstack(results)
	"""

# Some stuff to play with the output (in construnction) ...

"""
from enthought.mayavi import mlab

biga = pycs.gen.util.readpickle("biga.pkl").astype(float)

x,y,z = np.mgrid[-10:10:5j, -10:10:5j, -10:10:5j]


src = mlab.pipeline.scalar_field(biga)
mlab.pipeline.iso_surface(src)
#mlab.pipeline.iso_surface(src, contours=[s.max()-0.1*s.ptp(), ],)
mlab.pipeline.image_plane_widget(src,
                            plane_orientation='z_axes',
                            slice_index=1,
                        )

#mlab.contour3d(x,y,z,biga)

#mlab.show()


#print biga
sys.exit()

"""
"""

x,y,z = np.mgrid[-10:10:5j, -10:10:5j, -10:10:5j]


#mlab.contour3d(x, y, z, biga)


src = mlab.pipeline.scalar_field(biga)
mlab.pipeline.iso_surface(src)
#mlab.pipeline.iso_surface(src, contours=[s.max()-0.1*s.ptp(), ],)
mlab.pipeline.image_plane_widget(src,
                            plane_orientation='z_axes',
                            slice_index=1,
                        )

mlab.show()

sys.exit()
"""

"""
lcs = pycs.gen.util.readpickle("lcsopt.pkl")
s = pycs.gen.util.readpickle("source.pkl")

tss = np.mgrid[-10:10:5j, -10:10:5j, -10:10:5j]

#print tss.shape

biga = pycs.spl.multispec.explore(lcs, s, tss=tss)
pycs.gen.util.writepickle(biga, "biga.pkl")
#print biga.shape

sys.exit()
"""

"""

print biga

pycs.gen.util.writepickle(biga, "biga.pkl")

sys.exit()

biga = pycs.gen.util.readpickle("biga.pkl").transpose()

mlab.surf(biga[1], biga[2], biga[4])

"""
"""
mlab.clf()

src = mlab.pipeline.scalar_scatter(biga[1], biga[2], biga[3], biga[4])

mlab.pipeline.scalar_cut_plane(src)
#mlab.colorbar(title='r2', orientation='vertical')
"""

