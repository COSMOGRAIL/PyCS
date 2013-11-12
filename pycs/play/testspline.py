import pycs
import numpy as np

jds = np.linspace(0.0, 100.0, 50)
mags = np.sin(0.1*jds) + 0.1*np.random.randn(len(jds))
magerrs = 0.1*np.ones(len(jds))

datapoints = pycs.gen.spl.DataPoints(jds, mags, magerrs, splitup=False, deltat=0.000001, sort=True, stab=False, stabext=300.0, stabgap = 30.0, stabstep = 5.0, staberr = 0.1)

spline = pycs.gen.spl.Spline(datapoints, t = None, c = None, k = 3, plotcolour="black")

spline.uniknots(4, n=True)
spline.optc()

spline.display(showbounds = False, showdatapoints = True, showerrorbars=True)

print spline.getintt()

print spline.getc()