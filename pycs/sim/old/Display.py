"""
Display method for ploting the original curve and the simulated light curve
"""
	

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np 

def Display(curvelist):
	"""
	Take a list of curve and plot it according to its type
	@type	curvelist:	list of curve
	@param	curvelist:	list of curve	 for plotting

	@todo: improve the display with legend, color...
	"""
 
	plt.figure(figsize=(12,8))	# sets figure size
	axes = plt.gca()
	# Something for astronomers only : we invert the y axis direction !
	axes.set_ylim(axes.get_ylim()[::-1])
	# Astronomers like minor tick marks :
	minorxLocator = MultipleLocator(100)

	for curve in curvelist:
		if curve.type=="purelightcurve":
			plt.plot(linspace(0,curve.length,curve.length*curve.res),curve.data)	

		if curve.type=="simlightcurve":
			if curve.name=="A": plt.plot(np.linspace(0,curve.originalcurve.length,curve.originalcurve.length*curve.originalcurve.res),curve.originalcurve.data+curve.dmag,color="black",label="Original curve")	
			plt.errorbar(curve.datatime,curve.datamagoff(),curve.dataerr,ls='None',marker='.',ecolor="#BBBBBB", color=curve.plotcolor,label=str(curve.name)+str(curve.shift))
			plt.plot(np.linspace(0,curve.originalcurve.length,curve.originalcurve.length*curve.originalcurve.res),curve.mlcurve.data+curve.dmag,color=curve.plotcolor)#,label="microlensing for "+str(curve.name))	
		plt.xlabel('time [j]')
		plt.ylabel('Magnitude ')
		plt.legend( numpoints = 1, prop = fm.FontProperties(size = 10))
	plt.show()
