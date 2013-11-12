"""
Define a simulated light curve
"""

import numpy as np

from util import *
class SimLightCurve:
	"""
	This class has anly a few method because it's just a definition class for the simulated light curve. It contain the data of the simulated light curve, some information and the original light curve.
	The data will be created by Generator class which can produce N simulated light curve like this one.
	"""

	def __init__(self,name,plotcolor,length,period,originalcurve,dmag):
		"""
		Define all te data for a simulated light curve

		
		@type	name:	string
		@param	name:	define the name of the light curve (the name of the deformed quasar image)

		@type	length:	int
		@param	length:	length of the light curve ( in days)

		@type 	period:	float
		@param	period:	define the period of the mesure ( in days), the mesure is around this period 
		
		@type	originalcurve:	PureLightCurve
		@param	originalcurve:	the original data for the light curve ( for comparaison an display)	
		
		@type	dmag:	float
		@param	dmag:	define a offset for the magnitue( for a more realistic display)
	
		"""
	
		self.name=name
		self.plotcolor=plotcolor
		self.length=length
		self.period=period
		self.originalcurve=originalcurve
		self.dmag=dmag

		self.type="simlightcurve"
		"""
		@type:	string
		@ivar:	the type of the curve (for information and display)	
		"""


		self.datamag=np.empty(self.length/self.period)
		"""
		@type:	float array
		@ivar:	magnitude data of the light curve
		"""

		self.datatime=np.empty(self.length/self.period)
		"""
		@type:	float array
		@ivar:	time data equivqlent to the date of the mesurements	
		"""

		self.dataerr=np.empty(self.length/self.period)
		"""
		@type:	float array
		@ivar:	simulated error in the data
		"""

		self.shift=0
		"""
		@type:	float
		@ivar:	define the virtual time shift of the curve
		"""
		

	def save(self,file): #save the data
		"""
		Save the data in a pickle
		maybe obsolete
		@todo: complete it
		"""

		writepickle([self.datatime,self.datamag],file)

	def timelength(self):
		"""
		@return:	return the maximum time of the mesure
		@todo:		can be more intelligent

		"""

		return np.max(self.datatime)
	
	def datamagoff(self):
		"""
		@return: 	return the magnitude data with the magnitude offset dmag

		"""

		return self.datamag+self.dmag


