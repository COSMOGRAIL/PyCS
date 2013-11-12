"""
Define a microlensing curve
"""

import numpy as np
import pylab as py

class MlCurve:
	"""
	Create a typical microlensing curve with two different method: law noise generator or magnification map
	@todo: magnification map...or not...
	"""

	def __init__(self,length, res):
		"""
		Set some parameter
		
		@type 	length:	int
		@param	length:	lenght of the curve(in days)

		@type	res:	float
		@param	res:	resolution of the curve

		"""

		self.length=length
		self.res=res
		self.data=[]
	
	def LawNoise(self,beta,lawnoisestd):
		"""
		Create the data with a lawnoise generator( as for the Pure light Curve
		
		@type	beta:	float
		@param	beta:	beta parameter for the generator
	
		@type 	lawnoisestd:	float
		@param	lawnoisestd:	std of the generator (equivalent to amplitude)	
	
		"""

		omega = np.linspace(1,self.length*self.res,self.length*self.res )
		spec = lawnoisestd*np.random.randn(self.length*self.res)/np.sqrt(omega**beta)
		
		self.data=np.fft.irfft(spec)		
		self.data=np.delete(self.data,np.arange(self.length*self.res,2*self.length*self.res))
