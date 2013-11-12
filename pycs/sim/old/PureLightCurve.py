"""
Define a pure lightcurve
"""

from util import *
import numpy as np 
import pylab as py 
import scipy.ndimage as nd

class PureLightCurve:
	"""
	This class generate a discret hypothetic quasar light curve. At this time we can generate it by:
		1. a random walk
		2. a law noise generator

	(see the funtion description for further information)
	
	In order to control the curve, we can plot the psd the periodigram and the variability of the curve.
	
	This function generator remain basic. We have a special method to import some more realistic curve given by a fit of a real data. 
	
	"""		
	def __init__(self,length,res):
		
		"""
		Define all the main parameters


		@type	length:	int		
		@param	length: lenght of the light curve in days

		@type	res:	float	
		@param	res:	resolution of the light curve (number of point per day)

		"""

		self.type="purelightcurve"
		"""
		@type:	string
		@ivar:	type of the curve ( this is only for information and display)		
		"""
		self.length=length
		self.res=res
		self.data=[]
		"""
		@type:	float array	
		@ivar:	array containing all the light curve data
		"""

	def RandomWalk(self,walkstd):	
		"""
		create a basic random walk in 1D
		@type	walkstd:	float
		@param 	walkstd:	standart deviation for the normal random distribution
		"""

		randvec=np.random.normal(0,walkstd,size=self.length*self.res)
		self.data=np.cumsum(randvec)

	def smooth(self,smoothstd):
		"""
		method which smooth the curve with a gaussian filter(not used)
		@type	smoothstd:	float
		@param	smoothstd:	standart deviation for the gaussian filter
		"""

		self.data=nd.filters.gaussian_filter1d(np.abs(self.data),smoothstd)

	def normalize(self,norm):
		"""
		nomalize the data(not used)
		@type	norm:	float
		@param	norm:	norm of the data
		"""

		fact=norm/np.max(self.data)
		self.data=fact*self.data	

	def PowerLawNoise(self,beta,lawnoisestd):

		"""
		Law noise generator with controled power spectrum.
		This generator is supposed to be more realistic but it's quite hard to control the fft.
		
		@type	beta:	float
		@param	beta:	the order of the power spectrum
		
		@type	lawnoisestd:	float
		@param	lawnoisestd:	std of the generator (equivalent to amplitude of the curve)	
	
		The algorithm create a random spectrum with an order beta power spectrum and then transform it with a numpy fft into the light curve.
		
		@todo: Compared the curve with the real curve to find the good parameter as the frequency in the spectrum an beta(value used:beta=3.5 and amplitude depend of the quasar)

		"""
		omega = np.linspace(1,self.length*self.res,self.length*self.res )
		"""
		@type:	float array
		@ivar:	array containing all the frequencies for the generator
		"""

		spec = lawnoisestd*np.random.randn(self.length*self.res)/np.sqrt(omega**beta)
		"""
		@type:	float array
		@ivar:	the random spectrum of our light curve	
		"""
	
		self.data=np.fft.irfft(spec)
		self.data=np.delete(self.data,np.arange(self.length*self.res,2*self.length*self.res))

		#self.data=-2.5*np.log10(100+(self.data))	#convert magnitude data to flux data (just for testing)


	def Import(self,file):
		"""
		Import a pickle with a Pure curve data. Used for taking a pseudo real curve(HE0435) which we had fit with some fit program(gaussian process, spline...)
		Usefull if we always want the same original data which correspond to real data
		@type	file:	pickle
		@param	file:	pickle file with the data
		
		@todo: better gestion of the resolution and the length for future intensive use
		"""

		self.data=readpickle(file)
		self.res=2
		self.length=len(self.data)/self.res

	def PlotPsd(self):
		"""
		compute and plot the matplotlib psd

		"""

		py.psd(self.data,512,self.length)
		py.show()

	def ComputePer(self,beta):
		"""
		compute a periodigram for the data and compared it to a theorical line

		"""
			
		datafft=np.fft.fft(self.data)
		per=np.abs(datafft)**2
		py.plot(per)
		x=np.array([100.*i+1 for i in range(self.length*self.res)])
		y=per[0]*(1./x)**(beta/2)
		py.loglog(x,y,c='red')
		py.show()

	def variability(self):
		"""
		A test to compare variability of our curve with the experimental value given by Hook 1993. It doesn't match well but I don't know if it's very important
		@todo:finish the method and try to know if it could be use
		"""

		L=100	#maximum time difference
		l=9000	#length of the scan
		dt=np.arange(0,L,1/float(self.res))
		#print dt
		#print len(self.data)
		dmtmp=np.empty([l,L*self.res])
		dmmax=np.array(L*self.res*[0.])	
		dmmin=np.array(L*self.res*[0.])	
		#print dm
		for j in range(l): #loops for save all the delta magnitude
			for i in range(L*self.res):
				dmtmp[j,i]=np.abs(self.data[i+j]-self.data[j])
				if np.abs(self.data[i+j]-self.data[j])> dmmax[i]:
					dmmax[i]=np.abs(self.data[i+j]-self.data[j])

		dmmed=np.median(dmtmp,0)
		dmmean=np.mean(dmtmp,0)
		#print dm
		M=-26.0
		z=1
		dmth=abs(0.155+0.023*(M+25.7))*(dt/(1+z))**(0.18) #experimental values (Hook)
		dmth1=abs(0.155+0.023*(M+25.7))*(dt/(100+z))**(0.18)
		py.plot(dt,dmmed,'-b')
		py.plot(dt,dmmean,'-m')
		py.plot(dt,dmmax,'-y')
		py.plot(dt,dmth,'-r')
		py.plot(dt,dmth1,'-g')
		py.show()

	def __getitem__(self,key):
		"""
		Defining __getitem__ just for convenience
		@return:	return the data in the key position
		"""

		return self.data[key]	

	def save(self,file):
		"""
		Just a save method if needed
		"""

		writepickle(self.data,file)
