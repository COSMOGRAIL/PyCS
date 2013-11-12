"""
Generate simulated light curve with the PureLightCurve
"""

from PureLightCurve import *
from SimLightCurve import *
from MlCurve import *
from util import *
import sys
sys.path.append("../../")
from pycs.gen import lc


class Generator:
	"""
	This class is the main class of this program, it takes a pure light curve and create N simulated light curves in a list.
	This class is composed by some methods which create the differents parts of the curve( sampling, errorbar, seasons...)

	@todo: Maybe a general method like Construct(all the parameter) will be good	
	
	"""
	
	def __init__(self,OriginLightCurve,NLens=4,names=["A","B","C","D"],color=["red","blue","green","magenta"],length=20,period=1,dmag=[0]):
		"""
		Initialise some general variable

		@type	NLens:	int
		@param	NLens:	number of simulated light curve

		@type	names:	list of string
		@param	names: list of the simulated light curve names

		@type	length:	float
		@param	length: lenght of the simulated data in days

		@type 	period: float
		@param	period:	period of the simulated mesure
		
		@type	dmag:	array of float
		@param	dmag:	array containig the magnitude shift for each curve
		"""

		self.olcurve=OriginLightCurve

		"""
		@type: PureLightCurve
		@ivar: The pure light curve which generate the simulation
		"""

		self.curves=NLens*[SimLightCurve]

		"""
		@type:list of SimLightCurve
		@ivar: a list of NLens simulated light curve
		"""

		self.length=length
		self.period=period

		for i in range(NLens):	#construct the list with empty SimLightCurve
			self.curves[i]=SimLightCurve(names[i],plotcolor=color[i],length=self.length,period=self.period,originalcurve=self.olcurve,dmag=dmag[i])	

		
	def Shift(self,shift=[0, 1, -2,10]):
		"""
		Pseudo shifting method because it just save the shifting parameter in the each SimLightCurve object
		
		@type	shift:	array of float
		@param	shift:	array of NLens length which take all the time shift parameters

		"""
		self.shift=shift
		for i in range(len(shift)):	
			self.curves[i].shift=shift[i]	

	def create_timeline(self,jdsamplingstd):	
		"""
		Create an array with the value in day of each mesure. It is based on the period of mesurement but with some deviations according to the bad wether for a night or something else.
		@type	jdsamplingstd:	float
		@param	jdsamplingstd:	standart deviation for the normal distribution which simulate the unperiodicity of the mesures
		
		"""
		if jdsamplingstd != 0:	#we can remove this error for some test 
			std=np.random.normal(0,jdsamplingstd,self.length/self.period)	#create the vector of error
		else:
			std=np.array(self.length/self.period*[0.])

		for curve in self.curves:
			f = np.arange(self.length/self.period)
			curve.datatime=f*self.period+std

	def create_seasons(self,pos, size):
		"""
		This method create the seasons. It's a brut force method because it completly delete the data beween each season. Maybe we can change that...
		
		@type	pos:	array of float
		@param	pos: 	array of the position for each begining of inter-seasons

		@type	size:	array of foat
		@param	size:	array of the lenght of each inter-seasons

		"""

		seas=[0]

		for i in range(len(size)):
			seas=np.append(seas,np.linspace(pos[i],pos[i]+size[i],size[i]))	#define the inter seasins period

		for curve in self.curves:	#delete it from the data
			curve.datatime=np.delete(curve.datatime,seas)
			curve.datamag=np.delete(curve.datamag,seas)
			curve.dataerr=np.delete(curve.dataerr,seas)
			
	def create_errorbar(self,low,high): 
		"""
		Create an array of errorbar with a normal distribution(obsolete)
		Now this is a uniform distribution between low and high values.

		@type	low:	float
		@param	low:	minimum of the uniform distribution

		@type	high:	float
		@param	high:	maximum of the uniform distribution

		"""

		for curve in self.curves:
			#curve.dataerr=np.abs(np.random.normal(0,errorstd,self.length/self.period)) #normal generator(obsolete)
			curve.dataerr=np.abs(np.random.uniform(low,high,self.length/self.period))	#uniform generator

	def create_data(self,erramp,mlbeta,mlstd):
		"""
		This is the main method which create the sampling data. It shift and sample the data according to the shift parameter and the timeline.It add some micro lensing effect in each curve too.

		@type	erramp:	float
		@param	erramp: amplitude of the error of each point( the errorbar give the standart deviation), 0 or 1 in general

		@type	mlbeta:	float
		@param	mlbeta:	beta parameter for microlensing curve generator

		@type	mlstd:	float
		@param	mlstd:	std parameter for microlensing generator (equivalent to an amplitude parameter)

		"""

		for curve in self.curves:

			std=erramp*np.random.randn(len(curve.dataerr))*curve.dataerr	#error vector
			curve.mlcurve=MlCurve(self.length,self.olcurve.res)		#microlensing vector
			curve.mlcurve.LawNoise(mlbeta,mlstd)				#generate microlensing

			for (i, time) in enumerate(curve.datatime):
				if (time+curve.shift)*self.olcurve.res> self.olcurve.length*self.olcurve.res: 
					break
				if (time+curve.shift)>0:
					curve.datamag[i]=self.olcurve[(time+np.floor(curve.shift)-1)*self.olcurve.res]+std[i]+curve.mlcurve.data[(time+np.floor(curve.shift)-1)*self.olcurve.res]
			#curve.datamag=-2.5*np.log10(curve.datamag)	#convert magnitude to flux (just for testing purpose don't use it with micro lensing which must be multiplicative in this case)

	def save(self):
		"""
		Save the data in a dictionnary of lightcurve define by the shifting part of the program
		So we can directly take this pickle to test the main program

		@todo: maybe give different name for archive
		"""
	
		dict={}
		for curve in self.curves:
			dict[curve.name]=lc.factory(curve.datatime,curve.datamag,magerrs=curve.dataerr,telescopename="simulation",object=curve.name)

		writepickle(dict,'lcs.pkl')
