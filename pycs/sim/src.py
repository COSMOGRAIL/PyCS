"""
Stuff to represent a light curve as points on a regular grid, and power spectrum plots.
This allows e.g. to tweak its power spectrum by adding "correlated noise", resample it, etc.

"""
import sys
import numpy as np
import copy as pythoncopy
import matplotlib.pyplot as plt
from matplotlib.mlab import psd
import scipy.optimize as spopt
import scipy.interpolate as si
import pycs.gen.spl


class Source():
	"""
	A class representing a "source" lightcurve, i.e. any artificial signal whose powerspectrum we can modify,
	and still evaluate the curve at given irregular jds.
	
	To do this, we work with "internal regular arrays" (and do FFTs and so on on these), and then interpolate these
	fine arrays to get back values for given jds.
	
	We store the data in magnitudes.
	
	Ideas : use the mask carried around by self.inispline.datapoints to define regions that
	really matter for the power specturm and stuff.
	
	"""

	def __init__(self, spline = None, name = "Source", range=(0, 10000), sampling = 0.2):
		"""
		range are in jds, only used if you don't give me a sourcespline.
		sampling is in days.
		
		I will store the spline as inispline, but I will not modify it, so no need to pass me a copy.
		
		"""
		if spline != None:
			self.inispline = spline
		
			self.jdmin = spline.datapoints.jds[0]
			self.jdmax = spline.datapoints.jds[-1]
			self.plotcolour = spline.plotcolour
		else:
			
			# Some default spline to play with and for tests :
			self.jdmin = range[0]
			self.jdmax = range[-1]
			jds = np.linspace(self.jdmin, self.jdmax, 100)
			mags = np.zeros(jds.shape)
			magerrs = 0.1 * np.ones(jds.shape)
			dp = pycs.gen.spl.DataPoints(jds, mags, magerrs, splitup=False, sort=False)
			self.inispline = pycs.gen.spl.Spline(dp)
			self.inispline.uniknots(3) # Puts coeffs to 0, that's ok.
			#self.inispline.optc()
			self.plotcolour = "grey"
		
		self.name = name

		self.sampling = sampling # Number of days between points of the internal arrays.
		# Note that we want an even number of internal points.
		# Hence, this sampling will be adapted a bit so to make the number of points even.
		self.setup()
		
	def setup(self):
		"""
		Sets the internal regular arrays by sampling the spline.
		We update self.sampling, so to get an even number of points (for FFTs).
		"""
		
		nsteps = int(((self.jdmax - self.jdmin)/self.sampling)/2.0)*2
		self.sampling = float(self.jdmax - self.jdmin)/float(nsteps)
		self.ijds = np.linspace(self.jdmin, self.jdmax, nsteps)
		self.imags = self.inispline.eval(jds = self.ijds)
	
	
	def __str__(self):
		return "%s (%i)" % (str(self.name), len(self.ijds))
		#return "%s from %s [%f, %f], %i/%.2f" % (self.name, self.inispline, self.jdmin, self.jdmax, len(self.ijds), self.sampling)
		
		
		
	def copy(self):
		return pythoncopy.deepcopy(self)

	
	
	
	def setmean(self, mean=-12.0):
		"""
		Shifts the magnitudes so that their mean value are at the specified level.
		Don't do this, just for testing purposes !
		When doing flux PS, this scales the power, that's all.
		"""
		self.imags += mean - np.mean(self.imags)
	
	def addgn(self, sigma=0.1, seed = None):
		"""
		Don't do this, just for testing purposes ! There is no reason to add white noise to the source !
		"""
		rs = np.random.RandomState(seed) # we create a random state object, to control the seed.
		self.imags += rs.standard_normal(self.imags.size)*sigma
	
	def addrw(self, sigma=0.1, seed = None):
		"""
		Add a random walk, also for experiement (power law param is -2)
		"""
		rs = np.random.RandomState(seed) # we create a random state object, to control the seed.
		randvec = rs.standard_normal(self.imags.size)*sigma
		self.imags += np.cumsum(randvec)
	
# 	def addplaw(self, beta=-3.0, sigma=0.01, flux=False, fmin=None, fmax=None, seed=None):
# 		"""
# 		Version 0.1, full of errors ...
# 		Adds noise according to a power law PSD.
# 		See Timmer & Koenig 1995
# 		"""
# 		# To simplify, we will generate a symmetric curve and use only half of it. Hence this symmetric curve will be twice as long.
# 		
# 		n = self.imags.size # The real number of points in our curve. This is even by construction.
# 		n2 = 2*n # This is even more even.
# 		
# 		freqs2 = np.linspace(0, 0.5, n2/2 + 1) # The 0 frequency is included, but we will leave it's power at 0.
# 		
# 		# To help representing what we do, here are :
# 		freqs = np.linspace(0, 0.5/self.sampling, n2/2 + 1) # Our associated frequencies, in units of days^-1
# 		
# 		# Now we generate the random coefficents for those freqs.
# 		rs = np.random.RandomState(seed) # we create a random state object, to control the seed.
# 		
# 		specs = rs.standard_normal(n2/2 +1) * freqs2**(beta/2.0)
# 		specs[0] = 0.0 # To get 0 power for the 0 frequency.
# 		
# 		# We want to set all frequencies lower than cutf to 0 :
# 		if fmin:
# 			specs[freqs <= fmin] = 0.0
# 		if fmax:
# 			specs[freqs >= fmax] = 0.0
# 		# Hmm, this is quite brutal. Try a window funtion, to avoid getting regular oscillations here.
# 		
# 		# We add the points to the values
# 		if flux:
# 			
# 			iflux = 10.0**(-0.4 * self.imags)
# 			iflux -= sigma * np.fft.irfft(specs)[:n2/2]
# 			print np.min(iflux)
# 			self.imags = -2.5 * np.log10(iflux)
# 		else:
# 			self.imags += sigma * np.fft.irfft(specs)[:n2/2]
		
		
	def addplaw2(self, beta=-3.0, sigma=0.01, flux=False, fmin=None, fmax=None, hann=False, seed=None):
		"""
		Next version, better
		Adds noise according to a power law PSD.
		See Timmer & Koenig 1995
		
		power law would be -2
		
		if hann, we soften the window (Hann window instead of tophat).
		
		"""
		# To simplify, we will generate a symmetric curve and use only half of it. Hence this symmetric curve will be twice as long.
		
		n = self.imags.size # The real number of points in our curve. This is even by construction.
		n2 = 2*n # This is even more even.
		
		# The positive FFT frequences, in relative units :
		freqs2 = np.linspace(0, 0.5, n2/2 + 1) # The 0 frequency is included, but we will leave it's power at 0.
		
		#print freqs2
		# To help representing what we do, here are the same in days^-1 :
		freqs = np.linspace(0, 0.5/self.sampling, n2/2 + 1) # Our associated frequencies, in units of days^-1
		
		# Now we generate the random coefficents for those freqs.
		rs = np.random.RandomState(seed) # we create a random state object, to control the seed.
		
		# Complex and imaginary part
		rspecs = rs.standard_normal(n2/2 +1) # same length as freqs2, associated
		rspecs[1:] *= freqs2[1:]**(beta/2.0)
		ispecs = rs.standard_normal(n2/2 +1)
		ispecs[1:] *= freqs2[1:]**(beta/2.0)
		rspecs[0] = 0.0 # To get 0 power for the 0 frequency.
		ispecs[0] = 0.0 # To get 0 power for the 0 frequency.
		
		# As we work with even number of signals, the Nyquist frequency term is real :
		ispecs[-1] = 0
		specs = rspecs + 1j * ispecs
		
		# We now build a mask
		if fmin == None: # if fmin is None
			fmin = freqs[0] - 1.0
		if fmax == None:
			fmax = freqs[-1] + 1.0
		windowmask = np.logical_and(freqs <= fmax, freqs >= fmin) # We will later set everything outside of this to 0.0
		
		outofwindowmask = windowmask == False # we invert the mask
		specs[outofwindowmask] = 0.0
		
		if hann:	
			bell = np.zeros(freqs2.shape)
			bell[windowmask] = np.hanning(np.sum(windowmask))
			specs *= bell		
		
		#windowfct = np.hanning(np.sum(windowmask))
		#windowfct = np.ones(np.sum(windowmask))
		#rs = int(len(windowfct)/100)
		#rs = 10
		#windowfct[:rs] = np.linspace(0.0, 1.0, rs)
		#specs[windowmask] *= windowfct
		#print specs
		#print np.fft.irfft(specs)[:n2/2]
		
		
		
		# We add the points to the values
		if flux:
			print "Implement flux !"
			print "Hmm, change sigma somehow, can't be the same as for mags !"
			
			iflux = 10.0**(-0.4 * self.imags) # The fluxes
			iflux += sigma * np.fft.irfft(specs)[:n2/2] # We add our "noise"
			#print np.min(iflux)
			self.imags = -2.5 * np.log10(iflux) # and get back to fluxes
		else:
			
			"""
			iflux = 10.0**(-0.4 * self.imags) # The fluxes
			medflux = np.median(iflux)
			magsigma = 
			"""
			
			self.imags -= sigma * np.fft.irfft(specs)[:n2/2] # -, to make it coherent with the fluxes.
		
	
	def eval(self, jds):
		"""
		I interpolate my ijds/imags to give you an array of mags corresponding to your jds
		This could in principle be done using the spline object made by the function below, but this is safer and faster.
		"""
		if np.min(jds) < self.jdmin or np.max(jds) > self.jdmax:
			raise RuntimeError("Sorry, your jds are out of bound !")
		
		f = si.interp1d(self.ijds, self.imags, kind="linear", bounds_error=True)
		return f(jds)



	def spline(self):
		"""
		I return a new pycs.gen.spl.Spline object corresponding to the source.
		So this is a bit the inverse of the constructor.
		You can then put this spline object as ML of a lightcurve, as source spline, or whatever.
		
		Note that my output spline has LOTs of knots... it is an interpolating spline, not a regression spline !
		
		
		..note:: This spline is a priori for display purposes only. To do an interpolation, it might be safer (and faster) to use
			the above linear interpolation eval() function.
			But making a plot, you'll see that this spline seems well accurate.
		"""
		
		x = self.ijds.copy()
		y = self.imags.copy()
		magerrs = np.zeros(len(x))
		
		out = si.splrep(x, y, w=None, xb=None, xe=None, k=3, task=0, s=0.0, t=None, full_output=1, per=0, quiet=1)
		# s = 0.0 means interpolating spline !
		
		tck = out[0]
		
		# From this we want to build a real Spline object.
		datapoints = pycs.gen.spl.DataPoints(x, y, magerrs, splitup=False, sort=False, stab=False)
		outspline = pycs.gen.spl.Spline(datapoints, t = tck[0], c = tck[1], k = tck[2], plotcolour=self.plotcolour)
		
		# Conditions for the knots (no, everything is ok with them, no need to tweak anything).
		#intt = outspline.getintt()
		#outspline.setintt(intt)
		#print tck[2]
		#sys.exit()
		
		outspline.knottype = "MadeBySource"
		outspline.showknots = False # Usually we have a lot of them, slows down.
		return outspline
		


def sourceplot(sourcelist, filename=None, figsize=(12, 8), showlegend=True, showspline=True, marker=None):
	"""
	I show you a plot of a list of Source objects.
	"""
	
	plt.figure(figsize=figsize)
	for s in sourcelist:
		if marker == None:
			plt.plot(s.ijds, s.imags, marker="None", color=s.plotcolour, linestyle="-", label = "%s" % (s.name))
		else:
			plt.plot(s.ijds, s.imags, marker=marker, color=s.plotcolour, linestyle="none", label = "%s" % (s.name))
	
	
		if showspline:
			spline = s.spline()
			xs = np.arange(s.ijds[0], s.ijds[-1], 0.02)
			ys = spline.eval(jds = xs)
			plt.plot(xs, ys, "-", color="red", zorder=+20, label="%s.spline()" % (s.name))
				
			
			plt.plot()
	
	# Something for astronomers only : we invert the y axis direction !
	axes = plt.gca()
	axes.set_ylim(axes.get_ylim()[::-1])
	plt.xlabel("HJD - 2400000.5 [days]", fontsize=14)
	plt.ylabel("Magnitude (relative)", fontsize=14)	
		
	if showlegend:
		plt.legend()

	if filename:
		plt.save(filename)
	else:
		plt.show()
	
	



def window_hanning(x):
	"""return x times the hanning window of len(x)"""
	return np.hanning(len(x))*x

def detrend_mean(x):
	"""Return x minus the mean(x)"""
	return x - x.mean()

class PS():
	"""
	A class representing a power spectrum of a Source object.
	"""
	
	def __init__(self, source, flux=False):
		"""
		Constructor, simply takes a Source object and calculates the power spectrum
		
		The Source is always expressed in magnitudes, but you might want to get a powerspectrum in terms
		of flux or "counts".
		Put flux = True, and I will convert your magnitudes into fluxes before calculating the power spectrum.
		"""
		
		self.flux = flux
		self.n = len(source.ijds) # Even by construction.
		
		if flux == False:
			x = source.imags
		else:
			x = 10.0**(-0.4 * source.imags)
		
		#(self.p, self.f) = psd(x, NFFT=self.n, Fs=1.0/source.sampling, sides='onesided')
		# Note that these frequencies of self.f have physical units, they are directly in days^-1 !
		#print self.f
		
		# An alternative using only fft by myself, strongly "inspired" by the mlab psd above 
		# (It gives the same results, of course, but avoids the stupid figure)
		
		# The frequencies, and preparing a window :
		self.f = np.linspace(0, 0.5/source.sampling, self.n/2 + 1)
		windowVals = window_hanning(np.ones(self.n))
		
		# The FFT and power :
		fx = np.fft.rfft(windowVals * x)
		p = np.abs(fx)**2
		
		# Scale the spectrum by the norm of the window to compensate for
		# windowing loss; see Bendat & Piersol Sec 11.5.2.
		p *= 1 / (np.abs(windowVals)**2).sum()

		# Also include scaling factors for one-sided densities and dividing by the
		# sampling frequency, if desired. Scale everything, except the DC component
		# and the NFFT/2 component:
		p[1:-1] *= 2.0 * source.sampling
		#But do scale those components by Fs, if required
		p[[0,-1]] *= source.sampling
		
		self.p = p
	
		self.plotcolour = source.plotcolour
	
		self.slope = None
		self.name = source.name
	
	def __str__(self):
		if self.flux:
			return "%s(flux)" % str(self.name)	
		else:	
			return "%s(mag)" % str(self.name)
	
	def copy(self):
		return pythoncopy.deepcopy(self)
	
	def calcslope(self, fmin=1./1000.0, fmax=1./2.0):
		"""
		Measures the slope of the PS, between fmin and fmax.
		All info about this slope is sored into the dictionnary self.slope.
		
		Should this depend on self.flux ?? No, only the constructor depends on this !
		This is just fitting the powerspectrum as given by the constructor.
		
		"""
		if fmin == None:
			fmin = self.f[1]
		if fmax == None:
			fmax = self.f[-1]
		 
		reg = np.logical_and(self.f <= fmax, self.f >= fmin)
		
		fitx = np.log10(self.f[reg]).flatten()
		fity = np.log10(self.p[reg]).flatten()
		
		if not np.all(np.isfinite(fity)):
			print "Skipping calcsclope for flat function !"
			return
			
		self.slope = {}
		
		def func(x, m, h):
			return m*x + h
		popt = spopt.curve_fit(func, fitx, fity, p0=[0.0, 0.0])[0]
		sepfit = func(fitx, popt[0], popt[1]) # we evaluate the linear fit.
		
		self.slope["slope"] = popt[0]
		self.slope["f"] = 10.0**fitx
		self.slope["p"] = 10.0**sepfit
		self.slope["fmin"] = fmin
		self.slope["fmax"] = fmax
		

def psplot(pslist, nbins = 0, filename=None, figsize=(12, 8), showlegend=True):
		"""
		Plots a list of PS objects.
		If the PS has a slope, it is plotted as well.
		
		if nbins > 0, I bin the spectra.
		
		add option for linear plot ?
		"""
		
		plt.figure(figsize=figsize)
		for ps in pslist:
		
			
			if not np.all(np.isfinite(np.log10(ps.p))):
				print "No power to plot (probably flat curve !), skipping this one."
				continue
			# We bin the points
			
			if nbins > 0:
				logf = np.log10(ps.f[1:]) # we remove the first one
				logbins = np.linspace(np.min(logf), np.max(logf), nbins+1) # So nbins +1 numbers here.
				bins = 10**logbins
				bincenters = 0.5*(bins[:-1] + bins[1:]) # nbins centers
				logbins[0] -= 1.0
				logbins[-1] += 1.0
				binindexes = np.digitize(logf, logbins) # binindexes go from 1 to nbins+1
				binvals = []
				binstds = []
				for i in range(1, nbins+1):
					vals = ps.p[1:][binindexes == i]
					binvals.append(np.mean(vals))
					binstds.append(np.std(vals)/np.sqrt(vals.size))
			
				bincenters = np.array(bincenters)
				binvals = np.array(binvals)
				binstds = np.array(binstds)
			
				plt.loglog(bincenters, binvals, marker=".", linestyle="-", color=ps.plotcolour, label = "%s" % (ps))
			
			else:
				plt.loglog(ps.f, ps.p, marker=".", linestyle="None", color=ps.plotcolour, label = "%s" % (ps))
			if ps.slope != None:
				plt.loglog(ps.slope["f"], ps.slope["p"], marker="None", color=ps.plotcolour, label = "Slope %s = %.3f" % (ps, ps.slope["slope"]))
				plt.axvline(ps.slope["fmin"], color = ps.plotcolour, dashes = (5,5))
				plt.axvline(ps.slope["fmax"], color = ps.plotcolour, dashes = (5,5))
		
		plt.xlabel("Frequency [1/days]")
		plt.ylabel("Power")
		
		if showlegend:
			plt.legend()
		
		#plt.text(np.min(10**fitx), np.max(10**pfit), "Log slope : %.2f" % (popt[0]), color="red")
		
		
		if filename:
			plt.save(filename)
		else:
			plt.show()
