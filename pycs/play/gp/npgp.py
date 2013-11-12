"""
Playing with GP

"""
import numpy as np


def rbf_kernel(x, y, sigma = 100.0):
	"""
	
	RBF = radial basis function
	This is the covariance matrix, high value if i is close to j (strong corelation), 0 if i is far from j.
	
	function K = Krbf(trX,tsX,sig)
	trL = size(trX,1);
	tsL = size(tsX,1);

	>>>for i=1:trL
	>>> for j=1:tsL
	>>>  K(i,j)=exp(-norm(trX(i,:)-tsX(j,:))^2/(2*sig^2));
   	>>> end
	>>>end
	"""
	xl = x.size
	yl = y.size
	K = np.zeros((xl, yl))
	for i in range(xl):
		for j in range(yl):
			K[i, j] = np.exp( -(x[i] - y[j])**2.0 / (2*sigma*sigma))	
	return K

def regression(lc, n, sigma = 100.0):
	
	
	x = lc.getjds()
	y = lc.getmags() - np.mean(lc.getmags())

	K = rbf_kernel(x, x, sigma = sigma)
	print "K Done !"
	beta = 0.1

	smoothx = np.linspace(min(x), max(x), n)
	
	a = K + beta**2.0 * np.eye(x.size)
	s = np.linalg.solve(a,y)
	smoothy = np.dot(rbf_kernel(x, smoothx, sigma=sigma).transpose(), s)
	
	# axY = Krbf(trX,axX,sigma)' * ((K + beta^2 * eye(trL)) \ trY);
	
	return (smoothx, smoothy + np.mean(lc.getmags()))



