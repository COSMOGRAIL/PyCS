import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pycs

"""
Machinery to combine various probability distributions that are not necessarily Gaussian, nor analytical.

To combine various distributions, they need to be "linearized" on the same basis. For each delay, I define a interval range (see binlists below). This range will be used to draw numerical distributions matching the theoretical one, if you give me e.g. mean and std of an analytical Gaussian distrib. Make sure to have ranges large enough and binnings thin enough to properly reproduce the distributions.


.. note:: A good idea would be to rewrite this using a Class to store the estimates. They are nothing more than dicts but doing so might improve the readability of the code.


.. todo:: Use the terminology of the PG1115 paper: an estimate is a time-delay measurement + uncertainty, a Group is all the estimates between all possible pairs, and a series is a list of groups. Group should be the Class here.
"""

class CScontainer():
	"""
	Curve-shifting container: basically just a bunch of informations related to the PyCS structure of storing results.

	Cloud be a dictionnary as well, e.g.

	{"data": "C2", "knots": 25, "ml": "splml-150", "name": "C2_ks25_ml150", "drawopt": "spl1", "runopt": "spl1t5", "ncopy": 200, "nmocks": 1000, "truetsr": 3, "colour": "royalblue"}

	This is simply a nicer interface between a PyCS simulation output and a Group.
	"""

	def __init__(self, data, knots, ml, name, drawopt, runopt, ncopy, nmocks, truetsr, colour):

		self.data = data
		self.knots = knots
		self.ml = ml
		self.name = name
		self.drawopt = drawopt
		self.runopt = runopt
		self.ncopy = ncopy
		self.nmocks = nmocks
		self.truetsr = truetsr
		self.colour = colour


class Group():
	"""
	A class representing all the time-delay estimates from a single PyCS curve-shifting technique (1 choice of estimator, 1 set of estimator parameters, 1 data set) - i.e. all the delays between the various pairs of images.

	Groups should contain:

	  * a list of labels, i.e. ["AB", "AC", "BC"] to which all other lists should correspond
	  * a list of median delays (50th percentiles)
	  * a list of errors_up and errors_down on the delay, corresponding to the 16th and 84th percentiles
	  * an ugly but unique name

	Extra arguments computed inside the class:

	  * a list of "linearized" distributions
	  * a fancy name for display purposes

	"""

	def __init__(self, labels, medians, errors_up, errors_down, name, binslist=None, lins=None, nicename=None):

		# start with some assertion tests
		assert (len(labels) == len(medians) == len(errors_up) == len(errors_down))

		# mandatory params
		self.labels = labels
		self.medians = medians
		self.errors_up = errors_up
		self.errors_down = errors_down
		self.name = name


		# optional params
		self.binslist = binslist
		self.nicename = nicename
		self.lins = None

	def linearize(self, testmode=True, verbose=True):
		"""
		I create a Gaussian distribution with mean = my median and error = (errors_up + errors_down)/2 that I spread on my binslist

		.. warning:: Of course, if your mean and error are not expected to follow a Gaussian distribution, calling this function makes no sense.

		:param testmode: bool. Defines the number of samples drawn from the Normal distribution. More samples give more precise results but take longer to create and process. Setting testmode to True create 50k samples, otherwise 500k.
		"""
		assert (self.binslist is not None)

		if testmode:
			nsamples = 50000
		else:
			nsamples = 500000

		lins = []
		for ind, label in enumerate(self.labels):
			scale = (self.errors_up[ind] + self.errors_down[ind])/2.0
			xs = np.random.normal(loc=self.medians[ind], scale=scale, size=nsamples)

			digitized = np.digitize(xs, self.binslist[ind])
			bin_means = [float(len(xs[digitized == i])) for i in range(1, len(self.binslist[ind]))]
			bin_means = np.nan_to_num(bin_means)
			bin_means = bin_means / np.sum(bin_means)

			# To check we are not too far off
			ci = confinterval(self.binslist[ind][:-1], weights=bin_means, testmode=testmode)
			if verbose:
				print "="*45
				print "Delay %s:" % label, "%.2f" % self.medians[ind], "+/-", "%.2f" % scale," --> ", "%.2f" % ci[0], "+%.2f-%.2f" % (ci[2], ci[1])

			lins.append(bin_means)

		self.lins = lins


	def sort(self, sortlabels):
		"""
		Sort all my list according to sortlabel

		"""

		assert len(sortlabels) == len(self.labels), "There are more sorting labels than actual labels in your object %s" % self.name

		pass
		#@todo: to be written!!!


	def niceprint(self):
		"""
		Nicely print in the terminal the median and 1sigma values of the current group.
		"""
		if not self.nicename:
			name = self.name
		else:
			name = self.nicename

		print "="*5, name, "="*5
		for ind, l in enumerate(self.labels):
			print "%s:" % l, "%.2f +%.2f-%.2f" % (self.medians[ind], self.errors_up[ind], self.errors_down[ind])



def getcombweightslist(groups):
	"""
	Give me 2+ groups. If they have the same binning (return an error otherwise), I return their bin-by-bin linear combination

	@param groups: list of Group. must have a list of linearized distributions

	@return: weights of the combined distribution
	"""

	# assert that each group has linearized distributions
	for group in groups:
		assert hasattr(group, "lins"), "Group %s has no linearized distributions - you must compute them firts!" % group.name


	# assert that the binnings are the same and that labels are in the same order
	for group in groups[1:]:
		assert group.binslist == groups[0].binslist, "Group %s has different bins than group %s" % (group.name, groups[0].name)
		assert group.labels == groups[0].labels, "Group %s has different labels (or labels order) than group %s" % (group.name, groups[0].name)

	# we assume the groups are ordered according to their labels (assert test would have failed otherwise)
	name = "|".join([g.name for g in groups])
	weightslist = {"weights": [], "binslist": groups[0].binslist, "labels": groups[0].labels, "name": name}

	for ind, label in enumerate(groups[0].labels):
		indweights = [g.lins[ind] for g in groups]

		"""
		This is trickily written so let me comment
		indweights contains a list of various linearized corresponding normalised distributions binned the same way (i.e. AB pdf for spl20 and spl30). The values of distrib are added bin by bin, then reweighted so that the sum of the result is still normalised to 1 (so they can be combined later with other results from various techniques, instruments, etc...) Since we work with numpy arrays, this can be done in a single line.
		"""
		weights = sum(indweights)/float(len(indweights))
		weightslist["weights"].append(weights)
	return weightslist



def getresults(csc):
	"""
	Give me a CScontainer, I read the corresponding result folder and create a group from it.

	.. warning:: This works **only** if you already collected the runresults in delaycontainer, using e.g. the :py:func:`pycs.sim.plot.hists` and :py:func:`pycs.sim.plot.measvstrue` functions.

	@param cscontainer: CScontainer object
	@return: a Group object
	"""

	# get the correct results path
	cp = "%s_ks%i_%s/sims_copies_%s_n%i_opt_%s_delays.pkl" % (
	csc.drawopt, csc.knots, csc.ml, csc.data, csc.ncopy, csc.runopt)
	mp = "%s_ks%i_%s/sims_mocks_%s_n%it%i_opt_%s_errorbars.pkl" % (
	csc.drawopt, csc.knots, csc.ml, csc.data, csc.nmocks, csc.truetsr, csc.runopt)

	# collect the results
	result = (pycs.gen.util.readpickle(cp), pycs.gen.util.readpickle(mp))

	# read them
	labels = [d["label"] for d in result[0].data]
	means = []
	errors = []

	for label in labels:
		means.append([d["mean"] for d in result[0].data if d["label"] == label][0])
		errors.append([d["tot"] for d in result[1].data if d["label"] == label][0])

	# create a Group out of them
	return Group(labels=labels, medians=means, errors_up=errors, errors_down=errors, name=csc.name)


def asgetresults(weightslist, testmode=True):
	"""
	Transform any given weightlist into a Group.

	We keep getcombweightlist and asgetresults separated, so we can import any binned result as a Group

	@param weightslist: a dictionnary containing a list of labels, a list of bins and a list of weights, whose order matches. Typically the output of getcombweightlist.
	@param testmode:
	@return: a Group object
	"""
	medians = []
	errors_up = []
	errors_down = []
	lins = []
	for b, w in zip(weightslist["binslist"], weightslist["weights"]):
		ci = confinterval(b[:-1], w, testmode=testmode)
		medians.append(ci[0])
		errors_up.append(ci[2])
		errors_down.append(ci[1])
		lins.append(w)

	return Group(labels=weightslist["labels"], medians=medians, errors_up=errors_up, errors_down=errors_down, lins=lins, name=weightslist["name"])




def confinterval(xs, weights=None, conflevel=0.68, cumulative=False, testmode=True):
	"""
	Hand-made quantile/percentile computation tool

	@param xs: list, distribution
	@param weights: list, weights of the distribution
	@param conflevel: float, quantile (between 0 and 1) of the conf interval
	@param cumulative: bool. If set to True, return the "bound" (max) value of xs corresponding to the conflevel. For example, if conflevel=0.68, xs=nnu, return the value such that nnu<=val at 68% confidence.
	@param testmode: bool. Defines the finesse of the binning of the histogram used when computing the percentiles. A finer binning gives more precise results but takes longer to process. Setting it to True sets 500 bins, otherwise 10'000.

	@return: return 3 floats: the median value, minimal value and maximal value matching your confidence level
	"""

	if testmode:
		nbins = 500
	else:
		nbins = 10000

	hist, bin_edges = np.histogram(xs, weights=weights, bins=nbins)

	if not cumulative:
		frac = 0.
		for ind, val in enumerate(hist):
			frac += val
			if frac/sum(hist) >= 0.5:
				meanval = (bin_edges[ind]+bin_edges[ind-1])/2.0
				break

		lowerhalf = 0.
		upperhalf = 0.
		for ind, val in enumerate(bin_edges[:-1]):
			if val < meanval:
				lowerhalf += hist[ind]
			else:
				upperhalf += hist[ind]


		limfrac = (1.0 - conflevel) / 2.0
		lowerfrac, upperfrac = 0., 0.
		for ind, val in enumerate(hist):
			lowerfrac += val
			if lowerfrac/sum(hist) > limfrac:
				bottomlimit = (bin_edges[ind]+bin_edges[ind-1])/2.0
				break

		for ind, val in enumerate(hist):
			upperfrac += val
			if upperfrac/sum(hist) > 1.0 - limfrac:
				toplimit = (bin_edges[ind]+bin_edges[ind-1])/2.0
				break

		return meanval, (meanval-bottomlimit), -1.0*(meanval-toplimit)

	else:
		frac = 0.
		for ind, val in enumerate(hist):
			frac += val
			if frac/sum(hist) >= conflevel:
				return (bin_edges[ind]+bin_edges[ind-1])/2.0




def get_bestprec(groups, refimg=None, verbose=False):
	"""
	Return the most precise Group among a list of Groups.

	The most precise Group is computed in a straightforward, naive way: simply sum the relative precision over all the pair of delays that contain refimg in their label (or all, if no refimg given), then take the estimate that has the lowest sum.

	To avoid giving too much weight to the larger delays, we compute the precision relative to the "mean" delay from all estimates.

	@param groups: list of Groups
	@param refimg: reference image (e.g. "A"), if you want to discard all the delays that have not "A" in their label

	@return: index of the most precise Group in the list
	"""

	# compute average means:
	avgmeans = []
	for l in groups[0].labels:
		medians = [group.medians[group.labels.index(l)] for group in groups]
		avgmeans.append(np.mean(medians))


	precs = []
	for group in groups:
		thisgroupprecs = []

		if refimg is None:
			mylabels = group.labels
		else:
			mylabels = [l for l in group.labels if refimg in l]

		for l in mylabels:
			mean = avgmeans[groups[0].labels.index(l)]
			err = (group.errors_up[group.labels.index(l)] + group.errors_down[group.labels.index(l)])/2.0

			thisgroupprecs.append(err/np.abs(mean))
			if verbose:
				print "="*5
				print "%s: %.2f +- %.2f" % (l, np.abs(mean), err)
				print "P = %.2f" % (err/np.abs(mean))

		precs.append(sum(thisgroupprecs))

	if verbose:
		for ind, prec in enumerate(precs):
			print "Total prec of %i: %f" % (ind, prec)

	return precs.index(min(precs))



def compute_sigmas(refgroup, groups, verbose=False):

	"""
	Compute the tension (in sigma units) between a reference group and the rest of the series

	@param refgroup: Reference delays to which you want to compute the tension with. Dict must be the output of getresults function.
	@param groups: List of Groups on which the tension is computed.

	@return: a list (one elt per dict in dictslist, same order) of list (one val per label, same order than each dict) of sigmas.
	"""

	sigmaslist = []
	for group in groups:
		sigmas = []
		for l in group.labels:

			# the reference estimate
			refmedian = refgroup.medians[refgroup.labels.index(l)]
			referr_up = refgroup.errors_up[refgroup.labels.index(l)]
			referr_down = refgroup.errors_down[refgroup.labels.index(l)]

			# the estimate compared to the ref
			median = group.medians[group.labels.index(l)]
			err_up = group.errors_up[group.labels.index(l)]
			err_down = group.errors_down[group.labels.index(l)]
			if median < refmedian:
				sigma = np.abs(median-refmedian)/np.sqrt(err_up**2 + referr_down**2)
			else:
				sigma = np.abs(median-refmedian)/np.sqrt(err_down**2 + referr_up**2)

			sigmas.append(sigma)
			if verbose:
				print "="*15
				print "%s: " % l,
				# reference estimate
				print "ref: %.2f + %.2f - %.2f" % (refmedian, referr_up, referr_down)

				#  comparision estimate
				print "ref: %.2f + %.2f - %.2f" % (median, err_up, err_down)
				print "tension: ", sigma

		sigmaslist.append(sigmas)
	return sigmaslist


### REWRITING DONE UNTIL HERE
def niceprint_sigmas(sigmaslist, names, refname, labelslist):
	"""
	Nice terminal printing.

	:param sigmaslist: list of sigmas, from compute_sigmas
	:param names: names of all the estimates used for sigma computation. MUST correspond to the :dictslist: you passed to compute_sigmas.
	:param refname: Name of the :refdict: you passed to compute_sigmas.
	:param labels: list of delay pairs, strings. MUST correspond to the order in which sigmas were computed in compute_sigmas. Nasty working trick, you can use [d["labels"] for d in dictslist], with dictslist being the same you used in compute sigmas. The order is important here, hence this little trick instead of using only a global variable labels.
	:return: None, just print stuff
	"""

	print "="*5+"Tension comparison"+"="*5
	print "Reference is %s" % refname
	for name, sigmas, labels in zip(names, sigmaslist, labelslist):
		print "---%s---" % name
		for sigma, label in zip(sigmas, labels):
			print "  %s: %.2f" % (label, sigma)


def convolve_estimates(results, errors, binslist, labels, testmode=True):
	"""
	Convolve the estimates in results by the corresponding error distributions in errors. This is typically used to add microlensing time-delay bias and error to a time-delay distribution.

	The two distributions need to have the same binning! binslist is the distribution on which the results are draw. I rebin the errors distributions to binslist in order to perform the convolution

	:param results: dictionnary containing one estimate. Ouput of getresults, asgetresults or any combination of these. If dict has already been linearized and thus contains the "lin" field I do nothing, otherwise I linearize it using the bins in binslist.
	:param errors: roughly similar to asgetresults: dict with elements {"labels", "means", "errors", "lin"}. If there's no lins, then I assume the errors follow Gaussian distribution. Else, I use lins that I rebin on binslist so that they matches the binning of results. WARNING: lin should have a corresponding binning "bins" as extra attribute !!
	:param binslist: np.linspace, one per delay in results

	:return: convolved estimates
	"""

	# we linarize the results if not done yet
	if not "lin" in results:
		linresults = linearize([results], binslist, testmode=testmode)[0]
		results["lin"] = linresults

	# easy case: errors are given as an analytical Gaussian distrib, I directly linearize on binslist
	if not "lin" in errors and "means" in errors and "errors" in errors:
		# I need to put the zero of the error distribution at the center of the padding
		# --> I simply shift the binslist I use so that it is centered on 0
		newbinslist = []
		for b in binslist:
			bmean = np.mean(b)
			newbins = b-bmean
			newbinslist.append(newbins)

		linerrors = linearize([errors], newbinslist, testmode=testmode)[0]

		#errors["lin"] = linerrors

	# lins is defined, but not necessarily on the same binning
	else:
		linerrors = []
		for lin in errors["lin"]:
			oldvals = lin["lin"]
			oldbins = lin["bins"]

			mybins = binslist[results["labels"].index(lin["label"])]
			newbins = mybins - np.mean(mybins) # centering the padding
			newvals = np.interp(x=newbins, xp=oldbins, fp=oldvals, left=0, right=0)
			newvals = newvals / sum(newvals)

			import matplotlib.pyplot as plt

			#plt.figure()

			#plt.plot(oldbins, oldvals)
			#plt.plot(newbins, newvals)

			#plt.show()
			#sys.exit()
			linerrors.append({"label": lin["label"], "lin": newvals[:-1]})



	# and do the convolution

	convols = []
	for rl in results["lin"]:
		label = rl["label"]
		el = linerrors[[l["label"] for l in linerrors].index(label)]
		import matplotlib.pyplot as plt
		if 0:
			plt.figure()
			plt.plot(rl["lin"])

			plt.figure()
			plt.plot(el["lin"])

			plt.figure()
			plt.plot(np.convolve(rl["lin"], el["lin"], mode="same"))
			plt.show()
		convols.append(np.convolve(rl["lin"], el["lin"], mode="same"))



	convolved = asgetresults(labels, convols, binslist, testmode=testmode)
	convolved["name"] = results["name"]+"-convolved"
	return convolved



def combine_estimates(ests, sigmathresh, binslist, labels, fancynames=None, verbose=True, testmode=True):
	"""
	Top-level function doing everything in one line (or close to...)

	Give me a list of estimates, I pick the best one by myself, compute the tension with all the others. If the tension exceed a certain threshold, I combine (i.e. sum) the best estimate with the most precise estimate that it is in tension with. I iteratively compute the tension and combine with the remaining estimates if necessary.

	:param ests: list of dicts, each containing one estimate. Dicts are outputs of getresults or asgeresults
	:param sigmathresh: threshold above which two estimates are combined
	:param binslist: binning used for the combined estimate pdf sampling
	:param labels: delays expected to be found in the dicts (in case you want to discard some...)

	:return: combined estimate from asgetresults.
	"""

	if fancynames is not None:
		assert (len(fancynames) == len(ests))
		for r, n, in zip(ests, fancynames):
			r["name"] = n

	estsnames = [r["name"] for r in ests]

	# I first compute the tension between the best estimate and the others
	bestest = ests[get_bestprec(ests)]
	sigmaslist = compute_sigmas(bestest, ests)
	combname = bestest["name"]
	if verbose:
		niceprint_sigmas(sigmaslist, estsnames, combname, [d["labels"] for d in ests])

	# If the tension exceed a threshold in a least one of the delays, I identify the corresponding estimates.
	threshfolders = []
	threshnames = []
	for ind, sigmas in enumerate(sigmaslist):
		for sigma in sigmas:
			if sigma > sigmathresh:
				threshfolders.append(ests[ind])
				threshnames.append(estsnames[ind])
				break

	i = 1  # to keep track of the iterations
	# Iterative combination of estimates in tension
	while len(threshfolders) > 0:
		# when iteratively computing the tension, the combined estimates must weight more and more since it results from the combination of an increasing number of estimates. Yet, getcombweighslist return normalized pdf. To bypass this, I add bestest to the combination list one more time at each iteration.

		tocombine = []
		for ind in np.arange(i):
			tocombine.append(bestest)
		tocombine.append(threshfolders[get_bestprec(threshfolders)])
		linearized_tocombine = []

		for est in tocombine:
			if "lin" in est:
				linearized_tocombine.append(est["lin"])
			else:
				linest = linearize([est], binslist, testmode=testmode)[0]
				linearized_tocombine.append(linest)
				est["lin"] = linest

		weightslist = getcombweightslist(linearized_tocombine, labels)
		combined = asgetresults(labels, weightslist, binslist, testmode=testmode)
		combname += '+' + threshnames[get_bestprec(threshfolders)]

		# remove the already combined estimates from results and update resultsnames
		ests = [r for r in ests if r["name"] not in (bestest["name"], threshfolders[get_bestprec(threshfolders)]["name"])]
		estsnames = [r["name"] for r in ests]

		# recompute the tension
		sigmaslist = compute_sigmas(combined, ests)
		if verbose:
			niceprint_sigmas(sigmaslist, estsnames, combname, [d["labels"] for d in ests])

		threshfolders = []
		threshnames = []
		for ind, sigmas in enumerate(sigmaslist):
			for sigma in sigmas:
				if sigma > sigmathresh:
					threshfolders.append(ests[ind])
					threshnames.append(estsnames[ind])
					break

		# Loop again if threshfolders > 0, otherwise exit the loop and return the combined estimate
		bestest = combined  # this is the new ref we want to compare the remaining estimates to.
		bestest["name"] = combname
		i += 1

	return bestest

def mult_estimates(ests, labels, binslist, testmode=True, verbose=True):
	"""
	Give me a list of estimates, I return the combined distribution as if they were independent measurements of the same variable. I basically just multiply the pdfs

	:param ests:
	:return:
	"""
	linearized = []
	for est in ests:
		if "lin" in est:
			linearized.append(est["lin"])
		else:
			linest = linearize([est], binslist, testmode=testmode, verbose=verbose)[0]
			linearized.append(linest)
			est["lin"] = linest

	weightslist = []
	for ind, label in enumerate(labels):
		labellins = [[l for l in lins if l["label"]==label][0] for lins in linearized]
		indweights = [lin["lin"] for lin in labellins]
		multweights = indweights[0]
		for e in indweights[1:]:
			multweights = multweights*e
		weightslist.append(multweights)

	comb = asgetresults(labels, weightslist, binslist, testmode=testmode)
	comb["name"] = "-x-".join([est["name"] for est in ests])
	return comb


