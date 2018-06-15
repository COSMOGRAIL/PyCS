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


def confinterval(xs, weights=None, conflevel=0.68, cumulative=False, testmode=True):
	"""
	Hand-made quantile/percentile computation tool

	:param xs: list, distribution
	:param weights: list, weights of the distribution
	:param conflevel: float, quantile (between 0 and 1) of the conf interval
	:param cumulative: bool. If set to True, return the "bound" (max) value of xs corresponding to the conflevel. For example, if conflevel=0.68, xs=nnu, return the value such that nnu<=val at 68% confidence.
	:param testmode: bool. Defines the finesse of the binning of the histogram used when computing the percentiles. A finer binning gives more precise results but takes longer to process. Setting it to True sets 500 bins, otherwise 10'000.

	:return: return 3 floats: the median value, minimal value and maximal value matching your confidence level
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


def binlist(list, bins):
	"""
	Spread a list of values on a given binning and return the mean value of each bin, normalized to 1.

	:param list: list of values
	:param bins: list of bins on which you want to spread your values
	:return: mean value of the bins
	"""
	digitized = np.digitize(list, bins)
	bin_means = [float(len(list[digitized == i])) for i in range(1, len(bins))]
	bin_means = np.nan_to_num(bin_means)
	bin_means = bin_means/np.sum(bin_means)
	return bin_means


def getresults(allfolders, folder, dataname, methodkw, labels):
	"""
	Give me a list of folders in which your PyCS runresults are. I extract the measured time delay that correspond to a given folder and method key work you give me as argument.

	.. warning:: This works **only** if you already collected the runresults in delaycontainer, using e.g. the :py:func:`pycs.sim.plot.hists` and :py:func:`pycs.sim.plot.measvstrue` functions.

	:param allfolders: all possible folders that you want me to search into. Could be *
	:param folder: specific folder you want to get the results from. Name related to the way the mocks were drawn. E.g. spl1_ks20_splml
	:param dataname: name of the dataset. E.g. WFI-5
	:param methodkw: name related to the way the mocks were analysed. E.g. spl1t3
	:param labels: list of labels (e.g. "AB") corresponding to the delay pairs you want to collect.

	:return: dict of lists with corresponding indexes, containing the labels, means and errors of the delays.
	"""
	os.chdir(allfolders[allfolders.index(folder)])
	copies = [line.split('.pkl')[0][12:] for line in os.popen('ls -d *%s*' % dataname).readlines() if 'delays.pkl' in line]
	mocks = [line.split('.pkl')[0][11:] for line in os.popen('ls -d *%s*' % dataname).readlines() if 'errorbars' in line]

	mylabels, means, errors = [], [], []

	for copie, mock in zip(copies, mocks):
		if dataname in copie and dataname in mock and methodkw in copie and methodkw in mock:
			result = (pycs.gen.util.readpickle('sims_copies_%s.pkl' % copie), pycs.gen.util.readpickle('sims_mocks_%s.pkl' % mock))

			# making it robust against potentially mixed dicts...
			mylabels = [d["label"] for d in result[0].data]
			means = []
			errors = []
			for label in labels:
				means.append([d["mean"] for d in result[0].data if d["label"] == label][0])
				errors.append([d["tot"] for d in result[1].data if d["label"] == label][0])

	os.chdir("../")
	return {"labels": mylabels, "means": means, "errors":errors}



def linearize(resultslist, binslist, verbose=True, testmode=True):
	"""
	Give me a list of results (each one being a dict with means and errors values for various delays) and a list of corresponding bins. I draw samples from a Normal distribution according to your mean and error (1sigma) and spread them on the given bins.

	Simply said, I create Gaussian distributions on the bins of your choice.

	.. warning:: Of course, if your theoretical values are not expected to follow a Gaussian distribution, calling this function makes no sense.

	:param resultslist: list of dict, each dict is the output of getresults. Each dict should contain the fields "labels", "means" and "errors"
	:param binslist: np.linspace, one per delay in resultlist
	:param testmode: bool. Defines the number of samples drawn from the Normal distribution. More samples give more precise results but take longer to create and process. Setting testmode to True create 50k samples, otherwise 500k.

	"""

	if testmode:
		nsamples = 50000
	else:
		nsamples = 500000

	linslist = []
	for results in resultslist:
		lins = []
		for ind, label in enumerate(results["labels"]):
			xs = np.random.normal(loc=results["means"][ind], scale=results["errors"][ind], size=nsamples)
			xs = binlist(xs, binslist[ind])

			# To check we are not too far off
			ci = confinterval(binslist[ind][:-1], weights=xs, testmode=testmode)
			if verbose:
				print "="*45
				print "Delay %s:" % label, "%.2f" % results["means"][ind], "+/-", "%.2f" % results["errors"][ind] ," --> ", "%.2f" % ci[0], "+%.2f-%.2f" % (ci[2], ci[1])
			lins.append({"lin": xs, "label": label})
		linslist.append(lins)
	return linslist




def getcombweightslist(linslist, labels):
	"""
	Give me 2+ distributions defined on the same binning, and I return a linear combination (a sum, if you prefer) of them.


	:param linslist: list of linearized distributions you want to combine (linear combination, no choice)
	:param labels: labels of the distributions you want to combine
	:return: weights of the combined distribution
	"""
	weightslist = []
	for ind, label in enumerate(labels):
		labellins = [[l for l in lins if l["label"]==label][0] for lins in linslist]

		indweights = [lin["lin"] for lin in labellins]

		"""
		This is trickily written so let me comment
		indweights contains a list of various linearized corresponding normalised distributions binned the same way (i.e. AB pdf for spl20 and spl30). The values of distrib are added bin by bin, then reweighted so that the sum of the result is still normalised to 1 (so they can be combined later with other results from various techniques, instruments, etc...) Since we work with numpy arrays, this can be done in a single line.
		"""
		weights = sum(indweights)/float(len(indweights))
		weightslist.append(weights)
	return weightslist



def asgetresults(labels, weightslist, binslist, testmode=True):
	"""
	Transform the output of confinterval function into a dict, with the same field than what getresults would give.

	Small exception, since confinterval usually return asymmetric distribution, I save these info errors_up and errors_down dicts. I also save the linearized distribution of this estimate that actually corresponds to the input weightslist.
	"""

	means = []
	errors_up = []
	errors_down = []
	lins = []
	for l, w, b in zip(labels, weightslist, binslist):
		ci = confinterval(b[:-1], w, testmode=testmode)
		means.append(ci[0])
		errors_up.append(ci[2])
		errors_down.append(ci[1])
		lins.append({"lin": w, "label": l})

	return {"labels": labels, "means": means, "errors_up": errors_up, "errors_down": errors_down, "lin": lins}


def niceprint_values(name, mydict):
	"""
	:param name: string, give a name to your results for the display
	:param mydict: output of getresults or asgetresults
	:return: print nicely in the terminal the delays of your dict
	"""

	print "="*5, name, "="*5
	for ind, l in enumerate(mydict["labels"]):

		if "errors_up" in mydict and "errors_down" in mydict:
			print "%s:" % l, "%.2f +%.2f-%.2f" % (mydict["means"][ind], mydict["errors_up"][ind], mydict["errors_down"][ind])
		else:
			print "%s:" % l, "%.2f +-%.2f" % (mydict["means"][ind], mydict["errors"][ind])

def get_bestprec(dictslist, refimg=None, verbose=False):
	"""
	Return the most precise estimate among a list of estimates.

	The best estimate is computed in a straightforward, naive way: simply sum the relative precision over all the pair of delays that contain refimg in their label (or all, if no refimg given), then take the estimate that has the lowest sum.

	To avoid giving too much weight to the larger delays, we compute the precision relative to the "mean" delay from all estimates.

	:param dictslist: list of dic, output of getresult or ...
	:param refimg: reference image, if you want to discard some delays from the process

	:return: index of the most precise estimate in the list
	"""

	# compute average means:
	avgmeans = []
	for l in dictslist[0]["labels"]:
		means = [dict["means"][dict["labels"].index(l)] for dict in dictslist]
		avgmeans.append(np.mean(means))


	precs = []
	for dict in dictslist:
		thisdictprecs = []
		if refimg==None:
			for l in dict["labels"]:
				mean = avgmeans[dictslist[0]["labels"].index(l)]
				# if asymmetric distribution, we have errors_up and errors down instead of errors only
				if "errors_up" in dict and "errors_down" in dict:
					# for precision, we simply average over err_up and err_down
					err = (dict["errors_up"][dict["labels"].index(l)] + dict["errors_down"][dict["labels"].index(l)])/2.0
				else:
					err = dict["errors"][dict["labels"].index(l)]
				thisdictprecs.append(err/np.abs(mean))
				if verbose:
					print "="*5
					print "%s: %.2f +- %.2f" % (l, np.abs(mean), err)
					print "P = %.2f" % (err/np.abs(mean))
		else:
			for l in [l for l in dict["labels"] if refimg in l]:
				mean = avgmeans[dictslist[0]["labels"].index(l)]
				if "errors_up" in dict and "errors_down" in dict:
					err = (dict["errors_up"][dict["labels"].index(l)] + dict["errors_down"][dict["labels"].index(l)])/2.0
				else:
					err = dict["errors"][dict["labels"].index(l)]
				thisdictprecs.append(err/np.abs(mean))
				if verbose:
					print "="*5
					print "%s: %.2f +- %.2f" % (l, np.abs(mean), err)
					print "P = %.2f" % (err/np.abs(mean))

		precs.append(sum(thisdictprecs))

	if verbose:
		for ind, prec in enumerate(precs):
			print "Total prec of %i: %f" % (ind, prec)

	return precs.index(min(precs))


def compute_sigmas(refdict, dictslist, verbose=False):
	"""
	Compute the tension (in sigma units) between a reference group and the rest of the series


	:param refdict: Reference delays to which you want to compute the tension with. Dict must be the output of getresults function.
	:param dicslist: List of dicts on which the tension is computed. Each dict must be the output of getresults function
	:return: a list(one elt per dict in dictslist, same order) of list(one val per label, same order than each dict) of sigmas.
	"""
	sigmaslist = []
	for dict in dictslist:
		sigmas = []
		for l in dict["labels"]:

			# the reference estimate
			refmean = refdict["means"][refdict["labels"].index(l)]
			if "errors_up" in refdict and "errors_down" in refdict:  # this is an asymmetric pdf
				referr_up = refdict["errors_up"][refdict["labels"].index(l)]
				referr_down = refdict["errors_down"][refdict["labels"].index(l)]
			else:
				referr = refdict["errors"][refdict["labels"].index(l)]
				# "copy referr into referr_up and down to simplify the conditional expressions below"
				referr_up = referr
				referr_down = referr

			# the estimate compared to the ref
			mean = dict["means"][dict["labels"].index(l)]
			if "errors_up" in dict and "errors_down" in dict:  # this is an asymmetric pdf
				err_up = dict["errors_up"][dict["labels"].index(l)]
				err_down = dict["errors_down"][dict["labels"].index(l)]
				if mean < refmean :
					sigma = np.abs(mean-refmean)/np.sqrt(err_up**2 + referr_down**2)
				else:
					sigma = np.abs(mean-refmean)/np.sqrt(err_down**2 + referr_up**2)


			else:  # then it's a Gaussian estimate
				err = dict["errors"][dict["labels"].index(l)]
				sigma = np.abs(mean-refmean)/np.sqrt(err**2 + err**2)
			sigmas.append(sigma)
			if verbose:
				print "="*15
				print "%s: " % l,
				# reference estimate
				if "errors_up" in refdict and "errors_down" in refdict:
					print "ref: %.2f + %.2f - %.2f" % (refmean, referr_up, referr_down)
				else:
					print "ref: %.2f +- %.2f" % (refmean, referr)

				#  comparision estimate
				if "errors_up" in dict and "errors_down" in dict:
					print "ref: %.2f + %.2f - %.2f" % (mean, err_up, err_down)
				else:
					print "comp: %.2f +- %.2f" % (mean, err)
				print "tension: ", sigma

		sigmaslist.append(sigmas)
	return sigmaslist


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


