import numpy as np
import pycs

def combiconf1(estimates):
	"""
	Give me a group of estimates of the same TDC pair, I return the combined confidence.
	1 is for doubtless
	2 for plausible
	3 for multimodal
	4 for uninformative

	Simple version

	@param estimates: list of Estimate objects

	@return: dictionary that contains the combined confidence code and the estimate set, rung and pair.
	"""

	(tdcset, rung, pair) = pycs.tdc.est.checkallsame(estimates)
	idconfs = [(est.id,est.confidence) for est in estimates]
	
	
	### now we define criteria by criteria
	
	# doubtless, conflevels = 1 only
	if all(idconf[1]==1 for idconf in idconfs):
		combiconfcode = 1
	
	# plausible, conflevels = 1 or 2
	elif all(idconf[1]<=2 for idconf in idconfs):
		combiconfcode = 2
	
	# doubless+multimodal, need a manual check
	elif all(idconf[1]<=3 for idconf in idconfs) and 1 in (idconf[1] for idconf in idconfs):
		combiconfcode = 3
		
	# plausible+multimodal, no doubtless
	elif all(idconf[1]<=3 for idconf in idconfs):
		combiconfcode = 4
		
	# uninformative+doubless/plausible, need a manual check
	elif 1 in (idconf[1] for idconf in idconfs) or 2 in (idconf[1] for idconf in idconfs):
		combiconfcode = 5
	
	# uninformative and/or multimodal, to shoot
	elif 3 in (idconf[1] for idconf in idconfs) or 4 in (idconf[1] for idconf in idconfs):
		combiconfcode = 6
	
	else:
		raise RuntimeError("I didn't catch %s, should not happen !" % est.id)

	return {"code":combiconfcode, "set":tdcset, "rung":rung, "pair":pair}
	
	
def combiconf2(estimates):
	"""

	Give me a group of estimates of the same TDC pair, I return the combined confidence.
	1 is for doubtless
	2 for plausible
	3 for multimodal
	4 for uninformative

	Updated version of combiconf1, now with accordance criterias between the estimates, and weight given to Vivien/Malte estimates

	@param estimates: list of Estimate objects

	@return: dictionary that contains the combined confidence code and the estimate set, rung and pair.

	"""	
	(tdcset, rung, pair) = pycs.tdc.est.checkallsame(estimates)

	combiconfcode = 0 # the output value.
	
	idconfs = [(est.id,est.confidence) for est in estimates]	
	
	hugefact = 0.4
	sig = 2.0
	
	### First step, define "advanced" criterias
	#
	# Idea:	see if estimates agree (they stay in each other's errorbars)
	#	 	see if the errorbars are not too big compared to the delay
	#	if not, if Malte has an estimate, see if it agrees with Vivien
	#		if not, see if Vivien has a high confidence level (doubt/plaus)		  

	
	# doubtless, conflevels = 1 only.
	if all(idconf[1]==1 for idconf in idconfs):
		for est in estimates:
			if combiconfcode == 999:
				break
			subests = [subest for subest in estimates if subest.methodpar != est.methodpar]
			for subest in subests:
				if est.td > subest.td-sig*subest.tderr and est.td < subest.td+sig*subest.tderr:
					if min([subest.tderr/(subest.td+0.01) for subest in subests]) < hugefact and combiconfcode != 11:
						combiconfcode = 10 # then, estimates agree
					else:
						combiconfcode = 11 # they agree, but the errorbars are "huge"
				else:
					combiconfcode = 999
					break														
								
		if combiconfcode == 999: 	
			if "mtewes" in (est.methodpar for est in estimates): 
				mest = [est for est in estimates if est.methodpar == 'mtewes'][0]
				vest = [est for est in estimates if est.methodpar == 'Vivien'][0]
				if vest.td > mest.td-sig*mest.tderr and vest.td < mest.td+sig*mest.tderr and mest.td > vest.td-sig*vest.tderr and mest.td < vest.td+sig*vest.tderr:
					combiconfcode = 12 # Malte and Vivien agree, someone else disagree
				else:
					combiconfcode = 13 # Malte and Vivien disagree		
			else:
				combiconfcode = 14 # No Malte, only Vivien, disagree with someone else	


	# plausible, conflevels = 1 or 2. 
	elif all(idconf[1]<=2 for idconf in idconfs):
		for est in estimates:
			if combiconfcode == 999:
				break
			subests = [subest for subest in estimates if subest.methodpar != est.methodpar]
			for subest in subests:
				if est.td > subest.td-sig*subest.tderr and est.td < subest.td+sig*subest.tderr:
					if min([subest.tderr/(subest.td+0.01) for subest in subests]) < hugefact and combiconfcode != 21:
						combiconfcode = 20 # then, estimates agree
					else:
						combiconfcode = 21 # they agree, but the errorbars are "huge"

				else:
					combiconfcode = 999
					break
									
		if combiconfcode == 999: 
			if "mtewes" in (est.methodpar for est in estimates): 
				mest = [est for est in estimates if est.methodpar == 'mtewes'][0]
				vest = [est for est in estimates if est.methodpar == 'Vivien'][0]
				if vest.td > mest.td-sig*mest.tderr and vest.td < mest.td+sig*mest.tderr and mest.td > vest.td-sig*vest.tderr and mest.td < vest.td+sig*vest.tderr:
					combiconfcode = 22 # Malte and Vivien agree, someone else disagree
				else:
					combiconfcode = 23 # Malte and Vivien disagree		
			else:
				combiconfcode = 24 # No Malte, only Vivien, disagree with someone else							


	# doubless+multimodal
	elif all(idconf[1]<=3 for idconf in idconfs) and 1 in (idconf[1] for idconf in idconfs):
		for est in estimates:
			if combiconfcode == 999:
				break
			subests = [subest for subest in estimates if subest.methodpar != est.methodpar]
			for subest in subests:
				if est.td > subest.td-sig*subest.tderr and est.td < subest.td+sig*subest.tderr:
					if min([subest.tderr/(subest.td+0.01) for subest in subests]) < hugefact and combiconfcode != 31:
						combiconfcode = 30 # then, estimates agree
					else:
						combiconfcode = 31 # they agree, but the errorbars are "huge"
				else:
					combiconfcode = 999
					break														
						
		if combiconfcode == 999:
		 
			if "mtewes" in (est.methodpar for est in estimates): 
				mest = [est for est in estimates if est.methodpar == 'mtewes'][0]
				vest = [est for est in estimates if est.methodpar == 'Vivien'][0]
				if vest.td > mest.td-sig*mest.tderr and vest.td < mest.td+sig*mest.tderr and mest.td > vest.td-sig*vest.tderr and mest.td < vest.td+sig*vest.tderr:
					combiconfcode = 32 # Malte and Vivien agree, someone else disagree
				else:
					combiconfcode = 33 # Malte and Vivien disagree		
			else:
				combiconfcode = 34 # No Malte, only Vivien, disagree with someone else	
			
										
	# plausible+multimodal, no doubtless
	elif all(idconf[1]<=3 for idconf in idconfs):
		for est in estimates:
			if combiconfcode == 999:
				break
			subests = [subest for subest in estimates if subest.methodpar != est.methodpar]
			for subest in subests:
				if est.td > subest.td-sig*subest.tderr and est.td < subest.td+sig*subest.tderr:
					if min([subest.tderr/(subest.td+0.01) for subest in subests]) < hugefact and combiconfcode != 41:
						combiconfcode = 40 # then, estimates agree
					else:
						combiconfcode = 41 # they agree, but the errorbars are "huge"
				else:
					combiconfcode = 999
					break														
						
		if combiconfcode == 999: 
			if "mtewes" in (est.methodpar for est in estimates): 
				mest = [est for est in estimates if est.methodpar == 'mtewes'][0]
				vest = [est for est in estimates if est.methodpar == 'Vivien'][0]
				if vest.td > mest.td-sig*mest.tderr and vest.td < mest.td+sig*mest.tderr and mest.td > vest.td-sig*vest.tderr and mest.td < vest.td+sig*vest.tderr:
					combiconfcode = 42 # Malte and Vivien agree, someone else disagree
				else:
					combiconfcode = 43 # Malte and Vivien disagree		
			else:
				combiconfcode = 44 # No Malte, only Vivien, disagree with someone else

				
	# uninformative+doubless/plausible
	elif 1 in (idconf[1] for idconf in idconfs) or 2 in (idconf[1] for idconf in idconfs):
		for est in estimates:
			if combiconfcode == 999:
				break
			subests = [subest for subest in estimates if subest.methodpar != est.methodpar]
			for subest in subests:
				if est.td > subest.td-sig*subest.tderr and est.td < subest.td+sig*subest.tderr:
					if min([subest.tderr/(subest.td+0.01) for subest in subests]) < hugefact and combiconfcode != 51:
						combiconfcode = 50 # then, estimates agree
					else:
						combiconfcode = 51 # they agree, but the errorbars are "huge"
				else:
					combiconfcode = 999
					break														
						
		if combiconfcode == 999:
		 	vest = [est for est in estimates if est.methodpar == 'Vivien'][0]
			if "mtewes" in (est.methodpar for est in estimates): 
				mest = [est for est in estimates if est.methodpar == 'mtewes'][0]
				if vest.td > mest.td-sig*mest.tderr and vest.td < mest.td+sig*mest.tderr and mest.td > vest.td-sig*vest.tderr and mest.td < vest.td+sig*vest.tderr:
					if mest.confidence < 4 or vest.confidence < 4: 
						combiconfcode = 52 # Malte or Vivien tagged plausible/doubtless and both agree
					else:
						combiconfcode = 53 # Malte and Vivien tagged uninformative (but agreed...) 
				else:
					if mest.confidence < 4 or vest.confidence < 4:
						combiconfcode = 54 # Malte or Vivien tagged plausible/doubtless, but both disagree 
					else:
						combiconfcode = 55 # Malte and Vivien tagged uninformative and disagree				
				
				
										
			else:
				if vest.confidence < 4:
					combiconfcode = 56 # ests disagree, but Vivien flagged plausible/doubtless.
				else:
					combiconfcode = 57 # ests disagree, Vivien flagged uninformative						

					
	# uninformative and/or multimodal, to shoot
	elif 3 in (idconf[1] for idconf in idconfs) or 4 in (idconf[1] for idconf in idconfs):
		combiconfcode = 60					
										
	else:
		raise RuntimeError("WARNING -- I didn't catch %s, should not happen !" % est.id)

	if combiconfcode == 0:
		raise RuntimeError("WARNING -- combiconfcode for %s is 0, should not happen" % est.id)
		
	return {"code":combiconfcode, "set":tdcset, "rung":rung, "pair":pair}					
					

					
def reroll(estimates):
	"""
	Give me a list of estimates. I remove all estimates not from Malte or Vivien, then recompute the combiconfcode
	Return the new combiconfcode, along with the id of the estimate

	@param estimates: list of Estimate objects

	@return: dictionnary containing the new combiconfcode computed by combiconf2, the estimate set, rung and pair and the corresponding estimates from Malte and Vivien only
	"""					
					
	new_estimates = [est for est in estimates if est.methodpar == "Vivien" or est.methodpar == "mtewes"]
	
	# we return the same parameters as combiconf2, and the new estimates that are to be taken into account
	outcombiconf2 = combiconf2(new_estimates)				
	return {"code":outcombiconf2["code"], "set":outcombiconf2["set"], "rung":outcombiconf2["rung"], "pair":outcombiconf2["pair"], "ests":new_estimates}				
					
													

def combiquads(estimates):
	"""

	This function is based on a exploit discovered during the rungs generation: the quads with the same pair indices are the same amongst all rungs, thus have the same delay.

	Not used for TDC1 official submissions, just for fun.

	Give me a list of already combined estimates (i.e. one per pair only).
	I compute for each quad pair the best td and tderr among the corresponding rungs (which is an exploit !)
	I return the modified estimates list, with all the quads modified, and with only the doubtless double
	
	# BIG WARNING !!! This function (pycs.tdc.combiconf.combiquads) actually modify estimates, we DO NOT want that !

	@param estimates: list of Estimate objects. Each id must be present no more than once.

	@return: list of modified Estimate objects.
	"""

	print 'BIG WARNING !!! This function (pycs.tdc.combiconf.combiquads) actually modify estimates, we DO NOT want that ! '
	
	# Check unicity
	pycs.tdc.est.checkunique(estimates)
	
	# Split double and quads
	doubleests = pycs.tdc.est.select(estimates, pairs = np.arange(1037)[:720])
	doubleests = [est for est in doubleests if est.confidence == 1] # we keep only the doubless double !
	quadests = pycs.tdc.est.select(estimates, pairs = np.arange(1037)[720:])
	
	
	# Now, give all the quadests a fake id, to ease their combination
	for est in quadests:
		est.method = 'exploit'
		est.methodpar = 'rung_%s' % est.rung
		est.originalid = est.id # brutally added field, for later recomposition of the full sample
		est.id = 'tdc1_exploit_%i' % est.pair
	
	
	groupests = pycs.tdc.est.group(quadests) # group function group the estimates according to their id...which explain the id tweak above
	
	
	# Compute the mean td and tderr according to confidence criterias from the wiki page (here, a simple version for testing)
	for groupest in groupests:
		
		# two or more estimates have a doubtless/plausible confidence level
		if sum(est.confidence == 1 or est.confidence == 2 for est in groupest)	>= 2:
			td = np.median([est.td for est in groupest if est.confidence in [1,2]])
			tderr = np.median([est.tderr for est in groupest if est.confidence in [1,2]])/np.sqrt(sum(est.confidence == 1 or est.confidence == 2 for est in groupest))
			confidence = 1
	
		# one estimate at least have a doubtless or plausible confidence level
		elif 1 in (est.confidence for est in groupest) or 2 in (est.confidence for est in groupest):
			td = np.median([est.td for est in groupest if est.confidence in [1,2]]) 
			tderr = np.median([est.tderr for est in groupest if est.confidence in [1,2]])/np.sqrt(sum(est.confidence == 1 or est.confidence == 2 for est in groupest))
			confidence = 2
			
		# at least one multimodal...
		elif 3 in (est.confidence for est in groupest):
			td = np.median([est.td for est in groupest if est.confidence in [3]])
			tderr = np.median([est.tderr for est in groupest if est.confidence in [3]])/np.sqrt(sum(est.confidence == 3 for est in groupest))
			confidence = 3			
					
		elif 4 in (est.confidence for est in groupest):
			td = 99.00
			tderr = -99.00
			confidence = 4	
	
		else:
			print est.confidence
			raise RuntimeError("WARNING -- I didn't catch %s, should not happen !" % est.originalid)
			

		for est in groupest:
			est.td = td
			est.tderr = tderr
			est.confidence = confidence
			est.id = est.originalid # get back to normal...ouf !

	# Rebuild the list of estimates to return
	
	exploitests = []
	
	for est in doubleests:
		exploitests.append(est)
	for groupest in groupests:
		for est in groupest:
			exploitests.append(est)	
	
	return exploitests

