import numpy as np
import pycs

def combiconf1(estimates):
	"""
	Give me a group of ests of the same pair, I return the combined confidence.
	"""
	
	(tdcset, rung, pair) = pycs.tdc.est.checkallsame(estimates)
	
	identity = estimates[0].id
	combiconfcode = 0 # the output value.
	
	idconfs = [(est.id,est.confidence) for est in estimates]
	
	
	### now we define criteria by criteria
	
	# doubtless, conflevels = 0 only
	if all(idconf[1]==1 for idconf in idconfs):
		combiconfcode = 1
	
	# plausible, conflevels = 0 or 1
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
		raise RuntimeError("I didn't catch this one, should not happen !")

	return {"code":combiconfcode, "set":tdcset, "rung":rung, "pair":pair}
	
	
	
