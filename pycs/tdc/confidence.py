import pycs
import os,sys
import numpy as np
import matplotlib.pyplot as plt

"""
Determine which estimates is in which confidence level --- add this to pycs once it is cleaned.
"""


# Start with import of db...

iniests = pycs.tdc.est.importfromd3cs("web/d3cslog.txt")



"""
confidence:

doubtless : 2x single doubtless estimates where meanestimates +- 2sigma includes the two estimates
plausible : 1x doubtless and 1x plausible estimates that matches

These two makes 67% of the estimates...
... see other criterias in the code below


"""

groupests = pycs.tdc.est.group(iniests)

#pycs.gen.util.writepickle(groupests,'groupests.pkl')



basedir = '' # Where do you want the .pkl files grouping the estimates to be written...S 

#groupests = pycs.gen.util.readpickle('groupests.pkl')

selections = ['doubtless','plausible','doubtless_to_multimodal','multimodal','doubtless_to_uninformative','uninformative']


idconfuniques = []
doubtlesses = []
plausibles = []
doubtless_to_multimodals = []e
multimodals = []
doubtless_to_uninformatives = []
uninformatives = []

for groupest in groupests:
	idconfs = [(est.id,est.confidence) for est in groupest]
	identity = idconfs[0][0]
	if identity == 'tdc1_1_341': # skip that one, as it bug in the code (don't know why yet...)
		print "WARNING !!!  I skip %s" % identity 
		continue
	
	### now we define criteria by criteria
	
	# doubtless, conflevels = 0 only
	if all(idconf[1]==1 for idconf in idconfs):
		selection=selections[0]
		idconfuniques.append((identity,selection))
		doubtlesses.append(identity)
	
	# plausible, conflevels = 0 or 1
	elif all(idconf[1]<=2 for idconf in idconfs):
		selection=selections[1]
		idconfuniques.append((identity,selection))
		plausibles.append(identity)
	
	# doubless+multimodal, need a manual check
	elif all(idconf[1]<=3 for idconf in idconfs) and 1 in (idconf[1] for idconf in idconfs):
		selection=selections[2]	
		idconfuniques.append((identity,selection))
		doubtless_to_multimodals.append(identity)
		
	# plausible+multimodal, no doubtless
	elif all(idconf[1]<=3 for idconf in idconfs):
		selection=selections[3]
		idconfuniques.append((identity,selection))
		multimodals.append(identity)
		
	# uninformative+doubless/plausible, need a manual check
	elif 1 in (idconf[1] for idconf in idconfs) or 2 in (idconf[1] for idconf in idconfs):
		selection=selections[4]
		idconfuniques.append((identity,selection))
		doubtless_to_uninformatives.append(identity)
	
	# uninformative and/or multimodal, to shoot
	elif 3 in (idconf[1] for idconf in idconfs) or 4 in (idconf[1] for idconf in idconfs):
		selection=selections[5]			
		idconfuniques.append((identity,selection))
		uninformatives.append(identity)

# Check that we got all the estimates...

doubtless=[idconfunique[1] for idconfunique in idconfuniques].count(selections[0])
plausible=[idconfunique[1] for idconfunique in idconfuniques].count(selections[1])
doubtless_to_multimodal=[idconfunique[1] for idconfunique in idconfuniques].count(selections[2])
multimodal=[idconfunique[1] for idconfunique in idconfuniques].count(selections[3])
doubtless_to_uninformative=[idconfunique[1] for idconfunique in idconfuniques].count(selections[4])
uninformative=[idconfunique[1] for idconfunique in idconfuniques].count(selections[5])
		
if len(groupests)-1 != len(idconfuniques):
	print 'WARNING !!! -- confidence criterias do not define all the estimates !!!'
	sys.exit()


# write pickles

if not os.path.isdir(basedir):
	os.mkdir(basedir)

print 'removing old pkls...'	
os.system('rm %s/*' %basedir)

pycs.gen.util.writepickle(doubtlesses,os.path.join(basedir,'doubtlesses.pkl'))
pycs.gen.util.writepickle(plausibles,os.path.join(basedir,'plausibles.pkl'))
pycs.gen.util.writepickle(doubtless_to_multimodals,os.path.join(basedir,'doubtless_to_multimodals.pkl'))
pycs.gen.util.writepickle(multimodals,os.path.join(basedir,'multimodals.pkl'))
pycs.gen.util.writepickle(doubtless_to_uninformatives,os.path.join(basedir,'doubtless_to_uninformatives.pkl'))
pycs.gen.util.writepickle(uninformatives,os.path.join(basedir,'uninformatives.pkl'))

print '---------------'


sys.exit()


# Now, plot these values

plt.figure('counts')
counts = [doubtless,plausible,doubtless_to_multimodal,multimodal,doubtless_to_uninformative,uninformative]
percents = [count/5180.0*100 for count in counts]
labels=[]
for count,name in zip(percents,selections):
	labels.append(name+': '+'%.0f' %count +'%')
 
y_pos = np.arange(len(labels))
plt.barh(y_pos, counts, align='center', alpha=0.4)
plt.yticks(y_pos, labels)
plt.xlabel('Number of pairs')



plt.show()
