import sys
sys.path.append("../")			
from pycs.gen import *
from pycs.pelt import *
from numpy import *

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm


lcdict = util.readpickle("lcs+ml3.pkl")
	
lca = lcdict['A']
lcb = lcdict['B']
lcc = lcdict['C']
lcd = lcdict['D']
	
lcs = [lca, lcb, lcc, lcd]


#rawdispersionmethod = lambda lc1, lc2 : disps.linintnp(lc1, lc2, interpdist = 30.0)
#dispersionmethod = lambda lc1, lc2 : disps.symmetrize(lc1, lc2, rawdispersionmethod)


#lca.shifttime(8.4)
#lcc.shifttime(7.8)
#lcd.shifttime(-6.5)


#lc.display(lcs)

for l in lcs:
	if l.ml != None:
	
		plt.plot(l.getjds(), l.ml.calcmlmags(l), ".", color = l.plotcolour, label=str(l))
	
		stats = l.ml.stats(l)
		plt.axhline(stats["mean"], color=l.plotcolour)
		plt.axhspan(stats["mean"]-stats["std"], stats["mean"]+stats["std"], facecolor=(0.8, 0.8, 0.8))
		
		plt.text(mean(l.getjds()), stats["mean"], "%s : %5.2f +/- %5.2f" % (l.object, stats["mean"], stats["std"]), color="black", fontsize=20)
		
		stats["std"] = sqrt(stats["std"]*stats["std"] + 0.1*0.1)
		
		
		mup = stats["mean"] - stats["std"]
		mmean = stats["mean"]
		mlow = stats["mean"] + stats["std"]
		
		fup = 10**(-0.4 * mup)
		fmean = 10**(-0.4 * mmean)
		flow = 10**(-0.4 * mlow)
		
		fstd = ( (fup-fmean) + (fmean-flow) ) / 2.0
		
		#print fup, fmean, flow
		
		print "%s (ratio) : %5.3f +/- %5.3f" % (l.object, fmean, fstd)
		print "%s (mag)   : %5.3f +/- %5.3f" % (l.object, stats["mean"], stats["std"])

#axes = plt.gca()
#axes.set_ylim(axes.get_ylim()[::-1])

plt.legend(loc = "best", numpoints = 1, prop = fm.FontProperties(size = 10))

plt.xlabel("HJD - 2400000.5 [days]", fontsize=14)
plt.ylabel("Microlensing magnitude", fontsize=14)

plt.show()