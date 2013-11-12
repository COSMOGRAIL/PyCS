

def boot(lcs, rawdispersionmethod, n=1000, saven=100, filename="multiopt_boot.pkl", bootmode="mags", optmode="full", jdamplitude=1.0, tdrandom=20.0):
	"""
	
	abs : we save the "absolute" timeshifts of your n lcs into the pickle, not just the delays.
	(note that *nothing* else changes, no impact in the actual curve shifting !)
	
	First big bootstrapping to produce time delay histogram.
	Warning we change the lcs here... Even if they is no reason as we could make a copy.
	
	Think about how to handle a first optimization of the non-bootstrapped curve,
	then always work on copies of this guy ...
	
	initrandom is a range of days in which the initial time delays of each bootstrap are made to vary
	(around the initial time delays you give to the lightcurves before starting this bootstrapping.)
	You can put 0.0 and we will always start from the same location.
	
	Test if setting the delays to a random position does change something ...
	tdrandom = 0.0 suppresses time delay randomization, and we allways start the optimization from the same initial delays.
	
	bootmode = "mags" : only mags are bootstrapped
	bootmode = "jds"  : only jds are bootstrapped
	bootmode = "both" : both are boostrapped, simply one after the other
	bootmode = "none" : we do not bootstrap at all. Useful to test optimizations etc.
	
	optmode = "full"  : both
	optmode = "ml"    : only ml
	optmode = "ts"    : only time shifts
	optmode = "none"  : no opt at all, just return the current value of the dispersion.
	"""
	
	# This was a nice idea, but now that we use also jdbootstrapping, let's just make a plain copy 
	# of the full lightcurves ...
	#def backupmags(lcs):
	#	return [l.mags.copy() for l in lcs]
	#
	#def restoremags(lcs, bumg):
	#	for i in range(len(lcs)):
	#		lcs[i].mags = bumg[i].copy()
	#bumg = backupmags(lcs)
	
	
	initdelays = lc.multigettimedelays(lcs)
	
	tslist = []
	d2list = []
	
	def savefile():	
		tss =  np.column_stack(tslist)
		d2s = np.array(d2list)
		util.writepickle({"tss":tss, "d2s":d2s, "lcs":lcs, "n":b+1, "filename":filename,
		"bootmode":bootmode, "optmode":optmode, "tdrandom":tdrandom, "jdamplitude":jdamplitude}, filename)
		
	
	for b in range(n):
			
		#lc.display(lcs)
		print "Bootstrap %4i / %4i" % (b+1, n)
		
		lcsc = [l.copy() for l in lcs]
		
		if bootmode == "mags" or bootmode == "both":
			for l in lcsc:
				l.montecarlomags()
				#l.pseudobootstrap()
		if bootmode == "jds" or bootmode == "both":
			for l in lcsc:
				l.montecarlojds(amplitude = jdamplitude, keepml=True)
		
		#lc.display(lcs)
		# We set the time delays to random positions, just to be sure ...
		randinitdelays = tdrandom * (np.random.rand(len(initdelays)) - 0.5) + initdelays
		lc.multisettimedelays(lcsc, randinitdelays)
		
		if optmode == "full" :
			finald2 = opt_full(lcsc, rawdispersionmethod, verbose=False)
		if optmode == "ml" :
			finald2 = opt_ml(lcsc, rawdispersionmethod, verbose=True)
		if optmode == "ts" :
			finald2 = opt_ts_mix(lcsc, rawdispersionmethod, verbose=True)
		if optmode == "none":
			couplelist = [couple for couple in [[lc1, lc2] for lc1 in lcsc for lc2 in lcsc] if couple[0] != couple[1]]
			d2values = np.array([rawdispersionmethod(*couple)["d2"] for couple in couplelist])
			finald2 = np.mean(d2values)
		
		timeshifts = np.array([l.timeshift for l in lcsc])
		tslist.append(timeshifts)
		d2list.append(finald2)
		
		# Crucial ... we revert to the original mags :
		#restoremags(lcsc, bumg) # no, not needed anymore. But i keep it so that you do not forget about this...

		if (b+1) % saven == 0 or b+1 == n:
			savefile()



def bootmerge(infilenamelist, outfilename):
	"""
	We read some pickles of several realisations made by boot (of the SAME lightcurves and with
	the SAME parameters !!!) and save them into a single pickle file.
	"""

	data = util.readpickle(infilenamelist[0])
	lcs = data["lcs"]
	bootmode = data["bootmode"]
	optmode = data["optmode"]
	tdrandom = data["tdrandom"]
	jdamplitude = data["jdamplitude"]
	
	n = 0
	tsslist = []
	d2slist = []

	for infilename in infilenamelist:
	
		data = util.readpickle(infilename)
		tsslist.append(data["tss"])
		d2slist.append(data["d2s"])
		#print data["tss"]
	
	tss = np.concatenate(tsslist, axis=1)
	d2s = np.concatenate(d2slist)
	#print tss.shape
	#print d2s.shape

	util.writepickle({"tss":tss, "d2s":d2s, "lcs":lcs, "n":n, "filename":outfilename,
		"bootmode":bootmode, "optmode":optmode, "tdrandom":tdrandom, "jdamplitude":jdamplitude}, outfilename)
		



#def threadboot(lcs, rawdispersionmethod, n=1000, saven=100 filename="multiopt_boot.pkl", bootmode="mags", optmode="full", jdamplitude=1.0, tdrandom=20.0):
#	"""
#	"""
#	import threading
#		
#	initdelays = lc.multigettimedelays(lcs)
#	
#	tdlist = []
#	d2list = []
#	b = 1 # the boot count...
#	
#		
#	def savefile():	
#		tds =  np.column_stack(tdlist)
#		d2s = np.array(d2list)
#		util.writepickle({"tds":tds, "d2s":d2s, "lcs":lcs, "n":b+1, "filename":filename,
#		"bootmode":bootmode, "optmode":optmode, "tdrandom":tdrandom, "jdamplitude":jdamplitude}, filename)
#	
#	
#	class singleboot ( threading.Thread ):
#		
#		def run ( self ):
#			global b, tdlist, d2list
#			print "Bootstrap %i" %b
#			
#			lcsc = [l.copy() for l in lcs]
#		
#			if bootmode == "mags" or bootmode == "both":
#				for l in lcsc:
#					l.magbootstrap()
#			if bootmode == "jds" or bootmode == "both":
#				for l in lcsc:
#					l.jdbootstrap(amplitude = jdamplitude, keepml=True)
#			
#			#lc.display(lcs)
#			# We set the time delays to random positions, just to be sure ...
#			randinitdelays = tdrandom * (np.random.rand(len(initdelays)) - 0.5) + initdelays
#			lc.multisettimedelays(lcsc, randinitdelays)
#			
#			if optmode == "full" :
#				finald2 = opt_full(lcsc, rawdispersionmethod, verbose=False)
#			if optmode == "ml" :
#				finald2 = opt_ml(lcsc, rawdispersionmethod, verbose=True)
#			if optmode == "ts" :
#				finald2 = opt_ts_mix(lcsc, rawdispersionmethod, verbose=True)
#			if optmode == "none":
#				couplelist = [couple for couple in [[lc1, lc2] for lc1 in lcsc for lc2 in lcsc] if couple[0] != couple[1]]
#				d2values = np.array([rawdispersionmethod(*couple)["d2"] for couple in couplelist])
#				finald2 = np.mean(d2values)
#			
#			tds = lc.multigettimedelays(lcsc)
#			tdlist.append(tds)
#			d2list.append(finald2)
#			
#			# Crucial ... we revert to the original mags :
#			#restoremags(lcsc, bumg) # no, not needed anymore. But i keep it so that you do not forget about this...
#	
#			if (b+1) % saven == 0:
#				savefile()
#
#			
#	
	
	

	
# 
# def plot_boot(filename="multiopt_boot.pkl", title=None, tdmin=None, tdmax=None, dispmax=None):
# 	"""
# 	The basic kind of time-delay histogram that you can plot from your bootstrapping.
# 	This is mainly here to give an example ...
# 	"""
# 	
# 	#print "WARNING I AM TWEAKED !!!"
# 	
# 	data = util.readpickle(filename)
# 	lcs = data["lcs"]
# 	tds = data["tds"]
# 	d2s = data["d2s"]
# 	
# 	ntds = len(tds) # 3 if you have 4 lcs...
# 	nboot = len(d2s)
# 	
# 	
# 	description = "%i bootstraps (bootmode = %s, optmode=%s) around :\n"% (nboot, data["bootmode"], data["optmode"])
# 	description += "\n".join([str(l) for l in lcs])
# 	print description
# 	
# 	if tdmax == None:
# 		tdmax = max(tds.ravel())
# 	if tdmin == None:
# 		tdmin = min(tds.ravel())
# 	
# 	fig = plt.figure(figsize=(6,3*ntds))
# 	
# 	for i, delays in enumerate(tds):
# 		ax = plt.subplot(ntds, 1, i+1)
# 		
# 		# Write the name of the delay on the graph
# 		ax.annotate("%s%s"%(lcs[0].object, lcs[i+1].object), xy=(0.05, 0.8),  xycoords='axes fraction')
# 		
# 		if i +1 == ntds:
# 			plt.xlabel("Time delay")
# 		if i == 0:
# 			if title != None:
# 				plt.title(title)
# 			
# 		n, bins, patches = plt.hist(delays, 200, range=(tdmin, tdmax), histtype='stepfilled', facecolor='#7777BB')
# 		
# 		if dispmax != None:
# 			ax.set_ylim( (0, dispmax) )
# 		
# 		for tick in ax.yaxis.get_major_ticks(): # reduce the label fontsize.
# 			tick.label1.set_fontsize(8)
# 	
# 	
# 	plt.figtext(0.25, 0.82, description, fontsize=8)
# 	
# 
# 	plt.show()
# 	
# 
# def plot_bootc(filename="multiopt_boot.pkl", title=None, tdmin=None, tdmax=None, dispmax=None):
# 	"""
# 	like plot_boot, but with colors ...
# 	
# 	"""
# 	
# 	#print "WARNING I AM TWEAKED !!!"
# 	
# 	data = util.readpickle(filename)
# 	lcs = data["lcs"]
# 	tds = data["tds"]
# 	d2s = data["d2s"]
# 	
# 	ntds = len(tds) # 3 if you have 4 lcs...
# 	nboot = len(d2s)
# 	
# 	
# 	description = "%i bootstraps (bootmode = %s, optmode=%s) around :\n"% (nboot, data["bootmode"], data["optmode"])
# 	description += "\n".join([str(l) for l in lcs])
# 	print description
# 	
# 	if tdmax == None:
# 		tdmax = max(tds.ravel())
# 	if tdmin == None:
# 		tdmin = min(tds.ravel())
# 	
# 	fig = plt.figure(figsize=(10,3.4*ntds))
# 	
# 	norm = colors.normalize(min(d2s), max(d2s))
# 	
# 	for i, delays in enumerate(tds):
# 		ax = plt.subplot(ntds, 1, i+1)
# 		
# 		# Write the name of the delay on the graph
# 		namebla = "%s%s"%(lcs[0].object, lcs[i+1].object)
# 		ax.annotate(namebla, xy=(0.05, 0.8),  xycoords='axes fraction')
# 		
# 		meandelay = np.mean(delays)
# 		stddelay = np.std(delays)
# 		tdbla = "%8.1f +/- %3.1f"%(meandelay, stddelay)
# 		
# 		ax.annotate(tdbla, xy=(0.01, 0.40),  xycoords='axes fraction')
# 		print namebla, tdbla
# 		
# 		if i +1 == ntds:
# 			plt.xlabel("Time delay")
# 		if i == 0:
# 			if title != None:
# 				plt.title(title)
# 			
# 		#n, bins, patches = plt.hist(delays, 200, range=(tdmin, tdmax), histtype='stepfilled', facecolor='#7777BB')
# 		n, bins, patches = plt.hist(delays, 200, range=(tdmin, tdmax))
# 		
# 		
# 		# We construct the mean dispersions for each bin
# 		meandisps = []
# 		for i in range(len(bins)-1):
# 			a = bins[i]
# 			b = bins[i+1]
# 			d2sofbin = []
# 			for j in range(len(delays)):
# 				if delays[j] < b and delays[j] > a:
# 					d2sofbin.append(d2s[j])
# 			#d2sofbin.append(2.4)
# 			meandisps.append(np.mean(np.array(d2sofbin)))
# 		
# 		meandisps = np.array(meandisps)
# 		meandisps = np.clip(np.nan_to_num(meandisps), np.nanmin(meandisps), np.nanmax(meandisps))
# 		#print meandisps
# 		#print (min(meandisps), max(meandisps))
# 		
# 
# 		#
# 		
# 		#print norm(2.35)
# 		#for patch in patches:
# 		#	print patch
# 
# 		#print meandisps.min(), meandisps.max()
# 		
# 		for thismeandisp, thispatch in zip(meandisps, patches):
# 			color = cm.jet(norm(thismeandisp))
# 	  		thispatch.set_facecolor(color)
# 		
# 		if dispmax != None:
# 			ax.set_ylim( (0, dispmax) )
# 		
# 		for tick in ax.yaxis.get_major_ticks(): # reduce the label fontsize.
# 			tick.label1.set_fontsize(8)
# 	
# 	
# 	plt.figtext(0.25, 0.82, description, fontsize=8)
# 	
# 
# 	#plt.show()
# 	
# 	filename = os.path.basename(filename)
# 	filebase = os.path.splitext(filename)[0]
# 	
# 	print filebase
# 	plt.savefig(filebase + ".png")
# 	


def plot_bootshifthist(filename="multiopt_boot.pkl", ccurve=0, title=None, shiftmin=None, shiftmax=None, dispmax=None):
	"""
	Plots histograms of the time shifts resulting from boot ...
	
	ccurve = 0 : use mean dispersion per bin as color
	ccurve = n : use the nth lc's timeshift as color (1 = first light curve)
	
	"""
	
	data = util.readpickle(filename)
	lcs = data["lcs"]
	tss = data["tss"] # These are the time shifts of the 4 curves.
	d2s = data["d2s"]
	
	nlcs = len(lcs)
	nboot = len(d2s)
	
	# For this first histogram we want to plot the 4 timeshifts with respect to the mean of them.
	#print tss
	timeshifts = tss - np.mean(tss, axis=0)
	
	#sys.exit()
	
	description = "%i bootstraps (bootmode = %s, optmode=%s) around :\n"% (nboot, data["bootmode"], data["optmode"])
	description += "\n".join([str(l) for l in lcs])
	print description
	
	if shiftmax == None:
		shiftmax = np.max(timeshifts.ravel())
	if shiftmin == None:
		shiftmin = np.min(timeshifts.ravel())
	
	fig = plt.figure(figsize=(10,3.4*nlcs))
	
	if ccurve == 0:
		norm = colors.normalize(min(d2s), max(d2s))
	else:
		norm = colors.normalize(min(timeshifts[ccurve-1,:]), max(timeshifts[ccurve-1,:]))
	
	for i, curveshifts in enumerate(timeshifts):
		ax = plt.subplot(nlcs, 1, i+1)
		
		minorLocator = MultipleLocator(1) # so to have a tick every day
		#minorFormattor = FormatStrFormatter('%0.1f')
		ax.xaxis.set_minor_locator(minorLocator)
		
		# Write the name of the curve on the graph
		#namebla = "%s"%(lcs[i].object)
		namebla = "%s" % (str(lcs[i]))
		ax.annotate(namebla, xy=(0.02, 0.85),  xycoords='axes fraction')
		
		#meanshift = np.mean(curveshifts)
		#stdshift = np.std(curveshifts)
		#tdbla = "%8.1f +/- %3.1f"%(meandelay, stddelay)
		
		#ax.annotate(tdbla, xy=(0.01, 0.40),  xycoords='axes fraction')
		#print namebla, tdbla
		
		if i+1 == nlcs:
			plt.xlabel("Relative time shift [days]")
		if i == 0:
			if title != None:
				plt.title(title)
			
		#n, bins, patches = plt.hist(delays, 200, range=(tdmin, tdmax), histtype='stepfilled', facecolor='#7777BB')
		n, bins, patches = plt.hist(curveshifts, 200, range=(shiftmin, shiftmax))
		
		
		# We construct the mean dispersions for each bin
		meandisps = []
		meancshifts = []
		for i in range(len(bins)-1):
			a = bins[i]
			b = bins[i+1]
			d2sofbin = []
			cshiftsofbin = []
			for j in range(nboot):
				if curveshifts[j] < b and curveshifts[j] > a:
					d2sofbin.append(d2s[j])
					if ccurve > 0:
						cshiftsofbin.append(timeshifts[ccurve-1, j])
			meancshifts.append(np.mean(np.array(cshiftsofbin)))
			meandisps.append(np.mean(np.array(d2sofbin)))
		
		meandisps = np.array(meandisps)
		meandisps = np.clip(np.nan_to_num(meandisps), np.nanmin(meandisps), np.nanmax(meandisps))
		meancshifts = np.array(meancshifts)
		meancshifts = np.clip(np.nan_to_num(meancshifts), np.nanmin(meancshifts), np.nanmax(meancshifts))
		
		
		#print meandisps
		#print (min(meandisps), max(meandisps))
		

		#
		
		#print norm(2.35)
		#for patch in patches:
		#	print patch

		#print meandisps.min(), meandisps.max()
		
		if ccurve == 0:
			for thismeandisp, thispatch in zip(meandisps, patches):
				color = cm.jet(norm(thismeandisp))
	  			thispatch.set_facecolor(color)
		else :
			for thismeandisp, thispatch in zip(meancshifts, patches):
				color = cm.jet(norm(thismeandisp))
	  			thispatch.set_facecolor(color)
		
		
		if dispmax != None:
			ax.set_ylim( (0, dispmax) )
		
		for tick in ax.yaxis.get_major_ticks(): # reduce the label fontsize.
			tick.label1.set_fontsize(8)
	
	
	#plt.figtext(0.25, 0.82, description, fontsize=8)
	

	#plt.show()
	
	filename = os.path.basename(filename)
	filebase = os.path.splitext(filename)[0]
	
	#print filebase
	plt.savefig(filebase + ".png")
	
	



def plot_bootdelayhist(filename="multiopt_boot.pkl", tdmin=None, tdmax=None, nbins=200, showplot=False):
	"""
	Plots one single histogram for each delay ( = pair of lcs)
	
	ccurve = 0 : use mean dispersion per bin as color
	ccurve = n : use the nth lc's timeshift as color (1 = first light curve)
	
	"""
	
	data = util.readpickle(filename)
	lcs = data["lcs"]
	tss = data["tss"] # These are the time shifts of the 4 curves.
	d2s = data["d2s"]
	
	nlcs = len(lcs)
	nboot = len(d2s)
	
	print "We have %i lightcurves." % (nlcs)
	# We prepare a list of indexes to consider
	couples = [(i, j) for i in range(nlcs) for j in range(nlcs) if i<j]
	
	filename = os.path.basename(filename)
	filebase = os.path.splitext(filename)[0]
	
	description = "%i runs (bootmode = %s, optmode=%s) around :\n"% (nboot, data["bootmode"], data["optmode"])
	description += "\n".join([str(l) for l in lcs])
	
	for couple in couples:
	
		outfilebase = "%s_%s%s" % (filebase, lcs[couple[0]].object, lcs[couple[1]].object)
 		#ashifts = tss[couple[0]]
		#bshifts = tss[couple[1]]
		delays = tss[couple[1]] - tss[couple[0]]
		
		meandelay = np.mean(delays)
		stddelay = np.std(delays)
		tdbla = "%8.1f +/- %3.1f"%(meandelay, stddelay)
		#topline = "%s%s : %s" % (lcs[couple[0]].object, lcs[couple[1]].object, tdbla)
		topline = "%s%s" % (lcs[couple[0]].object, lcs[couple[1]].object)
		
		print topline
		
		fig = plt.figure(figsize=(6,3))
		ax = fig.add_axes((0.1, 0.15, 0.85, 0.75))
	
		minorLocator = MultipleLocator(1) # so to have a tick every day
		ax.xaxis.set_minor_locator(minorLocator)
		
		ax.annotate(topline, xy=(0.03, 0.88),  xycoords='axes fraction', fontsize=18, color="blue")
		#ax.annotate("%s" % (str(lcs[couple[0]])), xy=(0.03, 0.81),  xycoords='axes fraction', fontsize=10)
		#ax.annotate("%s" % (str(lcs[couple[1]])), xy=(0.03, 0.76),  xycoords='axes fraction', fontsize=10)
		ax.annotate("%s" % (description), xy=(0.03, 0.85),  xycoords='axes fraction', fontsize=10, verticalalignment="top", color="gray")
		
		
		
		ax.hist(delays, nbins, range=(tdmin, tdmax))
		ax.set_title(outfilebase)
		ax.set_xlabel("Delay [days]")
		plt.savefig(outfilebase + ".png")
		if showplot:
			plt.show()

	
	
	