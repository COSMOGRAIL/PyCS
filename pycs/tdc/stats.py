"""
Compute some statistics about TDC1
"""

import numpy as np

def progress(estimates, htmlcode=False, path='leaderboard.html'):

	"""
	Print the overall progress of TDC1 estimations through D3CS

	@param estimates: list of Estimate objects.
	@param htmlcode: boolean. If True, I will write the results in an html file
	@param path: path of the html output file.
	"""

	print 'TDC1 Progress: \n'
	
	users = list(set([est.methodpar for est in estimates]))
	print '=======================================|'
	print '    User \t estimations \t  %    |'
	print '=======================================|'
	
	stats=[]
	for user in users:
		
		userests = [est for est in estimates if est.methodpar == user]
		nests = len(list(set([est.id for est in userests])))
		frac = nests/5180.0*100.0
		stats.append((user,nests,frac))
		
		print '--  %s \t      %i \t %.2f   |'   % (user,nests,frac)

	estids       = [est.id for est in estimates]		
	nperestids   = [estids.count(estid) for estid in set(estids)] # number of estimations per pairs
	nperestids_corr = []
	for n in nperestids:
		if n > 5:
			nperestids_corr.append(5)
		else:
			nperestids_corr.append(n)
		

	nestims      = [nperestids_corr.count(nestim) for nestim in set(nperestids_corr)] # number of pairs with 1 estimation, 2, 3,...
	print nestims
	
	if min(nperestids_corr) > 1:
		nestims.insert(0,0)
		
	if min(nperestids_corr) > 2:
		nestims.insert(0,0)		
	
	nestimcumuls = [sum(nestims[ind:]) for ind in np.arange(len(nestims))]
	print nestimcumuls
	
	fracs	     = [nestimcumul / 5180.0 * 100 for nestimcumul in nestimcumuls]
	maxoverall   = 3
	print '=======================================|'
	
	for ind in np.arange(min(len(nestimcumuls),maxoverall)):
		print '-- Overall %i x \t %i/5180 \t %.2f |'   % (ind+1,nestimcumuls[ind],fracs[ind])
		print '=======================================|'

	if htmlcode:
		
		stats = sorted(stats, key=lambda stat: stat[1])
		stats = stats[::-1]
		f=open(path,'w')
		# File initialisation
		f.write('<!DOCTYPE html> \n')
		f.write('<html> \n')
		f.write('<body> \n')
		f.write('<h1>Leaderboard:</h1> \n')
		f.write('<table border="1"> \n')
		f.write('<tr> \n')
  		f.write('<th>User</th> \n')
  		f.write('<th>Estimates</th> \n')
  		f.write('<th>%</th> \n')
		f.write('</tr> \n')
		
		# Write one line per user
		
		for stat in stats:
			f.write('<tr> \n')
			f.write('<td>'+str(stat[0])+'</td> \n')
			f.write('<td>'+str(stat[1])+'</td> \n')
			f.write('<td>%.2f</td> \n'% float(stat[2])) 
			f.write('</tr> \n')
		
		# Overall estimation (buggy, need to investigate...)

		
		for ind in np.arange(min(len(nestimcumuls),maxoverall)):
			f.write('<tr> </tr>\n')
			f.write('<td>Overall '+str(ind+1)+'x</td> \n')
			f.write('<td>'+str(nestimcumuls[ind])+'</td> \n')
			f.write('<td>%.2f</td> \n' % float(fracs[ind]))

		# File closing
		f.write('</table> \n')
		f.write('</body> \n')
		f.write('</html> \n')
		f.close()



