"""
Compute some statistics about TDC1
"""

import os,sys
import numpy as np
import math
import pycs



def progress(estimates,htmlcode=False,path='leaderboard.html'):

	"""
	Print the overall progress of TDC1 estimations through D3CS
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
	nestims      = [nperestids.count(nestim) for nestim in set(nperestids)] # number of pairs with 1 estimation, 2, 3,...
	nestimcumuls = [sum(nestims[ind:]) for ind in np.arange(len(nestims))]
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
		
		# Overall estimation 

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


def achievements(estimates,user, htmlcode=False,path='leaderboard.html'):

	"""
	D3CS Achievements, fun stuff only (no science here !)
	Spoiler Alert : contains sensible informations (such the winning pairs)
	"""


	# Load estimates	
	userests = [est for est in estimates if est.methodpar == user]
	if len(userests) < 50:
		return 


	### Wall of Fame initialisation
	print '='*55
	print 'Wall of Fame of %s ' % user
 	if htmlcode:		
		f=open(path,'a')
		sep = '='*33
		
		f.write('<p>'+sep+'</p> \n')
		f.write('<h3>Wall of Fame of %s :</h3> \n' % user)
		
		
	# Number of estimations		
	nests = len(list(set([est.id for est in userests])))
	frac = nests/5180.0*100	
	print '-- %i unique estimations ' % (nests)
	if htmlcode:
		f.write('<p>-- %i unique estimations</p> \n'% (nests))


	# Medals

	if frac>5:
		print '-- Chocolate Medal awarded !  '
		if htmlcode:
			f.write('<p>-- Chocolate Medal awarded !  </p> \n')
	if frac>15:
		print '-- Bronze Medal awarded !  '
		if htmlcode:
			f.write('<p>-- Bronze Medal awarded !  </p> \n')
	if frac>30:
		print '-- Silver Medal awarded !  '
		if htmlcode:
			f.write('<p>-- Silver Medal awarded !  </p> \n')
	if frac>60:
		print '-- Gold Medal awarded !  '
		if htmlcode:
			f.write('<p>-- Gold Medal awarded !  </p> \n')
	if frac>100:
		print '-- Awesome Medal awarded !  '
		if htmlcode:
			f.write('<p>-- Awesome Medal awarded !  </p> \n')

	
	# Exploration
	
	if len(list(set([est.rung for est in userests]))) == 5:
		print '-- Curious'
		if htmlcode:
			f.write('<p>-- Curious </p> \n')
	 
	if len(list(set([est.rung for est in userests]))) == 1 and nests>100:
		print '-- I am a keeper '
		if htmlcode:
			f.write('<p>-- I am a keeper </p> \n')		


	# Time Related
	
	fastests = [est for est in userests if est.timetaken < 10 and est.confidence == 1]
	if len(fastests) > 0:
		print '-- Fast as the wind'
		if htmlcode:
			f.write('<p>-- Fast as the wind</p> \n')
			
	longests = [est for est in userests if est.timetaken > 180 and est.tderr < 5]
	if len(longests) > 0:
		print '-- Slowly but surely'
		if htmlcode:
			f.write('<p>-- Slowly but surely</p> \n')
			
						
	esttime  = [est.timetaken for est in estimates]
	usertime = [est.timetaken for est in userests]

	if max(usertime) == max(esttime):
		print '-- Background task'
		if htmlcode:
			f.write('<p>-- Background task</p>')
			
			
	if min(usertime) == min(esttime):
		print '-- Faster than light'
		if htmlcode:
			f.write('<p>-- Faster than light</p>')
			
	esttderr  = [est.tderr for est in estimates]
	usertderr = [est.tderr for est in userests]
	
	if max(usertderr) == max(esttderr):
		print '-- I will be back'
		if htmlcode:
			f.write('<p>-- I will be back</p>')				
						
	if min(usertderr) == min(esttderr):
		print '-- Surgical Strike'
		if htmlcode:
			f.write('<p>-- Surgical Strike</p>')
			
	# Confidence related
	
	userconf = [est.confidence for est in userests]
	if np.median(userconf) < 1.5:
		print '-- Self-confident'
		if htmlcode:
			f.write('<p>-- Self-confident</p>')
			
	if np.median(userconf) >2.2:
		print '-- Hesitating'
		if htmlcode:
			f.write('<p>-- Hesitating</p>')	
			
	

			
	### Special Achievements
		
	# The Big Lottery
	
	winning_combs = [(1, 971), (1, 365), (0, 612), (4, 887), (3, 107), (2, 89), (4, 92), (0, 230), (0, 236), (3, 693)] # randomly chosen

	
	ncomb = 0 
	for comb in winning_combs:
		combid = 'tdc1_%i_%i' % (comb[0],comb[1])

		if combid in [est.id for est in userests]:
			ncomb += 1
	
	print '-- Winning pairs estimated: %i/10' % ncomb
	print ''	
	if htmlcode:
		f.write('<p>-- Winning pairs estimated: %i/10</p> \n'% ncomb)


	# The B.E.a.R (for Belgian Ethanol As Reward)
	
	bearlevel = frac//10
	if bearlevel > 0: 
		print '-- B.E.a.R level %i reached' % bearlevel
		if htmlcode:
			f.write('<p>-- B.E.a.R level %i reached</p> \n' % bearlevel)
			
	
	
	
	# The sweet pairs (one chocolate per estimate)
	

	sweet_combs  = [(2, 987), (1, 287), (1, 679), (0, 13), (4, 899), (2, 971), (2, 977), (0, 884), (1, 802), (3, 859),(4, 1036), (4, 131), (2, 223), (3, 496), (1, 278), (3, 96), (0, 678), (0, 162), (1, 772),(2, 833), (3, 467), (2, 355), (1, 330), (3, 698), (1, 374), (1, 972), (2, 976), (2, 775)]
	
	
	for comb in sweet_combs:
		combid = 'tdc1_%i_%i' % (comb[0],comb[1])		
		ids     = [est.id for est in estimates]
		try:

			combpos = ids.index(combid)

			sweetest = estimates[combpos]
			if sweetest.methodpar == user:
				print '-- Sweet pair estimated ! %s ' % sweetest.niceid
				if htmlcode:
					f.write('<p>-- Sweet pair estimated ! %s</p>' % sweetest.niceid)
		
		except:
			pass		 		
	
	
	f.close()
	


