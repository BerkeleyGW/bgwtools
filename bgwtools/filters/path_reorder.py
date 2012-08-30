#!/usr/bin/env python

from __future__ import division

#maximum spacing between bands
max_delta=2.0 #eV
#maximum lapacian compared to the std
max_lap=1.0
#max_lap=0.5
#should we keep the Dirac cone alone?
#graphene_like=False
graphene_like=True
#E_Fermi=-1.54
E_Fermi=-1.0845


import sys
from numpy import *
from itertools import permutations

if len(sys.argv)<3:
	print 'Usage: %s input output'%(sys.argv[0])
	exit()

class Cfixer:
	def __init__(self, fname):
		self.fname=fname
		self.data=loadtxt(self.fname,dtype=float64)
		self.nbands=self.data.shape[1]-3
		self.npts=self.data.shape[0]
		#print nbands,npts
		self.bands_arange = arange(self.nbands)
		self.bands = self.data[:,3:] #this is a view!
		self.iteration=0

	def lap_cost(self, k, per):
		#bands: array of [perm_left, perm_center, perm_right]
		lap = self.bands[k-1, per[0]] - 2*self.bands[k, per[1]] + self.bands[k+1, per[2]]
		return sum(lap**2)

	def fix_clusters(self, k, crit_clusters):
		for cluster in crit_clusters:
			perm0 = tuple(cluster)
			best_perm = perm0
			no_swap = [perm0]*3
			cost_center = self.lap_cost(k, no_swap)
			cost0 = (self.lap_cost(k-1, no_swap) + cost_center,\
				cost_center + self.lap_cost(k+1, no_swap))
			best_gain = 0.0
			for perm in permutations(cluster):
				if perm==perm0:	continue
				perms=(perm0,perm0,perm,perm,perm)
				for shift in (0,1):
					#swap the bond (shift: 0=>left; 1=>right)
					cost = self.lap_cost(k-1+shift, perms[0:]) + self.lap_cost(k+shift, perms[1:])
					gain = cost - cost0[shift]
					if (gain < best_gain):
						best_gain = gain
						best_perm = perm
						best_shift = shift
					
			if best_perm!=perm0:
				kp = k+best_shift
				print '\t',perm0,'->',best_perm,'@',kp
				self.bands[kp:,perm0] = self.bands[kp:,best_perm]
				self.lap_mean[array(perm0)] = self.lap_mean[array(best_perm)]
				self.lap_lim [array(perm0)] = self.lap_lim [array(best_perm)]

	def iterate(self):
		self.iteration+=1
		#laplacian
		lap_static = self.bands[:-2,:] + self.bands[2:,:] - 2*self.bands[1:-1,:]
		self.lap_mean = mean(lap_static, axis=0)
		self.lap_lim  = std (lap_static, axis=0)*max_lap
	
		for k in range(1, self.npts-3):
			#check if the laplacian is bigger than threshold
			lap = self.bands[k-1,:] + self.bands[k+1,:] - 2*self.bands[k,:]
			is_crit = abs(lap-self.lap_mean) > self.lap_lim 
			if not any(is_crit): continue
			#graphene test:
			if graphene_like:
				for b in self.bands_arange[is_crit]:
					if abs(self.bands[k,b]-E_Fermi) < max_delta:
						#print 'Band',b,'is near the Fermi energy, so it will be removed.'
						is_crit[b] = False
			if not any(is_crit): continue
			crit_bands = self.bands_arange[is_crit]
			if (self.iteration==4)and(k==15):
				print self.bands[k,:]
			#print is_crit
			#print 'Critical bands:',crit_bands
			#group all the the critical bands into clusters
			crit_clusters = []

			def add_pt(pt, new_cluster=True):
				#create a new group with the band
				if not new_cluster:
					crit_clusters[-1].append(pt)
				else:
					crit_clusters.append([pt])
				#look for close bands
				diffs = self.bands[k,:] - self.bands[k, pt]
				is_near = abs(diffs)<max_delta
				near_bands = self.bands_arange[is_near]
				for b in near_bands:
					if not (b in crit_clusters[-1]):
						add_pt(b, False)

			for crit_b in crit_bands:
				#check if the band was already added to crit_clusters
				old_b = False
				for cluster in crit_clusters:
					if crit_b in cluster:
						old_b = True
						break
				if not old_b: add_pt(crit_b)
				if len(crit_clusters[-1])==1:
					#print 'Dropping single-populated cluster. There might be missing bands!'
					del crit_clusters[-1]

			#do the magic
			if len(crit_clusters):
				print ' k=%2d'%(k),'\tcritical clusters:',crit_clusters
				self.fix_clusters(k, crit_clusters)

fixer = Cfixer(sys.argv[1])

for i in range(5):
	print 
	print 'ITERATION',i+1
	print
	fixer.iterate()
	
savetxt(sys.argv[2], fixer.data, '%8.5f')
