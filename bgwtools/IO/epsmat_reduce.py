#!/usr/bin/env python

# Given an epsmat file and a WFN file, create a new epsmat file with
# only the q-points present in the WFN

from epsmat import epsmatIO
from wfn import wfnIO
import sys
from numpy import *

S=1e-8

def reduce(eps, kpts):
	del_cnt = 0
	for iq in range(eps.nq):
		kp = eps.qpt[iq]
		#print kp
		cond = all(fabs(kp-kpts)<S, axis=1)
		if not any(cond):
			ii = iq - del_cnt
			
			eps.nmtx = delete(eps.nmtx, ii, 0)
			eps.q = delete(eps.q, ii, 0)
			eps.nq -= 1
			del eps.isort[ii]
			del eps.isort_i[ii]
			del eps.ekin[ii]
			del eps.epsmat[ii]	
			del_cnt += 1
		#else:
		#	print kpt
	eps.qpt = eps.q.copy()

if __name__=='__main__':
	if len(sys.argv)<2:
		print 'Usage: %s epsmat WFN'%(sys.argv[0])
		sys.exit()
	wfn = wfnIO(sys.argv[2])
	kpts = wfn.kpt
	print kpts
	eps = epsmatIO(sys.argv[1])
	eps.grid = wfn.kgrid
	reduce(eps, kpts)
	print eps
	eps.to_file('epsmat_out')
	#print eps

