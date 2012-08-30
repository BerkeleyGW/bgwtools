#!/usr/bin/env python

from epsmat import epsmatIO
import sys
from numpy import *

def merge(files):
	eps = epsmatIO(files[0])

	for n in range(1,len(files)):
		fname = files[n]
		eps2 = epsmatIO(fname)

		#merging
		eps.nq += eps2.nq
		eps.qpt = row_stack((eps.qpt, eps2.qpt))
		eps.nmtx = append(eps.nmtx, eps2.nmtx)
		eps.isort += eps2.isort
		eps.isort_i += eps2.isort_i
		eps.ekin += eps2.ekin
		eps.epsmat += eps2.epsmat
		eps.q = row_stack((eps.q, eps2.q))

	return eps

if __name__=='__main__':
	if len(sys.argv)<2:
		print 'Usage: %s epsmat1 [...] [epsmatN]'%(sys.argv[0])
		sys.exit()
	eps = merge(sys.argv[1:])
	print eps

