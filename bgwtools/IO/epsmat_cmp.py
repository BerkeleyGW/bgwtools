#!/usr/bin/env python

from epsmat import epsmatIO
import sys
from numpy import *

def _cmp(x,y):
	return all(fabs(array(x)-array(y)) < 1e-15)
	#return all(fabs(array(x)-array(y)) < 1e-2)

def compare(file1, file2):
	eps1 = epsmatIO(file1)
	eps2 = epsmatIO(file2)
	ok = True
	ok &= eps1.nq==eps2.nq
	ok &= _cmp(eps1.isort, eps2.isort)
	ok &= _cmp(eps1.isort_i, eps2.isort_i)
	ok &= _cmp(eps1.ekin, eps2.ekin)
	for em1,em2 in zip(eps1.epsmat, eps2.epsmat):
		#ok &= _cmp(eps1.epsmat, eps2.epsmat)
		ok &= _cmp(em1, em2)
	ok &= _cmp(eps1.q, eps2.q)
	#ok &= _cmp(eps1., eps2.)
	#ok &= _cmp(eps1., eps2.)

	return ok

if __name__=='__main__':
	if len(sys.argv)!=3:
		print 'Usage: %s epsmat1 epsmat2'%(sys.argv[0])
		sys.exit()
	ok = compare(sys.argv[1],sys.argv[2])
	print 'Same files?', ok

