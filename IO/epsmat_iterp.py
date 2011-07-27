#!/usr/bin/env python

'''
Duplicates the k-grid of an epsmat file
'''

from __future__ import division
from epsmat import epsmatIO
from numpy import *

if __name__=='__main__':
	import sys
	if len(sys.argv)<2:
		print 'Expecting one argument: [epsmat/eps0mat file]'
		sys.exit()

	epsmat = epsmatIO(sys.argv[1])
	qpts = empty((epsmat.nq,3))
	qpts_cub = empty_like(qpts)
	shift = empty_like(qpts)
	for n in xrange(epsmat.nq):
		qpt = epsmat.q[n]
		qpts[n,:] = qpt
		qpts_cub[n,:] = qpt - around(qpt)
		shift[n,:] = around(qpt)

	print qpts
	print qpts_cub
	print shift

	
