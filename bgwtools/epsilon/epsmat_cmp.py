#!/usr/bin/env python

# This utility compares two epsmat files and prints the larger difference.
# Felipe H. da Jornada, Apr 2013

from bgwtools.IO.epsmat import epsmatIO
import sys
from numpy import amax

fname1 = sys.argv[1]
print('File #1: %s'%(fname1))
fname2 = sys.argv[2]
print('File #2: %s'%(fname2))

eps1 = epsmatIO(fname1)
eps2 = epsmatIO(fname2)

for iq in range(eps1.nq):
	print('iq = %d'%iq)
	diff = eps2.epsmat[iq] - eps1.epsmat[iq]
	max_diff = amax(abs(diff))
	print('||eps2 - eps1||_inf = %e'%max_diff)
