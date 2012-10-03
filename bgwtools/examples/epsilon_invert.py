#!/usr/bin/env python

# This example opens an epsilon file, inverts (of rather, uninverts) the matrix,
# and saves the result to a new epsilon file.

# Felipe Homrich da Jornada <jornada@civet.berkeley.edu> (2012)

from bgwtools.IO.epsmat import epsmatIO
import sys
from numpy import *

if len(sys.argv)!=3:
	print('Usage: %s epsmat_input epsmat_output'%(sys.argv[0]))
	exit(1)

f_in = sys.argv[1]
f_out = sys.argv[2]

# Load the epsmat file into the epsmat object
epsmat = epsmatIO(f_in)

print
print 'Inverting epsilon matrices'
for iq in range(epsmat.nq):
	nmtx = epsmat.nmtx[iq]
	print ' - inverting a %d x %d matrix (%d/%d)'%(nmtx, nmtx, iq+1, epsmat.nq)
	mat = epsmat.epsmat[iq]
	epsmat.epsmat[iq] = linalg.inv(mat)

print
print 'Saving your the epsmat file to',f_out

epsmat.to_file(f_out)
