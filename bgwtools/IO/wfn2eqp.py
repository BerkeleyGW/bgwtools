#!/usr/bin/env python

# Gets all the energies from a wavefunction and outputs a eqp-file for a more
#  convenient manipulation of the data.
ryd = 13.60569193

from numpy import *
from wfn import wfnIO
import sys

if len(sys.argv)<3:
	print 'Usage: %s WFN eqp_file'%(sys.argv[0])
	sys.exit()

wfn = wfnIO(sys.argv[1])

kpts = wfn.kpt
el = wfn.energies*ryd
nb = el.shape[0]
nk = el.shape[1]
ns = el.shape[2]

eqp = open(sys.argv[2], 'w')
for ik in range(nk):
	eqp.write( (3*' %12.9f'+' %7d\n')%tuple(kpts[ik].tolist() + [nb*ns]) )
	for ib in range(nb):
		for _is in range(ns):
			#eqp.write( (' %7d %7d'+2*' %14.9f')%(_is+1, 1+ib, \
			eqp.write( (' %7d %7d'+2*' %23.18f')%(_is+1, 1+ib, \
				el[ib,ik,_is], 0.0) + '\n')
eqp.close()

