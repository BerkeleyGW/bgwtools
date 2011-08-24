#!/usr/bin/env python

import cPickle
from expand import epsmat_intp, tprint
from numpy import arange
import sys

if len(sys.argv)!=3:
	print 'Usage: interp.py start_point end_point'
	print 'Eg: interp.py 1 2 computes the first and second points'
	sys.exit()

tprint('Loading File')
f = open('epsmat_interp.pkl')
epsi = cPickle.load(f)
f.close()

tprint('Interpolating')
pt_s = int(sys.argv[1])
pt_e = int(sys.argv[2])
epsi.interpolate(arange(pt_s-1, pt_e))

tprint('Done')
