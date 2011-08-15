#!/usr/bin/env python

import cPickle
from expand import epsmat_intp, tprint
from numpy import arange

tprint('Loading File')
f = open('epsmat_interp.pkl')
epsi = cPickle.load(f)
f.close()

tprint('Interpolating')

#420
epsi.interpolate(arange(3))

tprint('Done')
