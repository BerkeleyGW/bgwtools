#!/usr/bin/env python

# This utility compares two wfn files and prints the larger difference.
# Felipe H. da Jornada, Apr 2013

from bgwtools.IO.wfn import wfnIO
from bgwtools.common.common import get_numpy_flavor
import sys
from numpy import amax, empty, empty_like

fname1 = sys.argv[1]
print('File #1: %s'%(fname1))
fname2 = sys.argv[2]
print('File #2: %s'%(fname2))

wfn1 = wfnIO(fname1)
wfn2 = wfnIO(fname2)

ng_max = amax(wfn1.ngk) * wfn1.nspinor
gvec = empty((3, ng_max), dtype=int, order='F')
data1 = empty((ng_max, wfn1.ns), dtype=get_numpy_flavor(wfn1.flavor), order='F')
data2 = empty_like(data1)
for ik in range(wfn1.nk):
	print('ik = %d'%ik)
        wfn1.read_gvectors(gvec)
        wfn2.read_gvectors(gvec)
        ngk = wfn1.ngk[ik] * wfn1.nspinor
        for ib in range(wfn1.nbands):
            wfn1.read_data(data1)
            wfn2.read_data(data2)
	    diff = (data2 - data1)[:ngk]
            max_diff = amax(abs(diff))
            print('  band %d:'%(ib))
            print('    ||cg2 - cg1||_inf = %e'%max_diff)
