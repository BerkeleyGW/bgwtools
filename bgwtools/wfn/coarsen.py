#!/usr/bin/env python

# Reduces the number of k-points in each direction by a constant factor.

# Felipe Homrich da Jornada (Feb 2013)

from numpy import *
from bgwtools.IO.wfn import wfnIO
import sys
from bgwtools.common.common import get_numpy_flavor

if len(sys.argv)!=4:
    print('Usage: %s fact wfn_in wfn_out'%(sys.argv[0]))
    sys.exit(1)

fact = int(sys.argv[1])
fname_in = sys.argv[2]
fname_out = sys.argv[3]

wfn_in = wfnIO(fname_in)

tol = 1e-6
delta0 = 1.0/wfn_in.kgrid
shift = wfn_in.kshift * delta0 #doesn`t work!!
shift = wfn_in.kpt[:,0]
dks = wfn_in.kpt - shift[:,newaxis]
#print dks*wfn_in.kgrid[:,newaxis]
should_keep = all( around(dks*wfn_in.kgrid[:,newaxis]).astype(int)%fact == 0, axis=0)

print wfn_in.nk
print sum(should_keep)

if raw_input('Writing data to %s. Are you sure? [y/N] '%(fname_out))!='y':
    sys.exit(0)

wfn_out = wfnIO()
wfn_out.__dict__ = wfn_in.__dict__.copy()
wfn_out.nk = sum(should_keep)
wfn_out.ngk = wfn_out.ngk[should_keep]
wfn_out.ngkmax = amax(wfn_out.ngk)
wfn_out.kw = wfn_out.kw[should_keep]
wfn_out.kw[:] = 1.0/wfn_out.nk
wfn_out.kpt = wfn_out.kpt[:, should_keep]
wfn_out.ifmin = wfn_out.ifmin[should_keep, :]
wfn_out.ifmax = wfn_out.ifmax[should_keep, :]
wfn_out.energies = wfn_out.energies[:, should_keep, :]
wfn_out.occupations = wfn_out.occupations[:, should_keep, :]
wfn_out.f = None
wfn_out.fname = fname_out
wfn_out.write_header(full=True)

ng_max = amax(wfn_in.ngk)
gvec = empty((3, ng_max), dtype=int, order='F')
data = empty((ng_max, wfn_in.ns), dtype=get_numpy_flavor(wfn_in.flavor), order='F')
k_done = 0
for ik in range(wfn_in.nk):
    if should_keep[ik]:
        k_done += 1
        print k_done,'/',wfn_out.nk
        wfn_in.read_gvectors(gvec)
        wfn_out.write_gvectors(gvec[:,:wfn_in.ngk[ik]])
        for ib in range(wfn_in.nbands):
            wfn_in.read_data(data)
            wfn_out.write_data(data[:wfn_in.ngk[ik],:])
        if k_done==wfn_out.nk:
            break
    else:
        wfn_in.read_gvectors()
        for ib in range(wfn_in.nbands):
            wfn_in.read_data()
    
print wfn_out

print 'All done!'
