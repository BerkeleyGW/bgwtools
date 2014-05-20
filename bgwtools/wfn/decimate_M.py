#!/usr/bin/env python

# This script decimates many k-points from graphene and keeps only one wedge
# corresponding to one M point.

# Felipe Homrich da Jornada (Jun 2013)

from numpy import *
from bgwtools.IO.wfn import wfnIO
import sys
from bgwtools.common.common import get_numpy_flavor


if len(sys.argv)!=3:
    print('Usage: %s wfn_in wfn_out'%(sys.argv[0]))
    sys.exit(1)

fname_in = sys.argv[1]
fname_out = sys.argv[2]

wfn_in = wfnIO(fname_in)

def get_delta(kpts, idim):
    delta = fabs(kpts[idim,:] - wfn_in.kpt[idim,0])
    cond = delta > 1e-8
    return amin(delta[cond])
dx = get_delta(wfn_in.kpt, 0)
dy = get_delta(wfn_in.kpt, 1)
x0 = amin(wfn_in.kpt[0,:])
y0 = amin(wfn_in.kpt[1,:])

print x0, y0
print dx, dy
DELTA = 1e-8

kpts = wfn_in.kpt
condLB = kpts[0,:] < kpts[1,:] + DELTA
condLT = 2.*kpts[0,:] < (1. - kpts[1,:]) + DELTA
condRB = 2.*(1. - kpts[0,:]) < kpts[1,:] + DELTA
condRT = 1. - kpts[0,:] < (1. - kpts[1,:]) + DELTA

cond = condLB & condLT
cond |= condRB & condRT
print cond
should_keep = where(cond)[0]
print should_keep

print len(should_keep)

kpt = kpts[:,cond]
if 0:
    import matplotlib.pyplot as plt
    plt.scatter(kpt[0,:], kpt[1,:])
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.show()
    sys.exit(0)

if raw_input('Writing data to %s. Are you sure? [y/N] '%(fname_out))!='y':
    sys.exit(0)

wfn_out = wfnIO()
wfn_out.__dict__ = wfn_in.__dict__.copy()
wfn_out.nk = len(should_keep)
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
    if ik in should_keep:
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
