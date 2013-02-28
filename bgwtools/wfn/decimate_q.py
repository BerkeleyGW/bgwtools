#!/usr/bin/env python

# This script creates a WFNq_fi with a reduced number of k-points given a full
# WFN_fi, a reduced/decimated WFN_fi, and a q-shift.

# Felipe Homrich da Jornada (Feb 2013)

from numpy import *
from bgwtools.IO.wfn import wfnIO
import sys
from bgwtools.common.common import get_numpy_flavor

ryd = 13.60569253
TOL = 1e-6

if len(sys.argv)!=5:
    print('Usage: %s wfn_full wfn_reduced wfn_out qx,qy,qz'%(sys.argv[0]))
    sys.exit(1)

fname_in = sys.argv[1]
fname_red = sys.argv[2]
fname_out = sys.argv[3]
q_shift = array(map(float,sys.argv[4].split(',')))
print q_shift

wfn_in = wfnIO(fname_in)
kpts_full = wfn_in.kpt

wfn_red = wfnIO(fname_red)
kpts_red = wfn_red.kpt

kpts_q = kpts_red + q_shift[:,newaxis]

should_keep = []
for jk in range(wfn_in.nk):
    kpt = kpts_full[:,jk]
    if any(all(abs(kpt[:,newaxis]-kpts_q)<TOL, axis=0)):
        should_keep += [jk]

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
