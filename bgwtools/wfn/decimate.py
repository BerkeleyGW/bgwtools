#!/usr/bin/env python

# This script reduces the number of k-points in a WFN file. It only keeps
# k-points between the Fermi energy + Edop +/- Edelta, where Edop and Edelta
# are input parameters.
# Use this script to generate a decimated WFN_fi file to use with BSE with the
# patched_sampling option. Then, use decimate_q.py to generate WFNq_fi.

# Felipe Homrich da Jornada (Feb 2013)

from numpy import *
from bgwtools.IO.wfn import wfnIO
import sys
from bgwtools.common.common import get_numpy_flavor

ryd = 13.60569253

if len(sys.argv)!=5:
    print('Usage: %s Edop Edelta wfn_in wfn_out'%(sys.argv[0]))
    sys.exit(1)

Edop = float(sys.argv[1])/ryd
Edelta = float(sys.argv[2])/ryd
fname_in = sys.argv[3]
fname_out = sys.argv[4]

wfn_in = wfnIO(fname_in, full=True)
FE = -inf
for jk in range(wfn_in.nk):
    for js in range(wfn_in.ns):
        FE = max(FE, amax(wfn_in.energies[wfn_in.ifmax[jk,js]-1,jk,js]))

FE = FE + Edop
should_keep = []
for jk in range(wfn_in.nk):
    if any((wfn_in.energies[:,jk,:]<FE+Edelta)&(wfn_in.energies[:,jk,:]>FE-Edelta)):
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
