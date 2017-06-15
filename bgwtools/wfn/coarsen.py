#!/usr/bin/env python

# This script reduces the number of k-points in a WFN file.

# Felipe Homrich da Jornada (Feb 2013)

import numpy as np
from numpy import *
from bgwtools.IO.wfn import wfnIO
import sys
from bgwtools.common.common import get_numpy_flavor

ryd = 13.60569253

if len(sys.argv)!=5:
    print('Usage: %s kgrid nb wfn_in wfn_out'%(sys.argv[0]))
    sys.exit(1)

kgrid = np.array(map(float, sys.argv[1].split(',')))
nb = int(sys.argv[2])
fname_in = sys.argv[3]
fname_out = sys.argv[4]

wfn_in = wfnIO(fname_in, full=True)

kpts = wfn_in.kpt[:,:]
should_keep = np.all(np.abs(np.round(kpts*kgrid[:,np.newaxis]) - kpts*kgrid[:,np.newaxis]) < 1e-6, axis=0)
should_keep = np.arange(len(should_keep))[should_keep]

if raw_input('Writing data to %s. Are you sure? [y/N] '%(fname_out))!='y':
    sys.exit(0)

if nb<1:
    nb = wfn_in.nbands
wfn_out = wfnIO()
wfn_out.__dict__ = wfn_in.__dict__.copy()
wfn_out.kgrid[:] = kgrid[:]
wfn_out.nk = len(should_keep)
wfn_out.ngk = wfn_out.ngk[should_keep]
wfn_out.ngkmax = amax(wfn_out.ngk)
wfn_out.kw = wfn_out.kw[should_keep]
wfn_out.kw[:] = 1.0/wfn_out.nk
wfn_out.kpt = wfn_out.kpt[:, should_keep]
wfn_out.ifmin = wfn_out.ifmin[should_keep, :]
wfn_out.ifmax = wfn_out.ifmax[should_keep, :]
wfn_out.energies = wfn_out.energies[:nb, should_keep, :]
wfn_out.occupations = wfn_out.occupations[:nb, should_keep, :]
wfn_out.nbands = nb
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
        sys.stdout.flush()
        wfn_in.read_gvectors(gvec)
        wfn_out.write_gvectors(gvec[:,:wfn_in.ngk[ik]])
        rec_sz = None
        for ib in range(nb):
            pos1 = file.tell(wfn_in.f)
            wfn_in.read_data(data)
            pos2 = file.tell(wfn_in.f)
            assert (rec_sz is None) or (pos2-pos1 == rec_sz)
            rec_sz = pos2 - pos1
            wfn_out.write_data(data[:wfn_in.ngk[ik],:])
        file.seek(wfn_in.f, (wfn_in.nbands-nb)*rec_sz, 1)
        #for ib in range(nb, wfn_in.nbands):
        #        wfn_in.read_data()

        if k_done==wfn_out.nk:
            break
    else:
        wfn_in.read_gvectors()
        for ib in range(wfn_in.nbands):
            wfn_in.read_data()
    
print wfn_out

print 'All done!'
