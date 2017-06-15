#!/usr/bin/env python

# This script calculates the phase of the G=0 component of the WFNs
# for all kpts and bands. This is useful when analyzing results from
# BSE calculations.

from bgwtools.IO.wfn import wfnIO
import sys
import numpy as np

fname_wfn = sys.argv[1]
fname_out = sys.argv[2]
wfn = wfnIO(fname_wfn)

gvecs = wfn.get_gvectors_bufer()
data = wfn.get_data_buffer()
heads = np.empty((wfn.nk,wfn.nbands), dtype='complex128')

print wfn.gvec
#exit()

for ik in xrange(wfn.nk):
    wfn.read_gvectors(gvecs)
    cond = np.all(gvecs[:,:wfn.ngk[ik]]==0, axis=0)
    idx = np.nonzero(cond)[0][0]
    print ik, idx, gvecs[:,idx]
    for ib in xrange(wfn.nbands):
        wfn.read_data(data)
        heads[ik,ib] = data[idx,0]

np.save(fname_out, heads)
