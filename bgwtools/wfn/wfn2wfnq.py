#!/usr/bin/env python

# This script removes all bands from the WFN except for the ones needed
# for WFNq (i.e., up to max(ifmax) + 1)

# Felipe Homrich da Jornada (Feb 2013)

import numpy as np
from bgwtools.IO.wfn import wfnIO
import sys
from bgwtools.common.common import get_numpy_flavor

if len(sys.argv)!=3:
    print('Usage: %s wfn_in wfn_out'%(sys.argv[0]))
    sys.exit(1)

fname_in = sys.argv[1]
fname_out = sys.argv[2]

wfn_in = wfnIO(fname_in, full=True)
nb_out = np.amax(wfn_in.ifmax)+1

wfn_out = wfnIO()
wfn_out.__dict__ = wfn_in.__dict__.copy()
wfn_out.energies = wfn_out.energies[:nb_out, :, :]
wfn_out.occupations = wfn_out.occupations[:nb_out, :, :]
wfn_out.f = None
wfn_out.fname = fname_out
wfn_out.nbands = nb_out
wfn_out.write_header(full=True)

ng_max = np.amax(wfn_in.ngk)
gvec = np.empty((3, ng_max), dtype=int, order='F')
data = np.empty((ng_max, wfn_in.ns), dtype=get_numpy_flavor(wfn_in.flavor), order='F')
k_done = 0
for ik in range(wfn_in.nk):
        k_done += 1
        print k_done,'/',wfn_out.nk
        wfn_in.read_gvectors(gvec)
        wfn_out.write_gvectors(gvec[:,:wfn_in.ngk[ik]])
        for ib in range(nb_out):
            wfn_in.read_data(data)
            wfn_out.write_data(data[:wfn_in.ngk[ik],:])
        for ib in range(nb_out,wfn_in.nbands):
            wfn_in.read_data()
    
print wfn_out

print 'All done!'
