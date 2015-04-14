#!/usr/bin/env python

import numpy as np
from bgwtools.IO.wfn import wfnIO
from bgwtools.IO.eigenvectors import eigenvectorsIO

ryd = 13.60535

def read_evecs(fname_evecs, N=None, nv=None, nc=None):
    evecs = eigenvectorsIO()
    evecs.fname = fname_evecs
    evecs.read_header(0)
    #print evecs.evecs.shape
    if N is None:
        N = evecs.nk*evecs.ns*evecs.nc*evecs.nv
    if nv is None:
        nv = evecs.nv
    if nc is None:
        nc = evecs.nc
    evs = np.empty((N,evecs.nk,nc,nv), order='C', dtype='complex128')
    evals = np.empty((N))
    f = evecs.f
    for i in range(N):
        print i
        print evecs.nk, evecs.nc, evecs.nv
        evals[i] = f.read('d')
        evs[i,:,:,:] = f.read('d').view('complex128').reshape((evecs.nk,evecs.nc,evecs.nv), order='C')[:,:nc,:nv]
    return evecs.kpt, evs, evals

def read_eqp(fname_wfn, fname_eqp, nv, nc):
    wfn = wfnIO(fname_wfn)
    n_occ = wfn.ifmax[0,0]

    if fname_eqp is None:
        #nv, nk
        en_v = wfn.energies[n_occ-nv:n_occ,:,0][::-1]*ryd
        #nc, nk
        en_c = wfn.energies[n_occ:n_occ+nc,:,0]*ryd
    else:
        from bgwtools.converters.read_eqp import get_data_from_eqp
        bands,kx,ky,kz,e_lda,e_gw = get_data_from_eqp(fname_eqp)
        bands_range = np.arange(len(bands))
        val_bands = bands_range[(bands<=n_occ) & (bands>n_occ-nv)][::-1]
        cond_bands = bands_range[(bands>n_occ) & (bands<n_occ+nc+1)]
        en_v = e_gw[:,val_bands].T
        en_c = e_gw[:,cond_bands].T
    return en_v, en_c

