#!/usr/bin/env python

# Decimates a WFN file by keeping only a subset of k-points.

# Felipe Homrich da Jornada (Aug 2013)

import numpy as np
from bgwtools.IO.wfn import wfnIO
from bgwtools.common.common import get_numpy_flavor


class Decimator:

    def __init__(self, fname_in):
        self.wfn_in = wfnIO(fname_in)
        self.keep_k = np.ones(self.wfn_in.nk, dtype=bool)

    def setup_kpts(self, arg):
        arg = np.asanyarray(arg).ravel('A')
        if arg.dtype==bool:
            if len(arg) != self.wfn_in.nk:
                raise TypeError('Incorrect dimension for boolean array arg')
            self.keep_k[:] = arg
        elif arg.dtype==int:
            if len(arg) > self.wfn_in.nk:
                raise TypeError('Incorrect dimension for integer array arg')
            self.keep_k[:] = [ik in arg for ik in range(self.wfn_in.nk)]
        print('Original number of k-points: %d'%(self.wfn_in.nk))
        print('Target number of k-points: %d'%(np.sum(self.keep_k)))

    def decimate_to(self, fname_out, verbose=True):
        wfn_in = self.wfn_in
        keep_k = self.keep_k
        nk_out = np.sum(keep_k)

        wfn_out = wfnIO()
        wfn_out.__dict__ = wfn_in.__dict__.copy()
        wfn_out.nk = nk_out
        wfn_out.ngk = wfn_out.ngk[keep_k]
        wfn_out.ngkmax = np.amax(wfn_out.ngk)
        wfn_out.kw = wfn_out.kw[keep_k]
        wfn_out.kw[:] = 1.0/wfn_out.nk
        wfn_out.kpt = wfn_out.kpt[:, keep_k]
        wfn_out.ifmin = wfn_out.ifmin[keep_k, :]
        wfn_out.ifmax = wfn_out.ifmax[keep_k, :]
        wfn_out.energies = wfn_out.energies[:, keep_k, :]
        wfn_out.occupations = wfn_out.occupations[:, keep_k, :]
        wfn_out.f = None
        wfn_out.fname = fname_out
        wfn_out.write_header(full=True)

        ng_max = np.amax(wfn_in.ngk)
        gvec = np.empty((3, ng_max), dtype=int, order='F')
        data = np.empty((ng_max, wfn_in.ns), dtype=get_numpy_flavor(wfn_in.flavor), order='F')
        k_done = 0
        for ik in range(wfn_in.nk):
            if keep_k[ik]:
                k_done += 1
                if verbose:
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
