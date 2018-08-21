#!/usr/bin/env python

# Functions to read and write eqp.dat files, originally from read_eqp.py.
# Felipe H. da Jornada (2018).

from __future__ import division
import numpy as np


def read_eqp(fname):
    kpts  = []
    en_map_lda = []
    en_map_gw = []
    f = open(fname)
    imag_part = False
    for line in iter(f.readline, ''):
        # Read k-point header
        items = line.split()
        kpts.append(list(map(float, items[0:3])))
        items_cnt = int(items[3])
        en_map_lda.append({})
        en_map_gw.append({})

        # Read data associated with k-point
        klines = [f.readline() for i in range(items_cnt)]
        for kline in klines:
            items = kline.split()
            band_idx = int(items[1])
            en_lda = float(items[2])
            if len(items) > 4:
                en_gw = float(items[3]) + 1.0j*float(items[4])
                imag_part = True
            else:
                en_gw = float(items[3])
            en_map_lda[-1][band_idx] = en_lda
            en_map_gw[-1][band_idx] = en_gw

    nk = len(kpts)
    kpts = np.array(kpts)
    bands = np.sort(np.unique([list(en_k.keys()) for en_k in en_map_lda]))
    nb = len(bands)
    en_lda = np.ma.masked_all((nk,nb))
    dtype = np.complex128 if imag_part else np.float64
    en_gw = np.ma.masked_all((nk,nb), dtype=dtype)
    with np.warnings.catch_warnings():
        np.warnings.filterwarnings('ignore', category=np.ma.core.MaskedArrayFutureWarning)
        for en_out, map_in in [[en_lda, en_map_lda], [en_gw, en_map_gw]]:
            for ik, en_k in enumerate(map_in):
                indices = [np.flatnonzero(bands==idx)[0] for idx in en_k.keys()]
                en_out[ik,:][indices] = list(en_k.values())
                en_out.mask[ik,:][indices] = False

    return bands, kpts, en_lda, en_gw


def write_eqp(fname, bands, kpts, en_lda, en_gw):
    nk = kpts.shape[0]
    nb = len(bands)
    if np.any(en_lda.mask) or np.any(en_gw.mask):
        raise NotImplementedError('K-points with different occupations not currently supported.')
    with open(fname, 'w') as f:
        for ik, kk in enumerate(kpts):
            f.write(' {:12.9f} {:12.9f} {:12.9f} {:7d}\n'.format(kk[0], kk[1], kk[2], nb))
            for ib in range(nb):
                f.write(' {:7d} {:7d} {:14.9f} {:14.9f}\n'.format(
                    1, bands[ib], en_lda[ik,ib], en_gw[ik,ib]))


if __name__=='__main__':
    import sys

    bands, kpts, en_lda, en_gw = read_eqp(sys.argv[1])
    en = en_gw
    if len(sys.argv) > 2:
        if 'lda' in sys.argv[2].lower():
            en = en_lda

    #also calculate the distance between kpts
    dists = np.zeros((kpts.shape[0],))
    # Missing metric! oh well...
    dists[1:] = np.linalg.norm(kpts[1:] - kpts[:-1])
    dists = np.cumsum(dists)
    data = np.column_stack((dists, kpts[:,0], kpts[:,1], kpts[:,2], en))
    np.savetxt(sys.stdout, data, '%10.7f')
