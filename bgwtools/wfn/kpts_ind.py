#!/usr/bin/env python

# Given two wfns, this script gives the indices of the points (in C order) from
# WFN_full which are present in WFN_sub

from bgwtools.IO.wfn import wfnIO
import numpy as np
from scipy.spatial import cKDTree

def compare_kpts(fname_full, fname_sub, tol=1e-8):
    wfn_full = wfnIO(fname_full)
    wfn_sub = wfnIO(fname_sub)
    kpts_full = (wfn_full.kpt + tol) % 1
    kpts_sub = (wfn_sub.kpt + tol) % 1
    tree = cKDTree(kpts_full.T)
    d, ind = tree.query(kpts_sub.T)
    return ind, kpts_full.shape[1]

if __name__=='__main__':
    import sys
    if len(sys.argv) != 3:
        print('usage: %s WFN_full WFN_sub')
        sys.exit(1)
    fname_full, fname_sub = sys.argv[1:]
    print compare_kpts(fname_full, fname_sub)
