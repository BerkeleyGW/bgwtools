#!/usr/bin/env python
#coding: utf-8

from __future__ import print_function
import numpy as np
import xml.etree.cElementTree as ET
from bgwtools.espresso.rho_io import read_rho
import os
bohr = 0.52917721092


def find_min_cell(fname_datafile, frac, fname_chargedensity=None, verbose=False):
    tree = ET.parse(fname_datafile)
    # order is (vector index, cartesian component)
    avecs = np.array([np.fromstring(
            tree.find('CELL/DIRECT_LATTICE_VECTORS/a{}'.format(i)).text.strip(),
            sep=' ')*bohr for i in (1,2,3)])
    alens = np.linalg.norm(avecs, axis=1)
    if verbose:
        print('Length of lattice vectors:')
        for i in range(3):
            print('a{} = {:.3f} Å'.format(i+1, alens[i]))
        print()

    if fname_chargedensity is None:
        t = tree.find('CHARGE-DENSITY')
        fname_chargedensity = t.get('iotk_link')
        fname_chargedensity = os.path.join(os.path.dirname(fname_datafile),
            fname_chargedensity)

    def min_range(x, L):
        eps = 1e-6
        y = x/L + eps
        y = y - np.floor(y + 0.5)
        return (y - eps)*L

    rho = read_rho(fname_chargedensity, verbose=False)
    rho /= np.sum(rho)
    rho_flatten = rho.flatten(order='F')
    order = np.argsort(-rho_flatten)
    order_inv = np.empty_like(order)
    order_inv[np.arange(len(order))] = order
    rho_cum = np.cumsum(rho_flatten[order])
    grids = [min_range(np.arange(rho.shape[i])*alens[i]/rho.shape[i], alens[i])
            for i in (0,1,2)]

    def get_inside(frac):
        cond = rho_cum<=frac
        ind_inside = np.flatnonzero(cond)
        ind_inside = order_inv[ind_inside]
        return ind_inside

    ind_inside = get_inside(frac)
    if verbose:
        print('Fraction of points inside region: {}'.format(
            len(ind_inside)/float(len(rho_cum))))
        print('Bounding box:')
    ipts_dims = np.unravel_index(ind_inside, rho.shape, order='F')
    alens_frac = np.empty((3,))
    for i in range(3):
        pts = grids[i][ipts_dims[i]]
        # Add delta_alens to each delta, otherwise if we compute the bounding
        # box that contains 100% of RHO, we will get a smaller box than the
        # original box
        delta = np.amax(pts) - np.amin(pts) + alens[i]/rho.shape[i]
        if verbose:
            print('a{}_min, a{}_max, delta_a{} = {:6.3f}, {:6.3f}, {:6.3f} Å'.format(
                   i+1, i+1, i+1, np.amin(pts), np.amax(pts), delta))
        alens_frac[i] = delta*2.
    cell = avecs * alens_frac[:,None] / alens[:,None]
    if verbose:
        print()
    return cell


if __name__=="__main__":
    import os
    import argparse
    parser = argparse.ArgumentParser(description='Finite box calculator')
    parser.add_argument('datafile', help='data-file.xml file')
    parser.add_argument('--frac', type=float, default=0.999,
        help=('Fraction of the charge density that should be encompased in '
            '1/8th of the unit cell'))
    parser.add_argument('--chargedensity', default=None,
        help='charge-density.dat file, defaults to file linked by data-file.xml')
    parser.add_argument('--verbose', default=False, action='store_true')
    args = parser.parse_args()
    
    cell = find_min_cell(args.datafile, args.frac, args.chargedensity, args.verbose)
    print('CELL_PARAMETERS angstrom')
    for i in range(3):
        print(' {:10.6f} {:10.6f} {:10.6f}'.format(*cell[i]))
