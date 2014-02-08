#!/usr/bin/env python

# Defines kpoints and kgrid classes, with support to reading kgrid.log files.

# Felipe Homrich da Jornada <jornada@berkeley.edu> - May 2013

import numpy as np


class Kpoints(np.ndarray):
    '''Stores any (3xN) ndarray.'''

    def __new__(self, kpts):
        kpts = np.asfortranarray(kpts)
        if len(kpts.shape) != 2:
            raise ValueError('kpts must be a 3xN array')
        if kpts.shape[0] != 3:
            raise ValueError('kpts must be a 3xN array')
        return kpts


class DualKpoints:
    '''Stores two sets of k-points and a mapping between them.'''

    def __init__(self, kpts_irr, kpts_full, map_f2i):
        self.kpts_irr = Kpoints(kpts_irr)
        self.nk_irr = self.kpts_irr.shape[1]
        self.kpts_full = Kpoints(kpts_full)
        self.nk_full = self.kpts_full.shape[1]
        self.map_f2i = np.asfortranarray(map_f2i)
        if len(self.map_f2i.shape) != 1:
            raise ValueError('arrays map_f2i must be a 1D array')
        if len(self.map_f2i) != self.kpts_full.shape[1]:
            raise ValueError('arrays map_f2i and kpts_full are incompatible')


class KgridKpoints(DualKpoints):
    '''Reads a kgrid.log file and stores the mapping from irreducible wedge
    to full BZ zone.'''

    def __init__(self, fname):
        def read_block(f):
            nk = int(f.readline())
            kpts = np.empty((3,nk), order='F', dtype='d')
            maps = np.empty((nk,), order='F', dtype='i')
            for i in range(nk):
                items = f.readline().split()[1:6]
                kpts[:,i] = map(float, items[:3])
                maps[i] = int(items[4]) - 1
            return nk, kpts, maps

        f = open(fname)
        nk_full = 0
        nk_irr = 0
        for line in iter(f.readline, b''):
            if 'k-points in the original uniform grid' in line:
                nk_full, kpts_full, map_full = read_block(f)
            elif 'k-points in the irreducible wedge' in line:
                nk_irr, kpts_irr, map_irr = read_block(f)

        self.full_is_irr = map_full == -1
        if self.full_is_irr.sum() != nk_irr:
            raise ValueError('irreducible grid is inconsistent')
        # First caveat:
        # map_full is "-1" wherever a particular point of the full zone is
        # exactly the same in the irreducible wedge. But we want it to point
        # index in the irr wedge.
        for i in xrange(nk_irr):
            map_full[map_irr[i]] = i
        # Second caveat:
        # map_full maps from kpts_full to kpts_full. But we want to map
        # kpts_full to kpts_irr. We have to correct the pts that were not
        # originally "-1"
        map_full[~self.full_is_irr] = map_full[map_full[~self.full_is_irr]]
        if (map_full>=nk_irr).any():
            raise ValueError('mapping from kgrid is inconsistent')

        DualKpoints.__init__(self, kpts_irr, kpts_full, map_full)


def test():
    kg = KgridKpoints('kgrid.log')
    #print kg.kpts_full
    #print kg.kpts_irr
    print kg.map_f2i
    print kg.full_is_irr


if __name__=='__main__':
    test()

