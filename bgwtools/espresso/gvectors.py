#!/usr/bin/env python

# Script to read/write G-vector files (G+k) from QE.
#
# Felipe H. da Jornada, Feb 2015.

from __future__ import print_function
import numpy as np
from bgwtools.IO.qe_xml import IotkFile
from scipy.spatial import cKDTree
import re


class GvectorsIO(object):

    def __init__(self, fname=None, verbose=False):
        if fname is not None:
            self.read(fname, verbose)

    def read(self, fname, verbose=False):
        re_tag = re.compile(r'''([^\s]*)=['"]([^'"]+)['"]''')
        re_name = re.compile(r'''<([^\s]*)''')
        f = IotkFile(fname)
        print('Reading G-vectors from '+fname)

        def my_print(*args):
            if verbose:
                print(*args)

        first = True
        while True:
            # Read all tags
            #try:
            if first:
                obj = f.read_field(force_rec_len=4)
            elif name=='g':
                print('!')
                obj = f.read_field(force_rec_len=1446030*4 + 4)
                #obj = f.read_field(force_rec_len=4)
            else:
                obj = f.read_field()
            first = False
            #except IOError:
            #    break

            if obj[0]==0:
                # This is a text field
                #print obj[1], repr(obj[2])
                tag = obj[2]
                name = re_name.search(tag).group(1)
                attrs = dict(re_tag.findall(tag))
            my_print(obj[0], tag)

            if obj[0]==1:
                if name=='NUMBER_OF_GK-VECTORS':
                    self.ngk = obj[2][0]
                    my_print(self.ngk)
                if name=='MAX_NUMBER_OF_GK-VECTORS':
                    self.ngk_max = obj[2][0]
                    my_print(self.ngk_max)
                if name=='GAMMA_ONLY':
                    self.gamma_only = obj[2][0].astype(bool)
                    my_print(self.gamma_only)
                if name=='K-POINT_COORDS':
                    self.kpt = obj[2]
                    my_print(self.kpt)
                if name=='INDEX':
                    self.index = obj[2]
                    assert np.all(self.index==np.arange(1, len(self.index)+1))
                if name=='GRID':
                    self.gvecs = obj[2].reshape((3,-1), order='F')
                    my_print(self.gvecs.T)

        print('Done reading G-vectors from file '+fname)
        print('')

    def write(self, fname):
        assert self.ngk_max>=self.ngk
        assert self.ngk==len(self.index)
        assert self.ngk==self.gvecs.shape[1]
        file = IotkFile(fname, mode='wb')
        f = file.write_field
        f((0,5,'<?iotk version="1.2.0"?>'))
        f((0,5,'<?iotk file_version="1.0"?>'))
        f((0,5,'<?iotk binary="T"?>'))
        f((0,5,'<?iotk qe_syntax="F"?>'))
        f((0,1,'<GK-VECTORS>'))
        f((0,1,'<NUMBER_OF_GK-VECTORS type="integer" size="1" kind="4">'))
        f((1,'i',self.ngk))
        f((0,2,'</NUMBER_OF_GK-VECTORS>'))
        f((0,1,'<MAX_NUMBER_OF_GK-VECTORS type="integer" size="1" kind="4">'))
        f((1,'i',self.ngk_max))
        f((0,2,'</MAX_NUMBER_OF_GK-VECTORS>'))
        f((0,1,'<GAMMA_ONLY type="logical" size="1" kind="4">'))
        f((1,'i',int(self.gamma_only)))
        f((0,2,'</GAMMA_ONLY>'))
        f((0,1,'<K-POINT_COORDS type="real" size="3" kind="8" UNITS="2 pi / a">'))
        f((1,'d',self.kpt.ravel(order='K')))
        f((0,2,'</K-POINT_COORDS>'))
        f((0,1,'<INDEX type="integer" size="{}" kind="4">'.format(self.ngk)))
        f((1,'i',self.index))
        f((0,2,'</INDEX>'))
        f((0,1,'<GRID type="integer" size="{}" kind="4" columns="3">'.format(3*self.ngk)))
        f((1,'i',self.gvecs.ravel(order='K')))
        f((0,2,'</GRID>'))
        f((0,2,'</GK-VECTORS>'))

    def copy(self):
        gvectors = GvectorsIO()
        for key, val in self.__dict__.iteritems():
            if isinstance(val, np.ndarray):
                val = val.copy()
            setattr(gvectors, key, val)
        return gvectors

    def unfold_trs(self):
        assert self.gamma_only
        gvectors = self.copy()
        gvecs = gvectors.gvecs
        assert np.all(gvecs[:,0]==0), 'First k-point is not Gamma!'
        tree = cKDTree(gvecs[:,1:].T)
        dist, ind = tree.query(-gvecs[:,1:].T)
        assert np.all(dist>0.1), 'List of G-vectors was not folded with time-reversal symmetry!'
        gvecs_k = np.empty((3,2*gvectors.ngk-1), dtype=gvecs.dtype, order='F')
        gvecs_k[:,0] = 0
        gvecs_k[:,1::2] = gvecs[:,1:]
        gvecs_k[:,2::2] = -gvecs[:,1:]

        gvectors.ngk = 2*gvectors.ngk - 1
        gvectors.ngk_max = 2*gvectors.ngk_max - 1
        gvectors.index = np.arange(1, gvectors.ngk+1)
        gvectors.gvecs = gvecs_k
        gvectors.gamma_only = False
        return gvectors


if __name__=="__main__":
    import sys

    if len(sys.argv)<2:
        print('usage: {} gkvectors.dat'.format(sys.argv[0]))
        print()
        sys.exit(1)
    fname = sys.argv[1]

    gvectors = GvectorsIO(fname, verbose=True)
    gvectors.write('gkvectors_tmp.dat')
