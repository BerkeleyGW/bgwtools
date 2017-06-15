#!/usr/bin/env python

# Script to read/write eigenvectors files from QE.
#
# Felipe H. da Jornada, Feb 2015.

from __future__ import print_function
import numpy as np
from bgwtools.IO.qe_xml import IotkFile
import re

class EigenvectorsIO(object):

    def __init__(self, fname=None, verbose=False):
        if fname is not None:
            self.read(fname, verbose)

    def read(self, fname, verbose=False):
        re_tag = re.compile(r'''([^\s]*)=['"]([^'"]+)['"]''')
        re_name = re.compile(r'''<([^\s]*)''')
        f = IotkFile(fname)
        print('Reading eigenvectors from '+fname)

        def my_print(*args):
            if verbose:
                print(*args)

        while True:
            # Read all tags
            try:
                obj = f.read_field()
            except IOError:
                break

            if obj[0]==0:
                # This is a text field
                #print obj[1], repr(obj[2])
                tag = obj[2]
                name = re_name.search(tag).group(1)
                attrs = dict(re_tag.findall(tag))
            my_print(obj[0], tag)

            if name=='INFO':
                for tag in ('ngw', 'igwx', 'nbnd', 'ik', 'nk', 'ispin', 'nspin'):
                    val = int(attrs[tag])
                    setattr(self, tag, val)
                    my_print(tag, val)
                assert self.ngw==self.igwx
                self.gamma_only = {'F':False, 'T':True}[attrs['gamma_only']]
                self.scale_factor = float(attrs['scale_factor'])
                assert np.allclose(self.scale_factor, 1.)
                self.evcs = np.empty((self.ngw, self.nbnd), dtype=np.complex128, order='F')
            if obj[0]==1 and name[:4]=='evc.':
                    # This is a binary field containing a slab of the charge density
                    ib = int(name.split('.')[-1])
                    size = int(attrs['size'])
                    assert size==self.ngw
                    data = obj[2]
                    assert len(data)==size
                    self.evcs[:,ib-1] = data

        print('Done reading eigenvectors from file '+fname)
        print('')


    def write(self, fname):
        assert self.ngw==self.igwx
        assert self.evcs.shape[0] == self.ngw
        assert self.evcs.shape[1] == self.nbnd
        file = IotkFile(fname, mode='wb')
        f = file.write_field
        f((0,5,'<?iotk version="1.2.0"?>'))
        f((0,5,'<?iotk file_version="1.0"?>'))
        f((0,5,'<?iotk binary="T"?>'))
        f((0,5,'<?iotk qe_syntax="F"?>'))
        f((0,1,'<WFC>'))
        scale_factor = '1.000000000000000E+000'
        if np.fabs(self.scale_factor-1)>1e-12:
            scale_factor = '{:19.14E}'.format(self.scale_factor)
        f((0,3,'<INFO ngw="{}" igwx="{}" gamma_only="{}" nbnd="{}" '
               'ik="{}" nk="{}" ispin="{}" nspin="{}" scale_factor="{}"/>'.format(
               self.ngw, self.igwx, {False:'F', True:'T'}[self.gamma_only], self.nbnd,
               self.ik, self.nk, self.ispin, self.nspin, scale_factor)))
        for ib in range(self.nbnd):
            buf = self.evcs[...,ib].ravel(order='K').view(np.float64)
            f((0,1,'<evc.{} type="complex" size="{}" kind="8">'.format(ib+1,self.ngw)))
            f((1,'d',buf))
            f((0,2,'</evc.{}>'.format(ib+1)))
        f((0,2,'</WFC>'))

    def copy(self):
        evecs = EigenvectorsIO()
        for key, val in self.__dict__.iteritems():
            if isinstance(val, np.ndarray):
                val = val.copy()
            setattr(evecs, key, val)
        return evecs

    def unfold_trs(self):
        assert self.gamma_only
        evecs = self.copy()
        evcs = evecs.evcs
        evcs_k = np.empty((2*evecs.ngw-1,evecs.nbnd), dtype=evcs.dtype, order='F')
        evcs_k[0,:] = evcs[0]
        evcs_k[1::2,:] = evcs[1:]
        evcs_k[2::2,:] = evcs[1:].conj()
        norms_k = np.linalg.norm(evcs_k, axis=0)
        assert np.allclose(norms_k, 1)

        evecs.ngw = 2*evecs.ngw - 1
        evecs.igwx = 2*evecs.igwx - 1
        evecs.evcs = evcs_k
        evecs.gamma_only = False
        return evecs

if __name__=="__main__":
    import sys

    if len(sys.argv)<2:
        print('usage: {} evcs.dat'.format(sys.argv[0]))
        print()
        sys.exit(1)
    fname = sys.argv[1]

    evecs = EigenvectorsIO(fname, verbose=True)
    norms = np.linalg.norm(evecs.evcs, axis=0)
    print(norms)
    evecs.write('evc_tmp.dat')
