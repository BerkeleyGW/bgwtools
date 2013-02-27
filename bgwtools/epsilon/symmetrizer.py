#!/usr/bin/env python

# Make the dielectric matrix Hermitian

from bgwtools.IO.epsmat import epsmatIO
from bgwtools.IO.wfn import wfnIO
import sys
import numpy as np

class symmetrizer:
    #'''Get epsillon at q=0 Symmetrized eps_{G1,G2}(q)
    '''Symmetrizes the Coulomb interaction for an epsilon file'''
    def __init__(self, wfn, epsmat):
        self.epsmat = epsmat
        self.wfn = wfn
        
    def symmetrize(self):
        pi = np.pi
        for iq in xrange(len(self.epsmat.epsmat)):
            mat = self.epsmat.epsmat[iq]
            nmtx = len(mat)
	    gidx = self.epsmat.isort[iq][:nmtx] - 1
            ekin0 = np.array(self.epsmat.ekin)[iq][gidx]
            #print 1/ekin*(8*pi)
	    qpt = self.epsmat.qpt[:,iq]
	    gvec_q = self.epsmat.gvec_k[:,gidx] + qpt[:,np.newaxis]
            M = np.linalg.cholesky(self.wfn.bdot).T
	    ekin = np.sum(np.dot(M,gvec_q)**2, axis=0)
            ekin0[0] = ekin[0]

            print 1/ekin*(8*pi)
            print 1/ekin0*(8*pi)

            # TODO: don`t touch the wings!
            for ig in xrange(nmtx):
                mat[ig,ig] -= 1.0
                mat[:,ig] *= ekin
            #mat = 0.5*(mat + np.conj(mat.T))
            for ig in xrange(nmtx):
                mat[:,ig] /= ekin0
                mat[ig,ig] += 1.0

            self.epsmat.epsmat[iq] = mat

if __name__=='__main__':
    import sys

    if len(sys.argv)!=4:
        print 'Usage: %s wfn epsmat_in epsmat_out'%(sys.argv[0])
        exit(1)

    f_wfn = sys.argv[1]
    f_in = sys.argv[2]
    f_out = sys.argv[3]

    wfn = wfnIO(f_wfn)
    eps = epsmatIO(f_in)
    sym = symmetrizer(wfn, eps)
    sym.symmetrize()
    eps.to_file(f_out)
