#!/usr/bin/env python

from kgrid import Kpoints
import numpy as np

def _min_range_abs(x):
    '''Moves vector x to the [-0.5, 0.5) range, and apply abs().'''
    return np.fabs(x - np.floor(x + 0.5))
    #return np.fabs(x)

class Kgrids:

    def __init__(self, syms, kpts, auto_map=True):
        self.syms = syms
        self.kpts = Kpoints(kpts)
        self.nk = self.kpts.shape[1]
        self.nsyms = syms.shape[2]
        # This is the indices of S(i)*k, where S is a symmetry op. This maps
        # a k-point from the FBZ to another k-point.
        self.idx_sk = np.empty((self.nsyms, self.nk), dtype='i')
        # Index of the k-point in the irr wedge equivalent to each point in the full BZ
        self.map_f2i = np.empty(self.nk, dtype='i')
        self.sym_f2i = np.empty(self.nk, dtype='i')
        if (auto_map):
            self.get_mapping()

    def get_mapping(self):
        '''Apply symmetry ops. on all-kpoints and create irreducible wedge.'''

        # Apply symmetry ops, get idx_sk and irreducible wedge
        idx_irr = []
        for ik in xrange(self.nk):
            #print('ik %d/%d'%(ik+1,self.nk))
            tmp = np.tensordot(self.kpts[:,ik], self.syms, axes=([0], [1]))

            is_irr = True
            for iop in xrange(self.nsyms):
                dists = np.sum(_min_range_abs(self.kpts - tmp[:,iop][:,np.newaxis]), axis=0)
                idx_k = np.argmin(dists)
                self.idx_sk[iop, ik] = idx_k
                if is_irr and (idx_k in idx_irr):
                    is_irr = False
                    self.map_f2i[ik] = idx_k
                    self.sym_f2i[ik] = iop

            if is_irr: # Not found in irr wedge => append
                idx_irr.append(ik)
                self.map_f2i[ik] = ik
                self.sym_f2i[ik] = 0

        self.idx_irr = np.array(idx_irr, dtype='i')

    def __repr__(self):
        import fractions

        s = ''
        def ar_str(ar, limit=1000):
            f = [str(fractions.Fraction(l).limit_denominator(limit)).rjust(6) for l in ar]
            return '[' + ','.join(f) + ']'
            #return '[ %9.6f, %9.6f, %9.6f ]'%tuple(ar.tolist())

        for ik in xrange(len(self.map_f2i)):
            ikp = self.map_f2i[ik]
            s += '\n\t%4d = %s =[%d]=> %s = %4d'%(ik, ar_str(self.kpts[:,ik]),
                self.sym_f2i[ik], ar_str(self.kpts[:,ikp]), ikp)
        return '''<Kgrid>
size(IBZ): %d
size(FBZ): %d
Mapping: %s
</Kgrid>'''%(len(self.idx_irr), self.nk, s)


def test_bse(fname_wfn, fname_evec):
    from bgwtools.IO.wfn import wfnIO
    from bgwtools.IO.bsemat import bsematIO

    wfn = wfnIO(fname_wfn)
    bsemat = bsematIO(fname_evec)
    syms = wfn.mtrx
    cond = syms[2,2,:] == 1
    syms = syms[:,:,cond]

    bsemat.kpt[np.fabs(bsemat.kpt)<1e-15] = 0
    kg = Kgrids(syms, bsemat.kpt)

    #print [ sum(kg.idx_sk[:,ik]==ik) for ik in xrange(kg.nk) ]
    for ik in xrange(kg.nk):
        print kg.idx_sk[:,ik]
    print kg

if __name__=='__main__':
    import sys
    test_bse(*sys.argv[1:3])
