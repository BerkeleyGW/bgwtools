#!/usr/bin/env python

# Given a kgrid.log file output from kgrid.x, this script defines a class
# with information about the symmetry and full and irreducible k-grids.
# Diana Y. Qiu (2014)

import numpy as np

class kgridlogIO:
    def __init__(self, fname):
        
        f = open(fname)

        f.readline()
        lines = f.readline().split()
        self.kgrid = np.array([int(lines[2]),int(lines[3]),int(lines[4])])
        lines = f.readline().split()
        self.kshift = np.array([float(lines[2]),float(lines[3]),float(lines[4])])
        f.readline()
        f.readline()

        lines = [f.readline() for ii in range(3)]
        info = np.genfromtxt(lines, dtype=None)
        self.avec = info[['f2','f3','f4']].view(np.float64).reshape((3,3), order='F')
        f.readline()
        self.vol = np.float64(f.readline().split()[4])
        f.readline()
        f.readline()
        lines = [f.readline() for ii in range(3)]
        info = np.genfromtxt(lines, dtype=None)
        self.bvec = info[['f2','f3','f4']].view(np.float64).reshape((3,3), order='F')
        f.readline()
        f.readline()

        self.nsymm_bravais = int(f.readline())
        f.readline()
        lines = [f.readline() for ii in range(self.nsymm_bravais)]
        info = np.genfromtxt(lines,dtype=None)
        self.mtrx_bravais = info[['f2','f3','f4','f5','f6','f7','f8','f9','f10']].view(np.int).reshape((3,3,-1), order='F')
        f.readline()
        f.readline()
        

        self.nsymm = int(f.readline())
        f.readline()
        lines = [f.readline() for ii in range(self.nsymm)]
        info = np.genfromtxt(lines,dtype=None)
        self.mtrx = info[['f2','f3','f4','f5','f6','f7','f8','f9','f10']].view(np.int).reshape((3,3,-1), order='F')
        f.readline()

        line = f.readline().split()
        self.fft = np.array([int(line[2]),int(line[3]),int(line[4])])
        f.readline()
        f.readline()
        self.nsymm_fft = int(f.readline())
        f.readline()
        lines = [f.readline() for ii in range(self.nsymm)]
        info = np.genfromtxt(lines,dtype=None)
        self.mtrx_fft = info[['f2','f3','f4','f5','f6','f7','f8','f9','f10']].view(np.int).reshape((3,3,-1), order='F')
        f.readline()
        f.readline()

        # Symmetries that reduce the k-points
        line = f.readline().split()
        self.trans = []
        while line!=[]:            
            self.trans.append(int(line[0][1:]))
            #print self.trans
            line = f.readline().split()
        self.trans = np.array(self.trans)

        # k-points in the original uniform grid
        f.readline()
        self.nf = int(f.readline())
        lines = [f.readline() for ik in range(self.nf)]
        info = np.genfromtxt(lines, dtype=None)
        self.fk = info[['f1','f2','f3']].view(np.float64).reshape((3,-1), order='F')
        self.indr = info['f5']
        # Replace 0's with index
        idx_zero = np.where(self.indr==0)[0]
        idx_nonzero = np.where(self.indr!=0)[0]
        self.indr[idx_zero] = np.arange(len(idx_zero))+1
        self.indr[idx_nonzero] = self.indr[self.indr[idx_nonzero]-1]
        self.itran = info['f6']
        idx_zero = np.where(self.itran=='---')[0]
        self.itran = np.delete(self.itran,idx_zero)
        f.readline()
        f.readline()
        # k-points in the irreducible wedge
        self.nrk = int(f.readline())
        lines =[f.readline() for ik in range(self.nrk)]
        info = np.genfromtxt(lines, dtype=None)
        self.rk = info[['f1','f2','f3']].view(np.float64).reshape((3,-1), order='F')
        self.weight = info['f4']
        
        f.close()

    def __repr__(self):
        return '''
    nf: {nf}
    fk:
{fk.T}
    itran:
{itran}
    indr:
{indr}
    nsymm: {nsymm}
    mtrx:
{mtrx}
    nrk: {nrk}
    rk:
{rk.T}
        '''.format(**self.__dict__)

if __name__ == '__main__':
    import sys
    klog = kgridlogIO(sys.argv[1], sys.argv[2])

