#!/usr/bin/env python

# A class for reading the fullbz.dat file from analyzebz.x 
# This is pretty much the same as rotator.py from ipsa without reading the wave function.

import numpy as np

class fullbzIO:

    def __init__(self, fname_bz):

        f = open(fname_bz)

        self.nf = int(f.readline())
        print self.nf
        lines = [f.readline() for ik in range(self.nf)]
        info = np.genfromtxt(lines, dtype=None)
        if self.nf==1:
            self.fk= np.array([info[['f0','f1','f2']]]).view(np.float64).reshape((3,-1), order='F')
            self.itran = np.array([info['f3']])
            self.indr = np.array([info['f4']])
        else:
            self.fk = info[['f0','f1','f2']].view(np.float64).reshape((3,-1), order='F')
            self.itran = info['f3']
            self.indr = info['f4']

        self.ntran = int(f.readline())
        print self.ntran
        lines = [f.readline() for ik in range(self.ntran)]
        info = np.genfromtxt(lines, dtype=int)
        self.mtrx= info.reshape((self.ntran,3,3), order='F')
        self.mtrx_i= np.zeros(np.shape(self.mtrx),dtype=np.float64)
        for itran in range(self.ntran):
            self.mtrx_i[itran,:,:] = np.linalg.inv(self.mtrx[itran,:,:])
        line=f.readline()
        lines=[]
        while len(line.split())!=1:
            line = np.genfromtxt(line.split(),dtype=np.float64)
            #line=line.reshape((-1,3))
            lines.append(line)
            line=f.readline()
        lines=np.array(lines)
        self.tnp = lines
        assert np.all(np.abs(self.tnp)<1e-12)

        self.nrk = int(line)
        lines = [f.readline() for ik in range(self.nrk)]
        if self.nrk==1:
            self.rk = np.array([np.genfromtxt(lines, dtype=float)]).T
        else:
            self.rk = np.genfromtxt(lines, dtype=float).T
        
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
    ntran: {ntran}
    mtrx:
{mtrx}
    mtrx_i:
{mtrx_i}
    tnp:
{tnp}
    nrk: {nrk}
    rk:
{rk.T}
        '''.format(**self.__dict__)

if __name__ == '__main__':
    import sys
    fullbz = fullbzIO(sys.argv[1])

    print fullbz
