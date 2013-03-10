#!/usr/bin/env python


from numpy import *
from FortranIO import FortranFile

class bsematIO:
    def __init__(self, fname):
        self.fname = fname
        self.f = FortranFile(self.fname)
        self.read_header()

    def read_header(self):
        f = self.f
        buf = f.read('i')
        self.nkpt = buf[0]
        self.nc = buf[1]
        self.nv = buf[2]
        self.nspin = buf[3]
        self.kpt = empty((3,self.nkpt), dtype=double, order='F')
        for ik in range(self.nkpt):
            buf = f.read_record()
            buf.read('i', 1)
            self.kpt[0:3,ik] = buf.read('d',3)

    def read_vpcpkp(self, data=None):
        if not data is None:
            rec = self.f.read_record()
            ik,ic,iv = rec.read('i',3)
            buf = rec.read('d')
            if len(buf) == 2*(self.nspin**2*self.nv*self.nc*self.nkpt):
                data[:] = buf.view(complex128).copy().reshape(
                    (self.nspin, self.nspin, self.nv, self.nc, self.nkpt), order='F')
            else:
                data[:] = buf.copy().reshape(
                    (self.nspin, self.nspin, self.nv, self.nc, self.nkpt), order='F')
            del rec,buf
            return ik,ic,iv
        else:
            self.f.next()
            return 0,0,0

    def __repr__(self):
        return '''<bsemat %s>
    File name: %s
    Num. k-points:   %d
    Num. cond bands: %d
    Num. val  bands: %d
    Num. spins:      %d
    K-pts:
        %s
</bsemat>
'''%\
        (self.fname, self.fname, self.nkpt, self.nc, self.nv, self.nspin,
        array_str(self.kpt.T,50,6).replace('\n','\n\t\t'),
                )

if __name__=='__main__':
    import sys
    bsemat = bsematIO(sys.argv[1])
    print bsemat
