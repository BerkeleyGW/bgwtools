#!/usr/bin/env python


from numpy import *
from FortranIO import FortranFile

class bsematIO:
    def __init__(self, fname=None):
        self.fname = fname
        self.f = None
        self.nkpt = 0
        self.nc = 0
        self.nv = 0
        self.nspin = 0
        self.kpt = empty((3,0), order='F')
        if fname:
            self.read_header()

    def read_header(self):
        if not self.f:
            self.f = FortranFile(self.fname)
        else:
            self.f.seek(0)
        f = self.f

        #READ
        buf = f.read('i')
        self.nkpt = buf[0]
        self.nc = buf[1]
        self.nv = buf[2]
        self.nspin = buf[3]
        self.kpt = empty((3,self.nkpt), dtype=double, order='F')

        for ik in range(self.nkpt):
            #READ
            buf = f.read_record()
            buf.read('i', 1)
            self.kpt[0:3,ik] = buf.read('d',3)

    def write_header(self):
        if not self.f:
            self.f = FortranFile(self.fname,'=','i','wb')
        else:
            self.f.seek(0)
        f = self.f

        #WRITE
        fmt = 'i'*4
        data = [self.nkpt, self.nc, self.nv, self.nspin]
        f.write_vals(fmt, *data)

        for ik in range(self.nkpt):
            #WRITE
            fmt = 'iddd'
            data = [ik+1] + self.kpt[:,ik].tolist()
            f.write_vals(fmt, *data)

    def read_vpcpkp(self, data=None):
        f = self.f
        if not data is None:
            rec = f.read_record()
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

    def write_vpcpkp(self, ik, ic, iv, data):
        f = self.f

        buf = data.view(float64).copy().ravel(order='F')
        fmt = 'i'*3 + 'd'*len(buf)
        f.write_vals(fmt, *([ik, ic, iv] + buf.tolist()))

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
