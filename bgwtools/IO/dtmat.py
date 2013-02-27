#!/usr/bin/env python

import numpy as np
from FortranIO import FortranFile

class dtmatIO:
    def __init__(self, fname):
        self.fname = fname
        self.f = FortranFile(self.fname)
        self.read_header()
        self.read_coefs(self.dcn, self.nc_co, self.nc_fi)
        self.read_coefs(self.dvn, self.nv_co, self.nv_fi)

    def read_header(self):
        f = self.f
        buf = f.read('i')
        self.nk_co = buf[0]
        self.nc_co = buf[1]
        self.nv_co = buf[2]
        self.nk_fi = buf[3]
        self.nc_fi = buf[4]
        self.nv_fi = buf[5]
        self.nspin = buf[6]
        self.kpts = np.empty((3,self.nk_co), dtype='d', order='F')
        for ik in range(self.nk_co):
            self.kpts[:, ik] = f.read('d')
        self.dcn = np.empty((self.nspin, self.nc_co, self.nc_fi, self.nk_fi),
            dtype='d', order='F')
        self.dvn = np.empty((self.nspin, self.nv_co, self.nv_fi, self.nk_fi),
            dtype='d', order='F')

    def read_coefs(self, d12, n_co, n_fi):
        f = self.f
        for jk in range(self.nk_fi):
            for jfi in range(n_fi):
                for jco in range(n_co):
                    for js in range(self.nspin):
                        buf = f.read_record()
                        ik, ifi, imap, ico, is_ = buf.read('i', 5)
                        d12[is_-1, ico-1, ifi-1, ik-1] = buf.read('d', 1)
        #norms = np.sum(d12**2, axis=1)
        #print np.max(norms), np.min(norms)

    def __repr__(self):
        return '''<dtmatIO %s>
\tFile name: %s
\tNum. coarse k-points:   %d
\tNum. coarse cond bands: %d
\tNum. coarse val  bands: %d
\tNum. fine k-points:   %d
\tNum. fine cond bands: %d
\tNum. fine val  bands: %d
\tNum. spins:      %d
\tCoarse k-points:
\t\t%s
\tBest norm for dvn coefs: %f
\tBest norm for dcn coefs: %f
\tWorst norm for dvn coefs: %f
\tWorst norm for dcn coefs: %f
</dtmatIO>
'''% (self.fname, self.fname, self.nk_co, self.nc_co, self.nv_co,
        self.nk_fi, self.nc_fi, self.nv_fi, self.nspin,
	np.array_str(self.kpts.T,50,6).replace('\n','\n\t\t'),
        np.max(np.sum(self.dvn**2, axis=1)), np.max(np.sum(self.dcn**2, axis=1)),
        np.min(np.sum(self.dvn**2, axis=1)), np.min(np.sum(self.dcn**2, axis=1)),
        )

if __name__=='__main__':
    import sys
    dtmat = dtmatIO(sys.argv[1])
    print dtmat

