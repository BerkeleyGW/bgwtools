#!/usr/bin/env python

# An abstract class that reads/writes epsmat files.
# Handles REAL and COMPLEX data, but only static calculations.

# Felipe Homrich da Jornada <jornada@civet.berkeley.edu> (2012)

from __future__ import division
from bgwtools.common import common
import FortranIO
import numpy as np
from numpy import *
import numpy.core as core

class epsmatIO:
    def __init__(self, fname=None, auto_read=True, read_all=True, flavor=common.flavor.NONE):
        self.fname=fname
        self.f=None
        self.cur_q=-1 #index of current q_pt (0-based)

        self.name=''
        self.date=''
        self.freq_dep=0
        self.num_freq=1
        self.grid=array([1,1,1])
        self.ecuts=0.
        self.nq=0
        self.qpt=empty((3,0), dtype=float64)
        self.ng=0
        self.gvec_k=empty((3,0), dtype=int)
        self.nmtx=empty(0, dtype=int)
        self.isort = []
        self.isort_i = []
        self.ekin = []
        self.epsmat = []
        self.q = empty((3,0), dtype=float64)
        self.flavor = flavor

        if fname and auto_read:
            self.from_file(self.fname, read_all)

    def read_header(self):
        if not self.f:
            self.f = FortranIO.FortranFile(self.fname)
        else:
            self.f.seek(0)
        f = self.f

        #READ
        rec = f.readline()
        self.name=rec[:6]
        self.date=rec[6:6+11]

        #READ
        self.freq_dep, self.num_freq = f.read('i')

        #READ
        self.grid = f.read('i') #k-grid (3)

        #READ
        if self.freq_dep==2:
            tmp = f.read('d')
            assert(len(tmp)==3*self.num_freq)
            self.freq_grid = np.array(tmp[:self.num_freq])
            self.freq_broad = np.array(tmp[self.num_freq:]).view(np.complex128)
        else:
            f.read()

        #READ
        f.read() #Empty by GSM
        #READ
        f.read() #Empty by GSM

        #READ
        self.ecuts = f.read('d')[0] #maximum energy, or gmax_in

        #READ
        buf = f.read_record()
        self.nq = buf.read('i',1) #number of q-pts in this file
        self.qpt = buf.read('d') #q-pts
        self.qpt = self.qpt.reshape((len(self.qpt)//self.nq, self.nq), order='F')
        #nq is same as nq at sigma_main

        #READ
        buf = f.read_record()
        self.ng = buf.read('i',1) #number of g-vects
        self.gvec_k = buf.read('i') #g-vects
        self.gvec_k = self.gvec_k.reshape((len(self.gvec_k)//self.ng, self.ng), order='F')

        self.nmtx = empty( self.nq, dtype=int )
        self.isort = []
        self.isort_i = []
        self.ekin = []
        self.epsmat = [] #array of eps matrices
        self.q = empty( (3, self.nq), order='F', dtype=float64 ) #the points

    def write_header(self):
        if not self.f:
            self.f = FortranIO.FortranFile(self.fname,'=','i','wb')
        else:
            self.f.seek(0)
        f = self.f

        #WRITE
        data=self.name.ljust(6) + self.date.ljust(11)
        f.writeline(data)

        #WRITE
        f.write_vals('i', self.freq_dep, self.num_freq)

        #WRITE
        f.write_vals('i', *self.grid.ravel('F'))

        #WRITE
        if self.freq_dep==2:
            f.write_vals('d', *np.hstack((self.freq_grid, self.freq_broad.view(np.float64))))
        else:
            f.writeline()

        #WRITE
        f.writeline() #Empty by GSM
        #WRITE
        f.writeline() #Empty by GSM

        #WRITE
        f.write_vals('d',self.ecuts) #maximum energy, or gmax_in

        #WRITE
        fmt='i'+'d'*self.qpt.size
        data=[self.nq]+list(self.qpt.ravel('F'))
        f.write_vals(fmt, *data)

        #WRITE
        fmt='i'+'i'*self.gvec_k.size
        data=[self.ng]+list(self.gvec_k.ravel('F'))
        f.write_vals(fmt, *data)

    def read_qpt_header(self):
        #still not used!
        if not self.f:
            raise IOError('File not opened')
        if self.nq==0:
            raise IOError('Header not initialized')
        self.cur_q += 1
        f = self.f

        #READ
        buf = f.read_record()
        buf.read('i',1)
        #all ng's are the same, so we can discart this value
        nmtx = buf.read('i',1) #matrix elements, epsi(G,Gp)
        self.nmtx[self.cur_q] = nmtx #size of eps_inv matrix
        #at sigma_main and eps_cop neps = max(nmtx)

        tmp = buf.read('i')
        tmp = tmp.reshape( 2, len(tmp)//2, order='F')
        self.isort.append( tmp[0,:].copy() )
        self.isort_i.append( tmp[1,:].copy() )
        del tmp

        #READ
        self.ekin.append( f.read('d') )
        #note: len(ekin[self.cur_q])=ng
        #ekin is just G*BDOT*Gp (?)
        #ekin < ecutb = bare_coulomb_cutoff

        #READ
        self.q[:,self.cur_q] = f.read('d')

    def _read_qpt_matrix_static(self, in_place=True, ignore=False):
        if not self.f:
            raise IOError('File not opened')
        if self.nq==0:
            raise IOError('Header not initialized')
        f = self.f

        nmtx = self.nmtx[self.cur_q] #size of eps_inv matrix
        if ignore:
            for line in xrange(nmtx):
                f.next()
            return None

        if self.flavor == common.flavor.NONE:
            #autodetect flavor
            #READ
            tmp = f.read('d')
            nmat = len(tmp)
            if tmp.shape[0]==2*nmtx*self.num_freq:
                #complex
                self.flavor = common.flavor.COMPLEX
            elif tmp.shape[0]==nmtx*self.num_freq:
                #real
                self.flavor = common.flavor.REAL
            else:
                raise ValueError('Invalid dimension for qpt %d'%(self.cur_q))
            flavor_str = common.get_numpy_flavor(self.flavor)
            buf = empty((nmtx, nmtx), order='F', dtype=flavor_str)
            tmp = tmp.view(dtype=flavor_str)
            buf[:,0] = tmp
            row_start = 1
        else:
            #we know the flavor
            flavor_str = common.get_numpy_flavor(self.flavor)
            buf = empty((nmtx, nmtx), order='F', dtype=flavor_str)
            row_start = 0

        for line in xrange(row_start, nmtx):
            #READ
            tmp = f.read('d')
            buf[:,line] = tmp.view(dtype=flavor_str)
        if in_place:
            self.epsmat.append(buf)
        else:
            return buf

    def _read_qpt_matrix_ff(self, in_place=True, ignore=False):
        if not self.f:
            raise IOError('File not opened')
        if self.nq==0:
            raise IOError('Header not initialized')
        f = self.f

        nmtx = self.nmtx[self.cur_q] #size of eps_inv matrix
        if ignore:
            for line in xrange(nmtx**2):
                f.next()
            return None

        if self.flavor == common.flavor.NONE:
            raise ValueError('Cannot automatically identify flavor in ff epsmat files.')

        #we know the flavor
        flavor_str = common.get_numpy_flavor(self.flavor)
        buf = []
        buf.append(empty((self.num_freq, nmtx, nmtx), order='F', dtype=np.complex128))
        if self.flavor == common.flavor.COMPLEX:
            buf.append(empty((self.num_freq, nmtx, nmtx), order='F', dtype=np.complex128))

        for line_ in xrange(nmtx):
            for line in xrange(nmtx):
                #READ
                tmp = f.read('d')
                buf[0][:,line, line_] = tmp.view(np.complex128)
            if self.flavor==common.flavor.COMPLEX:
                for line in xrange(nmtx):
                    #READ
                    tmp = f.read('d')
                    buf[1][:,line, line_] = tmp.view(np.complex128)

        if in_place:
            self.epsmat.append(buf)
        else:
            return buf

    def read_qpt_matrix(self, in_place=True, ignore=False):
        if self.freq_dep==0:
            return self._read_qpt_matrix_static(in_place, ignore)
        else:
            return self._read_qpt_matrix_ff(in_place, ignore)

    def read_qpt(self, in_place=True, ignore=False):
        self.read_qpt_header()
        return self.read_qpt_matrix(in_place, ignore)

    def write_qpt(self, cnt):
        if not self.f:
            raise IOError('File not opened')
        if self.nq==0:
            raise IOError('Header not initialized')
        f = self.f

        nmtx = self.nmtx[cnt]
        flavor_str = common.get_numpy_flavor(self.flavor)

        #WRITE
        tmp = column_stack((self.isort[cnt],self.isort_i[cnt]))
        data = [self.ng,nmtx]+list(tmp.flatten())
        f.write_vals('i', *data)
        del tmp
        del data

        #WRITE
        f.write_vals('d', *self.ekin[cnt])

        #WRITE
        f.write_vals('d', *self.q[:,cnt])

        for line in xrange(nmtx):
            #WRITE
            f.write_vals('d', *self.epsmat[cnt][:,line].view(float64))

    def from_file(self, fname, read_all=True):
        self.read_header()

        if read_all:
            for n in xrange(self.nq):
                self.read_qpt()

    def to_file(self, fname, write_all=True):
        f_old = self.f
        fname_old = self.fname
        self.fname = fname
        self.f = None

        self.write_header()
        if write_all:
            for n in xrange(self.nq):
                self.write_qpt(n)

        self.f.close()
        self.fname = fname_old
        self.f = f_old

    def insert_element(self, q, isort, isort_i, ekin, epsmat, deb=False):
        data = row_stack((self.q, q.reshape((1,3))))
        V = core.records.fromarrays( data.T )
        order = argsort(V)
        sz = len(V)
        idx = arange(sz)[order==(sz-1)][0]

        if len(ekin)!=self.ng:
            raise ValueError('Size mismatch: ekin and ng')

        if deb:
            print idx
            print epsmat
        nmtx = epsmat.shape[0]
        self.nq += 1
        self.isort.insert(idx, isort)
        self.isort_i.insert(idx, isort_i)
        self.nmtx = insert(self.nmtx, idx ,nmtx, axis=0)
        self.q = insert(self.q, idx, q, axis=0)
        self.qpt = insert(self.qpt, idx, q, axis=0)
        self.ekin.insert(idx, ekin)
        self.epsmat.insert(idx, epsmat)
        if deb:
            print self.epsmat[idx]

    def remove_element(self, index):
        #assuming q==qpt!
        self.nq -= 1
        self.isort.pop(index)
        self.isort_i.pop(index)
        self.nmtx = delete(self.nmtx, index, 0)
        self.q = delete(self.q, index, 0)
        self.qpt = delete(self.qpt, index, 0)
        self.ekin.pop(index)
        self.epsmat.pop(index)

    def get_diag(self, qpt=0, iw=0):
        if (qpt<0)or(qpt>=self.nq):
            raise ValueError('qpt out of range')

        g_idx = []
        epsinv = []
        for n in xrange(self.nmtx[qpt]):
            srt = self.isort[qpt][n]-1
            #G = self.gvec_k[:, srt]
            #g_idx += [G]
            g_idx += [srt]
            if self.freq_dep==0:
                epsinv += [self.epsmat[qpt][n, n]]
            else:
                epsinv += [self.epsmat[qpt][0][iw, n, n]]
        return array(g_idx), array(epsinv)

    def semi_diag(self, qpt=0):
        if (qpt<0)or(qpt>=self.nq):
            raise ValueError('qpt out of range')

        for dx in [-1,0,1]:
         for dy in [-1,0,1]:
          for dz in [-1,0,1]:
            cond0 = all(self.gvec_k==[dx,dy,dz], axis=1)
            pos0 = where(cond0)[0][0]
            srt_i = self.isort_i[qpt][pos0]-1

            if srt_i < self.nmtx[qpt]:
                srt = self.isort[qpt][srt_i]-1
                G = self.gvec_k[:,srt]
                print G, self.epsmat[qpt][srt_i, srt_i]
            else:
                print array([dx,dy,dz]), 'out of range'

    def __repr__(self):
        bool_repr = { False:'False', True:'True' }
        if self.nq==0: return '<epsmatIO/>'
        return '''<epsmatIO %s>
    File name: %s
    File date: %s
    Flavor: %s
    Freq. dependent: %s
    Num. of freqs.:  %d
    Monk.-Pack grid: %s
    W cut-off:       %f
    Num. of G-vects: %d
    G-vects:
        %s
    Num. of q-pts:   %d
    q-pts:
        %s
    Num. of G-vecs per q-pt:
        %s
    Bare Coulomb energy:
        %s
    q-pts read:
        %s
    Matrices heads:
        %s
</epsmatIO>
'''%\
        (self.fname, self.fname, self.date, common.flavors[self.flavor],
        bool_repr[self.freq_dep], self.num_freq,
        self.grid.__str__(), self.ecuts, self.ng,
        array_str(self.gvec_k,50,6).replace('\n','\n\t\t'),
        self.nq,
        array_str(self.qpt,   50,6).replace('\n','\n\t\t'),
        array_str(self.nmtx,  50,6).replace('\n','\n\t\t'),
        array_str(array(self.ekin),  50,6).replace('\n','\n\t\t'),
        array_str(self.q,     50,6).replace('\n','\n\t\t'),
        array_str(\
            array([self.epsmat[i][0,0] for i in arange(len(self.epsmat))])\
         ,50,6).replace('\n','\n\t\t'),
        )

if __name__=='__main__':
    import sys
    if len(sys.argv)<2:
        print 'Usage: %0 epsmat|eps0mat [...]'%(sys.argv[0])
        sys.exit()

    for fname in sys.argv[1:]:
        epsmat = epsmatIO(sys.argv[1], read_all=False)
        print epsmat
