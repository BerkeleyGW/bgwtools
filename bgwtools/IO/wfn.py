#!/usr/bin/env python

# An abstract class that reads/writes WFN and RHO files.
# Handles REAL and COMPLEX data, but not spinors.

# Felipe Homrich da Jornada <jornada@civet.berkeley.edu> (2012)

from __future__ import division
from bgwtools.common import common
import FortranIO
from numpy import *

class wfnIO:
	ftypes=['UNK','WFN','RHO','VXC']
	def __init__(self, fname=None, full=True):
		self.fname = fname
		self.name = ''
		self.nat = 0
		self.f = None
		self.ftype = 0 #flavor
                self.nbands = 0
                self.ns = 0
                self.nk = 0
                self.ng = 0 
		self.gvec = empty((3,0), order='F')
		self.kpt = empty((3,0), order='F')
		self.ifmin = empty((0,0), order='F')
		self.ifmax = empty((0,0), order='F')
		self.energies = empty((0,0,0), order='F')
		self.occupations = empty((0,0,0), order='F')
		self.flavor = common.flavor.NONE
                self.ns = 0
                self.nspinor = 0

		if fname:
			self.from_file(self.fname, full)

	def read_header(self, full=True):
		if not self.f:
			self.f = FortranIO.FortranFile(self.fname)
		else:		
			self.f.seek(0)
		f = self.f
		
		#READ
		rec = f.readline()
		self.name=rec[:32]
		self.date=rec[32:32+32]
		self.time=rec[32+32:32+32+32]
		self.ftype=common.get_ftype(self.name[:3], die=True)
		self.flavor=common.get_flavor(self.name[4:11], die=True)

		#READ
		buf = f.read_record()

		self.ns = buf.read('i',1)
                self.nspinor = 1
                if (self.ns==4):
                    self.ns = 1
                    self.nspinor = 2
		self.ng = buf.read('i',1)
		self.ntran = buf.read('i',1)
		self.cell_symmetry = buf.read('i',1)
		self.nat = buf.read('i',1)
		self.ecutrho = buf.read('d',1)
		#only for wfn
		if self.ftype==1:
			self.nk = buf.read('i',1)
			self.nbands = buf.read('i',1)
			self.ngkmax = buf.read('i',1)
			self.ecutwfn = buf.read('d',1)

		#READ
		buf = f.read_record()
		self.FFTgrid = buf.read('i',3)
		#only for wfn
		if self.ftype==1:
			self.kgrid = buf.read('i',3)
			self.kshift = buf.read('d',3)

		start = file.tell(f)
		#print start

		#READ
		buf = f.read_record()
		self.celvol = buf.read('d',1)
		self.alat = buf.read('d',1)
		self.avec = buf.read('d',9).reshape((3,3), order='F')
		self.adot = buf.read('d',9).reshape((3,3), order='F')
		
		#READ
		buf = f.read_record()
		self.recvol = buf.read('d',1)
		self.blat = buf.read('d',1)
		self.bvec = buf.read('d',9).reshape((3,3), order='F')
		self.bdot = buf.read('d',9).reshape((3,3), order='F')
		
		#READ
		self.mtrx = f.read('i').reshape((3,3,self.ntran), order='F') #symmetry el
		
		#READ
		self.tnp = f.read('d').reshape((3,self.ntran), order='F') #frac. translation
		end = file.tell(f)
		#print end-start

		#READ
		buf = f.read_record()
		self.apos = empty((3,self.nat), order='F')
		self.atyp = empty(self.nat)
		for n in range(self.nat):
			self.apos[:,n] = buf.read('d',3)
			self.atyp[n] = buf.read('i',1)

		#only for wfn
		if self.ftype==1:
			#READ
			self.ngk = f.read('i')
			#READ
			self.kw = f.read('d')
			#READ
			self.kpt = f.read('d').reshape((3,self.nk), order='F')
			#READ
			self.ifmin = f.read('i').reshape((self.nk, self.ns), order='F')
			#READ
			self.ifmax = f.read('i').reshape((self.nk, self.ns), order='F')
			#same as fortran order
			#READ
			self.energies = f.read('d').reshape((self.nbands, self.nk, self.ns), order='F')
			#same as fortran order
			#READ
			self.occupations = f.read('d').reshape((self.nbands, self.nk, self.ns), order='F')

			if (not full): return
			self.gvec = empty((3,self.ng), order='F', dtype='i')
			#READ
			self.read_gvectors(self.gvec)

	def write_header(self, full=False):
		if not self.f:
			self.f = FortranIO.FortranFile(self.fname,'=','i','wb')
		else:		
			self.f.seek(0)
		f = self.f

		#WRITE
		data=self.name.ljust(32) + self.date.ljust(32) + self.time.ljust(32)
		f.writeline(data)
		
		#WRITE
		fmt = 'iiiiid'
		data = [self.ns, self.ng, self.ntran, self.cell_symmetry, self.nat, self.ecutrho]
		#only for wfn
		if self.ftype==1:
			fmt += 'iiid'
			data += [self.nk, self.nbands, self.ngkmax, self.ecutwfn]
		f.write_vals(fmt, *data)

		#WRITE
		fmt = 'iii'
		data = list(self.FFTgrid.ravel('F'))
		#only for wfn
		if self.ftype==1:
			fmt += 'iiiddd'
			data += list(self.kgrid.ravel('F')) + list(self.kshift.ravel('F'))
		f.write_vals(fmt, *data)

		#WRITE
		fmt = 'd'*20
		data = [self.celvol, self.alat] + list(self.avec.ravel('F')) + list(self.adot.ravel('F'))
		f.write_vals(fmt, *data)
		
		#WRITE
		fmt = 'd'*20
		data = [self.recvol, self.blat] + list(self.bvec.ravel('F')) + list(self.bdot.ravel('F'))
		f.write_vals(fmt, *data)
		
		#WRITE
		tmp = self.mtrx.ravel('F')
		f.write_vals('i'*len(tmp), *tmp )
		
		#WRITE
		tmp = self.tnp.ravel('F')
		f.write_vals('d'*len(tmp), *tmp )

		#WRITE
		fmt='dddi'*self.nat
		data = []
		for i in range(self.nat):
			data += list(self.apos[:,i].ravel('F')) + [self.atyp[i]]
		f.write_vals(fmt, *data)
		
		#only for wfn
		if self.ftype==1:
			nk = self.nk
			ns = self.ns
			nb = self.nbands
			#WRITE
			f.write_vals('i'*nk, *self.ngk)
			#WRITE
			f.write_vals('d'*nk, *self.kw.ravel('F'))
			#WRITE
			f.write_vals('d'*3*nk, *self.kpt.ravel('F'))
			#WRITE
			f.write_vals('i'*nk*ns, *self.ifmin.ravel('F'))
			#WRITE
			f.write_vals('i'*nk*ns, *self.ifmax.ravel('F'))
			#WRITE
			f.write_vals('d'*nk*ns*nb, *self.energies.ravel('F'))
			#WRITE
			f.write_vals('d'*nk*ns*nb, *self.occupations.ravel('F'))
			if (not full): return
			#WRITE
			self.write_gvectors(self.gvec)


	def read_gvectors(self, gvec=None):
		f = self.f
		#READ
		nrecord_internal = f.read('i')[0]
		ig = 0
		for i in xrange(nrecord_internal):
			#READ
			ng_irecord = f.read('i')[0]
			if gvec is None:
                                f.next()
                        else:
                                #READ
                                buf = f.read('i')
				gvec[:,ig:ig+ng_irecord] = buf.reshape((3,ng_irecord), order='F')
				del buf

                        
			ig += ng_irecord

	def write_gvectors(self, gvec=None):
		f = self.f
		#WRITE
		f.write_vals('i', 1)
		if gvec is None:
			gvec = self.gvec
		#WRITE
		f.write_vals('i', gvec.shape[1])
		#WRITE
		#f.write_vals2('i', gvec.ravel('F'))
		f.write_vals('i'*len(gvec.ravel()), *gvec.ravel('F'))

        def get_gvectors_buffer(self):
            '''Results an ndarray that serves as buffer for calls to read_gvectors'''

            return empty((3, self.ngkmax), dtype=int, order='F')

        def get_data_buffer(self):
            '''Results an ndarray that serves as buffer for calls to read_data'''

            if self.flavor==common.flavor.NONE:
                raise Exception('Flavor is unknown.')

            dtype = common.get_numpy_flavor(self.flavor)
            return empty((self.ngkmax, self.ns), dtype=dtype, order='F')

	def read_data(self, data=None):
		f = self.f
		#READ
		nrecord_internal = f.read('i')[0]
		ig = 0
		ns = self.ns
		for i in xrange(nrecord_internal):
			#READ
			ng_irecord = f.read('i')[0] * self.nspinor

			if data is None:
                            f.next()
                        else:
                            #READ
                            buf = f.read('d')
                            if self.flavor==common.flavor.COMPLEX:
				buf = buf.view(complex128)
			    data[ig:ig+ng_irecord] = buf.reshape((ng_irecord,ns), order='F')
                            del buf
			ig += ng_irecord
	
	def write_data(self, data):
		f = self.f
		#WRITE
		f.write_vals('i', 1)
		sz = data.size
		if self.flavor==common.flavor.COMPLEX and data.dtype=='d':
			sz //= 2
		#WRITE
		f.write_vals('i', sz)
		#WRITE
		data_view = data.ravel().view(float64)
		f.write_vals('d'*len(data_view), *data_view)
		#f.write_vals2('d', data_view)

	def from_file(self, fname, full=True):
		self.read_header(full)

	def to_file(self, fname, full=True):
		f_old = self.f
		fname_old = self.fname
		self.fname = fname
		self.f = None
                pos = 0
		try:
			self.write_header(full)
                        pos = file.tell(self.f)
			self.f.close()
		finally:
			self.fname = fname_old
			self.f = f_old
                return pos

	def __repr__(self):
		bool_repr = { False:'False', True:'True' }
		#if self.nat==0: return '<wfnIO/>'
		str='''<wfnIO %s>
	File name: %s
	Work name: %s
	Work date: %s
	Work time: %s
	Flavor: %s
	Density cut-off: %f
	Reciprocal vectors:
		%s
	Reciprocal tensor:
		%s
        Number of atoms: %d
        Elements:
                %s
        Atomic positions:
                %s
        Fractional translations:
                %s
	Symmetry els:
		%s'''%\
		(self.fname, self.fname, self.name, self.date, self.time,
		common.flavors[self.flavor], self.ecutrho,
		array_str(self.bvec,50,6).replace('\n','\n\t\t'),
		array_str(self.bdot,50,6).replace('\n','\n\t\t'),
                self.nat, self.atyp,
                array_str(self.apos.T,50,6).replace('\n','\n\t\t'),
                array_str(self.tnp.T,50,6).replace('\n','\n\t\t'),
		#array_str(transpose(self.mtrx,[2,1,0]),50,6).replace('\n','\n\t\t'),
		array_str(self.mtrx,50,6).replace('\n','\n\t\t'),
		)
		if self.ftype==1:
			str+='''
	Energy cutoff: %f
	Number of bands: %d
	Number of k-pts: %d
	Number of G-gvectors per k-point:
		%s
	K-grid:
		%s
	K-shifts:
		%s
	K-pts:
		%s
	Lowest  occ band (ifmin):
		%s
	Highest occ band (ifmax):
		%s
	Energies^T:
		%s
	Occupations^T:
		%s
	K-pts weights:
		%s
</wfnIO>
'''%\
			(self.ecutwfn, self.nbands, self.nk,
			array_str(self.ngk,50,6).replace('\n','\n\t\t'),
			array_str(self.kgrid,50,6).replace('\n','\n\t\t'),
			array_str(self.kshift,50,6).replace('\n','\n\t\t'),
			array_str(self.kpt.T,50,6).replace('\n','\n\t\t'),
			array_str(self.ifmin.T,50,6).replace('\n','\n\t\t'),
			array_str(self.ifmax.T,50,6).replace('\n','\n\t\t'),
			array_str(self.energies.T,50,6).replace('\n','\n\t\t'),
			array_str(self.occupations.T,50,6).replace('\n','\n\t\t'),
			array_str(self.kw,50,6).replace('\n','\n\t\t'),
			)
		return str

if __name__=='__main__':
	import sys
	if len(sys.argv)<2:
		print 'Usage: %0 WFN|RHO [...]'%(sys.argv[0])
		sys.exit()

        set_printoptions(suppress=True)
	for fname in sys.argv[1:]:
		wfn = wfnIO(fname, full=False)
		print wfn
