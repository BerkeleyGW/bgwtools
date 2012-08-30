#!/usr/bin/env python

# An abstract class that reads/writes WFN and RHO files.
# Handles REAL and COMPLEX data, but not spinors.

# Felipe Homrich da Jornada <jornada@civet.berkeley.edu> (2012)

from __future__ import division
from bgwtools.common import common
import FortranIO
from numpy import *

class wfnIO:
	ftypes=['UNK','WFN','RHO']
	def __init__(self, fname=None, full=True):
		self.fname=fname
		self.name=''
		self.nat=0
		self.f=None
		self.ftype=0 #flavor
		self.kpt=empty((3,1))
		self.ifmin=empty((3,3))
		self.ifmax=empty((3,3))
		self.energies=empty((3,3))
		self.occupations=empty((3,3))
		self.flavor=common.flavor.NONE

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
		self.flavor=common.get_flavor(self.name[4:8], die=True)

		#READ
		buf = f.read_record()

		self.ns = buf.read('i',1)
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
		self.kmax = buf.read('i',3)
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
			#READ
			self.gvec = empty((3,self.ng), order='F', dtype='i')

			if (not full): return
			self.read_gvectors(self.gvec)

	def write_header(self, full=True):
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
		data = list(self.kmax.ravel('F'))
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
		
		if (not full): return

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
			#WRITE
			f.write_vals('i'*3*self.ng, *self.gvec.ravel('F'))

	def read_gvectors(self, gvec=None):
		f = self.f
		nrecord_internal = f.read('i')
		ig = 0
		for i in xrange(nrecord_internal):
			ng_irecord = f.read('i')
			buf = f.read('i')
			if not gvec is None:
				gvec[:,ig:ig+ng_irecord] = buf.reshape((3,ng_irecord), order='F')
			del buf
			ig += ng_irecord

	def read_data(self, data=None):
		f = self.f
		nrecord_internal = f.read('i')
		ig = 0
		ns = 1
		for i in xrange(nrecord_internal):
			ng_irecord = f.read('i')
			buf = f.read('d')
			# do we have complex data?
			if len(buf) == 2*ng_irecord:
				buf = buf.view(complex128)
			if not data is None:
				data[ig:ig+ng_irecord] = buf.reshape((ng_irecord,ns), order='F')
			del buf
			ig += ng_irecord

	def from_file(self, fname, full=True):
		self.read_header(full)

	def to_file(self, fname, full=True):
		f_old = self.f
		fname_old = self.fname
		self.fname = fname
		self.f = None
		try:
			self.write_header(full)
			self.f.close()
		finally:
			self.fname = fname_old
			self.f = f_old

	def __repr__(self):
		bool_repr = { False:'False', True:'True' }
		if self.nat==0: return '<wfnIO/>'
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
	Symmetry els:
		%s'''%\
		(self.fname, self.fname, self.name, self.date, self.time,
		common.flavors[self.flavor], self.ecutrho,
		array_str(self.bvec,50,6).replace('\n','\n\t\t'),
		array_str(self.bdot,50,6).replace('\n','\n\t\t'),
		#array_str(transpose(self.mtrx,[2,1,0]),50,6).replace('\n','\n\t\t'),
		array_str(self.mtrx,50,6).replace('\n','\n\t\t'),
		)
		if self.ftype==1:
			str+='''
	Number of k-pts: %d
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
			(self.nk,
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

	for fname in sys.argv[1:]:
		wfn = wfnIO(fname, full=False)
		print wfn
