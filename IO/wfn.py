#!/usr/bin/env python

import FortranIO
from numpy import *

class wfnIO:
	ftypes=['UNK','WFN','RHO']
	def __init__(self, fname=None):
		self.fname=fname
		self.name=''
		self.nat=0
		self.f=None
		self.ftype=0

		if fname:
			self.from_file(self.fname)

	def read_header(self):
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
		try:
			self.ftype=self.ftypes.index(self.name[:3])
		except:
			self.ftype=0

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
		self.kmax = array([buf.read('i',1) for n in range(3)])
		#only for wfn
		if self.ftype==1:
			self.kgrid = array([buf.read('i',1) for n in range(3)])
			self.kshift = array([buf.read('d',1) for n in range(3)])

		#READ
		buf = f.read_record()
		self.celvol = buf.read('d',1)
		self.alat = buf.read('d',1)
		self.avec = array([buf.read('d',1) for n in range(9)]).reshape((3,3))
		self.adot = array([buf.read('d',1) for n in range(9)]).reshape((3,3))
		
		#READ
		buf = f.read_record()
		self.recvol = buf.read('d',1)
		self.blat = buf.read('d',1)
		self.bvec = array([buf.read('d',1) for n in range(9)]).reshape((3,3))
		self.bdot = array([buf.read('d',1) for n in range(9)]).reshape((3,3))
		
		#READ
		self.mtrx = f.read('i').reshape((self.ntran,3,3)).T #symmetry el
		
		#READ
		self.tnp = f.read('d').reshape((self.ntran,3)).T #frac. translation

		#READ
		buf = f.read_record()
		self.apos = empty((3,self.nat))
		self.atyp = empty(self.nat)
		for n in range(self.nat):
			self.apos[:,n] = buf.read('d',3)
			self.atyp[n] = buf.read('i',1)
		#only for wfn
		if self.ftype==1:
			self.ngk = f.read('i')
			self.kw = f.read('d')
			self.kpt = f.read('d').reshape((self.nk,3)).T

	def from_file(self, fname):
		self.read_header()

	def __repr__(self):
		bool_repr = { False:'False', True:'True' }
		if self.nat==0: return 'Empty eps_class'
		str='''wfnIO
	File name: %s
	Work name: %s
	Work date: %s
	Work time: %s
	Density cut-off: %f
	Reciprocal vectors:
		%s
	Reciprocal tensor:
		%s
	Symmetry els:
		%s'''%\
		(self.fname, self.name, self.date, self.time, self.ecutrho,
		array_str(self.bvec,50,6).replace('\n','\n\t\t'),
		array_str(self.bdot,50,6).replace('\n','\n\t\t'),
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
		%s'''%\
			(self.nk,
			array_str(self.kgrid,50,6).replace('\n','\n\t\t'),
			array_str(self.kshift,50,6).replace('\n','\n\t\t'),
			array_str(self.kpt.T,50,6).replace('\n','\n\t\t'),
			)
		return str

		
	
if __name__=='__main__':
	import sys
	if len(sys.argv)<2:
		print 'Expecting one argument: [RHO file]'
		sys.exit()

	wfn = wfnIO(sys.argv[1])
	print wfn
	print

