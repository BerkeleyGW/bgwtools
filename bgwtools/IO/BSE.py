#!/usr/bin/env python


from numpy import *
from FortranIO import FortranFile

class eigenvectorsIO:
	def __init__(self, fname=None, auto_read=False):
		self.fname = fname
		self.f = None
		self.nk = 0
		self.cur_evec = 0
		self.num_evecs = 0

		if fname:
			self.f = FortranFile(self.fname)
			if auto_read:
				self.read_header()

	def read_header(self):
		f = self.f
		self.nspin = f.read('i')[0]
		self.nvband = f.read('i')[0]
		self.ncband = f.read('i')[0]
		self.nkpt = f.read('i')[0]
		self.kpts = f.read('d').reshape((3,self.nkpt), order='F')
		self.neig = self.nspin*self.nvband*self.ncband*self.nkpt

		#self.evals = empty(num_evecs ,dtype='d')
		#self.evecs = empty((num_evecs,self.num_evecs), dtype='d')

	def next_evec(self):
		self.f.next(); self.f.next()
		self.cur_evec += 1

	def prev_evec(self):
		self.f.prev(); self.f.prev()
		self.cur_evec -= 1
		
	def get_evec(self):
		eval = self.f.read('d')[0]
		evec = self.f.read('d')
		self.cur_evec += 1
		return (eval,evec)

	def goto_evec(self, idx, whence=0):
		record_size = 4*self.f._header_length + 8 + 8*self.neig
		if whence==0:
			file.seek(self.f, 0)
			self.read_header()
			file.seek(self.f, record_size*idx, 1)
			self.cur_evec = idx
		elif whence==1:
			file.seek(self.f, record_size*(idx-self.cur_evec), 1)
			self.cur_evec += idx
		else:
			raise NotImplementedError('whence=2 is not implemented')

class vmtxelIO:
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
		self.opr = buf[4]

	def read_mtxel(self):
		return self.f.read('d')

	def __repr__(self):
		opr = ('velocity','momentum')[self.opr]
		return '''<vmtxelIO %s>
	File name: %s
	Num. k-points:   %d
	Num. cond bands: %d
	Num. val  bands: %d
	Num. spins:      %d
	Operator:        %s
</vmtxelIO>
'''%\
		(self.fname, self.fname, self.nkpt, self.nc, self.nv, self.nspin, opr)
