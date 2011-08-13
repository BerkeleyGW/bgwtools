#!/usr/bin/env python

from __future__ import division

SAVE=0

import FortranIO
from numpy import *
import numpy.core as core

class epsmatIO:
	def __init__(self, fname=None, auto_read=True):
		self.fname=fname
		self.f=None
		self.cur_q=0 #next q_pt to be read

		self.name=''
		self.date=''
		self.freq_dep=0
		self.num_freq=1
		self.grid=array([1,1,1])
		self.ecuts=0.
		self.nq=0
		self.qpt=empty((0,3), dtype=float64)
		self.ng=0
		self.gvec_k=empty((0,3), dtype=int)
		self.nmtx=empty(0, dtype=int)
		self.isort = []
		self.isort_i = []
		self.ekin = []
		self.epsmat = []
		self.q = empty((0,3), dtype=float64)

		if fname and auto_read:
			self.from_file(self.fname)

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
		f.read() #only important if freq dependent
		#READ
		f.read() #??
		#READ
		f.read() #??

		#READ
		self.ecuts = f.read('d')[0] #maximum energy, or gmax_in

		#READ
		buf = f.read_record()
		self.nq = buf.read('i',1) #number of q-pts in this file
		self.qpt = buf.read('d') #q-pts
		self.qpt = self.qpt.reshape(self.nq, len(self.qpt)//self.nq)
		#nq is same as nq at sigma_main

		#READ
		buf = f.read_record()
		self.ng = buf.read('i',1) #number of g-vects
		self.gvec_k = buf.read('i') #g-vects
		self.gvec_k = self.gvec_k.reshape(self.ng, len(self.gvec_k)//self.ng)

		self.nmtx = empty( self.nq, dtype=int )
		self.isort = []
		self.isort_i = []
		self.ekin = []
		self.epsmat = [] #array of eps matrices
		self.q = empty( (self.nq,3), dtype=float64 ) #the points

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
		f.write_vals('i', *self.grid.flatten())

		#WRITE
		f.writeline() #only important if freq dependent
		#WRITE
		f.writeline() #??
		#WRITE
		f.writeline() #??

		#WRITE
		f.write_vals('d',self.ecuts) #maximum energy, or gmax_in

		#WRITE
		fmt='i'+'d'*self.qpt.size
		data=[self.nq]+list(self.qpt.flatten())
		f.write_vals(fmt, *data)

		#WRITE
		fmt='i'+'i'*self.gvec_k.size
		data=[self.ng]+list(self.gvec_k.flatten())
		f.write_vals(fmt, *data)

	def read_qpt(self):
		if not self.f:
			raise IOError('File not opened')
		if self.nq==0:
			raise IOError('Header not initialized')
		self.cur_q += 1
		f = self.f
		
		cnt = len(self.epsmat)
		
		#READ
		buf = f.read_record()
		buf.read('i',1)
		#all ng's are the same, so we can discart this value
		nmtx = buf.read('i',1) #matrix elements, epsi(G,Gp)
		self.nmtx[cnt] = nmtx #size of eps_inv matrix
		#at sigma_main and eps_cop neps = max(nmtx)

		tmp = buf.read('i')
		tmp = tmp.reshape(len(tmp)//2, 2)
		self.isort.append( tmp[:,0].copy() )
		self.isort_i.append( tmp[:,1].copy() )
		del tmp

		#READ
		self.ekin.append( f.read('d') )
		#note: len(ekin[cnt])=ng
		#ekin is just G*BDOT*Gp (?)
		#ekin < ecutb = bare_coulomb_cutoff

		#READ
		self.q[cnt,:] = f.read('d')

		self.epsmat.append( empty((nmtx, nmtx)) )
		for line in xrange(nmtx):
			#READ
			tmp = f.read('d')
			self.epsmat[-1][:,line] = tmp

	def write_qpt(self, cnt):
		if not self.f:
			raise IOError('File not opened')
		if self.nq==0:
			raise IOError('Header not initialized')
		f = self.f

		nmtx = self.nmtx[cnt]		

		#WRITE
		tmp = column_stack((self.isort[cnt],self.isort_i[cnt]))
		data = [self.ng,nmtx]+list(tmp.flatten())
		f.write_vals('i', *data)
		del tmp
		del data

		#WRITE
		f.write_vals('d', *self.ekin[cnt])

		#WRITE
		f.write_vals('d', *self.q[cnt])

		for line in xrange(nmtx):
			#WRITE
			f.write_vals('d', *self.epsmat[cnt][:,line])


	def from_file(self, fname):
		self.read_header()

		for n in xrange(self.nq):
			self.read_qpt()

	def to_file(self, fname):
		f_old = self.f
		fname_old = self.fname
		self.fname = fname
		self.f = None

		self.write_header()
		for n in xrange(self.nq):
			self.write_qpt(n)

		self.f.close()
		self.fname = fname_old
		self.f = f_old

	def insert_element(self, q, isort, isort_i, ekin, epsmat):
		data = row_stack((self.q, q.reshape((1,3))))
		V = core.records.fromarrays( data.T )
		order = argsort(V)
		sz = len(V)
		idx = arange(sz)[order==(sz-1)][0]

		if len(ekin)!=self.ng:
			raise ValueError('Size mismatch: ekin and ng')

		nmtx = epsmat.shape[0]
		self.nq += 1
		self.isort.insert(idx, isort)
		self.isort_i.insert(idx, isort_i)
		self.nmtx = insert(self.nmtx, idx ,nmtx, axis=0)
		self.q = insert(self.q, idx, q, axis=0)
		self.qpt = insert(self.qpt, idx, q, axis=0)
		self.ekin.insert(idx, ekin)
		self.epsmat.insert(idx,epsmat)

	def get_diag(self, qpt=0):
		if (qpt<0)or(qpt>=self.nq):
			raise ValueError('qpt out of range')

		for n in xrange(self.nmtx[qpt]):
			srt = self.isort[qpt][n]-1
			G = self.gvec_k[srt, :]
			print G, self.epsmat[qpt][n, n]
		print	

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
				G = self.gvec_k[srt, :]
				print G, self.epsmat[qpt][srt_i, srt_i]
			else:
				print array([dx,dy,dz]), 'out of range' 

	def __repr__(self):
		bool_repr = { False:'False', True:'True' }
		if self.nq==0: return 'Empty eps_class'
		return '''eps_class
	File name: %s
	File date: %s
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
	Rank of eps matrix per q-pt:
		%s
	Bare Coulomb energy:
		%s
	q-pts read:
		%s'''%\
		(self.fname, self.date, bool_repr[self.freq_dep], self.num_freq,
		self.grid.__str__(), self.ecuts, self.ng,
		array_str(self.gvec_k,50,6).replace('\n','\n\t\t'),
		self.nq,
		array_str(self.qpt,   50,6).replace('\n','\n\t\t'),
		array_str(self.nmtx,  50,6).replace('\n','\n\t\t'),
		array_str(array(self.ekin),  50,6).replace('\n','\n\t\t'),
		array_str(self.q,     50,6).replace('\n','\n\t\t'),
		) 
	
if __name__=='__main__':
	import sys
	if len(sys.argv)<2:
		print 'Expecting one argument: [epsmat/eps0mat file]'
		sys.exit()

	epsmat = epsmatIO(sys.argv[1])
	print epsmat
	print
	if SAVE:
		epsmat.to_file('/tmp/epsout')
	#print epsmat.epsmat
	
