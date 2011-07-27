#!/usr/bin/env python

import FortranIO
from numpy import *

class epsmatIO:
	def __init__(self, fname=None):
		self.fname=fname
		self.name=''
		self.nq=0
		self.f=None

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
		self.nq = buf.read('i',1)[0] #number of q-pts in this file
		self.qpt = buf.read('d') #q-pts
		self.qpt = self.qpt.reshape(self.nq, len(self.qpt)//self.nq)

		#READ
		buf = f.read_record()
		self.ng = buf.read('i',1)[0] #number of g-vects
		self.gvec_k = buf.read('i') #g-vects
		self.gvec_k = self.gvec_k.reshape(self.ng, len(self.gvec_k)//self.ng)

	def read_qpt(self):
		if not self.f:
			raise IOError('File not opened')
		if self.nq==0:
			raise IOError('Header not initialized')
		f = self.f
		
		cnt = len(self.epsmat)
		
		#READ
		buf = f.read_record()
		self.ng_q[cnt] = buf.read('i',1)[0] #same as ng?
		nmtx = buf.read('i',1)[0] #matrix elements, epsi(G,Gp)
		self.nmtx[cnt] = nmtx #size of eps_inv matrix

		#READ
		tmp = buf.read('i')
		tmp = tmp.reshape(len(tmp)//2, 2)
		self.isort.append( tmp[:,0].copy() )
		self.isort_i.append( tmp[:,1].copy() )
		del tmp

		#READ
		self.ekin.append( f.read('d') )

		#READ
		self.q[cnt,:] = f.read('d')

		self.epsmat.append( empty((nmtx, nmtx)) )
		for line in xrange(nmtx):
			#READ
			tmp = f.read('d')
			self.epsmat[-1][line,:] = tmp

	def from_file(self, fname):
		self.read_header()

		self.ng_q = empty( self.nq, dtype=int )
		self.nmtx = empty( self.nq, dtype=int )
		self.isort = []
		self.isort_i = []
		self.ekin = []
		self.epsmat = [] #array os eps matrices
		self.q = empty( (self.nq,3), dtype=float64 ) #the points
		for n in xrange(self.nq):
			self.read_qpt()

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
	Cut-off value:   %f
	Num. of G-vects: %d
	G-vects:
		%s
	Num. of q-pts:   %d
	q-pts:
		%s
	Matrix rank per q-pt:
		%s
	q-pts read:
		%s'''%\
		(self.fname, self.date, bool_repr[self.freq_dep], self.num_freq,
		self.grid.__str__(), self.ecuts, self.ng,
		array_str(self.gvec_k,50,6).replace('\n','\n\t\t'),
		self.nq,
		array_str(self.qpt,   50,6).replace('\n','\n\t\t'),
		array_str(self.nmtx,  50,6).replace('\n','\n\t\t'),
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
	#epsmat.get_diag()
	#for pt in range(epsmat.nq):
	#	print '%.4f %.4f %.4f  %.8f'%\
	#		tuple(epsmat.q[pt].tolist()+[epsmat.epsmat[pt][0,0]])
	#epsmat.f.seek()
	#data = epsmat.f.read(12)
	#epsmat.get_diag()
	#for n in range(epsmat.nq):
	#	q = epsmat.q[n]
	#	print
	#	print q
	#	print
	#	epsmat.semi_diag(n)

