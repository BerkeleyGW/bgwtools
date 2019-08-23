#!/usr/bin/env python

import FortranIO
import numpy as np
#from numpy import *

class eigenvectorsIO:
	def __init__(self, fname=None, num_evecs=None, full=False):
		self.fname=fname
		self.f=None
		self.nk=0
		self.cur_evec=0
		self.num_evecs = 0
                self.full=full

		if fname:
			self.from_file(self.fname, num_evecs)

	def read_header(self, num_evecs=None):
		if not self.f:
			self.f = FortranIO.FortranFile(self.fname)
		else:		
			self.f.seek(0)
		f = self.f
		
		#READ
		self.ns = f.read('i')[0]
		#READ
		self.nv = f.read('i')[0]
		#READ
		self.nc = f.read('i')[0]
		#READ
		self.nk = f.read('i')[0]

		self.num_evecs = self.ns * self.nv * self.nc * self.nk
		if not (isinstance(num_evecs,int)):
			num_evecs = self.num_evecs

		#READ
		self.kpt = f.read('d').reshape((3,self.nk), order='F')

		self.evals = np.empty(num_evecs, dtype=np.float64)
                fact = 1
                if self.full:
                    fact = 2
		self.evecs = np.empty((num_evecs,fact*self.num_evecs), dtype=np.complex128)

	def read_evecs(self, num_evecs=None):
		if not self.f:
			raise IOError('File not opened')
		if self.num_evecs==0:
			raise IOError('Header not initialized')
		f = self.f
		
		if not (isinstance(num_evecs,int)):
			num_evecs = self.num_evecs
		for i in range(num_evecs):
			#READ
			self.evals[i] = f.read('d')
			#READ
			self.evecs[i,:] = f.read('d').view(np.complex128)

	def from_file(self, fname, num_evecs):
		self.read_header(num_evecs)
		self.read_evecs(num_evecs)

	def __repr__(self):
		bool_repr = { False:'False', True:'True' }
		if self.nk==0: return 'Empty eps_class'
		str='''wigenvectorsIO
	File name: %s
	ns: %d
	nv: %d
	nc: %d
	nk: %d'''%\
		(self.fname, self.ns, self.nv, self.nc, self.nk,)
		str+='''
	K-grid:
		%s
	evals:
		%s
	evecs:
		%s
	'''%\
			(\
			array_str(self.kpt.T,50,6).replace('\n','\n\t\t'),
			array_str(self.evals,50,6).replace('\n','\n\t\t'),
			array_str(self.evecs,50,6).replace('\n','\n\t\t'),
			)
		return str

		
	
if __name__=='__main__':
	import sys
	if len(sys.argv)<2:
		print 'Expecting one argument: [eigenvectors file]'
		sys.exit()

	evecs = eigenvectorsIO(sys.argv[1])
	print evecs

	import matplotlib.pyplot as plt
	#plt.plot(evecs.evals[:100])
	#plt.spy(evecs.evecs[:100,:])
	print evecs.evals[5000-1]
	plt.plot(evecs.evecs[5000-1,:])
	plt.plot(evecs.evecs[5000-2,:])
	plt.show()


