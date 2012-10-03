#!/usr/bin/env python

#Truncate/untruncate an epsmat file

from bgwtools.IO.epsmat import epsmatIO
from bgwtools.IO.wfn import wfnIO
import sys
import numpy as np
from numpy import *


class truncator:
	def __init__(self, bdot):
		self.bdot = bdot
		self.Mdot = np.linalg.cholesky(bdot).T

	def get_factor(self, gvec_q):
		return ones_like(gvec_q)

class no_truncator(truncator):
	pass

class slab_truncator(truncator):
	def get_factor(self, gvec_q):
		qk = gvec_q
		zc = pi/sqrt(self.bdot[2,2])
		qkxy = qk.copy()
		qkxy[2,:] = 0
		qkz = qk.copy()
		qkz[0:2,:] = 0
		kxy=sqrt( sum( dot(self.Mdot,qkxy)**2, axis=0) )
		kz=sqrt( sum( dot(self.Mdot,qkz)**2, axis=0) )

		return 1.0 - exp(-kxy*zc)*cos(kz*zc)


class epsilon_truncator:
	def __init__(self, epsmat, wfn, truncation):
		self.epsmat = epsmat
		self.wfn = wfn
		#self.ekin2 = sum(dot(self.Mdot,qk)**2, axis=0)
		if issubclass(truncation, truncator):
			self.truncation = truncation(wfn.bdot)
		elif isinstance(truncation, truncator):
			self.truncation = truncation
		else:
			raise TypeError('Unknown data type for truncation scheme.')
		
	def truncate(self):
		for iq in xrange(len(self.epsmat.epsmat)):
			nmtx = self.epsmat.nmtx[iq]  #number of G vectors = rows=cols of epsinv
			qpt = self.epsmat.qpt[:,iq] 
			mat = self.epsmat.epsmat[iq] #epsinv
			#gidx, epsinv = epsmat.get_diag(i)
			gidx = self.epsmat.isort[iq][:nmtx] - 1
			#This is the array of |G+k|^2
			#ekin = array(self.epsmat.ekin)[0]
			#ekin = ekin[gidx]
			#equivalent to ekin!
			#ekin2 = sum(dot(M,gvec_q)**2, axis=0)

			gvec_q = self.epsmat.gvec_k[:,gidx] + qpt[:,np.newaxis]
			g = self.truncation.get_factor(gvec_q)
			g[abs(g)<1e-10] = 1e-10

			#uninvert, do operations, invert back
			#print mat.diagonal()
			mat = linalg.inv(mat)
			for ig in xrange(nmtx): mat[ig,ig] -= 1.0
			mat[:,:] = mat[:,:] / g[:,np.newaxis]
			for ig in xrange(nmtx): mat[ig,ig] += 1.0
			mat = linalg.inv(mat)
			#print mat.diagonal()
			self.epsmat.epsmat[iq] = mat

if __name__=='__main__':
	import sys

	if len(sys.argv)!=4:
		print 'Usage: %s wfn epsmat_in epsmat_out'%(sys.argv[0])
		exit(1)

	wfn = wfnIO(sys.argv[1])
	f_in = sys.argv[2]
	f_out = sys.argv[3]

	eps = epsmatIO(f_in)
	conv = epsilon_truncator(eps, wfn, slab_truncator)
	conv.truncate()
	eps.to_file(f_out)
