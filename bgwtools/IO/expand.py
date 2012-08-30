#!/usr/bin/env python

from __future__ import division
import sys
from numpy import *
import numpy.core as core
import time
import cPickle
import scipy.interpolate
import itertools 

eps_small=1e-3
TOL_Small=1e-6
max_pts=8
verb=False

def tprint(str_):
	print time.strftime("[%H:%M:%S]", time.localtime()), str_

def triang_interp(x_,y_,eps_,q, dim):
	lx_ = len(x_)
	idx=tuple(arange(lx_))
	if lx_>2:
		for p_ in itertools.combinations(idx,3):
			#print p_
			p=asarray(p_)
			x=x_[p]
			y=y_[p]
			eps=eps_[p]
			A=ones((3,3))
			A[:,0] = x
			A[:,1] = y
			if linalg.det(A)<TOL_Small: continue
			try:
				a=linalg.solve(A,eps)
				return a[0]*q[0] + a[1]*q[1] + a[2]
			except:
				#raise Warning('Singular matrix')
				continue

	#print 'INFO: using linear interpolation/nn'
	if dim<0 or lx_<2:
		return eps_[0]

	if dim==0: x=x_
	else: x=y_

	#center points around q
	x = x - q[dim]
	srt = argsort(x)
	x=x[srt]

	if x[1]==x[0]:
		#points are useless, return nn
		return eps_[0]

	y=eps_[srt]
	if lx_==2:
		#no fit is necessary
		#a=(y[1]-y[0])/(x[1]-x[0])
		b=(y[0]*x[1] - y[1]*x[0])/(x[1]-x[0])
	else:
		a,b=polyfit(x,y,1)
	return b

class epsmat_intp:
	def __init__(self):
		pass

	def prepare(self, wfn, epsmat, kpts_fi, verb=True):
		self.wfn = wfn
		self.epsmat = epsmat
		self.kpts_fi = kpts_fi
		self.kpts_fold = epsmat.q
		kpts_fold = self.kpts_fold

		max_k = wfn.kgrid[0]*wfn.kgrid[1]*wfn.kgrid[2]
		kpts_unfold = empty((max_k,3))
		self.kpts_unfold = kpts_unfold
		self.n_co = 0

		if verb: tprint('Unfolding BZ')
		kpts_unfold[0] = kpts_fold[0]

		for ir in xrange(self.wfn.nk):
			for it in xrange(self.wfn.ntran):
				tmpf = dot(self.wfn.mtrx[:,:,it], self.wfn.kpt[ir,:])
				#put tmpf in box BZ
				tmpf = tmpf - floor(tmpf)
				if not any( all(abs(self.kpts_unfold[:self.n_co,:]-tmpf)<TOL_Small,axis=1) ):
					self.kpts_unfold[self.n_co,:] = tmpf
					#map_co[n_co] = ir
					self.n_co += 1
		
		if verb: tprint('Full box BZ has %s points'%(self.n_co))

		pts_ = core.records.fromarrays(self.kpts_unfold[:self.n_co].T)
		order = argsort(pts_)
		self.kpts_unfold = self.kpts_unfold[order]
		#map_co = map_co[order]

		if verb: tprint('Looking for points to interpolate')
		self.pts_intp = []
		n_fi = len(kpts_fi)
		for n in xrange(n_fi):
			q = self.kpts_fi[n]
			q_box = self.kpts_fi[n] - floor(self.kpts_fi[n])

			diff = self.kpts_unfold - q_box
			diff = around(diff) - diff
			cond = abs(diff) < TOL_Small
			is_coarse = any( all(cond,axis=1) )
			if not is_coarse:
				self.pts_intp.append(q)

		self.pts_intp = array(self.pts_intp)
		self.dq_fi = (self.kpts_fi[2]-self.kpts_fi[1])[1]
		self.dq_co = (self.kpts_fold[2]-self.kpts_fold[1])[1]
		
		if verb:
			print ' %s points will be interp., %d will be copyied'\
			%(self.pts_intp.shape[0], n_fi-self.pts_intp.shape[0])
			print ' dq_fi=', self.dq_fi
			print ' dq_co=', self.dq_co

	#saves at fname
	def interpolate(self, pts, fname=None):
		epsmat = self.epsmat
		wfn = self.wfn
		pts_intp = self.pts_intp
		kpts_fi = self.kpts_fi
		kpts_fold = self.kpts_fold
		dq_co = self.dq_co
		dq_fi = self.dq_fi

		#dummy loop
		for idx_inpt in pts:
			q = pts_intp[idx_inpt]
			tprint('Dealing with point #%d: %s'%(idx_inpt,str(q)))

			#calculating ekin
			gk = epsmat.gvec_k + q.reshape((1,3))
			u = expand_dims(gk,1)
			v = expand_dims(gk,2)
			uv= u*v
			KE = tensordot(uv, wfn.bdot, axes=([1,2],[0,1]))
			G = epsmat.gvec_k.T
			V = core.records.fromarrays([KE,G[0],G[1],G[2]])
			
			#eps2G[n] = (index of) G vector associated with nth element of eps
			#note: isort = eps2G+1, isort_i = G2eps+1
			eps2G = argsort(V)
			G2eps = zeros_like(eps2G)
			G2eps[eps2G] = arange(len(eps2G))
			
			nmtx = sum(KE<epsmat.ecuts)
			print '\t\tNeed %d,%d G,Gp vectors'%(nmtx,nmtx)
			
			delta = (q+TOL_Small) % dq_co - TOL_Small
			dim=-1	
			if (abs(delta[0])<TOL_Small):
				if (abs(delta[1])<TOL_Small):
					raise Exception('This is a coarse point!')
				else:
					#Interpolate on Y (X,Z fixed)
					dim=1
			else:
				if (abs(delta[1])<TOL_Small):
					#Interpolate on X (Y,Z fixed)
					dim=0
				else:
					#Different X & Y: vert. and horiz. interp.
					dim=-1
			print '\t\tInterpolation type: %d'%(dim)

			#Finding closets points, calculating all distances simultaneously
			diff = kpts_fold - q

			#0=cartesian, 1=wfn.bdot
			metric=0
			if metric!=0:
				u = expand_dims(diff,1)
				v = expand_dims(diff,2)
				uv = u*v
				dist = tensordot(uv, wfn.bdot, axes=([1,2],[1,0]))
			else:
				dist = sum(diff**2, axis=1)

			order = argsort(dist)[:10]
			range_order = arange(len(order))
			
			#eps2eps_co[i,n] = (index of) matrix element of eps (of a coarse
			# q-point i) associated with eps[n] (of the interpolated q-point)
			# Note: eps2eps is not isort_i!
			eps2eps_co = []
			#is_avail[i,n] = True if the G vector associated with the nth element
			# of eps was required for the coarse point i
			is_avail = []
			for close_co in order:
				qpt_co = close_co #assume WFN kpts == epsmat qpts!
				#isort_i could also be called "G2eps"
				map_ = epsmat.isort_i[qpt_co][eps2G]-1
				eps2eps_co.append(map_)
				is_avail.append(map_ < epsmat.nmtx[qpt_co])
			eps2eps_co = asarray(eps2eps_co)
			is_avail = asarray(is_avail)

			eps=empty((nmtx,nmtx))
			#eps_i, eps_j: loop over row/columns of eps
			#eps_co_i, eps_co_j: row/column of eps of a set of coarse pts
			for eps_i in range(nmtx):
				cond_i = is_avail[:,eps_i]
				#print 'G=',epsmat.gvec_k[indx_G]
				for eps_j in range(nmtx):
					cond_j = is_avail[:,eps_j]
					#print '\tGp=',epsmat.gvec_k[indx_Gp]
					cond = cond_i & cond_j
					if not(cond.any()):
						print 'WARNING: no point available for interpolation (%d,%d)'%(eps_i,eps_j)
						eps[eps_i,eps_j] = eps_small
					else:
						#selection: which pts of "order" should we take?
						# apply condition "cond" and limit number of points
						selection = range_order[cond][:max_pts]
						
						eps_co_i = eps2eps_co[:,eps_i]
						eps_co_j = eps2eps_co[:,eps_j]

						eps_list = asarray([epsmat.epsmat[order[n]]\
							[eps_co_i[n], eps_co_j[n]] for n in selection])
						x = asarray(kpts_fold[order[selection],0])
						y = asarray(kpts_fold[order[selection],1])

						eps[eps_i, eps_j] = triang_interp(x,y,eps_list, q, dim)

						if verb:
							print
							print eps_list
							print x
							print y
							print eps[eps_i, eps_j]

			#calculate ekin	
			g = epsmat.gvec_k + q.reshape((1,3))# + epsmat.q[0]
			u = expand_dims(g,1)
			v = expand_dims(g,2)
			uv = u*v
			ekin = tensordot(uv, wfn.bdot, axes=([1,2],[0,1]))

			isort = eps2G+1
			isort_i = G2eps+1
			
			tprint('\tDumping File')
			fname='eps_%d.pkl'%(idx_inpt)
			f=open(fname,'wb')
			cPickle.dump(q, f, cPickle.HIGHEST_PROTOCOL)
			cPickle.dump(isort, f, cPickle.HIGHEST_PROTOCOL)
			cPickle.dump(isort_i, f, cPickle.HIGHEST_PROTOCOL)
			cPickle.dump(ekin, f, cPickle.HIGHEST_PROTOCOL)
			cPickle.dump(eps, f, cPickle.HIGHEST_PROTOCOL)
			f.close()

			#epsmat.insert_element(q, isort, isort_i, ekin, eps)
			#epsmat.to_file('eps_%d'%(idx_inpt))

if __name__=='__main__':
	from wfn import wfnIO
	from epsmat import epsmatIO

	if len(sys.argv)!=3:
		print 'Usage: %s [WFN] [SCF_KPOINTS]'%(sys.argv[0])
		print 'Note: epsmat.pkl (produced by dump_eps.py) must be present' 
		sys.exit()

	print time.strftime("%a, %d %b %Y, %H:%M:%S", time.localtime())

	fname = sys.argv[2]
	tprint('Loading list of fine points from %s'%(fname))
	fkpts = open(fname)
	n_fi = int(fkpts.readline())
	kpts_fi = empty((n_fi,3))
	for n in xrange(n_fi):
		line = fkpts.readline()
		kpts_fi[n,:] = array(map(float,line.split()[5:8]))
	fkpts.close()

	tprint('Reading WFN header from %s'%(sys.argv[1]))
	wfn = wfnIO(sys.argv[1])

	tprint('Loading epsmat file from %s'%(sys.argv[2]))
	f=open('epsmat.pkl','rb')
	epsmat = cPickle.load(f)
	f.close()
	#epsmat = epsmatIO(sys.argv[2])

	tprint('Initializing interpolation instance')
	epsi = epsmat_intp()
	epsi.prepare(wfn, epsmat, kpts_fi)

	epsi.wfn.f = None
	epsi.epsmat.f = None	
	f=open('epsmat_interp.pkl','wb')
	cPickle.dump(epsi, f, cPickle.HIGHEST_PROTOCOL)
	f.close()
	
	#tprint('Interpolating points')
	#epsi.interpolate([22])

	tprint('Done')

