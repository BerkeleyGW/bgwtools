#!/usr/bin/env python

from __future__ import division

import sys
from wfn import wfnIO
from epsmat import epsmatIO
from numpy import *
import numpy.core as core
import time
import cPickle
import scipy.interpolate

eps_small=1e-3
TOL_Small=1e-6

def tprint(str_):
	print time.strftime("[%H:%M:%S]", time.localtime()), str_

if len(sys.argv)!=3:
	print 'Usage: %s [WFN] [SCF_KPOINTS]'%(sys.argv[0])
	sys.exit()

print time.strftime("%a, %d %b %Y, %H:%M:%S", time.localtime())

tprint('Reading WFN header')

wfn = wfnIO(sys.argv[1])

kpts_fold = wfn.kpt.T

max_k = wfn.kgrid[0]*wfn.kgrid[1]*wfn.kgrid[2]
kpts_co = empty((max_k,3))
map_co = empty(max_k)
n_co=0

tprint('Unfolding BZ')

kpts_co[0,:] =  wfn.kpt[:,0]
for ir in xrange(wfn.nk):
	#print ir, wfn.kpt[:,ir]
	for it in xrange(wfn.ntran):
		tmpf = dot(wfn.mtrx[:,:,it], wfn.kpt[:,ir])
		#put tmpf in box BZ
		tmpf = tmpf - floor(tmpf)
		if not any( all(abs(kpts_co[:n_co,:]-tmpf)<TOL_Small,axis=1) ):
			kpts_co[n_co,:] = tmpf
			map_co[n_co] = ir
			n_co+=1

pts_ = core.records.fromarrays(kpts_co[:n_co].T)
order = argsort(pts_)
kpts_co = kpts_co[order]
map_co = map_co[order]

tprint('Full box BZ has %s points'%(n_co))

tprint('Populating cells')
cells = ones(wfn.kgrid, dtype=int)*(-1)
for n in xrange(n_co):
	idx = tuple(floor((kpts_co[n]+TOL_Small)*wfn.kgrid))
	if (cells[idx]!=-1):
		raise ValueError('Cell already populated')
	cells[idx] = n

ncell=2
tprint('Constructing 1st BZ')

#Now, put the points in the 1st BZ
fq = empty((n_co,3))
tmpfm = empty((n_co,3))
bz_kpts = empty((n_co,3))
lmin = ones(n_co)*inf

for i1 in range(-ncell,ncell):
	fq[:,0] = kpts_co[:,0] - i1
	for i2 in range(-ncell,ncell):
		fq[:,1] = kpts_co[:,1] - i2
		for i3 in range(-ncell,ncell):
			fq[:,2] = kpts_co[:,2] - i3
			u = expand_dims(fq,1)
			v = expand_dims(fq,2)
			uv = u*v
			length = tensordot(uv, wfn.bdot, axes=([1,2],[1,0]))
			cond=length<lmin
			if cond.any():
				lmin[cond] = length[cond]
				bz_kpts[cond] = fq[cond]

bz_kpts_ = core.records.fromarrays(bz_kpts.T)
order = argsort(bz_kpts_)
bz_kpts = bz_kpts[order]
bz_map = map_co[order]

tprint('Full BZ build')

tprint('Loading SCF_KPOINTS')
fkpts = open(sys.argv[2])
n_fi = int(fkpts.readline())
kpts_fi = empty((n_fi,3))
for n in xrange(n_fi):
	line = fkpts.readline()
	kpts_fi[n,:] = array(map(float,line.split()[5:8]))
fkpts.close()

tprint('Looking for points to interpolate')
pts_intp = []
for n in xrange(n_fi):
	q = kpts_fi[n]
	q_box = kpts_fi[n] - floor(kpts_fi[n])

	#diff = bz_kpts - q
	diff = kpts_co - q_box
	diff = around(diff) - diff
	cond = abs(diff) < TOL_Small
	is_coarse = any( all(cond,axis=1) )
	if not is_coarse:
		pts_intp.append(q)

pts_intp = array(pts_intp)
print ' %s points will be interpolated, %d will be copyied'\
	%(pts_intp.shape[0], len(kpts_fi)-pts_intp.shape[0])
dq_fi = (kpts_fi[1]-kpts_fi[0])[1]
dq_co = (kpts_fold[1]-kpts_fold[0])[1]
print ' dq_fi=',dq_fi
print ' dq_co=',dq_co

tprint('Loading epsmat file')
f=open('/tmp/epsmat.pkl','rb')
epsmat = cPickle.load(f)
f.close()

tprint('Interpolating points')
max_pts=5

import itertools, random
def triang_interp(x_,y_,eps_,q, dim):
	lx_ = len(x_)
	indx=tuple(arange(lx_))
	if lx_>2:
		for p_ in itertools.combinations(indx,3):
			#print p_
			p=array(p_)
			x=x_[p]
			y=y_[p]
			eps=eps_[p]
			A=ones((3,3))
			A[:,0] = x
			A[:,1] = y
			#if linalg.det(A)<TOL_Small: continue
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
	srt = argsort(x)
	x=x[srt]
	y=eps_[srt]
	coeff=polyfit(x,y,1)
	return coeff[0]*q[dim]+coeff[1]

	#raise Exception('Interpolation Failed')

idx_inpt = 22
##dummy loop
if 1:
#for ind_inpt in range(pts_intp.shape[0]):
	q = pts_intp[idx_inpt]
	#q=epsmat.q[4]
	q_box = kpts_fi[n] - floor(kpts_fi[n])
	tprint('Dealing with point #%d: %s'%(idx_inpt,str(q)))

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
	u = expand_dims(diff,1)
	v = expand_dims(diff,2)
	uv = u*v
	dist = tensordot(uv, wfn.bdot, axes=([1,2],[1,0]))
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
	eps2eps_co = array(eps2eps_co)
	is_avail = array(is_avail)

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
				print 'WARNING: no point available for interpolation'
				eps[eps_i,eps_j] = eps_small
			else:
				#selection: which pts of "order" should we take?
				# apply condition "cond" and limit number of points
				selection = range_order[cond][:max_pts]
				
				eps_co_i = eps2eps_co[:,eps_i]
				eps_co_j = eps2eps_co[:,eps_j]

				eps_list = array([epsmat.epsmat[order[n]]\
					[eps_co_i[n], eps_co_j[n]] for n in selection])
				x = array(kpts_fold[order[selection],0])
				y = array(kpts_fold[order[selection],1])

				#print
				#print eps_list
				#print x
				#print y
				eps[eps_i, eps_j] = triang_interp(x,y,eps_list, q, dim)
				#print eps[eps_i, eps_j]

	#calculate ekin	
	g = epsmat.gvec_k + q.reshape((1,3))# + epsmat.q[0]
	u = expand_dims(g,1)
	v = expand_dims(g,2)
	uv = u*v
	ekin = tensordot(uv, wfn.bdot, axes=([1,2],[0,1]))

	isort = eps2G+1
	isort_i = G2eps+1
	
	tprint('\tDumping File')
	f=open('eps_%d.pkl'%(idx_inpt),'wb')
	cPickle.dump(G_req, f, cPickle.HIGHEST_PROTOCOL)
	cPickle.dump(eps, f, cPickle.HIGHEST_PROTOCOL)
	f.close()

	#epsmat.insert_element(q, isort, isort_i, ekin, eps)
	#epsmat.to_file('eps_%d'%(idx_inpt))
				
tprint('Done')
