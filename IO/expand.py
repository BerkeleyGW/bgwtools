#!/usr/bin/env python

from __future__ import division

import sys
from wfn import wfnIO
from numpy import *
import numpy.core as core
import time

TOL_Small=1e-6

def tprint(str_):
	print time.strftime("[%H:%M:%S]", time.localtime()), str_

if len(sys.argv)!=3:
	print 'Usage: %s [WFN] [SCF_KPOINTS]'%(sys.argv[0])
	sys.exit()

print time.strftime("%a, %d %b %Y, %H:%M:%S", time.localtime())

tprint('Reading WFN header')

wfn = wfnIO(sys.argv[1])

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
#order = sort(pts_).view(float64).reshape((n_co,3))

tprint('Full box BZ has %s points'%(n_co))
#print 'Points:'
#print pts
#print 'Mapping:'
#print map

tprint('Populating cells')
cells = ones(wfn.kgrid, dtype=int)*(-1)
#one_delta = 1/wfn.kgrid
for n in xrange(n_co):
	idx = tuple(floor((kpts_co[n]+TOL_Small)*wfn.kgrid))
	if (cells[idx]!=-1):
		raise ValueError('Cell already populated')
	cells[idx] = n
#print cells

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
#bz_map = map[order]
#order = sort(pts_).view(float64).reshape((n_co,3))

tprint('Full BZ build')
#print 'BZ Points:'
#print bz_kpts

tprint('Loading kpoints')
fkpts = open(sys.argv[2])
n_fi = int(fkpts.readline())
kpts_fi = empty((n_fi,3))
for n in xrange(n_fi):
	line = fkpts.readline()
	kpts_fi[n,:] = array(map(float,line.split()[5:8]))
fkpts.close()

g = mgrid[0:wfn.kgrid[0], 0:wfn.kgrid[1], 0:wfn.kgrid[2]]
gx,gy,gz = g

#get kpts that are in the neighboring cells
def get_local_indices(idx):
	cell_dist=1

	'''
	deltas = array([g[i]-idx[i] for i in range(3)])
	cond_fix = deltas
	deltas = deltas
	cond=((gx-idx[0])<cell_dist) & ((gy-idx[1])<cell_dist) & ((gz-idx[2])<cell_dist)
	'''
	cell_min = array(idx-cell_dist, dtype=int)
	cell_max = array(idx+cell_dist+1, dtype=int)
	for dim in xrange(3):
		if 2*(cell_dist+1)>wfn.kgrid[dim]:
			cell_min[dim] = 0
			cell_max[dim] = wfn.kgrid[dim]

	indices=[]
	for i in xrange(cell_min[0],cell_max[0]):
		while i<0: i+=wfn.kgrid[0]
		while i>=wfn.kgrid[0]: i-=wfn.kgrid[0]
		for j in xrange(cell_min[1],cell_max[1]):
			while j<0: j+=wfn.kgrid[1]
			while j>=wfn.kgrid[1]: j-=wfn.kgrid[1]
			for k in xrange(cell_min[2],cell_max[2]):
				while k<0: i+=wfn.kgrid[2]
				while k>=wfn.kgrid[2]: k-=wfn.kgrid[2]
				indices.append(array([i,j,k]))
				#kpts.append(kpts_co[cells[i,j,k]])
	return indices
				
	

#tprint('Fine Points:')
for n in xrange(n_fi):
#for n in xrange(3):
	q = kpts_fi[n] - floor(kpts_fi[n])
	#print 'Fine point n:', q_box

	idx_arr = floor((q+TOL_Small)*wfn.kgrid)
	idx = tuple(idx)
	if (cells[idx]==-1):
		raise ValueError('Empty cell!')

	local_indices = get_local_indices(idx_arr)
	sys.exit()

	delta = (kpts_co - q_box)
	delta = delta - around(delta)
	is_small = all(abs(delta)<TOL_Small, axis=1)
	if any(is_small):
		idx=where(is_small)[0][0]
		#print 'Coarse point:', kpts_co[idx]
	else:
		#print 'No coarse point found.'
		pass


