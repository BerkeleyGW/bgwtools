#!/usr/bin/env python

from __future__ import division
from numpy import *

eps_small=1e-3
TOL_Small=1e-6


pts_ = core.records.fromarrays(kpts_co[:n_co].T)
order = argsort(pts_)
kpts_co = kpts_co[order]
map_co = map_co[order]

tprint('Full box BZ has %s points'%(n_co))

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

