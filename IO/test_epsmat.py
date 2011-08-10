#!/usr/bin/env python

from __future__ import division

TOL_Small = 1e-6

import sys
from wfn import wfnIO
from epsmat import epsmatIO
from numpy import *
import numpy.core as core
import time

if len(sys.argv)!=3:
	print 'Usage: %s [WFN] [epsmat]'%(sys.argv[0])
	sys.exit()

print time.strftime("%a, %d %b %Y, %H:%M:%S", time.localtime())

print('Reading WFN header')
wfn = wfnIO(sys.argv[1])

print('Reading epsmat file')
epsmat = epsmatIO(sys.argv[2])

qpt=0

'''
#Check if computed g-vectores have KE<ecuts
for n in xrange(epsmat.nmtx[qpt]):
	srt = epsmat.isort[qpt][n]-1
	G = epsmat.gvec_k[srt, :]
	gk = G + epsmat.q[qpt]
	m = dot(wfn.bdot,gk)
	KE = dot(gk,m)
	print G, KE, epsmat.epsmat[qpt][n, n]
	#print G, self.epsmat[0][n, n]
'''	

print('Searching vectors')

vects=[]
vects_i=[]

#for G_ in arange(-3,4):
for G_ in arange(-1,2):
	G = array([0,0,G_])
	cond = all(abs(epsmat.gvec_k[:,:] - G) < TOL_Small, axis=1)
	if not any(cond): continue
	srt = where(cond)[0][0]
	#print idx,G,epsmat.gvec_k[srt]
	#break
	gk = G + epsmat.q[qpt]
	m = dot(wfn.bdot,gk)
	KE = dot(gk,m)
	n = epsmat.isort_i[qpt][srt]-1
	if n>=epsmat.nmtx[qpt]:
		print 'Vector',G,'out of bonds'
		continue
	print epsmat.gvec_k[srt], KE, epsmat.epsmat[qpt][n,n]
	vects.append(srt)
	vects_i.append(n)

print

for i in range(len(vects)):
	isrt = vects[i]
	isrt_i = vects_i[i]
	Gi = epsmat.gvec_k[isrt]
	qi = Gi + epsmat.q[qpt]
	for j in range(len(vects)):
		jsrt = vects[j]
		jsrt_i = vects_i[j]
		Gj = epsmat.gvec_k[jsrt]
		qj = Gj + epsmat.q[qpt]

		print Gi,Gj,epsmat.epsmat[qpt][isrt_i][jsrt_i]



