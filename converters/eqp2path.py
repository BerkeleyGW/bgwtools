#!/usr/bin/env python

from __future__ import division
from optparse import OptionParser
from scipy.spatial import KDTree
from numpy import *

usage = "usage: %prog [options] input output"
parser = OptionParser(usage)
parser.add_option('-c','--cubic',dest='cubic',default=False,action='store_true',
	help='move the data to the primitive cubic cell')
parser.add_option('-d','--data',dest='data',choices=['LDA','GW', 'diff'],default='GW',
	help='use LDA or GW energies, or the difference bewteen them')
parser.add_option('-k','--kpts',dest='kpts',type='str',default='0:0.5',
	help='set of kpts that should be plotted. Enter a list of kpts separated by columns. Eg: '+
	'0:0,0.5:1/3,1/3:0 will create the following path:'+
	'(0, 0, 0)->(0.5, 0, 0)->(1/3, 1/3, 0)->(0, 0, 0)'+
	'Points are assumed to be 3D, and omitted components are assumed to be null.')
parser.add_option('-m',dest='max',type='float',default='-1',
	help='only points closer than MAX w.r.t. the path are included')

(opts,args) = parser.parse_args()
if len(args)<2:
	parser.error('incorrect number of arguments')

def get_data_from_eqp(fname):
	X  = []
	Y  = []
	Z  = []
	Vl = []
	#read from eqp file
	f=open(args[0])
	item_num=-1
	for line in f.readlines():
		items=line.split()	
		if (item_num<0):
			kpts=map(float,items[0:3])
			#Kpts+=kpts
			X.append(kpts[0])
			Y.append(kpts[1])
			Z.append(kpts[2])
			items_cnt=int(items[3])
			Vl.append({})
			item_num=0
		else:
			band_idx=int(items[1])
			if opts.data=='GW':
				energ=float(items[3])
			elif opts.data=='LDA':
				energ=float(items[2])
			else:
				energ=float(items[3])-float(items[2])
			Vl[-1][band_idx]=energ
			item_num+=1
			if (item_num==items_cnt):
				item_num=-1
	bands=[v.keys() for v in Vl]
	bands=array(sorted(unique(array(bands))))
	V = ma.masked_all( (len(X),len(bands)) )
	i=0
	for v_ in Vl:
		indices = [where(bands==idx)[0][0] for idx in v_.keys()]
		V[i,:][indices] = v_.values()
		V.mask[i,:][indices] = False
		i+=1

	return array(X),array(Y),array(Z), V

def create_path(X,Y,Z,V):
	tree = KDTree(column_stack((X,Y,Z)))
	N = len(kpts)
	data_x = []
	data_y = []
	data_z = []
	bands_cnt = V.shape[1]
	data_v = [[] for i in range(bands_cnt)]
	last_idx = -1
	for n in range(N-1):
		N_fine=1000
		kpt_i = kpts[n]
		kpt_f = kpts[n+1]
		xi = linspace(kpt_i[0], kpt_f[0], N_fine)
		yi = linspace(kpt_i[1], kpt_f[1], N_fine)
		zi = linspace(kpt_i[2], kpt_f[2], N_fine)
		dist0,idx0 = [tree.query(column_stack((xi,yi,zi)))[i] for i in [0,1]]
		for m in range(N_fine):
			if opts.max>0:
				if dist0[m]>opts.max: continue
			if idx0[m] == last_idx: continue
			last_idx=idx0[m]
			data_x.append(X[last_idx])
			data_y.append(Y[last_idx])
			data_z.append(Z[last_idx])
			for i in range(bands_cnt):
				if not V.mask[last_idx][i]:
					data_v[i].append(V[last_idx][i])
				else:	
					data_v[i].append(NaN)
	data_x = array(data_x)
	data_y = array(data_y)
	data_z = array(data_z)
	data_v = ma.array(data_v, mask=isnan(data_v)).T
	return data_x, data_y, data_z, data_v

kpts=[]
labels=[]
kpts_str = opts.kpts
for items in kpts_str.split(':'):
	vals = items.split(',')
	if len(vals)>3:
		vals_pts = vals[:3]
	else:
		vals_pts = vals
	try:
		kpt = map(eval,vals_pts)
	except:
		kpt = map(float,vals_pts)
	while len(kpt)<3: kpt.append(0.0)
	kpts.append(kpt)
	if len(vals)>3:
		labels.append(vals[3])
	else:
		labels.append(' ')

if len(kpts)<2:
	kpts=[[0,0,0],[0.5,0.5,0]]
	labels=[1,2]

print kpts
print labels

X,Y,Z,V = get_data_from_eqp(args[0])
x,y,z,v = create_path(X,Y,Z,V)
save_data=column_stack([x,y,z,v])
savetxt(args[1], save_data, '%10.5f')
