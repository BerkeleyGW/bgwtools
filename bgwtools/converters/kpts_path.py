#!/usr/bin/env python

from __future__ import division
from optparse import OptionParser
from scipy.spatial import KDTree
from numpy import *

usage = "usage: %prog [options] input output"
parser = OptionParser(usage)
parser.add_option('-c','--cubic',dest='cubic',default=False,action='store_true',
	help='move the data to the primitive cubic cell')
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

#get the raw data
def get_data_from_file(fname):
	data = loadtxt(fname, dtype='f8,f8,f8,i4,f8')
	bands = unique(data['f3'])
	#note: we assume the second band is fully occupied!
	idx = data['f3']==bands[0]
	X = data['f0'][idx]
	Y = data['f1'][idx]
	Z = data['f2'][idx]
	if opts.cubic:
		X = X-X.round()
		Y = Y-Y.round()
		Z = Z-Z.round()
	Vl=[]
	for b in bands:
		idx = data['f3']==b
		Vl.append(data['f4'][idx])
	return X,Y,Z,Vl

def create_path(X,Y,Z,Vl):
	tree = KDTree(column_stack((X,Y,Z)))
	N = len(kpts)
	data_x = []
	data_y = []
	data_z = []
	Vl_cnt = len(Vl)
	data_vl = [[] for i in range(Vl_cnt)]
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
			for i in range(Vl_cnt):
				data_vl[i].append(Vl[i][last_idx])

	return data_x, data_y, data_z, data_vl

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

X,Y,Z,Vl = get_data_from_file(args[0])
x,y,z,vl = create_path(X,Y,Z,Vl)
save_data = column_stack([array(x),array(y),array(z)]+[array(v) for v in vl])
#print save_data
savetxt(args[1], save_data, '%8.5f')
