#!/usr/bin/env python

#silicon
#BAND=5
#kz=0.03100

from optparse import OptionParser
import matplotlib.pyplot as p
from scipy.spatial import KDTree
#from scipy.interpolate import *
from matplotlib.mlab import griddata
from numpy import *

imshow_interp={'nn':'nearest','linear':'bilinear'}

usage = "usage: %prog [options] file1 [file2] [...] [file N]"
parser = OptionParser(usage)
parser.add_option('-v','--val',dest='val',default=0.0,type='float',
	help='value that the fixed column must have. Defaults to VAL=0')
parser.add_option('-b','--band',dest='band',type='int',default=2,
	help='band to plot. default=2 (graphene)')
parser.add_option('-d','--diff',dest='diff',default=False,action='store_true',
	help='plot the difference of the first two files')
parser.add_option('-i','--interp',dest='interp',type='choice',choices=['nn','linear'], default='nn',
	help='set the interpolation (nn,linear)')
parser.add_option('-c','--cubic',dest='cubic',default=False,action='store_true',
	help='move the data to the primitive cubic cell')
parser.add_option('--min',dest='min',type='float',
	help='set minimum value for the plot')
parser.add_option('--max',dest='max',type='float',
	help='set maximum value for the plot')
parser.add_option('-k','--kpts',dest='kpts',type='str',default='0,0:0.5,0.5',
	help='set of kpts that should be plotted. Enter a list of kpts separated by columns. Eg: '+
	'0,0:0.5,0.5 will connect point (0,0) to (0.5,0.5). You can put the number of kpts that you want')

(opts,args) = parser.parse_args()
if len(args)<1:
	parser.error('incorrect number of arguments')

#get the raw data
def get_data_from_file(fname):
	data = loadtxt(fname)
	idx = data[:,3]==opts.band
	data=take(data, [0,1,2,4], axis=1)
	Band = data[idx,:]

	cond = abs(Band[:,2]-opts.val)<1e-6
	Band = Band[cond,:]

	X = Band[:,0]
	Y = Band[:,1]
	if opts.cubic:
		X = X-X.round()
		Y = Y-Y.round()
	V = Band[:,3]

	return X,Y,V

def local_interp(X,Y,V, xi,yi):
	TOL = 1e-3
	local_X=[]
	local_Y=[]
	local_V=[]
	
	for pt_X,pt_Y,pt_V in zip(X,Y,V):
		if ((abs(xi-pt_X)<TOL)&(abs(yi-pt_Y)<TOL)).any():
			local_X.append(pt_X)
			local_Y.append(pt_Y)
			local_V.append(pt_V)
	f = interp2d(array(local_X), array(local_Y), array(local_V), 'linear')
	return f(xi,yi)

#perform interpolation
def plot_data(X,Y,V):
	tree = KDTree(column_stack((X,Y)))
	last_pt=0
	N = len(kpts)
	for n in range(N-1):
		print 'n=',n
		cnt=150
		kpt_i = kpts[n]
		kpt_f = kpts[n+1]
		xi = linspace(kpt_i[0], kpt_f[0], cnt)
		yi = linspace(kpt_i[1], kpt_f[1], cnt)
		vi = zeros_like(xi)
		dist0,idx0 = [tree.query(column_stack((xi,yi)))[i] for i in [0,1]]
		idx = tree.query_ball_point(column_stack((xi,yi)),2e-2)
		for m in range(cnt):
			#print 'm=',m
			#print idx[m]
			if dist0[m]<1e-6:
				vi[m]=V[idx0[m]]
				continue
			if len(idx[m])>1:
				tmp=0
				#rbf=Rbf(X[idx[m]], Y[idx[m]], V[idx[m]] )
				#vi[m]=rbf(xi[m],yi[m])
				try:
					vi[m]=griddata(X[idx[m]],Y[idx[m]],V[idx[m]],
					 array([xi[m]]*2),array([yi[m]]*2), 'linear')[0,0]
				except:
					vi[m]=V[idx0[m]]
			else:
				vi[m]=V[idx0[m]]
				
		p.plot(last_pt+arange(cnt),vi)
		last_pt += cnt-1

kpts=[]
kpts_str = opts.kpts
for items in kpts_str.split(':'):
	kpt = map(float,items.split(','))
	kpts.append(kpt)
if len(kpts)<2:
	kpts=[[0,0],[0.5,0.5]]
print kpts

#p.subplots_adjust(0.05,0.01,0.98,0.98, 0.02,0.02)
plots_cnt = len(args)
if (opts.diff): plots_cnt += 1
for n in range(len(args)):
	p.subplot(1, plots_cnt, n+1)
	fname=args[n]
	p.title(fname)
	X,Y,V = get_data_from_file(fname)
	plot_data(X,Y,V)
	#if n: p.yticks(())

if (opts.diff):
	p.subplot(1, plots_cnt, plots_cnt)
	fname1=args[0]
	X,Y,V1 = get_data_from_file(fname1)
	fname2=args[1]
	X,Y,V2 = get_data_from_file(fname2)
	p.title('%s - %s'%(fname2,fname1))
	plot_data(X,Y,V2-V1)
	#p.yticks(())

p.show()
