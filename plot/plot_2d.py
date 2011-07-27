#!/usr/bin/env python

#silicon
#BAND=5
#kz=0.03100

from optparse import OptionParser
import matplotlib.pyplot as p

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

(opts,args) = parser.parse_args()
if len(args)<1:
	parser.error('incorrect number of arguments')

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

def plot_data(X,Y,V):
	xmin,xmax = min(X),max(X)
	ymin,ymax = min(Y),max(Y)
	print 'Min Value:', min(V)
	print 'Max Value:', max(V)
	if not (opts.min is None):
		V[V<opts.min] = opts.min
	if not (opts.max is None):
		V[V>opts.max] = opts.max

	xi=linspace(xmin, xmax, 250)
	yi=linspace(ymin, ymax, 250)
	vi=griddata(X,Y,V,xi,yi,interp=opts.interp)

	p.imshow(vi, origin='lower', interpolation=imshow_interp[opts.interp], \
	#	vmin=-1, vmax=1,\
		extent=(xmin,xmax,ymin,ymax)
	)
	p.colorbar(orientation='horizontal')

p.subplots_adjust(0.05,0.01,0.98,0.98, 0.02,0.02)
plots_cnt = len(args)
if (opts.diff): plots_cnt += 1
for n in range(len(args)):
	p.subplot(1, plots_cnt, n+1)
	fname=args[n]
	p.title(fname)
	X,Y,V = get_data_from_file(fname)
	plot_data(X,Y,V)
	if n: p.yticks(())

if (opts.diff):
	p.subplot(1, plots_cnt, plots_cnt)
	fname1=args[0]
	X,Y,V1 = get_data_from_file(fname1)
	fname2=args[1]
	X,Y,V2 = get_data_from_file(fname2)
	p.title('%s - %s'%(fname2,fname1))
	plot_data(X,Y,V2-V1)
	p.yticks(())

p.show()
