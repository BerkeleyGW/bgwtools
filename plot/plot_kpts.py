#!/usr/bin/env python

'''
Simple utility that show all the kpts that meet some criterium
'''

TOL_Small=1e-6

from numpy import *
import sys
import matplotlib.pyplot as plt
from optparse import OptionParser

usage = "usage: %prog [options] file"
parser = OptionParser(usage)
parser.add_option('-x',dest='x',default=1,type='int',
	help='index of the x column in your file. Defaults to x=1')
parser.add_option('-y',dest='y',default=0,type='int',
	help='index of the y column in your file. Defaults to y=x+1')
parser.add_option('-z',dest='z',default=0,type='int',
	help='index of the z column in your file. Defaults to z=y+1')
parser.add_option('-f','--fix',dest='fix',type='choice',choices=['x','y','z'], default='z',
	help='column that you want to hold fixed (x,y,z). Defaults to the z column')
parser.add_option('-v','--val',dest='val',default=0.0,type='float',
	help='value that the fixed column must have. Defaults to VAL=0')
parser.add_option('-p','--print',dest='_print',default=False,action='store_true',
	help='print the kpts that meet the condition')
parser.add_option('-c','--cubic',dest='cubic',default=False,action='store_true',
	help='move the data to the primitive cubic cell')

(opts,args) = parser.parse_args()
if len(args)!=1:
	parser.error('incorrect number of arguments')

x,y,z = opts.x,opts.y,opts.z
if (y<1):
	y=x+1
if (z<1):
	z=y+1
col_map={'x':x, 'y':y, 'z':z}

fname = args[0]
try:
	data_raw = loadtxt(fname)
except:
	print 'Cannot load file',fname
	raise

fix_col = col_map[opts.fix]-1
cond = abs(data_raw[:,fix_col]-opts.val)<TOL_Small
if not(cond.any()):
	raise Exception('No value to plot.\nMake sure your condition is'+
			'correct and your columns are well specified.')

data = data_raw[cond,:]
if opts._print:
	print data

axis=[]
keys=col_map.keys()
keys.sort()
for label in keys:
	if label!=opts.fix:
		axis.append(label)

plot_1 = data[:,col_map[axis[0]]-1]
plot_2 = data[:,col_map[axis[1]]-1]

if opts.cubic:
	plot_1 -= plot_1.round()
	plot_2 -= plot_2.round()

plt.plot(plot_1, plot_2, 'o')
plt.xlabel(axis[0])
plt.ylabel(axis[1])
plt.show()

