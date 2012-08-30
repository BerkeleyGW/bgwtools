#!/usr/bin/env python

from __future__ import division

from numpy import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from optparse import OptionParser, OptionGroup

############
#  Parser  #
############
usage = "usage: %prog [options] file_1 [file_2 ... file_N]"
parser = OptionParser(usage)

group_range = OptionGroup(parser, 'Data Range')
group_range.add_option('--ymin',dest='ymin', default=None, type='float',
	help='minimum value of y to plot')
group_range.add_option('--ymax',dest='ymax', default=None, type='float',
	help='maximum value of y to plot')
group_range.add_option('--xmin',dest='xmin', default=None, type='float',
	help='miminum value of x to plot')
group_range.add_option('--xmax',dest='xmax', default=None, type='float',
	help='maximum value of x to plot')
parser.add_option_group(group_range)

group_energy = OptionGroup(parser, 'Centering and Middle Energy')
group_energy.add_option('--has_order', default=False, action='store_true',
	help='use this flag if the data has an extra column informing the '
	'coordinate of each point. This flag is off by default.')
group_energy.add_option('--auto_FE_bands', default=0, type='int',
	help='use this flag if you want the script to automatically find the '
	'Fermi energy. If the value of auto_FE_bands is greater than zero, '
	'then then Fermi energy will be calculated by counting auto_FE_bands '
	'bands. For a typical graphene calculation, use 4.')
group_energy.add_option('--fermi', dest='FEs', type='str', default=None,
	help='unless you are working with graphene, you might want to enter '
	'the Fermi energies here as a list of comma-separated values. This way, '
	'all the graphs are centeres around zero. Eg: -1.5,-1.4')
parser.add_option_group(group_energy)

group_visual = OptionGroup(parser, 'Aesthetics')
group_visual.add_option('-t','--title',type='str',
	help='title')
group_visual.add_option('-l','--labels',type='str',
	help='comma-separated list of labels for the plots')
group_visual.add_option('--alphas', type='str', default='',
	help='comma-separated list of alphas (transparencies). Eg: 1.0,0.9')
group_visual.add_option('--colors', type='str', default='',
	help='comma-separated list of colors. Eg: red, blue')
group_visual.add_option('--cmaps', type='str', default='Reds, Greens, Blues',
	help='comma-separated list of cmaps. Each cmap must be a valid '
	'pyplot directive. This directive (turned on by default) overrides '
	'COLROS. Eg: Reds, Greens, Blues')
group_visual.add_option('--cmap_frac', type='float', default=0.7,
	help='fraction of the cmap to use. Useful if the cmap gets too bright.')
group_visual.add_option('--linestyles', type='str', default='',
	help='comma-separated list of line styles. Eg: -,:,;')
group_visual.add_option('--no_points', action='store_true', default=False,
	help='don\'t print points')
parser.add_option_group(group_visual)

parser.add_option('-o',dest='output',type='str',metavar='FILE',
	help='output the graph to FILE')
parser.add_option('--load', type='str',metavar='FILE',
	help='load external python file with configurations')

###################
#  Configuration  #
###################
(opts,args) = parser.parse_args()
if opts.load:
	execfile(opts.load)

if len(args)<1:
	parser.error('incorrect number of arguments')

labels=[]
if opts.labels:
	if len(opts.labels):
		labels=opts.labels.split(',')
colors=[]
if opts.colors:
	if len(opts.colors):
		colors=[s.strip() for s in opts.colors.split(',')]
cmaps=[]
if opts.cmaps:
	if len(opts.cmaps):
		cmaps=[plt.get_cmap(s.strip()) for s in opts.cmaps.split(',')]
linestyles=[]
if opts.linestyles:
	if len(opts.linestyles):
		linestyles=[s.strip() for s in opts.linestyles.split(',')]
alphas=[]
if opts.alphas:
	if len(opts.alphas):
		alphas=[s.strip() for s in opts.alphas.split(',')]
FEs=[]
if opts.FEs:
	if len(opts.FEs):
		FEs=map(eval, opts.FEs.split(','))

###########
#  Plot!  #
###########
def plot_file(fname, cnt, fact=1.0):
	print 'Plotting:', fname

	band_idx=3
	if opts.has_order: band_idx+=1
	data=loadtxt(fname)
	bands=data[:,band_idx:]
	num_pts=bands.shape[0]
	num_bands=bands.shape[1]

	if opts.has_order:
		X=data[:,0]
	else:
		X=arange(num_pts)

	if opts.auto_FE_bands>0:
		order = argsort(bands[num_pts//2,:]).tolist()
		idx1 = opts.auto_FE_bands - 1
		vb = order[idx1]
		cb = order[idx1+1]
		#if opts.has_order:
		#	FE = (bands[27,vb]+bands[27,cb])/2
		#else:
		FE = (min(bands[:,cb]) + max(bands[:,vb]))/2
		print '  Gap:', min(bands[:,cb]) - max(bands[:,vb])
	else:
		order = argsort(bands[0,:]).tolist()
		if cnt<len(FEs):
			FE = FEs[cnt]
		else:	
			FE = None
	if not (FE is None):
		print '  Fermi Energy:', FE
		bands[:,:] -= FE

	#c_min = min(bands[bands>0])
	#v_max = max(bands[bands<0])
	#print '  Gap arround E=0:',c_min-v_max

	label=fname
	if cnt<len(labels):
		label=labels[cnt]

	ls=None
	if cnt<len(linestyles):
		ls=linestyles[cnt]
	cmap=None
	if cnt<len(cmaps):
		cmap=cmaps[cnt]
	if cmap is None:
		color=None
		if cnt<len(colors):
			color=colors[cnt]
	alpha=None
	if cnt<len(alphas):
		alpha=alphas[cnt]

	for b in range(num_bands):
		if not (cmap is None):
			col=cmap(1-(order.index(b))/(num_bands))
			color=[min(1,c*opts.cmap_frac) for c in col]
		plot_data = plt.plot(X, bands[:,b], '-o', lw=2, ms=4, mew=0)
		if not (color is None):
			plot_data[0].set_color(color)
			plot_data[0].set_markerfacecolor(color)
			if opts.no_points:
				plot_data[0].set_markersize(0)
		if not (alpha is None):
			plot_data[0].set_alpha(color)
		if not (ls is None):
			plot_data[0].set_linestyle(ls)
		if order.index(b)==num_bands//2:
			plot_data[0].set_label(label)

	if cnt==0:
		plt.axhline(0, ls=':', color='black')
	print 

fcnt=0	
for fname in args:
	plot_file(fname, fcnt)
	fcnt+=1

plt.legend()

F = plt.gcf()
DPI = F.get_dpi()
DefaultSize = F.get_size_inches()
F.set_size_inches( (DefaultSize[0], DefaultSize[1]*0.75) )
plt.ylabel('Energy (eV)')

#plt.xlim(16,32)
#plt.ylim(-4,3)
#plt.savefig('bandstruc2.pdf')
plt.show()
