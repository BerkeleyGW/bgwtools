#!/usr/bin/env python

from __future__ import division
from numpy import *
import pylab as p
import scipy
import scipy.signal
from optparse import OptionParser, OptionGroup

plot_dict={'abs':0, 'e2':1, 'e1':2, 'dos':3}
ylabels=[r'$\mathrm{Absorbance}$', r'$\varepsilon_2$',r'$\varepsilon_1$',r'$\mathrm{DOS}$']

############
#  Parser  #
############
usage = "usage: %prog [options] file"
parser = OptionParser(usage)

group_range = OptionGroup(parser, 'Data Range')
group_range.add_option('--ymin',dest='ymin', default=0.0, type='float',
	help='minimum value of y to plot')
group_range.add_option('--ymax',dest='ymax', type='float',
	help='maximum value of y to plot')
group_range.add_option('--xmin',dest='xmin', default=0.0, type='float',
	help='miminum value of x to plot')
group_range.add_option('--xmax',dest='xmax', default=7.0, type='float',
	help='maximum value of x to plot')
opts=plot_dict.keys()
opts.sort()
group_range.add_option('-p',dest='plot',default='abs',type='choice',choices=opts,
	help='specify what is to be plotted. Defaults to abs. Choices are: '+
	', '.join(opts))
parser.add_option_group(group_range)

group_smooth = OptionGroup(parser, 'Data Smoothing and Asymmetry Analysis')
group_smooth.add_option('-s', dest='smooth', default='none',
	choices=['none','wiener','gaussian','lorentzian'],
	help='smooth the data. Options are none, wiener, gaussian and lorentzian')
group_smooth.add_option('--window', dest='window', type='int', default=31,
	help='size of the window used in the Wiener filter')
group_smooth.add_option('--sigma', type='float', default=0.5,
	help='stardard deviation of the gaussian/lorentzian filter')
group_smooth.add_option('--asymmetry', dest='asymmetry',  default=False, action='store_true',
	help='calculate and plot asymmetry of the peak')
parser.add_option_group(group_smooth)

group_visual = OptionGroup(parser, 'Aesthetics')
group_visual.add_option('-t','--title',type='str',
	help='title')
group_visual.add_option('-l','--labels',type='str',
	help='semi-column-separated list of labels for the plots')
group_visual.add_option('--colors', type='str', default='',
	help='comma-separated list of colors. Eg: red, blue')
group_visual.add_option('--linestyles', type='str', default='',
	help='comma-separated list of line styles. Eg: -,:,;')
parser.add_option_group(group_visual)

parser.add_option('-o',dest='output',type='str',metavar='FILE',
	help='output the graph to FILE')
parser.add_option('--load', type='str',metavar='FILE',
	help='load external python file with configurations')

###################
#  Configuration  #
###################
(opts,args) = parser.parse_args()
delta_E_l = None #used for asymmetry analysis
if opts.load:
	execfile(opts.load)

if len(args)<1:
	parser.error('incorrect number of arguments')

#labels
labels=[]
if opts.labels:
	if len(opts.labels):
		labels=opts.labels.split(';')
#colors
colors=[]
if opts.colors:
	if len(opts.colors):
		colors=[s.strip() for s in opts.colors.split(',')]
#linestyles
linestyles=[]
if opts.linestyles:
	if len(opts.linestyles):
		linestyles=[s.strip() for s in opts.linestyles.split(',')]
if opts.output:
	fig_width_pt = 246.0*2.0	# Get this from LaTeX using \showthe\columnwidth
	inches_per_pt = 1.0/72.27               # Convert pt to inch
	#golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
	golden_mean = 0.5
	fig_width = fig_width_pt*inches_per_pt  # width in inches
	fig_height = fig_width*golden_mean      # height in inches
	fig_size =  [fig_width,fig_height]
	params = {'backend': 'eps',
		  'axes.titlesize': 18,
		  'axes.labelsize': 14,
		  'text.fontsize': 12,
		  'legend.fontsize': 12,
		  'xtick.labelsize': 10,
		  'ytick.labelsize': 10,
		  'text.usetex': True,
		  'font.family': 'serif',
		  'font.serif': 'Times' ,
		  'figure.figsize': fig_size}
	#p.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)
else:
	#p.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)
	params = {'figure.figsize': [10,5]}

#Initilize Plot
p.rcParams.update(params)
if opts.output:
	p.axes([0.1,0.16,0.95-0.1,0.95-0.2])

plot_idx = plot_dict[opts.plot]

###########
#  Plot!  #
###########
cnt=0
data=[]
data_max=0
for fname in args:
	d = genfromtxt(fname)
	x=d[:,0]
	dx = x[1]-x[0]

	if opts.smooth=='wiener':
		for n in range(1,d.shape[1]):
			d[:,n] = scipy.signal.wiener(d[:,n],opts.window)
	elif opts.smooth=='gaussian':
		import scipy.ndimage.filters
		sigma = opts.sigma/(d[1,0]-d[0,0])
		for n in range(1,d.shape[1]):
			d[:,n] = scipy.ndimage.filters.gaussian_filter1d(d[:,n],sigma)
	elif opts.smooth=='lorentzian':
		n_window = int( ceil(4*opts.sigma/dx) )
		if (n_window%2==0):
			n_window+=1
		half_window = (n_window-1)/2
		lor_x = (arange(n_window)-half_window)*dx
		lor_y = 1/pi * (opts.sigma/(lor_x**2 + opts.sigma**2))
		for n in range(1,d.shape[1]):
			d[:,n] = scipy.signal.convolve(d[:,n], lor_y, mode='same') *dx

	if plot_idx==0:
		y=d[:,0]*d[:,1]
	else:
		y=d[:,plot_idx]


	data.append(d)
	idx=argmax(y[(x<opts.xmax)&(x>opts.xmin)])+sum(x<=opts.xmin)
	print 'File:',fname
	print '  peak: %.3f, %.3f'%(x[idx],y[idx])
	data_max = max(data_max, max(y[(x>opts.xmin)&(x<opts.xmax)]))

	if opts.asymmetry:
		#calculate the ration of the right and left slopes
		if not delta_E_l: 
			delta_E_l = 0.08
			delta_E_r = 0.15
			window_E_l = 0.45
			window_E_r = 0.25

		delta_idx_l = int(delta_E_l/(x[1]-x[0]))
		delta_idx_r = int(delta_E_r/(x[1]-x[0]))
		window_idx_l = int(window_E_l/(x[1]-x[0]))
		window_idx_r = int(window_E_r/(x[1]-x[0]))
		x_r = x[idx + delta_idx_r : idx + delta_idx_r + window_idx_r]
		y_r = y[idx + delta_idx_r : idx + delta_idx_r + window_idx_r]
		x_l = x[idx - delta_idx_l - window_idx_l : idx - delta_idx_l]
		y_l = y[idx - delta_idx_l - window_idx_l : idx - delta_idx_l]
		#linear regression
		a_r,b_r = scipy.polyfit(x_r, y_r, 1)
		a_l,b_l = scipy.polyfit(x_l, y_l, 1)
		y_reg_r = scipy.polyval([a_r,b_r],x_r)
		y_reg_l = scipy.polyval([a_l,b_l],x_l)
		print ' asymmetry:',-a_r/a_l

	if cnt>=len(labels):
		label=fname
	else:
		label=labels[cnt]

	#plot
	plot_data = p.plot(x,y, label=label, lw=1.5)
	if cnt<len(colors):
		plot_data[0].set_color(colors[cnt])
	if cnt<len(linestyles):
		plot_data[0].set_linestyle(linestyles[cnt])
	if opts.asymmetry:
		p.plot(x[idx],y[idx], 'o' ,ms=8, mew=0, alpha=0.5,
			color=plot_data[0].get_color())
		#p.axvline(x[idx], ls=':', lw=1.5, alpha=0.75)
		p.plot([x_l[0],x_l[-1]],[y_reg_l[0],y_reg_l[-1]],color='black',lw=3,ls='--')
		p.plot([x_r[0],x_r[-1]],[y_reg_r[0],y_reg_r[-1]],color='black',lw=3,ls='--')

	cnt+=1

p.legend(loc='upper left')

if opts.title:
	p.title(opts.title)
else:
	p.title(r'$\mathrm{Absorbance}/\varepsilon$')

p.xlabel(r'$E\;(eV)$')
p.ylabel(ylabels[plot_idx])

p.xlim(opts.xmin, opts.xmax)
if not (opts.ymax is None):
	p.ylim(opts.ymin, opts.ymax)
else:
	p.ylim(opts.ymin, data_max*1.2)


if opts.output:
	p.savefig(opts.output, transparent=True)
else:
	p.show()
