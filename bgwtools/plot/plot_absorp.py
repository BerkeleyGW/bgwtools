#!/usr/bin/env python

from __future__ import division
from numpy import *
import scipy
import scipy.signal
from optparse import OptionParser, OptionGroup

wp=0.937985
#wp=0.942002

plot_dict={'abs':0, 'e2':1, 'e1':2, 'dos':3, 'abs2':-1}
ylabels=[r'Absorbance', r'$\varepsilon_2',r'$\varepsilon_1$',r'DOS']

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
group_range.add_option('-p',dest='plot_what',default='abs',type='choice',choices=opts,
	help='specify what is to be plotted. Defaults to abs. Choices are: '+
	', '.join(opts))
parser.add_option_group(group_range)

parser.add_option('-L', default=1.0, type='float',
	help='Length of the unit cell (Bohrs)')

group_smooth = OptionGroup(parser, 'Data Smoothing and Asymmetry Analysis')
group_smooth.add_option('-s', '--smooth', dest='smooth', default='none',
	choices=['none','wiener','gaussian','lorentzian'],
	help='smooth the data. Options are none, wiener, gaussian and lorentzian')
group_smooth.add_option('--window', dest='window', type='int', default=31,
	help='size of the window used in the Wiener filter')
group_smooth.add_option('--sigma', type='float', default=0.5,
	help='stardard deviation of the gaussian/lorentzian filter')
group_smooth.add_option('--asymmetry', dest='asymmetry',  default=False, action='store_true',
	help='calculate and plot asymmetry of the peak')
parser.add_option('--export', action='store_true',
	help='Export smoothed data')
parser.add_option_group(group_smooth)

group_visual = OptionGroup(parser, 'Aesthetics')
group_visual.add_option('-t','--title',type='str',
	help='title')
group_visual.add_option('-l','--labels',type='str',
	help='semi-column-separated list of labels for the plots')
group_visual.add_option('--colors', type='str', default='',
	help='comma-separated list of colors. Eg: red, blue')
group_visual.add_option('--ls', '--linestyles', type='str', default='', dest='linestyles',
	help='comma-separated list of line styles. Eg: -,:,;')
group_visual.add_option('--lw', '--linewidths', type='str', default='', dest='linewidths',
	help='comma-separated list of line styles. Eg: 1.5,1.5,2')
parser.add_option_group(group_visual)

parser.add_option('--noplot', dest='plot', default=True, action='store_false',
	help='don\'t plot anything (but print all the analysis)')
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

#linewidths
linewidths=[]
if opts.linewidths:
	if len(opts.linewidths):
		linewidths=[s.strip() for s in opts.linewidths.split(',')]

if len(args)<1:
	parser.error('incorrect number of arguments')

if opts.output or 1:
	fig_width_pt = 246.0*3.0	# Get this from LaTeX using \showthe\columnwidth
	inches_per_pt = 1.0/72.27               # Convert pt to inch
	#golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
	golden_mean = 0.5
	fig_width = fig_width_pt*inches_per_pt  # width in inches
	fig_height = fig_width*golden_mean      # height in inches
	fig_size =  [fig_width,fig_height]
	params = {'backend': 'eps',
		  'axes.titlesize': 20,
		  'axes.labelsize': 18,
		  'text.fontsize': 16,
		  'legend.fontsize': 16,
		  'xtick.labelsize': 16,
		  'ytick.labelsize': 16,
		  'text.usetex': True,
		  #'text.usetex': False,
		  #'font.family': 'serif',
		  #'font.serif': 'Times' ,
		  'font.family': 'sans-serif',
		  #'font.sans-serif': ['Helvetica'] ,
		  'figure.figsize': fig_size}
	#p.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)
else:
	#p.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)
	params = {'figure.figsize': [10,5]}

#Initilize Plot
if opts.plot_what:
	import pylab as p
	p.rcParams.update(params)
	#if opts.output:
	p.axes([0.1,0.14,0.95-0.1,0.95-0.2])

plot_idx = plot_dict[opts.plot_what]

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

	L=opts.L
	if isinstance(L,list):
		try:
			L=L[cnt]
		except:
			L=L[-1]
	
	if plot_idx==-1:
		#plot absorbance2, the extra column that I added, if available
		if d.shape[1]>3:
			y=d[:,3] * 0.0002681729 * L * 100
		else:
			y=d[:,0]*d[:,1] * 0.0002681729 * L * 100
	elif plot_idx==0:
		#plot absorbance
		y=d[:,0]*d[:,1] * 0.0002681729 * L * 100
	else:
		y=d[:,plot_idx]

	try:
		sig = opts.sigma[cnt]
	except:
		sig = opts.sigma

	#int_window = int( ceil(2*1/dx) )
	int_window=1
	i1 = 0
	#i2 = floor(opts.xmax/dx)
	#i2 = floor(14.5/dx)
	i2=-1
	sum_rule = sum(y[i1:i2])
	
	print 'x1=',x[i1],'x2=',x[i2],'sum=',sum_rule*(2.0/(pi*wp**2))
	
	if opts.smooth=='wiener':
		#for n in range(1,d.shape[1]):
		#	d[:,n] = scipy.signal.wiener(d[:,n],opts.window)
		y = scipy.signal.wiener(y,opts.window)
		y = y * sum_rule/sum(y[int_window:-int_window])

	elif opts.smooth=='gaussian':
		import scipy.ndimage.filters
		sigma = sig/dx
		y = scipy.ndimage.filters.gaussian_filter1d(y,sigma)
		#for n in range(1,d.shape[1]):
		#	d[:,n] = scipy.ndimage.filters.gaussian_filter1d(d[:,n],sigma)
		y = y * sum_rule/sum(y[int_window:-int_window])

	elif opts.smooth=='lorentzian':
		n_window = int( ceil(100*sig/dx) )
		if (n_window%2==0):
			n_window+=1
		half_window = (n_window-1)/2
		lor_x = (arange(n_window)-half_window)*dx
		#lor_x = arange(n_pts)-
		lor_y = 1./pi * (sig/(lor_x**2 + sig**2))
		#lor_y = sqrt(1/(2*pi*sig**2)) * exp(-lor_x**2/(2*sig**2))
		#for n in range(1,d.shape[1]):
		#	d[:,n] = scipy.signal.convolve(d[:,n], lor_y, mode='same') *dx
		y = scipy.signal.convolve(y, lor_y, mode='same') * dx
		y = y * sum_rule/sum(y[int_window:-int_window])

	#simple peak
	data.append(d)
	offset=sum(x<=opts.xmin)
	idx=argmax(y[(x<opts.xmax)&(x>opts.xmin)])+offset
	#y *= 36.0/sum_rule

	print 'File:',fname
	print '  Peak: %.5f, %.5f'%(x[idx],y[idx])
	DX = 2.0
	cond = (x>x[idx]-DX)&(x<x[idx]+DX)
	x_cond = x[cond]
	y_cond = y[cond]
	x_mean = sum(x_cond*y_cond)/sum(y_cond)
	print '  Mean: %.5f'%(x_mean)

	#now, fitting polynomial
	from scipy.interpolate import UnivariateSpline as USpl
	DX = 0.15
	cond = (x>x[idx]-DX)&(x<x[idx]+DX)
	x_cond = x[cond]
	y_cond = y[cond]
	p_=polyfit(x_cond, log(y_cond), 4)
	
	#y0=y[idx]*0.85
	#DX = 2.0
	#cond = (x>x[idx]-DX)&(x<x[idx]+DX)
	#cond = cond&(y>y0)
	#x_cond = x[cond]
	#y_cond = y[cond]
	#print x_cond

	#spl = USpl(x_cond, (y_cond), s=1e3)
	#print p_
	#print -p_[1]/(2*p_[0])
	#y2 = (polyval(p_, x_cond))
	#y2 = (spl(x_cond))
	#p.plot(x_cond,y2)
	

	data_max = max(data_max, max(y[(x>opts.xmin)&(x<opts.xmax)]))

	if opts.asymmetry:
		#calculate the ration of the right and left slopes
		if not delta_E_l: 
			delta_E_l = 0.08
			delta_E_r = 0.15
			window_E_l = 0.45
			window_E_r = 0.25

		if (1):
			#find largest slope
			window = int(1/dx)
			step=5
			#left
			rng = arange(idx-window,idx,step)
			slope = y[rng] - y[rng-1]
			x_l = x[idx-window + argmax(slope)*step]
			y_l = y[idx-window + argmax(slope)*step]
			#print x_l, y_l
			a_l = amax(slope)/(dx)
			#right
			rng = arange(idx,idx+window,step)
			slope = y[rng+1] - y[rng]
			x_r = x[idx + argmin(slope)*step]
			y_r = y[idx + argmin(slope)*step]
			#print x_r,y_r
			a_r = amin(slope)/(dx)
		else:	
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

	#plot
	if opts.plot:
		if cnt>=len(labels):
			import re
			label=re.sub(r'([^0-9a-zA-Z\-\.\/])', r'\\\1', fname)
			#label=fname.replace(r'_',r'\_')
		else:
			label=labels[cnt]

		#plot_data = p.plot(x,y, '-',label=label, lw=2.5)
		cutl = 30
		#cutl = 0
		print x.shape
		print y.shape
		plot_data = p.plot(x[cutl:],y[cutl:], '-',label=label, lw=2.0)
		if cnt<len(colors):
			plot_data[0].set_color(colors[cnt])
		if cnt<len(linestyles):
			plot_data[0].set_linestyle(linestyles[cnt])
		if cnt<len(linewidths):
			plot_data[0].set_linewidth(linewidths[cnt])
		if opts.asymmetry:
			p.plot(x[idx],y[idx], 'o' ,ms=8, mew=0, alpha=0.5,
				color=plot_data[0].get_color())
			#p.axvline(x[idx], ls=':', lw=1.5, alpha=0.75)
			delta=0.2
			p.plot([x_l-delta, x_l+delta], [y_l-delta*a_l, y_l+delta*a_l], '--', color='black', lw=3 )
			p.plot([x_r-delta, x_r+delta], [y_r-delta*a_r, y_r+delta*a_r], '--', color='black', lw=3 )
			#p.plot([x_l[0],x_l[-1]],[y_reg_l[0],y_reg_l[-1]],color='black',lw=3,ls='--')
			#p.plot([x_r[0],x_r[-1]],[y_reg_r[0],y_reg_r[-1]],color='black',lw=3,ls='--')

	if opts.export:
		print 'Exporting'
		f_new = 'smoothed-'+fname
		savetxt(f_new, column_stack((x,y)) )

	cnt+=1
	print

if opts.plot:
	#p.legend(loc='upper left', prop={'size':20})
	p.legend(loc='upper left', prop={'size':16})

	if opts.title:
		p.title(opts.title, size=20)
	else:
		#p.title(r'$\mathrm{Absorbance}/\varepsilon$')
		#p.title(r'$\mathrm{Absorbance\;(7\%\,hole\,doping)}$',size=20)
		p.title(r'Absorbance',size=20)

	p.xlabel(r'E (eV)',size=18)
	if plot_idx!=-1:
		p.ylabel(ylabels[plot_idx], size=18)
	else:
		p.ylabel(ylabels[0], size=18)

	p.xlim(opts.xmin, opts.xmax)
	if not (opts.ymax is None):
		p.ylim(opts.ymin, opts.ymax)
	else:
		p.ylim(opts.ymin, data_max*1.2)


	if opts.output:
		#p.savefig(opts.output, transparent=True)
		p.savefig(opts.output, dpi=150)
	else:
		p.show()
