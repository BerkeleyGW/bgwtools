#!/usr/bin/env python

# Felipe Homrich da Jornada <jornada@berkeley.edu>
# Script to plot the correction (GW-LDA) to the DFT band structure. It plots
# both the real and imaginary part of Sigma.

#THIS VERSION SHOULD BE USED ONLY FOR GRAPHENE, AND THE M POINT BUST BE PRESENT

from __future__ import division

from numpy import *
import matplotlib.pyplot as plt
import sys,os
from scipy.interpolate import *
from optparse import OptionParser, OptionGroup

############
#  Parser  #
############
usage = "usage: %prog [options] file_1 [file_2 ... file_N]"
parser = OptionParser(usage)

group_data = OptionGroup(parser, 'Data Related')
group_data.add_option('--FE0', default=0.0, type='float',
	help='DFT Fermi Energy')
group_data.add_option('--shift-FE', dest='shift_FE', default=0.0, type='float',
	help='Shift to be applied to the DFT Fermi Energy')
group_data.add_option('--show-imag', dest='show_imag', type='choice', choices=('never', 'auto', 'always'), default='auto',
	help='Determine whether imaginary part of Sigma should be plotted.')
parser.add_option_group(group_data)

group_visual = OptionGroup(parser, 'Aesthetics')
group_visual.add_option('-t','--title',type='str',
	help='title')
group_visual.add_option('-l','--labels',type='str',
	help='comma-separated list of labels for the plots')
group_visual.add_option('--no-plot', dest='no_plot', action='store_true', default=False,
	help='don\'t plot anything, only output useful info')
parser.add_option_group(group_visual)

parser.add_option('--load', type='str',metavar='FILE',
	help='load external python file with configurations')

###################
#  Configuration  #
###################
(opts,args) = parser.parse_args()

#auto load GW_cor_conf.py
if os.path.exists('GW_cor_conf.py'):
	execfile('GW_cor_conf.py')
	
if opts.load:
	execfile(opts.load)

if len(args)<1:
	parser.error('incorrect number of arguments')

labels=[]
if opts.labels:
	if len(opts.labels):
		labels=opts.labels.split(',')
#labels=['$\mathrm{GPP:}\; G_0W_0$', '$\mathrm{GPP:}\; GW_0$', '$\mathrm{FF:}\; GW_0$']
#labels=['$\mathrm{GPP:}\; GW_0$', '$\mathrm{FF:}\; GW_0$']

shift_FE = opts.shift_FE
FE0 = opts.FE0

import bgwtools.converters.read_eqp as read_eqp

#output some scissors info onto params.py
f_params=open('params_out.py', 'w')
f_params.write('params.FE0=%f\n'%(FE0))
f_params.write('params.shift_FE=%f\n'%(shift_FE))

def params_start(fname):
	f_params.write('\nlast_file=%s\n'%(repr(fname)))
	f_params.write('params.add_file(last_file)\n')

def params_write(prop, value):
	f_params.write('params.files[last_file].%s=%s\n'%(prop,repr(value)))


#wrapper for splrep with custom values for smoothing and spline order
def get_tck(X, Y):
	global k
	k=1
	tck = splrep(X,Y,k=k,s=0)
	return tck

#plot region whose transitions are blocked by Pauli principle
def plot_blocking(color):
	plt.axhline(FE0_cor, color='green', ls=':')
	plt.axvline(FE0, color='green', ls=':')
	if not (shift_FE is None):
		#plt.axvline(FE0+shift_FE, color='red', ls=':')
		#plt.axhline(shift_FE2-shift_FE, color='red', ls=':')
		plt.axvspan(FE0-shift_FE, FE0+shift_FE, color=color, alpha=0.2)
		if plot_imag and has_imag:
			#blocking = (DFT blocking) + (correction)
			blocking = 2*shift_FE + (splev(FE0 + shift_FE, tck) - splev(FE0 - shift_FE, tck))
			print '\tBlocking: %.3f'%(blocking)
			plt.subplot(212)
			plt.axvspan(0, blocking, color=color, alpha=0.2)
			plt.subplot(211)

def on_click(event):
	print
	thisline = event.artist
	#print dir(thisline)
	xdata, ydata = thisline.get_data()
	for ind in event.ind:
		#print ind
		idx = thisline.idx[ind]
		print 'Point #',idx+1
		print '  E_DFT',xdata[ind]
		print '  kpt=',thisline.kpts[ind]
		

#plt.connect('button_press_event', on_click)
plt.connect('pick_event', on_click)

cnt=0
for fname in args:
	params_start(fname)
	if cnt<len(labels):
		label=labels[cnt]
	else:
		label=fname
	print 'File %s\n'%(fname)

	#should we plot the imaginary part of this GW band structure?
	if opts.show_imag == 'always':
		plot_imag = True
	else:
		plot_imag = False

	#read eqp file, organize stuff
	bands, x,y,z, e_dft,e_gw = read_eqp.get_data_from_eqp(fname)
	idx_v = arange(len(bands))[bands==4][0]
	idx_c = arange(len(bands))[bands==5][0]
	E_DFT = e_dft[:,idx_v:idx_c+1]
	E_GW = e_gw[:,idx_v:idx_c+1]
	dft = ravel(E_DFT)
	gw  = ravel(E_GW)
	#store the index of each kpoint
	idx = arange(len(x))
	idx = column_stack((idx,idx)).ravel()
	order = argsort(dft)
	X=dft[order]          #this is what will be plotted in the
	Y=real(gw-dft)[order] # top panel
	idx=idx[order]

	#find M point
	idx_M = arange(len(y))[abs(y-0.5)<1e-3][0]
	xl = E_DFT[idx_M,0] #xl = DFT energy of val band @ M
	xr = E_DFT[idx_M,1] #xr = DFT energy of cond band @ M

	#print info about M point
	def get_vals(A):
		return tuple([A[0], A[1], A[1]-A[0]])
	print '\tM point:       E_v     E_c     GAP'
	print '\t  DFT     '+(3*' %7.4f')%get_vals(E_DFT[idx_M,:])
	print '\t  GW      '+(3*' %7.4f')%get_vals(real(E_GW [idx_M,:]))
	print

	#if requested, autodetect if we have a complex eqp file
	has_imag = isinstance(gw[0], complex)
	if opts.show_imag == 'auto':
		plot_imag = has_imag

	if plot_imag:
		plt.subplots_adjust(hspace=0.3)
		plt.subplot(211)

	#interpolate GW corrections to band structure
	tck = get_tck(X, Y)
	params_write('tck', tck)

	#print info about FE
	FE0_cor = splev(FE0, tck) #correction for Fermi energy (:= FE0_GW - FE0)
	FE0_GW = FE0 + FE0_cor
	print '\tFermi Energy:       DFT      GW     cor'
	print '\t  Undoped FE   '+(3*' %7.4f')%(FE0, FE0_GW, FE0_cor)
	if not (shift_FE is None):
		#final FE
		FE1 = FE0 + shift_FE
		FE1_cor = splev(FE1, tck)
		FE1_GW = FE1 + FE1_cor
		print '\t  Doped FE     '+(3*' %7.4f')%(FE1, FE1_GW, FE1_cor)
		#this is what BerkeleyGW wants!!
		efermi = FE1_GW - FE0_GW
		print
		print '\t> Corrected undoped FE (fe0_cor): %7.4f'%(FE0_GW)
		print '\t  DFT shift in FE:                %7.4f'%(shift_FE)
		print '\t> QP shift in FE (efermi):        %7.4f'%(efermi)

		params_write('FE1_GW', FE1_GW)
		params_write('efermi', efermi)
	print

	#right: points higher in energy than the cond band @ M
	cond_r=(X > xr + 1e-3)
	ar,br = polyfit(X[cond_r],Y[cond_r],deg=1)
	x0r=-br/ar
	pr = poly1d((ar,br))
	#left: points lower in energy than the val band @ M
	cond_l=(X < xl - 1e-3)
	al,bl = polyfit(X[cond_l],Y[cond_l],deg=1)
	x0l=-bl/al
	pl = poly1d((al,bl))

	print '\tSplines Data (n,t,c,k):'
	print len(tck[0]) #n -- number of knots
	for i in range(len(tck[0])):
		print tck[0][i], #t -- pos of knots
	print
	for i in range(len(tck[1])):
		print tck[1][i], #c -- bsplines coeffs
	print
	print k #k -- degree

	params_write('cvfit', array([0.0,x0l,al,0.0,x0r,ar]) )
	print
	print ('\tcvfit      '+ (' %10.6f')*6)%(0.0, x0l, al, 0.0, x0r, ar)
	print ('\tcvfit_outer'+ (' %10.6f')*6)%(0.0, x0l, al, 0.0, x0r, ar)
	print 

	#plot (linear) scissors fit
	pts=plt.plot([X[0],xl], [pl(X[0]), pl(xl)], ls='-', lw=5, alpha=0.5)
	color=pts[0].get_color()
	plt.plot([xr,X[-1]], [pr(xr), pr(X[-1])], ls='-', lw=5, alpha=0.5, color=color)

	lines=plt.plot(X, Y, 'o', color=color, label=label, ms=5, picker=5)
	xx = column_stack((x,x)); yy = column_stack((y,y)); zz = column_stack((z,z));
	xx = ravel(xx); yy=ravel(yy); zz=ravel(zz);

	kpts = column_stack((xx, yy, zz))[order]
	#kpts = column_stack((x.tolist()*2, y.tolist()*2, z.tolist()*2))[order]
	lines[0].kpts=kpts
	lines[0].idx=idx
	x_intp=linspace(X[0],X[-1],1000)

	plt.plot(x_intp, splev(x_intp, tck), color=color)

	#Plot the imaginary part
	if plot_imag and has_imag:
		plt.subplot(212)
		#only plot if cond states are above original FE + shift
		unblocked = real(E_DFT[:,1]) >= (FE1)
		E_trans_DFT = real(E_DFT[:,1] - E_DFT[:,0])
		E_trans = real(E_GW[:,1] - E_GW[:,0])
		EM = E_trans[idx_M]
		Decay = -imag(E_GW[:,1] + E_GW[:,0])
		DecayM = Decay[idx_M]
		print '\tTotal Im(Sigma) @ M:', DecayM

		#FIXME
		print repr(real(dft)[order])
		print repr(imag(gw)[order])


		order_DFT = argsort(array(E_trans_DFT))
		E_trans_DFT = E_trans_DFT[order_DFT]
		Decay_DFT = Decay[order_DFT]

		order_GW = argsort(array(E_trans))
		E_trans = E_trans[order_GW]
		Decay = Decay[order_GW]

		#we output Decay vs. E_trans_DFT bc its user to use later
		print E_trans_DFT
		print Decay_DFT
		tck_en = splrep(E_trans_DFT, Decay_DFT, k=1, s=0)
		print '\tIm Splines (total Im vs. LDA excitation energy):'
		print repr(tck_en)

		#blocked transitions
		cond_b = (~unblocked)[order_GW]
		if cond_b.any():
			E_trans_b = E_trans[cond_b]
			Decay_b = Decay[cond_b]
			plt.plot(E_trans_b, Decay_b, 'o', color=color, ms=2, picker=5)
			plt.plot(E_trans, Decay, ':', color=color, lw=0.5)

		#unblocked transitions
		cond_u = unblocked[order_GW]
		if cond_u.any():
			E_trans_u = E_trans[cond_u]
			Decay_u = Decay[cond_u]
			plt.plot(E_trans_u, Decay_u, 'o-', color=color, ms=6, lw=1.25, picker=5)

		plt.xlabel('GW EH Excitation Energy', size=16)
		#plt.ylabel('$-\mathrm{Im}(\Sigma)$', size=16)
		plt.ylabel('$-\mathrm{Im}(\Sigma_v + \Sigma_c)$', size=16)
		#plt.axvline(2*shift_FE2, color=color, ls='--')
		plt.axvline(EM, color=color, ls='--')
		lims = plt.xlim()
		plt.xlim(0, lims[-1])
		lims = plt.ylim()
		plt.ylim(0, lims[-1])
		plt.subplot(211)

	plot_blocking(color)

	cnt+=1

f_params.close()

def on_pick(event):
	print '!'
	thisline = event.artist
	xdata, ydata = thisline.get_data()
	ind = event.ind
	print 'on pick line:', zip(xdata[ind], ydata[ind])

#print M points. Should only work for graphene
plt.axvline(xl, color='black', ls='--')
plt.axvline(xr, color='black', ls='--')

#labels
plt.xlabel('$E_{DFT}$', size=14)
plt.ylabel('$E_{GW}\;-\;E_{DFT}$', size=18)
plt.title('GW Correction to DFT Band Structure', size=20)
plt.legend(loc='upper left', prop={'size':12})

if opts.no_plot:
	sys.exit()

#F = plt.gcf()
#DPI = F.get_dpi()
#DefaultSize = F.get_size_inches()
#F.set_size_inches( (DefaultSize[0], DefaultSize[1]*0.75) )
plt.show()
#plt.savefig('scissors.svg')
