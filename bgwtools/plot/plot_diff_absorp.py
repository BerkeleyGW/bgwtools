#!/usr/bin/env python

from __future__ import division

delta=1e-7

from numpy import *
import pylab as p

'''
fig_width_pt = 246.0*1.5	# Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
'''
params = {'backend': 'eps',
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'text.usetex': True,
          'font.family': 'serif',
          'font.serif': 'Times'} #,
#          'figure.figsize': fig_size}
p.rcParams.update(params)

import sys

eps1 = genfromtxt(sys.argv[1])
eps2 = genfromtxt(sys.argv[2])

p.plot(eps1[:,0],eps2[:,1]-eps1[:,1])
p.plot(eps1[:,0],(eps2[:,1]-eps1[:,1])/(eps1[:,1]+delta))

#p.legend(loc='lower right')i
p.title('Absorptance - Difference')
p.xlabel(r'$E(eV)$')
p.ylabel(r'$\varepsilon$', rotation='horizontal')
p.show()
