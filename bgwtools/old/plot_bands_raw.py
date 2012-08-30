#!/usr/bin/env python
#coding:utf-8

from numpy import *
import matplotlib.pyplot as p
import sys

if len(sys.argv)<2:
	sys.exit()

data_raw = loadtxt(sys.argv[1])
data = take(data_raw, [0,1,4], axis=1)

def sel_band(idx):
	return data[data_raw[:,3]==idx,:]

kpts = take(sel_band(1), [0,1], axis=1)
indices = arange(len(kpts))

locs=[]
labels=[]
colors=('black', 'gray', 'blue', 'green', 'red', 'yellow', 'orange', 'purple')
#colors = ['black']*18
def plot_band(band_num, kpts_idx, color=None,
	 start_mark=None, end_mark=None, plot_2d=True):
	global last_kpt, locs, labels
	band = sel_band(band_num)[kpts_idx,:]

	if plot_2d:
		#left panel: kpts in a 2D surface
		p.subplot(121)
		p.plot(band[:,0],band[:,1], 'o', color=colors[band_num-1])

	#right panel: Energy vs. kpts
	p.subplot(122)
	idx = last_kpt + arange(len(band[:,0]))
	last_kpt = idx[-1]
	if start_mark:
		locs.append(idx[0])
		labels.append(start_mark)
		p.axvline(idx[0], ls='--', color='gray')
	if end_mark:
		locs.append(last_kpt)
		labels.append(end_mark)
		p.axvline(idx[-1], ls='--', color='gray')
	if (idx[0]==0):	leg = 'band %d'%(band_num)
	else: leg=None
	p.plot(idx, band[:,2] , '-o', lw=1.5, color=colors[band_num-1], \
		label=leg)

idx_GK = indices[kpts[:,0]==0] #Γ to K
idx_KM = indices[append(kpts[:-1,1]>kpts[1:,1], [True])] #K to M
idx_MG = indices[append([True], kpts[:-1,0]<kpts[1:,0])[::1]] [::-1] #M to Γ (+ reserve order)

cond=True
for i in [1,2,3,4,5,6,7,8]:
#for i in [1,2]:
#for i in [2]:
	last_kpt=0
	plot_band(i, idx_GK, plot_2d=cond, start_mark='$\Gamma$', end_mark='M')
	plot_band(i, idx_KM, plot_2d=cond)
	plot_band(i, idx_MG, start_mark='K', end_mark=r'$\Gamma$', plot_2d=cond)
	cond=False

p.subplot(122)
FE=-1.527932
delta=2
p.axhline(FE)
p.axhline(FE-delta)
p.axhline(FE+delta)
p.xticks(array(locs),labels)
p.legend(prop={'size':'small'})
p.show()
