#!/usr/bin/env python

# Plots the classical S(w) graph. Requires four parameters:
# - the eigenvalues file
# - the vmtxel file
# - the eigenvalues_noeh.dat file
# - the excitonic state to be plotted.

# NOT WORKING!

# If true, we'll plot S(w) and the integrated quantity.
# Else, we will plot... WHAT?
plot_S = False

gui=False

import numpy as np

def get_S(evecs, iS, Eqp, vmtxel):
    """Given an excitonic state Avck and the QP energies Eqp, returns the """



import matplotlib.pyplot as plt
import sys

evals = loadtxt(sys.argv[1])
f_comp = open(sys.argv[2])
noeh = loadtxt(sys.argv[3])
vmtxel = loadtxt(sys.argv[4])
vmtxel = vmtxel.ravel()

# Which state to plot? 0=brightest state, 1=second brightest, etc.
if len(sys.argv)>5:
	S_idx = int(sys.argv[5])
else:
	S_idx = 0

f_comp.readline()
sub_idx = []
sub_comp = []
n_elems = []
head=False
for line in f_comp.readlines():
	head = not head
	if head:
		try:
			#the index of each eigenvector
			sub_idx.append(int(line)-1)
		except:
			break
	else:
		#the components that make up that eigenvector
		# idx-> the index (wrt to vertical transition) of a part. component
		# A -> the coefficient
		comp = fromiter(map(float,line.split()), dtype=float)
		comp = comp.reshape( (len(comp)/2, 2) )
		comp = rec.fromarrays( [array(comp[:,0],dtype=int)-1,comp[:,1]], names='idx,A' )
		sub_comp.append(comp)
		n_elems.append(comp.shape[0])

sub_idx = array(sub_idx)
n_elems = array(n_elems)

#print sub_idx
#print sub_comp
#print n_elems

## Select the most delocalized mode
#idx_ = argmax(n_elems)
#idx_ = 500

# Select brightest mode
f = evals[sub_idx,1]
indices = argsort(-f)
#local_idx = argmax(f) #wrt to subset of evals that we have output
local_idx = indices[S_idx]
eig_idx = sub_idx[local_idx] #real eval/evec index
print evals[eig_idx, 1]

evec = sub_comp[local_idx]
freqs = noeh[evec.idx, 0]

print 'Selected mode:', evals[eig_idx]
print 'Importing %d components'%(len(freqs))

order = argsort(freqs)
freqs = freqs[order]

if plot_S:
	vmtx = vmtxel[ evec.idx[order] ]
	data = evec.A[order] * vmtx
else:
	data = evec.A[order]**2

N=50
df=1e-6
bins = linspace(freqs[0]-df, freqs[-1]+df, N)
inds = digitize(freqs, bins)
def mean_(x):
	if len(x):
		return mean(x)
	else:
		return 0
if plot_S:
	means = array([mean_(data[inds==i]) for i in xrange(1, len(bins))])
else:
	means = array([sum(data[inds==i]) for i in xrange(1, len(bins))])

dx = bins[1]-bins[0]
x = (bins[1:] + bins[:-1])*0.5
cum = cumsum(means)

means *= sign(cum[-1])
cum   *= sign(cum[-1])

plt.bar(x-dx, means, width=dx, color='r')

plt.plot(x, cum, 'b-', lw=2)
#plt.bar(x, means, 0.5)
#plt.plot(freqs, sub_comp[idx_max].A, 'o' )
plt.axvline(evals[eig_idx,0], color='gray')
#order = argsort(freqs)
#print freqs[order]
#plt.hist(freqs[order], sub_comp[idx_max].A[order])

if not plot_S:
	#normalize weights
	data /= sum(data)
	#calculate weighted averaged
	_avg = sum(data * freqs)
	_std = sqrt( sum(data * (freqs - _avg)**2) )
	print 'Average:', _avg
	print 'Std:', _std
	f=open('states.dat','a')
	f.write('%d %d %f %f %f\n'%(S_idx, eig_idx, evals[eig_idx][0], _avg, _std))
	f.close()

if gui:
	plt.show()

