#!/usr/bin/env python

# Gets the output from eval_components, which has the most important contributions
# for a particular exciton mode, and plot the S(w) graph

from numpy import *
import matplotlib.pyplot as plt
import sys
from BSE import *
			

evals = loadtxt(sys.argv[1]) #eigenvalues.dat
fevecs = eigenvectorsIO(sys.argv[2]) #eigenvectors

noeh = loadtxt(sys.argv[3]) #eig_noeh/eigenvalues_noeh.dat
freqs = noeh[:, 0]
fvmtxel = vmtxelIO(sys.argv[4]) #vmtxel
vmtxel = fvmtxel.read_mtxel().ravel()
	
order = argsort(freqs)
freqs = freqs[order]
print vmtxel[:5]
vmtxel = vmtxel[order]
#vmtxel *= freqs #multiply back the denominator...

def weighted_mean(X,W,do_print=True):
	mean = sum(W*X)/sum(W)
	if do_print:
		print '  mean:', mean
	std  = sqrt(sum(W*(X-mean)**2)/sum(W))
	if do_print:
		print '   std:', std
	return std

#only print index, eigenvalue and std
def print_std(idx):
	fevecs.goto_evec(idx, 0)
	eval,evec = fevecs.get_evec()
	evec = evec[order]
	data = evec * vmtxel

	# This is the only quantity that really makes sense!
	std = weighted_mean(freqs, abs(data), False)
	print "%d\t%.5f\t%.2f\t%.5f"%(idx, eval, evals[idx,1], std)
	return std

def print_evec(idx, do_plot=False):
	fevecs.goto_evec(idx, 0)
	eval,evec = fevecs.get_evec()

	print 'Mode #',idx
	print '  eval:', eval
	print '    cs:', evals[idx,1]

	evec = evec[order]
	data = evec * vmtxel

	# This is the only quantity that really makes sense!
	print '  DATA: |data|'
	weighted_mean(freqs, abs(data))
	#print '  DATA: |data|^2'
	#weighted_mean(freqs, data**2)
	#print '  DATA: |evec|^2'
	#weighted_mean(freqs, evec**2)

	N=250
	bins = linspace(0, 10, N)
	inds = digitize(freqs, bins)
	def mean_(x):
		if len(x):
			return sum(x)
		else:
			return 0

	means = array([mean_(data[inds==i]) for i in xrange(1, len(bins))])

	dx = bins[1]-bins[0]
	x = (bins[1:] + bins[:-1])*0.5
	cum_hp = cumsum(data)
	#cum = cumsum(means)

	means *= sign(cum_hp[-1])
	#cum   *= sign(cum_hp[-1])
	cum_hp*= sign(cum_hp[-1])

	print '  calculated cs:',cum_hp[-1]**2

	if do_plot:
		plt.bar(x-dx, means, width=dx, color='r', edgecolor='r')
		#plt.plot(x, cum, 'b-', lw=2)
		plt.plot( insert(freqs,0,0), insert(cum_hp,0,0), 'b-', lw=2)
		plt.axvline(eval, color='gray')

		plt.xlim(0,10)
		plt.show()

# Select brightest mode
sub = (evals[:,0]>4) & (evals[:,0]<5.0)
offset = arange(len(evals[:,0]))[sub][0]
f = -evals[:,1][sub]
order2 = argsort(f); del f
if 0: #some nice plots
	idx = order2[0] + offset
	print_evec(idx, True)
	idx = order2[1] + offset
	print_evec(idx, True)
	idx = order2[3000] + offset
	print_evec(idx, True)
if 1: #only show statistics
	print '#idx\teval\tcs*W\tstd'
	std = []
	for i in range(100):
		idx = order2[i] + offset
		std += [print_std(idx)]
	std = array(std)
	print '#mean(std):\t%.5f'%(std.mean())
	print '#std(std):\t%.5f'%(std.std())
	plt.hist(std)
	plt.show()

