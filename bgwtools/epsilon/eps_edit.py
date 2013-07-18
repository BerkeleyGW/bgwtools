#!/usr/bin/env python

# A tool to edit an epsmat file

# Felipe Homrich da Jornada <jornada@berkeley.edu> (2013)

from bgwtools.IO.epsmat import epsmatIO
import numpy as np

# This is the function that modifies the epsilon matrix.
# Currently, it's setting the head to 1.0 for all q-pts.
def eps_edit(eps):
    '''Do something interesting with the epsilon matrix'''

    nq = eps.nq
    # We will set the head of epsilon to these values (1's).
    heads = np.ones(nq)

    #self.epsmat.epsmat
    for iq in xrange(nq):
            # "nmtx" = number of G vectors = rows = cols of epsinv
            nmtx = eps.nmtx[iq]
            # "qpt" is the vector for the current q-point
            qpt = eps.qpt[:,iq]
            # "mat" is the epsinv matrix
            mat = eps.epsmat[iq]
            # "gvec_q" is G + qpt, where "G" is the Gvector in the order 
            # they appear in the row/columns of epsinv matrix ("mat")
            gidx = eps.isort[iq][:nmtx] - 1
            gvec_q = eps.gvec_k[:,gidx] + qpt[:,np.newaxis]

            # set the head to some value..
            mat[0,0] = heads[iq]

            eps.epsmat[iq] = mat

if __name__=='__main__':
	import sys

	if len(sys.argv)!=3:
		print 'Usage: %s epsmat_in epsmat_out'%(sys.argv[0])
		exit(1)

	f_in = sys.argv[1]
	f_out = sys.argv[2]

	eps = epsmatIO(f_in)
        eps_edit(eps)
	eps.to_file(f_out)
