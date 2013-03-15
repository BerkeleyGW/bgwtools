#!/usr/bin/env python

# Given a decimated grid, this script creates a pair of WFN and WFNq files so
# that all q transitions are present
# User should provide the Fermi energy, N, and the q-shift

# Felipe Homrich da Jornada (Feb 2013)

from numpy import *
from bgwtools.IO.wfn import wfnIO
import sys
from bgwtools.common.common import get_numpy_flavor

ryd = 13.60569253

if len(sys.argv)!=5:
    print('Usage: %s qshift wfn_in wfn_out wfnq_out'%(sys.argv[0]))
    sys.exit(1)

qshift = array(map(float, sys.argv[1].split(',')))
fname_in = sys.argv[2]
fname_out = sys.argv[3]
fnameq_out = sys.argv[4]

print 'WFN:', fname_in
wfn_in = wfnIO(fname_in)
print 'q-shift', qshift

#c -> k
#v -> k + q
tol=1e-8
#kpts_qmq0 = (wfnq_in.kpt - qshift[:,newaxis] + 1) % 1
# calculate k + q
#kpts_pq = wfn_in.kpt + qshift[:,newaxis]
iks = []
ikqs = []
for ik in range(wfn_in.nk):
    ks_pq = wfn_in.kpt[:,ik] + qshift
    delta = (wfn_in.kpt - ks_pq[:,newaxis] + tol) % 1
    #delta = (kpts_qmq0 - wfn_in.kpt[:,ik][:,newaxis])
    #kk_q = wfn_in.kpt[:,ik] + qshift
    #delta = kpts_q - kk_q[:,newaxis]
    #put in [-0.5,0.5) range
    #delta = delta - floor(delta + 0.5)
    #delta = delta - floor(delta + 0.5)
    ikq = argmin( sum(fabs(delta), axis=0 ) )
    #print ikq
    if all( fabs(wfn_in.kpt[:,ikq] - wfn_in.kpt[:,ik] - qshift) < tol):
        iks.append(ik)
        ikqs.append(ikq)

print len(iks), len(ikqs), wfn_in.nk
iks = array(iks)
ikqs = array(ikqs)

if raw_input('Writing data to %s and %s. Are you sure? [y/N] '%(fname_out, fnameq_out))!='y':
    sys.exit(0)

def truncate(wfn_in, fname_out, ik_keep):
    print
    print 'Truncating',wfn_in.fname,'->',fname_out
    wfn_out = wfnIO()
    wfn_out.__dict__ = wfn_in.__dict__.copy()
    wfn_out.nk = len(ik_keep)
    wfn_out.ngk = wfn_out.ngk[ik_keep]
    wfn_out.ngkmax = amax(wfn_out.ngk)
    wfn_out.kw = wfn_out.kw[ik_keep]
    wfn_out.kw[:] = 1.0/wfn_out.nk
    wfn_out.kpt = wfn_out.kpt[:, ik_keep]
    #HACK!
    #wfn_out.ifmin[:] = 1
    #wfn_out.ifmax[:] = 4
    wfn_out.ifmin = wfn_out.ifmin[ik_keep, :]
    wfn_out.ifmax = wfn_out.ifmax[ik_keep, :]
    wfn_out.energies = wfn_out.energies[:, ik_keep, :]
    wfn_out.occupations = wfn_out.occupations[:, ik_keep, :]
    wfn_out.f = None
    wfn_out.fname = fname_out
    wfn_out.write_header(full=True)

    ng_max = amax(wfn_in.ngk)
    gvec = empty((3, ng_max), dtype=int, order='F')
    data = empty((ng_max, wfn_in.ns), dtype=get_numpy_flavor(wfn_in.flavor), order='F')
    k_done = 0
    for ik in range(wfn_in.nk):
        if ik in ik_keep:
            k_done += 1
            print k_done,'/',wfn_out.nk
            wfn_in.read_gvectors(gvec)
            wfn_out.write_gvectors(gvec[:,:wfn_in.ngk[ik]])
            for ib in range(wfn_in.nbands):
                wfn_in.read_data(data)
                wfn_out.write_data(data[:wfn_in.ngk[ik],:])
            if k_done==wfn_out.nk:
                break
        else:
            wfn_in.read_gvectors()
            for ib in range(wfn_in.nbands):
                wfn_in.read_data()
    print wfn_out

truncate(wfn_in, fname_out, iks)
wfn_in = wfnIO(fname_in)
truncate(wfn_in, fnameq_out, ikqs)

print 'All done!'
