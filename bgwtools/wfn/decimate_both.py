#!/usr/bin/env python

# This script reduces the number of k-points in a WFN file.
# It only keeps the N smallest v->c transitions.
# User should provide the Fermi energy, N, and the q-shift

# Felipe Homrich da Jornada (Feb 2013)

from numpy import *
from bgwtools.IO.wfn import wfnIO
import sys
from bgwtools.common.common import get_numpy_flavor

ryd = 13.60569253

if len(sys.argv)!=8:
    print('Usage: %s Efermi N qshift wfn_in wfnq_in wfn_ou wfnq_out'%(sys.argv[0]))
    sys.exit(1)

Efermi = float(sys.argv[1])/ryd + 1e-9
N = int(sys.argv[2])
qshift = array(map(float, sys.argv[3].split(',')))
fname_in = sys.argv[4]
fnameq_in = sys.argv[5]
fname_out = sys.argv[6]
fnameq_out = sys.argv[7]

print 'WFN:', fname_in
wfn_in = wfnIO(fname_in)
print 'WFNq:', fnameq_in
wfnq_in = wfnIO(fnameq_in)
print 'q-shift', qshift

#c -> k
#v -> k + q
tol=1e-6
#kpts_qmq0 = (wfnq_in.kpt - qshift[:,newaxis] + 1) % 1
kpts_qmq0 = wfnq_in.kpt - qshift[:,newaxis]
iks = []
ikqs = []
evs = []
ecs = []
for ik in range(wfn_in.nk):
    delta = (kpts_qmq0 - wfn_in.kpt[:,ik][:,newaxis] + tol) % 1
    #delta = (kpts_qmq0 - wfn_in.kpt[:,ik][:,newaxis])
    #kk_q = wfn_in.kpt[:,ik] + qshift
    #delta = kpts_q - kk_q[:,newaxis]
    #put in [-0.5,0.5) range
    #delta = delta - floor(delta + 0.5)
    #delta = delta - floor(delta + 0.5)
    ikq = argmin( sum(fabs(delta), axis=0 ) )
    #print ikq
    if all( fabs(wfnq_in.kpt[:,ikq] - wfn_in.kpt[:,ik] - qshift) < tol):
        iks.append(ik)
        ikqs.append(ikq)
        ec = wfn_in.energies[:,ik,0]
        ecs.append(ec[ec>Efermi][0])
        ev = wfnq_in.energies[:,ikq,0]
        evs.append(ev[ev<Efermi][-1])

print len(iks), wfn_in.nk
iks = array(iks)
ikqs = array(ikqs)
ecs = array(ecs)
evs = array(evs)
dEs = (ecs - evs)*ryd

ind = argsort(dEs)[:N]

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

truncate(wfn_in, fname_out, iks[ind])
truncate(wfnq_in, fnameq_out, ikqs[ind])

print 'All done!'
