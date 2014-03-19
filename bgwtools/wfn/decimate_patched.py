#!/usr/bin/env python

# This script reduces the number of k-points in a WFN file to a small patch
# centered on a given point. The radius of the patch, qmax, is given in cartesian coordinates,
# while the center of the patch is given in crystal coordinates.
# Set print_q to True to print out all the q-points needed for Epsilon for the decimated
# wave function.
#
# Adapted from FHJ's decimate.py

from numpy import *
from bgwtools.IO.wfn import wfnIO
import sys
from bgwtools.common.common import get_numpy_flavor
from bgwtools.IO.fullbz import fullbzIO

ryd = 13.60569253
print_q=True
plot_patch=False
keep_Kprime=True
TOL=1e-6

if len(sys.argv)!=7:
    print('Usage: %s qmax kx ky kz wfn_in wfn_out'%(sys.argv[0]))
    sys.exit(1)

qmax=float(sys.argv[1])
k_center=array([float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4])])
fname_in = sys.argv[5]
fname_out = sys.argv[6]

wfn_in = wfnIO(fname_in, full=True)
# Convert to cartesian coordinates
rk_cart=tensordot(wfn_in.bvec,wfn_in.kpt, axes=([1],[0]))
kc_cart=dot(wfn_in.bvec,k_center)

fullbz=fullbzIO('fullbz.dat')
print shape(fullbz.fk)
fk_cart=tensordot(wfn_in.bvec,fullbz.fk,axes=([1],[0]))
print shape(fk_cart)

if plot_patch:
    from matplotlib import pyplot as plt
    plt.scatter(fk_cart[0,:],fk_cart[1,:],s=100,c='b',marker='o')

should_keep = []
keep_full =[]
for ifk in range(fullbz.nf):
    qvec = fk_cart[:,ifk] - kc_cart
    dist = sqrt(dot(qvec,qvec))
    # K' coordinate hard-coded in for now
    if keep_Kprime:
        qvecp = fk_cart[:,ifk] - kc_cart*2
        distp = sqrt(dot(qvecp,qvecp))
    else:
        distp = qmax+1
    if abs(dist)<qmax or abs(distp)<qmax:
        should_keep += [fullbz.indr[ifk]-1]
        keep_full += [ifk]
        if plot_patch:
            plt.scatter(fk_cart[0,ifk],fk_cart[1,ifk],s=100,c='r',marker='o')
if plot_patch:
    plt.show()
           
should_keep=unique(should_keep)
print "Keeping %i k-points."%(len(should_keep))

# Calculating allowed q for decimated grid
if print_q:
    f_qpt=open('qpoints.dat','w')
    nkf=len(keep_full)
    print nkf, shape(fullbz.fk)
    for ik in range(nkf):
        print "%i/%i q points complete"%(ik,nkf)
        for ikp in range(nkf):
            qvec=fullbz.fk[:,keep_full[ik]]-fullbz.fk[:,keep_full[ikp]]
            if qvec[0]<-TOL:
                qvec[0]=1+qvec[0]
            if qvec[1]<-TOL:
                qvec[1]=1+qvec[1]
            if ik==0 and ikp==0:
                qpoints=array([qvec])
            else:
                qpoints=append(qpoints,[qvec],axis=0)
                
    # Throw out repeats
    q_temp=ascontiguousarray(qpoints).view(dtype((void, qpoints.dtype.itemsize * qpoints.shape[1])))
    _, idx = unique(q_temp, return_index=True)
    qpoints = qpoints[idx]
    print "Calculated %i qpoints in full BZ"%(len(qpoints))

    # Fold q-points to irreducible wedge
    # Note: For now, it only works for  Gamma centered k-grid
    keep_q=[]
    print len(qpoints)
    for iq in range(len(qpoints)):
        if iq%500==0:
            print "%f percent Done with q-vecs \n"%(iq*100./len(qpoints))
        for ik in range(fullbz.nf):
            if all(abs(qpoints[iq,:]-fullbz.fk[:,ik])<TOL):
                keep_q += [fullbz.indr[ik]-1]
                break
            if ik==fullbz.nf-1:
                print "Error: q-point not found",qpoints[iq,0],qpoints[iq,1]
    keep_q=unique(keep_q)
    print "Number of q-points = ",len(keep_q)
    for iq in range(len(keep_q)):
        f_qpt.write('%13.10f %13.10f %13.10f 1.0 \n'%(fullbz.rk[0,iq],fullbz.rk[1,iq],fullbz.rk[2,iq]))

f_qpt.close()
if raw_input('Writing data to %s. Are you sure? [y/N] '%(fname_out))!='y':
    sys.exit(0)

wfn_out = wfnIO()
wfn_out.__dict__ = wfn_in.__dict__.copy()
wfn_out.nk = len(should_keep)
wfn_out.ngk = wfn_out.ngk[should_keep]
wfn_out.ngkmax = amax(wfn_out.ngk)
wfn_out.kw = wfn_out.kw[should_keep]
wfn_out.kw[:] = 1.0/wfn_out.nk
wfn_out.kpt = wfn_out.kpt[:, should_keep]
wfn_out.ifmin = wfn_out.ifmin[should_keep, :]
wfn_out.ifmax = wfn_out.ifmax[should_keep, :]
wfn_out.energies = wfn_out.energies[:, should_keep, :]
wfn_out.occupations = wfn_out.occupations[:, should_keep, :]
wfn_out.f = None
wfn_out.fname = fname_out
wfn_out.write_header(full=True)

ng_max = amax(wfn_in.ngk)
gvec = empty((3, ng_max), dtype=int, order='F')
data = empty((ng_max, wfn_in.ns), dtype=get_numpy_flavor(wfn_in.flavor), order='F')
k_done = 0
for ik in range(wfn_in.nk):
    if ik in should_keep:
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

print 'All done!'
