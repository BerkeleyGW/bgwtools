#!/usr/bin/env python

# Print out epsinv(q) for Gx=Gy=0, Gz=/=0
# Diana Y. Qiu <dqiu@berkeley.edu> (2013)

import sys
import os
from numpy import *
from bgwtools.IO.epsmat import epsmatIO
from bgwtools.IO.wfn import wfnIO

if len(sys.argv)<5:
    print "Usage: %s eps0mat epsmat WFN output"%(sys.argv[0])
    sys.exit(1)

f_eps0mat=sys.argv[1]
f_epsmat=sys.argv[2]
f_wfn=sys.argv[3]
f_out=open(sys.argv[4],'w')
f_out.write("%10s %10s %10s %5s %15s %10s \n"%('# qx','qy','|q|','Gz','eps','ekin'))

eps0mat = epsmatIO(f_eps0mat)
wfn = wfnIO(f_wfn)

print "Reading eps0mat"
print "   Number of G-vectors: ",eps0mat.ng

for ig in range(eps0mat.ng):
    g_indx = eps0mat.isort[0][ig]
    g_indx_i = eps0mat.isort_i[0][g_indx-1]
    if eps0mat.gvec_k[0,g_indx-1]==0 and eps0mat.gvec_k[1,g_indx-1]==0:
        print "   Found head"
        print "   Gvec: ",eps0mat.gvec_k[:,g_indx-1]
        print "   isort,isort_i: ",g_indx,g_indx_i
        print "   ig", ig
        Gz = eps0mat.gvec_k[2,g_indx-1]
        qG = [eps0mat.qpt[0,0], eps0mat.qpt[1,0], Gz]
        ekin = dot(qG,dot(wfn.bdot,qG))
        qq=sqrt(dot(eps0mat.qpt[:,0],dot(wfn.bdot,eps0mat.qpt[:,0])))
        if ekin>eps0mat.ecuts:
            break
        f_out.write("%10.5f %10.5f %10.5f %5i %15.10f %15.10f \n"%(eps0mat.qpt[0,0], \
                                                                       eps0mat.qpt[1,0],qq,Gz,eps0mat.epsmat[0][g_indx_i-1,g_indx_i-1],ekin))

eps = epsmatIO(f_epsmat)
print "\n"
print "Reading epsmat"
print "   Number of G-vectors: ",eps.ng
print "   Number of q-vectors: ",eps.nq
for iq in range(eps.nq):
    for ig in range(eps.ng):
        g_indx = eps.isort[iq][ig] - 1
        g_indx_i = eps.isort_i[iq][g_indx] - 1
        if eps.gvec_k[0,g_indx]==0 and eps.gvec_k[1,g_indx]==0:
            print "   Found head"
            print "   Gvec: ",eps.gvec_k[:,g_indx]
            print "   isort,isort_i: ",g_indx,g_indx_i
            print "   ig", ig
            Gz = eps.gvec_k[2,g_indx]
            qG = [eps.qpt[0,iq], eps.qpt[1,iq], Gz]
            ekin = dot(qG,dot(wfn.bdot,qG))
            qq=sqrt(dot(eps.qpt[:,iq],dot(wfn.bdot,eps.qpt[:,iq])))
            if ekin>eps.ecuts:
                break
            f_out.write("%10.5f %10.5f %10.5f %5i %15.10f %15.10f \n"%(eps.qpt[0,iq], \
                                                                           eps.qpt[1,iq],qq,Gz,eps.epsmat[iq][g_indx_i,g_indx_i],ekin))

print "Done!"
