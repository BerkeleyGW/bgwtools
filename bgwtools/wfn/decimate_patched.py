#!/usr/bin/env python

# This script takes a uniform k-grid and calculates k-points on a patch of given radius.
# Given a fullbz.dat or kgrid.x log file, the script finds k-points and q-points on a 
# patch in the brillouin zone.
# If given an input wave function, it can write out a smaller wave function containing only points in the patch.
# centered on a given point. The radius of the patch, qmax, is given in cartesian coordinates,
# while the center of the patch is given in crystal coordinates.
#
# Diana Y. Qiu (March 2014)
# Adapted from FHJ's decimate.py

from numpy import *
import sys

# Options
ryd = 13.60569253
TOL=1e-6
print_q=True         # Print a file with all the allowed q-vectors
print_k=True         # Print decimated k-points to file
plot_patch=False       # Plot the patch
keep_Kprime=False      # If the patch is centered on K, also include the K' point
read_wfn = False       # Read a wave function file. Requires fullbz.dat file if true. Requires kgrid.log if false.
write_wfn= False        # Write a wave function with the decimated k-points
use_symmetries=False  

def move_bz(k,bvec,return_length=False):
    # Moves vector k to the 1st Bz
    k2_min = inf
    for dGx in [-1,0,1]:
        for dGy in [-1,0,1]:
            k_tmp = k + [dGx,dGy,0]
            k2_tmp = sum(dot(bvec,k_tmp)**2)
            if (k2_tmp<k2_min):
                k_min = k_tmp
                k2_min = k2_tmp
    if return_length:
        return k_min, sqrt(k2_min)
    else:
        return k_min

def write_decimated(wfn,keep,bz,use_symmetries):
    if use_symmetries:
        keep_reduced = keep
    else:
        keep_reduced = bz.indr[keep]-1
        keep = keep[argsort(keep_reduced)] # sort k's by reduced index
    wfn_out = wfnIO()
    wfn_out.__dict__ = wfn.__dict__.copy()
    wfn_out.nk = len(keep)
    if use_symmetries:
        wfn_out.ngk = wfn_out.ngk[keep]
        wfn_out.kw = wfn_out.kw[keep]
    else:
        wfn_out.ngk = wfn_out.ngk[bz.indr[keep]-1]
        wfn_out.kw = wfn_out.kw[bz.indr[keep]-1]
    wfn_out.ngkmax = amax(wfn_out.ngk)
    wfn_out.kw[:] = 1.0/wfn_out.nk
    if use_symmetries:
        wfn_out.kpt = wfn_out.kpt[:, keep]
        wfn_out.ifmin = wfn_out.ifmin[keep, :]
        wfn_out.ifmax = wfn_out.ifmax[keep, :]
        wfn_out.energies = wfn_out.energies[:, keep, :]
        wfn_out.occupations = wfn_out.occupations[:, keep, :]
    else:
        wfn_out.kpt = bz.fk[:, keep]
        wfn_out.ifmin = wfn_out.ifmin[bz.indr[keep]-1, :]
        wfn_out.ifmax = wfn_out.ifmax[bz.indr[keep]-1, :]
        wfn_out.energies = wfn_out.energies[:, bz.indr[keep]-1, :]
        wfn_out.occupations = wfn_out.occupations[:, bz.indr[keep]-1, :]
    wfn_out.f = None
    wfn_out.fname = fname_out
    wfn_out.write_header(full=True)

    ng_max = amax(wfn_in.ngk)
    gvec = empty((3, ng_max), dtype=int, order='F')
    if use_symmetries:
        data = empty((ng_max, wfn_in.ns), dtype=get_numpy_flavor(wfn_in.flavor), order='F')
    else:
        data = empty((wfn_in.nbands,ng_max, wfn_in.ns), dtype=get_numpy_flavor(wfn_in.flavor), order='F')
    k_done = 0

    for ik in range(wfn_in.nk):
        if ik in keep_reduced:
            k_done += 1
            print k_done,'/',wfn_out.nk
            wfn_in.read_gvectors(gvec)
            if use_symmetries:
                wfn_out.write_gvectors(gvec[:,:wfn_in.ngk[ik]])
                for ib in range(wfn_in.nbands):
                    wfn_in.read_data(data)
                    wfn_out.write_data(data[:wfn_in.ngk[ik],:])
                if k_done==wfn_out.nk:
                    break
            else:
                # ik is the reduced index
                # find places where ik occurs in indr
                idx_k = where(keep_reduced==ik)[0]
                for kk in range(len(idx_k)):
                    wfn_out.write_gvectors(gvec[:,:wfn_in.ngk[ik]])
                    for ib in range(wfn_in.nbands):
                        if kk==0:
                            wfn_in.read_data(data[ib,:,:])
                        wfn_out.write_data(data[ib,:wfn_in.ngk[ik],:])
            if k_done==wfn_out.nk:
                break                    
        else:
            wfn_in.read_gvectors()
            for ib in range(wfn_in.nbands):
                wfn_in.read_data()
    print wfn_out


if __name__ == '__main__':

    if len(sys.argv)==5 or (len(sys.argv)==6 and write_wfn):
        print('Usage: %s qmax kx ky kz wfn_in|kgrid.log [wfn_out]'%(sys.argv[0]))
        print('       Don\'t forget to check options before running')
        sys.exit(1)

    # Reading arguments and input files
    qmax=float(sys.argv[1])
    k_center=array([float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4])])
    fname_in = sys.argv[5]
    if write_wfn:
        if not read_wfn:
            print "Error: cannot write decimated wave function without input wave function"
            sys.exit(1)
        fname_out = sys.argv[6]
    if read_wfn:
        from bgwtools.IO.wfn import wfnIO
        from bgwtools.common.common import get_numpy_flavor
        from bgwtools.IO.fullbz import fullbzIO
        wfn_in = wfnIO(fname_in, full=True)
        M = wfn_in.bvec
        # Convert to cartesian coordinates
        rk_cart=tensordot(wfn_in.bvec,wfn_in.kpt, axes=([1],[0]))
        fullbz=fullbzIO('fullbz.dat')
        print "Dimensions of input k-grid: %i %i %i"%(wfn_in.kgrid[0],wfn_in.kgrid[1],wfn_in.kgrid[2])
    else:
        from bgwtools.IO.kgrid_log import kgridlogIO
        fullbz = kgridlogIO(sys.argv[5])
        norm = sqrt(fullbz.bvec[0,0]**2 +fullbz.bvec[1,0]**2)
        fullbz.bvec =fullbz.bvec /norm
        M = fullbz.bvec
        print "Dimensions of input k-grid: %i %i %i"%(fullbz.kgrid[0],fullbz.kgrid[1],fullbz.kgrid[2])
    fk_cart=tensordot(M,fullbz.fk, axes=([1],[0]))
    nf = fullbz.nf
    nrk = fullbz.nrk
    print "%i k-points on full grid "%(nf)
    print "%i k-points on reduced grid"%(nrk)

    #Find k-points within the patch
    if plot_patch:
        from matplotlib import pyplot as plt
        plt.scatter(fk_cart[0,:],fk_cart[1,:],s=20,c='b',marker='o')
    keep_reduced = []
    keep_full =[]
    for gx in range(-1,2):
        for gy in range(-1,2):
            for ifk in range(nf):
                kc_cart=dot(M,k_center + array([gx,gy,0.0]))
                qvec = fk_cart[:,ifk] - kc_cart
                dist = sqrt(dot(qvec,qvec))
                # K' coordinate hard-coded in for now
                if keep_Kprime:
                    kc_cartp=dot(M,k_center+array([1./3+gx,1./3+gy,0.0]))
                    qvecp = fk_cart[:,ifk] - kc_cartp
                    distp = sqrt(dot(qvecp,qvecp))
                else:
                    distp = qmax+1
                if abs(dist)<qmax or abs(distp)<qmax:
                    keep_reduced += [fullbz.indr[ifk]-1]
                    keep_full += [ifk]
                    if plot_patch:
                        plt.scatter(fk_cart[0,ifk],fk_cart[1,ifk],s=20,c='r',marker='o')
    keep_full = unique(keep_full)
    print "Found %i kpts"%(len(keep_full))
    keep_reduced=unique(keep_reduced)
    print "Found %i reduced k-points."%(len(keep_reduced))
    if print_k:
        f_k=open('kpoint.dat','w')
        if use_symmetries:
            f_k.write('%i \n'%(len(keep_reduced)))
            for ik in range(len(keep_reduced)):
                f_k.write('%13.10f %13.10f %13.10f 1.0\n'%(fullbz.rk[0,keep_reduced[ik]],fullbz.rk[1,keep_reduced[ik]],fullbz.rk[2,keep_reduced[ik]]))
        else:
            f_k.write('%i \n'%(len(keep_full)))
            for ik in range(len(keep_full)):
                f_k.write('%13.10f %13.10f %13.10f 1.0\n'%(fullbz.fk[0,keep_full[ik]],fullbz.fk[1,keep_full[ik]],fullbz.fk[2,keep_full[ik]]))
        f_k.close()
    if plot_patch:
        plt.show()

    # Write decimated WFN
    if write_wfn:
        if raw_input('Writing data to %s. Are you sure? [y/N] '%(fname_out))!='y':
            sys.exit(0)
        if use_symmetries:
            write_decimated(wfn_in,keep_reduced,fullbz,use_symmetries)
        else:
            write_decimated(wfn_in,keep_full,fullbz,use_symmetries)            
        print "Finished writing decimated WFN!"

    # Find q-vecs in patch
    if not print_q:
        sys.exit(1)

    if not use_symmetries:
        qmax = qmax*2*1.005
        keep_reduced = []
        keep_full =[]
        for ifk in range(nf):
            qmin,dist = move_bz(fullbz.fk[:,ifk],M,return_length=True)
            if abs(dist)<qmax:
                keep_reduced += [fullbz.indr[ifk]-1]
                keep_full += [ifk]
        print "Found %i qpts"%(len(keep_full))
        keep_reduced=unique(keep_reduced)
        print "Keeping %i reduced q-points."%(len(keep_reduced))
        f_k=open('qpoints.dat','w')
        f_k.write('%i \n'%(len(keep_reduced)))
        for ik in range(len(keep_reduced)):
            f_k.write('%13.10f %13.10f %13.10f 1.0\n'%(fullbz.rk[0,keep_reduced[ik]],fullbz.rk[1,keep_reduced[ik]],fullbz.rk[2,keep_reduced[ik]]))
        f_k.close()
    else:
        nkf = len(keep_full)
        keep_q = []
        import scipy.spatial
        tree  = scipy.spatial.cKDTree(fullbz.fk.T)
        for ik in range(nkf):
            for ikp in range(nkf):
                for gx in range(-1,2):
                    for gy in range(-1,2):
                        q = fullbz.fk[:,keep_full[ik]] - fullbz.fk[:,keep_full[ikp]]  + array([gx,gy,0])
                        qmin = move_bz(q,M) + 1.e-6
                        qmin = qmin -  floor(qmin)
                        d,iq = tree.query(qmin)
                        assert d<5.e-6,"min dist iq %f,  %i for %f %f "%(d,iq,qmin[0],qmin[1])
                        keep_q += [fullbz.indr[iq]-1]
            print "%i / %i completed"%(ik,nkf)
        print "Found %i qpts"%(len(keep_q))
        keep_q = unique(keep_q)
        print "Keeping %i reduced q-points."%(len(keep_q))
        f_k=open('qpoints.dat','w')
        f_k.write('%i \n'%(len(keep_q)))
        for ik in range(len(keep_q)):
            f_k.write('%13.10f %13.10f %13.10f 1.0\n'%(fullbz.rk[0,keep_q[ik]],fullbz.rk[1,keep_q[ik]],fullbz.rk[2,keep_q[ik]]))
        f_k.close()

    print 'All done!'
