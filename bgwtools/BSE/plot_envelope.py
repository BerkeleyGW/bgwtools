#!/usr/bin/env python
# Given an eigenvectors file, plots the envelope of the exciton wave function, including the phase, 
# in k-space and real-space. Only works in 2D. Please read the various options before using.
# Diana Qiu (July 2014)

import numpy as np
import scipy.spatial
import sys
from matplotlib import pyplot as plt
from bgwtools.IO.wfn import wfnIO
#from bgwtools.IO.eigenvectors import eigenvectorsIO
import evecs_common

# Options
TOL_ZERO = 1.e-6
plot_realspace = True # If true, take FFT and plot real-space envelope
find_slice = True # Plot a slice of the real-space wave-function along the line between the center point and maximum
square_amplitude = True # Plot the sum of the squares of the wavefunction
plot_phasepath = False # Plot the path over which the phase is calculated


def get_phase(wfn,cg,nv,nc,ik,ikp,phase_old):
    # calculate phase between two bloch waves
    mtxel = np.zeros(nv+nc,dtype=np.complex128)
    #for ig in range(wfn.ngkmax):
    #    for igp in range(wfn.ngkmax):
    #        mtxel += np.conj(cg[:,ik,ig])*cg[:,ikp,igp]
    for ib in range(nv+nc):
        mtxel[ib] = np.tensordot( np.conj(cg[ib,ik,:]),cg[ib,ikp,:],axes=1)
    #print np.shape(mtxel)
    phase = np.conj(phase_old) * mtxel / np.absolute(mtxel)
    #print np.shape(phase)
    return phase

def get_phase_A(evecs,NS,nv,nc,ik,ikp,phase_old):
    # calculate phase between two Asvck's
    #mtxel = np.zeros((NS,nv,nc),dtype=np.complex128)
    mtxel = np.conj(evecs[0,ik,0,0])*evecs[0,ikp,0,0]
    #print mtxel
    phase = phase_old * mtxel/np.absolute(mtxel)
    #print phase_old*np.conj(evecs[0,ik,0,0])
    #print np.conj(phase)*evecs[0,ikp,0,0]

    assert np.all(np.imag(phase_old*np.conj(evecs[0,ik,0,0])) <1.e-14)
    assert np.all(np.imag(np.conj(phase)*evecs[0,ikp,0,0]) <1e-14)
    assert np.all(np.imag( np.conj(phase)*phase_old*mtxel/np.absolute(mtxel))<1.e-14)
    assert (np.all(np.real( np.conj(phase)*phase_old*mtxel/np.absolute(mtxel))< 1. + 1e-12) and 
        np.all(np.real( np.conj(phase)*phase_old*mtxel/np.absolute(mtxel))> 1. - 1e-12 ))

    return phase
    

def DFT_matrix(N):
    i, j = np.meshgrid(np.arange(N), np.arange(N))
    omega = np.exp( - 2 * np.pi * 1j / N )
    W = np.power( omega, i * j ) / np.sqrt(N)
    return W

def read_wfn(wfn,nv=1,nc=1,need_gvec=False):
    nk = wfn.nk
    vb = len(np.nonzero(wfn.occupations[:,0])[0]) - 1
    nb = wfn.nbands
    nvc = nv+nc
    #print wfn.ngkmax
    #print np.amax(wfn.ngk)
    gvecs_k = np.zeros((nk,3,wfn.ngkmax))
    data = np.zeros((nvc,nk,wfn.ngkmax),dtype=np.complex128)
    print "WFN valence band is: %i"%(vb+1)
    print "Reading %i bands (%i val, %i cond) from WFN"%(nvc,nv,nc)
    for ik in range(nk):
        print "Read %i/%i kpts from WFN"%(ik +1,nk)
        wfn.read_header()
        ngk = wfn.ngk[ik]
        if need_gvec:
            wfn.read_gvectors(gvec=gvecs_k[ik,:,:])
        else:
            wfn.read_gvectors()
        ib_vc = 0
        for ib in range(nb):
            if ib>vb-nv and ib<=vb+nc:
                wfn.read_data(data=data[ib_vc,ik,:,np.newaxis])
                ib_vc +=1
            else:
                wfn.read_data()
    print "Finished reading WFN!"
    if need_gvec:
        return gvecs_k,data
    else:
        return data

if __name__ == '__main__':

    if len(sys.argv)!= 5 and len(sys.argv)!=7:
        print "Usage: %s WFN_fi eigenvectors <first state> <last state> [nv nc]"
        print('Note: exciton indices follow Fortran convention (first index=1)')
        sys.exit(1)

    fname_wfn = sys.argv[1]
    fname_evecs = sys.argv[2]
    iS0 = int(sys.argv[3])
    iSf = int(sys.argv[4])
    if (iS0<1):
        print('Must have <first state> > 0')
        sys.exit(1)
    if (iSf<iS0):
        print('Must have <last state> >= <first state>')
        sys.exit(1)

    NS = iSf
    if len(sys.argv)==7:
        nv = int(sys.argv[5])
        nc = int(sys.argv[6])
    else:
        nv = 1
        nc = 1

    # Read eigenvectors
    wfn=wfnIO(fname_wfn)
    print 'Reading Evecs'
    k_crys, evecs, evals = evecs_common.read_evecs(fname_evecs, NS, nv, nc)
    print np.shape(evecs)
    # Convert kpts to cartesian coordinates
    Bv = wfn.bvec[:2,:][:,:2]
    Bv /= Bv[0,0]
    k_crys = k_crys[:2,:]
    k_crys = k_crys - np.floor(k_crys)

###########################################################################################
# Calculate Phase difference between k-points along a path starting from the central k-point
##########################################################################################

    # Set-up grid of k-points 
    # -1 = buffer zone, 0 = phase not yet calculated, 1 = phase already calculated
    #cg = read_wfn(wfn) #cg[nv+nc,nk,ngkmax]
    ind_x = np.around(k_crys[0,:] * wfn.kgrid[0])
    ind_y = np.round(k_crys[1,:] * wfn.kgrid[0])    
    nx = int(np.amax(ind_x) - np.amin(ind_x))
    nx_full = nx+3
    ny = int(np.amax(ind_y) - np.amin(ind_y))
    ny_full = ny+3
    print "Dimensions of grid: ",nx,ny
    calculated = -np.ones((nx_full,ny_full))
    ind_x_map = map(int,ind_x - np.amin(ind_x) + 1)
    ind_y_map = map(int,ind_y - np.amin(ind_y) + 1)
    calculated[ind_x_map,ind_y_map] = 0
    if plot_phasepath:
        plt.imshow(calculated)
        plt.show()
    #tree_wfn_kpt = scipy.spatial.cKDTree(wfn.kpt.T)
    tree_wfn_kpt = scipy.spatial.cKDTree(k_crys.T)

    # Cycle over points in path starting from the center and reacinh every k-point
    import operator
    new_pts = [(nx_full/2,ny_full/2)]
    phases = np.ones(1,dtype=np.complex128)
    ik_map = []# maps index of phases and new_pts list to index of wfn.kpts
    kk = np.array([0.,0.])
    kkp = np.array([0.,0.])
    ii=0
    for pt in new_pts:
        pt = new_pts[ii]
        calculated[pt[0],pt[1]]=1
        # Set the phase of the first point to 1
        if len(new_pts)==1:
            # find index of kpt corrsponding to pt
            kk[0] = np.float(pt[0] - 1 + np.amin(ind_x))/wfn.kgrid[0]
            kk[1] = np.float(pt[1] - 1 + np.amin(ind_y))/wfn.kgrid[1]
            dist, ik = tree_wfn_kpt.query(kk)
            assert dist<TOL_ZERO
            ik_map.append(ik)
            #phases[ii,:] = np.exp(1j*np.angle(evecs[iS0-1:iSf,ik,0,0]))
            phases [ii] = np.exp(1j*np.angle(evecs[0,ik,0,0]))
        for delta in [(-1,0),(1,0),(0,-1),(0,1)]:
            working_pt = tuple(map(operator.add,pt,delta))
            if calculated[working_pt[0],working_pt[1]]==0:
                # find index of kpt corresponding to working_pt
                kkp[0] = np.float(working_pt[0] - 1 + np.amin(ind_x))/wfn.kgrid[0]
                kkp[1] = np.float(working_pt[1] - 1 + np.amin(ind_y))/wfn.kgrid[1]
                dist, ikp = tree_wfn_kpt.query(kkp)
                assert dist<TOL_ZERO
                ik_map.append(ikp)
                calculated[working_pt[0],working_pt[1]]=1
                #new_phase = get_phase(wfn,cg,nv,nc,ik_map[ii],ikp,phases[ii])
                #new_phase = get_phase_A(evecs,NS-iS0+1,nv,nc,ik_map[ii],ikp,phases[ii,:])
                new_phase = get_phase_A(evecs,NS,nv,nc,ik_map[ii],ikp,phases[ii])
                phases=np.append(phases,[new_phase],axis=0)
                #print np.shape(phases)
                new_pts.append(working_pt)
        if plot_phasepath:
            plt.imshow(calculated)
            plt.savefig('kpt_phase_path-'+str(len(new_pts)).zfill(4)+'.png')
        print "%i / %i kpts completed"%(ii,wfn.nk)
        ii += 1
        #plt.show()

    assert len(phases[:])==wfn.nk
    assert len(ik_map)==wfn.nk

################################################################
# Multiply eigenvectors by phase, interpolate, and plot in k-space
################################################################

    # multiply by phase(ib,ik) sum over valence and conduction bands
    ik_map_i = np.argsort(ik_map)
    evecs[:iSf,:,0,0] = evecs[:iSf,:,0,0] * np.conj(phases[ik_map_i[:]]).T #*phases[ik_map_i[:],iv] #(iS,ik,ic,iv)
    #evecs = abs(evecs)**2
    evecs = np.real(np.sum(np.sum(evecs, axis=3), axis=2))
    if square_amplitude:
        v0 = np.sum(evecs[iS0-1:iSf,:], axis=0)**2
    else:
        v0 = evecs[iS0-1,:]
        if iSf-iS0>0: print "Warning: Summing wave functions of multiple states may not make sense!"
    print "Max and Min before norm: ",np.amax(v0), np.amin(v0)
    norm = np.amax(np.abs(v0))
    print "Norm: ", norm
    v0 = v0 / norm
    print  "Max and Min after norm: ",np.amax(v0), np.amin(v0)
    
    k_cart = np.tensordot(Bv, k_crys, axes=([1],[0]))
    x0 = k_cart[0,:]
    x0 += np.random.rand( len(x0))*1.e-5
    y0 = k_cart[1,:]
    y0 += np.random.rand( len(y0))*1.e-5
    N = 300
    from matplotlib.mlab import griddata
    import matplotlib.cm as cm
    xi = np.linspace(np.amin(x0), np.amax(x0), N)
    print "x-range: ",np.amin(x0), np.amax(x0)
    Lx = np.amax(x0)-np.amin(x0)
    yi = np.linspace(np.amin(y0),np.amax(y0), N)
    print "y-range: ",np.amin(y0),np.amax(y0)
    Xi,Yi = np.meshgrid(xi,yi)
    from scipy.interpolate import griddata as sci_griddata
    #V = griddata(x0, y0, v0, xi, yi, interp='nn')
    V_k =  sci_griddata(np.array([x0,y0]).T,v0,np.array([Xi,Yi]).T,fill_value=0.0)
    print V_k
    print np.shape(V_k),np.shape(x0),np.shape(v0)
    im = plt.imshow(-V_k, origin='lower',cmap=cm.RdBu)
    plt.colorbar()
    plt.clim(-.9,.9)
    #plt.axis('off')
    barlen = N*0.45 # bar is 0.1 Angstroms ^-1
    plt.plot([50,50+barlen],[50,50],lw=2,c='k') 
    plt.plot([50,50],[40,60],lw=2,c='k') # bar length
    plt.plot([50+barlen,50+barlen],[40,60],lw=2,c='k') # bar length
    plt.savefig('xcton_%d-%d.pdf'%(iS0,iSf), pad_inches=0)
    plt.show()

    if not plot_realspace:
        sys.exit()

####################################################################
# Fourier transform the envelope of the exciton wave function to real-space
######################################################################

    print "Starting FFT"
    
    nS = iSf-iS0+1
    # Put exciton  in box from [0,1)
    # Check k_crys = wfn.kpt ???
    xct_k = np.zeros((nS,wfn.kgrid[0],wfn.kgrid[1]),dtype=np.float64)
    ind_x = ind_x - 1./3.*300
    ind_y = ind_y - 1./3.*300
    ind_x = np.int32(np.round(ind_x))
    ind_y = np.int32(np.round(ind_y))
    print ind_x[409],ind_y[409]
    xct_k[:,ind_x,ind_y] = np.real(evecs[iS0-1:iSf,:])
    #im = plt.imshow(xct_k,origin='lower')
    #plt.show()

    from scipy import fftpack
    #xct_k = fftpack.ifftshift(xct_k)
    xct_r = np.zeros(np.shape(xct_k),dtype=np.float64)
    for iS in range(nS):
        xct_r[iS] = fftpack.fft2(xct_k[iS])
        xct_r[iS] = fftpack.fftshift( xct_r[iS] )
        print "done with FFT"
        xct_r[iS] = np.real(xct_r[iS])
        xct_flatten=np.ndarray.flatten(xct_r[iS])
        # Convert real-space coordinates from cyrstal to cartesian
        coord_crys = np.meshgrid(np.arange(wfn.kgrid[0]),np.arange(wfn.kgrid[0]))[0]
        xx_crys = np.ndarray.flatten(coord_crys)
        yy_crys = np.ndarray.flatten(coord_crys.T)
        coord_crys = np.array([xx_crys,yy_crys])
        Av = wfn.avec[:2,:][:,:2]
        print Av
        coord = np.dot(Av,coord_crys)
        print np.shape(coord)
        xx = coord[0,:]
        yy = coord[1,:]
        xi = np.linspace(np.amin(xx), np.amax(xx), N)
        yi = np.linspace(np.amin(yy), np.amax(yy), N)
        Xi,Yi = np.meshgrid(xi,yi)
        if iS==0:
            V = sci_griddata(coord.T,xct_flatten,np.array([Xi,Yi]).T,fill_value=0.0)
            if square_amplitude or nS>1:
                V = V**2
        else:
            V += sci_griddata(coord.T,xct_flatten,np.array([Xi,Yi]).T,fill_value=0.0)**2
            print np.shape(xx_crys),np.shape(xct_flatten),np.shape(Xi),np.shape(V)
    #V=V.reshape(wfn.kgrid[0],wfn.kgrid[1])
    #xct_r /= np.amax(xct_r)
    #im = plt.imshow(V[100:-100,100:-100],cmap=cm.RdBu)
    #print np.unravel_index(np.argmax(V),np.shape(V))
        
    # Print a slice of the real space wave function
    if find_slice:
        center = (149.9,149.9)#((np.amax(Xi)-np.amin(Xi))/2.,(np.amax(Xi)-np.amin(Xi))/2.)
        find_node = True
        if find_node:
            max_point = np.unravel_index(np.argmax(V.T),np.shape(V.T))
            #max_point = (Xi[max_point[0],max_point[1]],Yi[max_point[0],max_point[1]])
        else:
            max_point = (200,200)
        #Draw a line between the center and max point
        slope = np.float(max_point[1]-center[1])/np.float(max_point[0]-center[0])
        print "For real-space slice: "
        print "slope = ",slope
        print "center = ",center
        print "max_point = ",max_point
        b = center[1] - center[0]*slope
        print " y-intercept = ",b
        x = np.arange(100,center[0],0.1)
        line = b+slope*x
    im = plt.imshow(V,cmap=cm.jet)
    plt.colorbar(im)
    if find_slice:
        plt.plot(x,line,'--k',lw=2)
    plt.savefig('r_xcton_%d-%d.pdf'%(iS0,iSf), pad_inches=0)
    plt.show()
    if not find_slice:
        sys.exit()

    #print np.amax(Xi),np.amin(Xi),np.amax(Yi),np.amin(Yi)
    from scipy import ndimage
    print "Beginning real-space interpolation of exciton cross-section"
    zi = ndimage.map_coordinates(V, np.vstack((x,line)))
    print "Finished Interpolation"

    #xct_line = sci_griddata(np.array([np.ndarray.flatten(Xi),np.ndarray.flatten(Yi)]).T,np.ndarray.flatten(V),np.array([x,line]).T,fill_value=0.0)
    #print np.shape(np.array([np.ndarray.flatten(Xi),np.ndarray.flatten(Yi)])),np.shape(np.ndarray.flatten(V)),np.shape(xct_line)
    #print np.amax(xct_line),np.argmax(xct_line)
    r=np.sqrt((x-center[0])**2+(line-center[1])**2)
    nm = 0.0529177
    r *=wfn.alat * nm*np.pi*2.
    norm = np.sum(r[:]*zi[:])
    print np.sum(r*zi)
    zi = zi/norm
    print np.sum(r*zi)
    print "Writing to file"
    f = open('xct_slice_%d-%d.dat'%(iS0,iSf),'w')
    f.write("# r (nm), |Psi(r)|**2 \n")
    for ii in range(200):
       f.write('%15.9f %15.9f\n'%(r[-ii],zi[-ii]))
    plt.plot(r[-200:],zi[-200:])
    plt.show()
