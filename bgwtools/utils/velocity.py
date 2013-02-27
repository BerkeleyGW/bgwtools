#!/usr/bin/env python

# Calculates the velocity vector for each k-point.
# Requires a WFN file with a ramdomly-shifted k-grid.

from numpy import *
from bgwtools.IO.wfn import wfnIO
TOL_SMALL = 1e-6

def get_velocity(wfn, same_order=True, downsample=0):
    if not isinstance(wfn, wfnIO):
        raise TypeError('Expected wfn type, got %s.'%(type(wfn)))
    nb = wfn.nbands
    ns = wfn.ns
    if downsample>1:
        kgrid = wfn.kgrid
        wfn.energies = wfn.energies.reshape(nb, kgrid[0], kgrid[1], kgrid[2], ns, order='F')
        wfn.kpt = wfn.kpt.reshape(3, kgrid[0], kgrid[1], kgrid[2], order='F')
        ds = downsample
        wfn.kgrid = wfn.kgrid / ds
        wfn.kgrid[wfn.kgrid<1]=1
        wfn.nk = product(wfn.kgrid)
        wfn.energies = wfn.energies[:,::downsample,::ds,::ds,:].reshape(nb, wfn.nk, ns, order='F')
        wfn.kpt = wfn.kpt[:,::ds,::ds,::ds].reshape(3, wfn.nk, order='F')

    kgrid = wfn.kgrid
    kgrid[kgrid==0] = 1
    kpts = wfn.kpt
    nb = wfn.nbands
    nk = wfn.nk
    ns = wfn.ns
    en = wfn.energies
    #print en.flags['F_CONTIGUOUS']
    #print kpts.flags['F_CONTIGUOUS']
    if nk != product(kgrid):
        raise ValueError('kgrid is not uniform.')

    # Consistency check: dk == 1/kgrid
    all_delta_k = kpts[:,1:] - kpts[:,:-1]
    dk = ones(3)
    for idim in range(3):
        delta_k = all_delta_k[idim]
        cond = delta_k>0
        if any(cond):
            dk[idim] = amin(delta_k[cond])
    if any(fabs(dk - 1./kgrid) > TOL_SMALL):
        print wfn.nk
        print wfn.kgrid
        print dk
        print 1./kgrid
        raise ValueError('kgrid not consistent.')
    #dk = 1./kgrid

    # It's easier to roll the axis and calculate the velocities if the 
    # energies are indexed like this:
    en = en.reshape(nb, kgrid[2], kgrid[1], kgrid[0], ns, order='F')

    # Calcualte band velocity
    vel = zeros((3, nb, kgrid[2], kgrid[1], kgrid[0], ns), dtype=float, order='F')
    for idim in range(3):
        if kgrid[2-idim]<2: continue
        # In principle we could calculate both the forward and 
        # backward derivatives...
        #vel[2-idim] = 0.5*(roll(en, -1, idim+1) - roll(en, 1, idim+1))
        vel[2-idim] = roll(en, -1, idim+1) - en

    if same_order:
        return vel.reshape(3, nb, nk, ns, order='F')
    else:
        return en, vel

def get_velocity_kp(wfn, same_order=True):
    if not isinstance(wfn, wfnIO):
        raise TypeError('Expected wfn type, got %s.'%(type(wfn)))
    nb = wfn.nbands
    ns = wfn.ns
    kpts = wfn.kpt
    nk = wfn.nk
    ns = wfn.ns
    #print en.flags['F_CONTIGUOUS']
    #print kpts.flags['F_CONTIGUOUS']

    for ik in xrange(len(kpts)):
        

    # Consistency check: dk == 1/kgrid
    for idim in range(3):
        delta_k = all_delta_k[idim]
        cond = delta_k>0
        if any(cond):
            dk[idim] = amin(delta_k[cond])
    if any(fabs(dk - 1./kgrid) > TOL_SMALL):
        print wfn.nk
        print wfn.kgrid
        print dk
        print 1./kgrid
        raise ValueError('kgrid not consistent.')
    #dk = 1./kgrid

    # It's easier to roll the axis and calculate the velocities if the 
    # energies are indexed like this:
    en = en.reshape(nb, kgrid[2], kgrid[1], kgrid[0], ns, order='F')

    # Calcualte band velocity
    vel = zeros((3, nb, kgrid[2], kgrid[1], kgrid[0], ns), dtype=float, order='F')
    for idim in range(3):
        if kgrid[2-idim]<2: continue
        # In principle we could calculate both the forward and 
        # backward derivatives...
        #vel[2-idim] = 0.5*(roll(en, -1, idim+1) - roll(en, 1, idim+1))
        vel[2-idim] = roll(en, -1, idim+1) - en

    if same_order:
        return vel.reshape(3, nb, nk, ns, order='F')
    else:
        return en, vel

if __name__=='__main__':
    # A simple test, great for graphene, but should work for
    # any system.

    import sys
    fname = sys.argv[1]
    ib = 4
    if len(sys.argv)>2:
        ib = int(sys.argv[2])
    plot_abs = True
    plot_vel = True
    plot_en = True

    wfn = wfnIO(fname)
    en, vel = get_velocity(wfn, same_order=False)

    import matplotlib.pyplot as plt
    import scipy.linalg

    dk = 1./wfn.kgrid
    shift = dk*wfn.kshift
    x = arange(shift[0], 1.0, dk[0])
    y = arange(shift[1], 1.0, dk[1])
    X, Y = meshgrid(x, y)
    Z = zeros_like(X)
    XYZ = row_stack((X[newaxis],Y[newaxis],Z[newaxis]))
    M = linalg.cholesky(wfn.bdot).T
    #M = eye(3)
    XYZ = tensordot(M, XYZ, (1,0))
    M = linalg.inv(M)
    vel = vel[:,ib,0,:,:,0]
    vel = tensordot(M, vel, (0,0))

    if plot_abs:
        v = vel * wfn.kgrid[:,newaxis,newaxis]
        v = sqrt(sum(v**2, axis=0))

        from matplotlib.mlab import griddata
        Ni = 1000
        xi = linspace(amin(XYZ[0]), amax(XYZ[0]), Ni)
        yi = linspace(amin(XYZ[1]), amax(XYZ[1]), Ni)
        X_ = XYZ[0].flatten(order='F')
        Y_ = XYZ[1].flatten(order='F')
        V_ = v.flatten(order='F')

        vi = griddata(X_,Y_,V_,xi,yi,interp='linear')
        plt.imshow(vi, origin='lower', \
            extent=(amin(xi), amax(xi), amin(yi), amax(yi)))
        
        #plt.imshow(v, origin='lower', extent=(0,1,0,1))
        plt.colorbar()
    if plot_vel:
        ds = 1
        v = vel[:,::ds,::ds]*ds
        scale = 2
        plt.quiver(XYZ[0][::ds,::ds], XYZ[1][::ds,::ds], \
            v[0], v[1], angles='xy', scale=scale, zorder=10,\
            headwidth=5, headlength=5, width=1e-3)
    if plot_en:
        e = en[ib,0,:,:,0]
        plt.contour(XYZ[0], XYZ[1], e, linewidths=2)
    
    plt.show()
    
