#!/usr/bin/env python

# Calculates the velocity vector using kp perturbation theory.

from numpy import *
from bgwtools.IO.wfn import wfnIO

def get_velocity_kp(wfn, ib, verbose=True):
    '''Calculates the velocities for band ib and all kpts using:
    v_k = k + \sum_G G |u(G)|^2
    '''

    if not isinstance(wfn, wfnIO):
        raise TypeError('Expected wfn type, got %s.'%(type(wfn)))

    vel = zeros((3,wfn.nk), dtype=float)
    for ik in xrange(wfn.nk):
        if verbose:
            print 'k-point: %d/%d'%(ik+1, wfn.nk)
        kpt = wfn.kpt[:,ik]
        ng = wfn.ngk[ik]
        gvec = empty((3,ng), dtype=int)
        wfn.read_gvectors(gvec)
        cg = empty((wfn.ngk[ik], wfn.ns), dtype=float)
        for ibp in xrange(wfn.nbands):
            if ibp==ib:
                wfn.read_data(cg)
           	cg2 = abs(cg[:,0])**2
                vel[:,ik] = kpt + dot(gvec, cg2)
            else:
                wfn.read_data()
        del cg2
        del cg
        del gvec

    return vel

if __name__=='__main__':
    # A simple test, great for graphene, but should work for
    # any system.

    import sys

    if len(sys.argv)!=3:
        print('usage: %s WFN band')
        sys.exit(1)

    fname = sys.argv[1]
    ib = int(sys.argv[2])

    plot_vel = True
    plot_en = False
    use_cartesian = False

    wfn = wfnIO(fname)
    vel = get_velocity_kp(wfn, ib)
    kpts = wfn.kpt

    import matplotlib.pyplot as plt
    import scipy.linalg

    if use_cartesian:
        M = linalg.cholesky(wfn.bdot).T
        kpts = dot(M, kpts)
        M = linalg.inv(M)
        vel = tensordot(M, vel, (0,0))

    if plot_vel:
        scale = 10
        head = 25
        plt.quiver(kpts[0], kpts[1],
            vel[0], vel[1], angles='xy', scale=scale, zorder=10,
            headwidth=head, headlength=head, width=2e-3)

    if plot_en:
        en = wfn.energies[ib,:,0]
        plt.contour(kpts[0], kpts[1], en, linewidths=2)
    
    plt.show()
    
