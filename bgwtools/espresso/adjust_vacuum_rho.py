#!/usr/bin/env python

# TODO
#
# Felipe H. da Jornada, Apr 2015.

import numpy as np
from bgwtools.IO.qe_xml import IotkFile
import re
from rho_io import read_rho, write_rho
import scipy.ndimage


def get_int3_vector(s):
    obj = re.findall(r'[^ ,;]+', s)
    qq = np.array(map(int, obj))
    assert qq.shape==(3,)
    return qq


def get_float3_vector(s):
    obj = re.findall(r'[^ ,;]+', s)
    qq = np.array(map(float, obj))
    assert qq.shape==(3,)
    return qq


def wrap_array(A):
    '''Extends all axes by unity by appending the first hyperplanes to the end
    of each axis.
    '''
    ndim = len(A.shape)
    shape = np.array(A.shape) + 1
    data = np.empty(shape, dtype=A.dtype)
    ind = [slice(None, A.shape[jdim]) for jdim in range(ndim)]
    data[ind] = A
    for idim in range(len(data.shape)):
        ind_from = [slice(None) if jdim!=idim else 0 for jdim in range(ndim)]
        ind_to = [slice(None) if jdim!=idim else -1 for jdim in range(ndim)]
        data[ind_to] = data[ind_from]
    return data


def transfer_rho(fname_in, fname_out, shape_out, scale, plot):
    rho_in = read_rho(fname_in)
    shape_in = rho_in.shape

    print('Input FFTgrid: {} x {} x {}'.format(*shape_in))
    print('Output FFTgrid: {} x {} x {}'.format(*shape_out))

    rho_in_ext = wrap_array(rho_in)
    mode = 'constant'
    order = 1
    mat = np.array(shape_in, dtype='d')/np.array(shape_out, dtype='d') * scale
    rho_out = scipy.ndimage.interpolation.affine_transform(rho_in_ext, mat,
        output_shape=shape_out, mode='constant', order=order, cval=0)

    # Charge conservation:
    #charge_in = sum(rho_in) * vol_in / product(shape_in)
    #charge_out = sum(rho_out) * vol_out / product(shape_out) = charge_in
    rho_out *= rho_in.sum()/rho_out.sum() / np.product(mat)

    write_rho(fname_out, rho_out)

    if plot:
        import seaborn as sns
        plt = sns.plt
        for i, (title, rho) in enumerate(zip(('Input','Output'), (rho_in, rho_out))):
            nr3 = rho.shape[2]
            plt.subplot(2,1,i+1)
            plt.title(title)
            #plt.imshow(rho[:,:,int(nr3*8/10)])
            y = rho.sum(axis=0).sum(axis=0)/np.product(rho.shape[:2])
            plt.plot(np.arange(nr3), y)
            plt.xlim(0, nr3-1)
            plt.ylim(0, np.amax(y))
            #plt.colorbar()
        plt.tight_layout()
        plt.show()


if __name__=="__main__":
    import argparse
    import sys

    desc = 'Resamples the charge density from QE from one FFT grid to another.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('cd_in', help='Input charge-density.dat file.')
    parser.add_argument('cd_out', help='Output charge-density.dat file.')
    parser.add_argument('FFTgrid', type=get_int3_vector,
        help='New FFT grid. Format is nx,ny,nz.')
    parser.add_argument('scale', type=get_float3_vector, default=np.ones((3,)), nargs='?',
        help=('By how much are we scaling each real-space lattice vectors? '
        'We assume that we don\'t change the angles of the unit cell, and that '
        'we simply add/remove vacuum as we change the unit cell.'))
    parser.add_argument('--plot', action='store_true', default=False,
        help='Plot initial/final charge densities.')
    args = parser.parse_args()

    transfer_rho(args.cd_in, args.cd_out, args.FFTgrid, args.scale, args.plot)

