#!/usr/bin/env python

# This script takes the charge density from a file and resamples it on another
# FFT grid specified by the user.
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


def resample_rho(fname_in, fname_out, shape_out, plot):
    rho_in = read_rho(fname_in)
    shape_in = rho_in.shape

    print('Input FFTgrid: {} x {} x {}'.format(*shape_in))
    print('Output FFTgrid: {} x {} x {}'.format(*shape_out))

    rho_in_ext = wrap_array(rho_in)
    mode = 'constant'
    order = 1
    mat = np.array(shape_in, dtype='d')/np.array(shape_out, dtype='d')
    rho_out = scipy.ndimage.interpolation.affine_transform(rho_in_ext, mat, output_shape=shape_out, mode='constant', order=order, cval=np.nan)

    write_rho(fname_out, rho_out)

    if plot:
        import sns
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
    parser.add_argument('--plot', action='store_true', default=False,
        help='Plot initial/final charge densities.')
    args = parser.parse_args()

    resample_rho(args.cd_in, args.cd_out, args.FFTgrid, args.plot)

