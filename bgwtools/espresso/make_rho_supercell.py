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


def make_rho_supercell(fname_in, fname_out, supercell, plot):
    rho_in = read_rho(fname_in)
    shape_in = rho_in.shape
    shape_out = np.array(shape_in) * np.array(supercell)
    rho_out = np.empty(shape_out, dtype=rho_in.dtype)
    print('Input FFTgrid: {} x {} x {}'.format(*shape_in))
    print('Output FFTgrid: {} x {} x {}'.format(*shape_out))

    for i1 in range(supercell[0]):
        start1 = shape_in[0] * i1
        end1 = shape_in[0] * (i1 + 1)
        for i2 in range(supercell[1]):
            start2 = shape_in[1] * i2
            end2 = shape_in[1] * (i2 + 1)
            for i3 in range(supercell[2]):
                start3 = shape_in[2] * i3
                end3 = shape_in[2] * (i3 + 1)
                rho_out[start1:end1,start2:end2,start3:end3] = rho_in

    write_rho(fname_out, rho_out)

    if plot:
        try:
            import sns
            plt = sns.plt
        except:
            import matplotlib.pyplot as plt

        for i, (title, rho) in enumerate(zip(('Input','Output'), (rho_in, rho_out))):
            nr1, nr2, nr3 = rho.shape
            plt.subplot(2,1,i+1)
            plt.title(title)
            #plt.imshow(rho[:,:,int(nr3*8/10)])
            if False:
                y = rho.sum(axis=0).sum(axis=0)/np.product(rho.shape[:2])
                plt.plot(np.arange(nr3), y)
                plt.xlim(0, nr3-1)
                plt.ylim(0, np.amax(y))
            else:
                y = rho.sum(axis=2).sum(axis=1)/np.product(rho.shape[1:])
                plt.plot(np.arange(nr1), y)
                plt.xlim(0, nr1-1)
                plt.ylim(0, np.amax(y))
            #plt.colorbar()
        plt.tight_layout()
        plt.show()


if __name__=="__main__":
    import argparse
    import sys

    desc = 'Resamples the charge density from QE from one FFT grid to another.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('fname_in', help='Input charge-density.dat file.')
    parser.add_argument('fname_out', help='Output charge-density.dat file.')
    parser.add_argument('supercell', type=get_int3_vector,
        help='Supercell size. Format is nx,ny,nz.')
    parser.add_argument('--plot', action='store_true', default=False,
        help='Plot initial/final charge densities.')
    args = parser.parse_args()

    make_rho_supercell(**vars(args))
