#!/usr/bin/env python

# This script combines multiple charge density files from QE into a single
# charge-density.dat file. This is useful if you have a weakly interacting
# system, such as substrate + adsorbate, where the individual charge densities
# converge faster than the total one.
#
# This script also support FFT grids that are not the same. When this happens,
# we assume that you want to create a supercell out of the two materials, and
# we created repeated copies of each rho(r) until we get a comensureate cell.
#
# Felipe H. da Jornada, Feb 2015.

import numpy as np
from bgwtools.IO.qe_xml import IotkFile
import re
from get_rho import get_rho
from fractions import gcd

do_plot = False

def lcm(a,b):
    return abs(a)*abs(b)/gcd(a,b) if a and b else 0


def combine_rho(fname_in1, fname_in2, fname_out):
    rho1 = get_rho(fname_in1)
    rho2 = get_rho(fname_in2)
    shape1 = rho1.shape
    shape2 = rho2.shape
    assert len(shape1)==3 and len(shape2)==3
    assert rho1.dtype==rho2.dtype

    shape_out = tuple([lcm(shape1[i],shape2[i]) for i in range(3)])
    fact1 = np.array(shape_out)/np.array(shape1)
    fact2 = np.array(shape_out)/np.array(shape2)
    print fact1, fact2
    print('Output FFTgrid: {} x {} x {}'.format(*shape_out))
    rho_out = np.zeros(shape_out, dtype=rho1.dtype, order='F')
    shapes = (shape1, shape2)
    facts = (fact1, fact2)
    rhos = (rho1, rho2)
    for i in range(2):
        fact = facts[i]
        nx, ny, nz = shapes[i]
        for iz in range(fact[2]):
            for iy in range(fact[1]):
                for ix in range(fact[0]):
                    rho_out[ix*nx:(ix+1)*nx, iy*ny:(iy+1)*ny, iz*nz:(iz+1)*nz] += rhos[i]

    frho_out = IotkFile(fname_out, mode='wb')
    f = frho_out.write_field
    f((0,5,'<?iotk version="1.2.0"?>'))
    f((0,5,'<?iotk file_version="1.0"?>'))
    f((0,5,'<?iotk binary="T"?>'))
    f((0,5,'<?iotk qe_syntax="F"?>'))
    f((0,1,'<Root>'))
    f((0,1,'<CHARGE-DENSITY>'))
    f((0,3,'<INFO nr1="{}" nr2="{}" nr3="{}"/>'.format(*rho_out.shape)))
    for iz in range(rho_out.shape[2]):
        buf = rho_out[...,iz].ravel(order='F')
        size = len(buf)
        f((0,1,'<z.{} type="real" size="{}" kind="8">'.format(iz+1,size)))
        f((1,'d',buf))
        f((0,2,'</z.{}>'.format(iz+1)))
    f((0,2,'</CHARGE-DENSITY>'))
    f((0,2,'</Root>'))

    if do_plot:
        nr3 = rho_out.shape[2]
        import seaborn as sns
        plt = sns.plt
        #plt.imshow(rho[:,:,int(nr3*8/10)])
        plt.plot(np.arange(nr3), rho_out.sum(axis=0).sum(axis=0))
        #plt.colorbar()
        plt.show()


if __name__=="__main__":
    import sys

    if len(sys.argv)<4:
        print 'usage: {} charge-density_in1.dat charge-density_in2.dat charge-density_out.dat'.format(sys.argv[0])
        print
        sys.exit(1)
    fname_in1, fname_in2, fname_out = sys.argv[1:4]
    combine_rho(fname_in1, fname_in2, fname_out)

