#!/usr/bin/env python

# This script extracts the charge density from a QE calculation.
#
# Felipe H. da Jornada, Feb 2015.

import numpy as np
from bgwtools.IO.qe_xml import IotkFile
import re


def get_rho(fname):
    re_tag = re.compile(r'''([^\s]*)=['"]([^'"]+)['"]''')
    re_name = re.compile(r'''<([^\s]*)''')
    type_dict = {'integer':np.int32, 'real':np.float64, 'complex':np.complex128, 'logical':np.int32}
    f = IotkFile(fname)
    print('Reading charge density file '+fname)
    while True:
        # Read all tags
        try:
            obj = f.read_field()
        except IOError:
            break
        if obj[0]==0:
            # This is a text field
            #print obj[1], repr(obj[2])
            tag = obj[2]
            name = re_name.search(tag).group(1)
            attrs = dict(re_tag.findall(tag))
        #print obj[0], tag

        if name=='INFO':
            nr1, nr2, nr3 = [int(attrs[nr]) for nr in ('nr1','nr2','nr3')]
            print('Setting up FFTgrid: {} x {} x {}'.format(nr1, nr2, nr3))
            rho = np.empty((nr1,nr2,nr3), dtype=np.float64, order='F')

        if obj[0]==1 and name[:2]=='z.':
            # This is a binary field containing a slab of the charge density
            iz = int(name.split('.')[-1])
            size = int(attrs['size'])
            #print('Reading slab iz={} with {} entries'.format(iz,size))
            data = obj[2]
            rho[:,:,iz-1] = data.reshape((nr1,nr2), order='F')

    print('Done reading charge density file '+fname)
    print('Sum of charge density: {}'.format(rho.sum()))
    print('')
    return rho


if __name__=="__main__":
    import sys

    if len(sys.argv)<2:
        print 'usage: {} charge-density.dat'.format(sys.argv[0])
        print
        sys.exit(1)
    fname = sys.argv[1]
    rho = get_rho(fname)

    #exit()
    nr3 = rho.shape[2]
    import seaborn as sns
    plt = sns.plt
    interp = 'bicubic'
    #interp = 'catrom'
    #interp = 'sinc'
    iz = nr3//2
    iz = 1045
    #plt.imshow(rho[:,:,iz], cmap=plt.cm.hot, interpolation=interp)
    #plt.colorbar()

    plt.plot(np.arange(nr3), rho.sum(axis=0).sum(axis=0))

    plt.show()
