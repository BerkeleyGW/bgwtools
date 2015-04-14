#!/usr/bin/env python

# Converts hdf5 format to old WFN format
# Diana Qiu (Sep 2014)

from numpy import *
import sys
from bgwtools.IO.wfn import wfnIO
from bgwtools.common.common import get_numpy_flavor
from bgwtools.common.common import get_ftype
import numpy as np

if len(sys.argv)!=4:
    print('Usage: %s WFN.h5 WFN_out flavor \n'%(sys.argv[0]))
    print(' Flavor: 1=real, 2=complex')
    sys.exit(1)

fname_in = sys.argv[1]
fname_out = sys.argv[2]
flavor=int(sys.argv[3])

# Read WFN.h5
import h5py
wfn=h5py.File(fname_in,'r')

# Copying Header

# File info
wfn_out = wfnIO()
wfn_out.fname = fname_out
if flavor==2:
    wfn_out.name = 'WFN-Complex'
else:
    wfn_out.name = 'WFN-Real'
wfn_out.ftype = get_ftype('WFN')                                                                                                                                                                    
wfn_out.flavor = flavor
wfn_out.date = ''
wfn_out.time = ''

# crystal
wfn_out.nat = int(wfn['/mf_header/crystal/nat'][...])
wfn_out.adot = wfn['/mf_header/crystal/adot'][...].T
wfn_out.alat = np.float(wfn['/mf_header/crystal/alat'][...])
wfn_out.apos = wfn['/mf_header/crystal/apos'][...].T
wfn_out.atyp = wfn['/mf_header/crystal/atyp'][...]
wfn_out.avec = wfn['/mf_header/crystal/avec'][...].T
wfn_out.bdot = wfn['/mf_header/crystal/bdot'][...].T
wfn_out.blat = np.float(wfn['/mf_header/crystal/blat'][...])
wfn_out.bvec = wfn['/mf_header/crystal/bvec'][...].T
wfn_out.celvol = np.float64(wfn['/mf_header/crystal/celvol'][...])
wfn_out.recvol = np.float64(wfn['/mf_header/crystal/recvol'][...])
wfn_out.cell_symmetry = int(wfn['/mf_header/symmetry/cell_symmetry'][...])
wfn_out.ntran = int(wfn['/mf_header/symmetry/ntran'][...])
wfn_out.mtrx = wfn['/mf_header/symmetry/mtrx'][...].T
wfn_out.tnp = wfn['/mf_header/symmetry/tnp'][...].T

# kpoints
wfn_out.ecutwfn = np.float(wfn['/mf_header/kpoints/ecutwfc'][...])
wfn_out.nbands = int(wfn['/mf_header/kpoints/mnband'][...])
wfn_out.ns = int(wfn['/mf_header/kpoints/nspin'][...])
wfn_out.nk = int(wfn['/mf_header/kpoints/nrk'][...])
wfn_out.kpt = wfn['/mf_header/kpoints/rk'][...].T
wfn_out.ifmin = wfn['/mf_header/kpoints/ifmin'][...].T # dimensions from (ns,nk) to (nk,ns)
wfn_out.ifmax = wfn['/mf_header/kpoints/ifmax'][...].T
wfn_out.energies = wfn['/mf_header/kpoints/el'][...].T # dimensions from (ns,nk,nb) to (nb,nk,ns)
wfn_out.occupations = wfn['/mf_header/kpoints/occ'][...].T # (nb,nk,ns)
wfn_out.nspinor = 1 # No spinors only!
wfn_out.kgrid = wfn['/mf_header/kpoints/kgrid'][...]
wfn_out.kshift = wfn['/mf_header/kpoints/shift'][...]
wfn_out.kw = wfn['/mf_header/kpoints/w'][...]
wfn_out.ngk = wfn['/mf_header/kpoints/ngk'][...]
wfn_out.ngkmax = int(wfn['/mf_header/kpoints/ngkmax'][...])

# gspace
wfn_out.ecutrho = np.float(wfn['/mf_header/gspace/ecutrho'][...])
wfn_out.ng = int(wfn['/mf_header/gspace/ng'][...])
wfn_out.gvec = wfn['/mf_header/gspace/components'][...].T
wfn_out.FFTgrid = wfn['/mf_header/gspace/FFTgrid'][...]


if raw_input('Writing data to %s. Are you sure? [y/N] '%(fname_out))!='y':
    sys.exit(0)

wfn_out.write_header(full=True)

# Writing Data

start_slice = 0
end_slice = 0
for ik in range(wfn_out.nk):
    end_slice += wfn_out.ngk[ik]
    print "kpt: ",ik
    gvec = wfn['/wfns/gvecs'][start_slice:end_slice].T
    wfn_out.write_gvectors(gvec)
    for ib in range(wfn_out.nbands):
        print "band: ",ib
        data = wfn['/wfns/coeffs'][ib,start_slice:end_slice,:].view(np.complex128) # dimensions are (nbands,ngk)
        wfn_out.write_data(data)
    start_slice = end_slice

print wfn_out

