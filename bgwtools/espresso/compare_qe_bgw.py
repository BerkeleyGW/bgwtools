#!/usr/bin/env python

# Script to read/write charge density files from QE.
#
# Felipe H. da Jornada, Feb 2015.

import numpy as np
from bgwtools.espresso.gvectors import GvectorsIO
from bgwtools.espresso.evcs import EigenvectorsIO
import h5py

if __name__=="__main__":
    import sys

    if len(sys.argv)<2:
        print 'usage: {} source_dir.save WFN.h5'.format(sys.argv[0])
        print
        sys.exit(1)
    src_dir = sys.argv[1]
    f_wfn = h5py.File(sys.argv[2])

    gvectors = GvectorsIO(src_dir+'/K00001/gkvectors.dat')
    print gvectors.gvecs.T
    print f_wfn['wfns/gvecs'][()]
    print f_wfn['mf_header/gspace/components'][()]
    #evecs = EigenvectorsIO(src_dir+'/K00001/evc.dat')

