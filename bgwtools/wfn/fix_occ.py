#!/usr/bin/env python

# An utility to fix the occupations in a WFN file
# Felipe H. da Jornada (Feb 2013)

import numpy as np
import sys
if len(sys.argv)!=3:
    print('Usage: %s wfn_in wfn_out'%(sys.argv[0]))
    sys.exit(1)
fname_in = sys.argv[1]
fname_out = sys.argv[2]

#copy *all* wfn file
import shutil
shutil.copyfile(fname_in, fname_out)

#open wfn file, write fixed header to file fname_tmp
from bgwtools.IO.wfn import wfnIO
fname_tmp = '_wfn_tmp_'
wfn = wfnIO(fname_in, full=False)
nocc = 1022
wfn.occupations = np.zeros_like(wfn.occupations, order='F')
wfn.occupations[:nocc] = 1.0
wfn.ifmin = np.zeros_like(wfn.ifmin, order='F')
wfn.ifmin[0,0] = 1
wfn.ifmax = np.zeros_like(wfn.ifmax, order='F')
wfn.ifmax[0,0] = nocc
#exit()

wfn.to_file(fname_tmp, full=False)

#copy fixed header to fname_out
f_tmp=open(fname_tmp, 'rb')
f_out=open(fname_out, 'r+b')
f_out.write(f_tmp.read())
f_out.close()
f_tmp.close()

#clean-up
import os
os.unlink(fname_tmp)

print('All done')
