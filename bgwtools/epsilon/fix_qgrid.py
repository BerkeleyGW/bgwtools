#!/usr/bin/env python

# A small utility to fix the q-grid from an epsmat file
# Felipe H. da Jornada (Dec 2012)

import numpy as np
import sys
if len(sys.argv)!=4:
    print('Usage: %s epsmat_input epsmat_output kx,ky,kz'%(sys.argv[0]))
    sys.exit(1)
fname_in = sys.argv[1]
fname_out = sys.argv[2]
grid = sys.argv[3].split(',')
if len(grid)!=3:
    print('Invalid q-grid entered: %s'%(sys.argv[3]))
    sys.exit(1)
grid = np.array(map(int, grid))

#copy the epsmat file, including q-points
import shutil
shutil.copyfile(fname_in, fname_out)

#open epsmat file, write fixed header to file fname_tmp1
from bgwtools.IO.epsmat import epsmatIO
fname_tmp = '_epsmat_tmp_'
epsmat = epsmatIO(fname_in, read_all=False)
print 'Original q-grid:', epsmat.grid
epsmat.grid = grid
print 'New q-grid:', epsmat.grid
epsmat.to_file(fname_tmp, write_all=False)

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
