#!/usr/bin/env python

# An utility to change the symmetries in a WFN file
# Felipe H. da Jornada (Feb 2013)

import numpy as np
import sys
if len(sys.argv)!=4:
    print('Usage: %s wfn_in wfn_ref wfn_out'%(sys.argv[0]))
    sys.exit(1)
fname_in = sys.argv[1]
fname_ref = sys.argv[2]
fname_out = sys.argv[3]

#open wfn file, write fixed header to file fname_tmp
from bgwtools.IO.wfn import wfnIO
wfn_in = wfnIO(fname_in, full=False)
wfn_ref = wfnIO(fname_ref, full=False)
print 'Original symmetries:'
print wfn_in.ntran
print wfn_in.mtrx
print wfn_in.tnp
print 'New symmetries:'
wfn_in.ntran = wfn_ref.ntran
wfn_in.mtrx = wfn_ref.mtrx
wfn_in.tnp = wfn_ref.tnp
print wfn_in.ntran
print wfn_in.mtrx
print wfn_in.tnp
pos_in = file.tell(wfn_in.f)
pos_out = wfn_in.to_file(fname_out, full=False)
print pos_in
print pos_out

#copy fixed header to fname_out
buf_sz = 128*1024*1024
f_in=open(fname_in, 'rb')
f_in.seek(pos_in)
f_out=open(fname_out, 'r+b')
f_out.seek(pos_out)
while True:
    buf = f_in.read(buf_sz)
    if buf=='':
        break
    f_out.write(buf)
