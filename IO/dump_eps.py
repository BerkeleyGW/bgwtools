#!/usr/bin/env python

import cPickle
from epsmat import epsmatIO
import epsmat_merge
import sys
import time

def tprint(str_):
	print time.strftime("[%H:%M:%S]", time.localtime()), str_

tprint('Reading and merging files')
epsmat = epsmat_merge.merge(sys.argv[1:])
tprint('Dumping /tmp/epsmat.pkl')
f=open('/tmp/epsmat.pkl','wb')
epsmat.f=None
cPickle.dump(epsmat, f, cPickle.HIGHEST_PROTOCOL)
f.close()
tprint('Done')
