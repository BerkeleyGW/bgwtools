#!/usr/bin/env python

# A small tool to dump espresso binary files created with the iotk library.
#
# Felipe H. da Jornada, Feb 2015.

import numpy as np
#from FortranIO import FortranFile, cast_data, cast_array
from bgwtools.IO.FortranIO import FortranFile, cast_data, cast_array
import re


class IOTKFile(FortranFile):

    def __init__(self, fname):
        FortranFile.__init__(self, fname)

    def read_header(self):
        h1 = FortranFile.read(self, 'i')[0]
        c1, l1 = h1%256, h1//256
        assert c1 in (1,2,3,4,5)
        return c1, l1

    def read_tag(self, l1):
        s = FortranFile.read(self)
        h2 = cast_data(s[0:4], 'i')
        c2, l2 = h2%256, h2//256
        assert c2==128
        assert l1==l2
        s = s[4:]
        i1, i2 = s.index('<'), s.rindex('>')
        return s[i1:i2+1]

    def read(self, dtype=None):
        c1, l1 = self.read_header()
        return self.read_tag(l1)

    def read_data(self, dtype):
        s = FortranFile.read(self)
        s = s[4:]
        return np.frombuffer(s, dtype)


def dump_file(fname):
    re_tag = re.compile(r'''([^\s]*)=['"]([^'"]+)['"]''')
    re_name = re.compile(r'''<([^\s]*)''')
    type_dict = {'integer':np.int32, 'real':np.float64, 'complex':np.complex128, 'logical':np.int32}
    f = IOTKFile(fname)

    while True:
        # Read all tags
        try:
            tag = f.read()
        except:
            break
        print tag
        name = re_name.search(tag).group(1)
        attrs = dict(re_tag.findall(tag))
        if len(attrs):
            print attrs

        # If the tag containts a "tag" argument, there will be a child data
        # node. In principle, you might want to skip data fields from tags
        # you don`t care about..
        if 'kind' in attrs:
            # As an example, we`ll skip the datasets with tags matching
            # evc*, unless the name is evc.6 (the 6th band)
            should_skip = name.startswith('evc') and name!='evc.6'
            if should_skip:
                print '(skipping data field)'
                f.next()
            else:
                kind = attrs['kind']
                type = attrs['type']
                size = int(attrs['size'])
                try:
                    cols = int(attrs['columns'])
                except:
                    cols = None
                dtype = type_dict[type]
                data = f.read_data(dtype)
                if cols is not None:
                    data = data.reshape((-1,cols))
                # The array data contains the resulting data structure
                print data
                print '  norm(data) =', np.linalg.norm(data)


if __name__=="__main__":
    import sys

    if len(sys.argv)<2:
        print 'usage: %0 file'
        print
        print 'where "file" is a QE binary xml file written with the iotk library'
        sys.exit(1)
    fname = sys.argv[1]
    dump_file(fname)
