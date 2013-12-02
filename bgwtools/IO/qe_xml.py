#!/usr/bin/env python

import numpy as np
from FortranIO import FortranFile
import re

class IotkFile(FortranFile):
    def __init__(self, fname, endian='=', header_prec='i', *args, **kwargs):
        FortranFile.__init__(self, fname, endian, header_prec, *args, **kwargs)
        self.next_format = ''
        self.re_obj = re.compile(r'type="([a-zA-Z]+)"')

    def read_xml_field(self):
        header1 = FortranFile.read(self, self.HEADER_PREC)[0]
        control = header1 % 256
        taglen = header1 // 256
        rec = FortranFile.read_record(self)
        header2 = rec.read('i', 1)
        assert(header2//256 == taglen)
        data = rec.read(None)
        assert(len(data)==taglen)
        return control, data

    def write_xml_field(self, control, data_str):
        header1 = control + 256*len(data_str)
        header2 = 128 + 256*len(data_str)
        FortranFile.write_vals(self, 'i', header1)
        buf = self.pack('i', header2)
        buf += data_str
        FortranFile.writeline(self, buf)

    def read_field(self):
        if self.next_format=='':
            obj = self.read_xml_field()
            match = self.re_obj.search(obj[1])
            if match is None:
                self.next_format = ''
            else:
                fmt = match.group(1)
                self.next_format = {'integer':'i', 'real':'d'}[fmt]
            return 0, obj[0], obj[1]
        else:
            fmt = self.next_format
            self.next_format = ''
            rec = self.read_record()
            rec.read('i', 1)
            #return 1, fmt, rec.read(fmt)
            return 1, fmt, np.frombuffer(rec.read(None), self.ENDIAN+fmt)

    def write_field(self, obj):
        assert(len(obj)==3)
        if obj[0]==0:
            self.write_xml_field(obj[1], obj[2])
        else:
            self.write_vals('i' + obj[1]*len(obj[2]), 0, *obj[2].tolist())


if __name__ == '__main__':
    import sys

    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    f_in = IotkFile(fname_in, endian='<', mode='rb')
    f_out = IotkFile(fname_out, endian='>', mode='wb')
    while True:
        try:
            obj = f_in.read_field()
        except IOError:
            break
        f_out.write_field(obj)
