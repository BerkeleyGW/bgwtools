#!/usr/bin/env python

# Another tool to dump espresso binary files created with the iotk library.
# This is almost the same as iotk_parse.py
# Yay to code duplication!
#
# Felipe H. da Jornada, Feb 2015.


import numpy as np
from FortranIO import FortranFile
import re

class IotkFile(FortranFile):
    def __init__(self, fname, endian='=', header_prec='i', *args, **kwargs):
        FortranFile.__init__(self, fname, endian, header_prec, *args, **kwargs)
        self.last_fmt = None
        self.last_tag = None
        self.cur_level = 0
        self.re_obj = re.compile(r'type="([a-zA-Z]+)"')

    def read_xml_field(self):
        header1 = FortranFile.read(self, self.HEADER_PREC)[0]
        control = header1 % 256
        taglen = header1 // 256
        rec = FortranFile.read_record(self)
        header2 = rec.read('i', 1)
        assert(header2//256 == taglen)
        data = rec.read(None)
        i1, i2 = data.index('<'), data.rindex('>')
        assert(len(data)==taglen)
        return control, data[i1:i2+1]

    def write_xml_field(self, control, data_str):
        if control==2:
            self.cur_level -= 1
        data_str = '\n' + '  '*self.cur_level + data_str + '\n'
        if control==1:
            self.cur_level += 1
        header1 = control + 256*len(data_str)
        header2 = 128 + 256*len(data_str)
        FortranFile.write_vals(self, 'i', header1)
        buf = self.pack('i', header2)
        buf += data_str
        FortranFile.writeline(self, buf)

    def read_field(self):
        if self.last_fmt==None:
            obj = self.read_xml_field()
            self.last_tag = obj[1]
            match = self.re_obj.search(obj[1])
            if match is None:
                self.last_fmt = None
            else:
                fmt = match.group(1)
                #self.last_fmt = {'integer':np.int32, 'real':np.float64,
                #    'complex':np.complex128}[fmt]
                self.last_fmt = {'integer':np.int32, 'real':np.float64,
                    'complex':np.complex128}[fmt]
            return 0, obj[0], obj[1]
        else:
            fmt = self.last_fmt
            self.last_fmt = None
            self.last_tag = None
            rec = self.read_record()
            rec.read('i', 1)
            #return 1, fmt, np.frombuffer(rec.read(None), self.ENDIAN+fmt)
            return 1, fmt, np.frombuffer(rec.read(None), fmt)

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
