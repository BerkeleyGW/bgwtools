#!/usr/bin/env python

# Script to read/write charge density files from QE.
#
# Felipe H. da Jornada, Feb 2015.

import numpy as np
from bgwtools.espresso.gvectors import GvectorsIO
from bgwtools.espresso.evcs import EigenvectorsIO
import os
import shutil
import fnmatch, re
try:
    from lxml import etree
except:
    import xml.etree.cElementTree as etree


if __name__=="__main__":
    import sys

    if len(sys.argv)<2:
        print 'usage: {} source_dir.save dest_dir.save'.format(sys.argv[0])
        print
        sys.exit(1)
    src_dir = sys.argv[1]
    dst_dir = sys.argv[2]

    if os.path.isdir(dst_dir):
        raise IOError('Directory {} already exists'.format(dst_dir))
    if os.path.exists(dst_dir):
        raise IOError('Cannot write to {}'.format(dst_dir))
    os.mkdir(dst_dir)
    os.mkdir(dst_dir+'/K00001')
    
    patterns = ('charge-density.dat', '*.upf')
    rule = re.compile('|'.join(('('+fnmatch.translate(p)+')' for p in patterns)), re.IGNORECASE)
    files = sorted([name for name in os.listdir(src_dir) if rule.match(name)] + 
                   ['K00001/eigenval.xml'])
    for fname in files:
        shutil.copyfile(os.path.join(src_dir,fname), os.path.join(dst_dir,fname))

    gvectors = GvectorsIO(src_dir+'/K00001/gkvectors.dat')
    gvectors_k = gvectors.unfold_trs()
    gvectors_k.write(dst_dir+'/K00001/gkvectors.dat')

    evecs = EigenvectorsIO(src_dir+'/K00001/evc.dat')
    assert gvectors.ngk==evecs.ngw
    evecs_k = evecs.unfold_trs()
    assert gvectors_k.ngk==evecs_k.ngw
    evecs_k.write(dst_dir+'/K00001/evc.dat')

    tree = etree.parse(src_dir+'/data-file.xml')
    def xml_replace(tag, orig, new):
        el = tree.find(tag)
        assert el.text.strip()==str(orig)
        spaces = el.text.split('\n')[-1]
        if isinstance(new, int):
            el.text = '\n{:10}\n{}'.format(new, spaces)
        else:
            el.text = '\n{}\n{}'.format(new, spaces)
    def xml_assert(tag, value):
        el = tree.find(tag)
        assert el.text.strip()==str(value)
    def xml_get(tag):
        el = tree.find(tag)
        return el.text.strip()
    xml_replace('PLANE_WAVES/MAX_NUMBER_OF_GK-VECTORS', gvectors.ngk_max, 2*gvectors.ngk_max - 1)
    xml_replace('PLANE_WAVES/GAMMA_ONLY', 'T', 'F')
    ng_rho = int(xml_get('PLANE_WAVES/GVECT_NUMBER'))
    xml_replace('PLANE_WAVES/GVECT_NUMBER', ng_rho, 2*ng_rho - 1)
    xml_replace('PLANE_WAVES/SMOOTH_GVECT_NUMBER', ng_rho, 2*ng_rho - 1)
    xml_replace('EIGENVECTORS/MAX_NUMBER_OF_GK-VECTORS', gvectors.ngk_max, 2*gvectors.ngk_max - 1)
    xml_replace('EIGENVECTORS/K-POINT.1/NUMBER_OF_GK-VECTORS', gvectors.ngk, 2*gvectors.ngk - 1)
    el = tree.find('BRILLOUIN_ZONE/MONKHORST_PACK_GRID')
    el.attrib['nk1'] = el.attrib['nk2'] = el.attrib['nk3'] = '1'

    '''
    el = tree.find('BRILLOUIN_ZONE')
    builder = etree.TreeBuilder()
    builder.start('STARTING_K-POINTS', dict(type='integer', size='1'))
    builder.data('\n{:10}\n    ')
    subel = builder.end('STARTING_K-POINTS')
    el.append(subel)

    builder = etree.TreeBuilder()
    builder.start('K-POINT_START.1', dict(
                  XYZ="0.000000000000000E+000 0.000000000000000E+000 0.000000000000000E+000",
                  WEIGHT="1.000000000000000E+000"))
    subel = builder.end('K-POINT_START.1')
    el.append(subel)
    '''

    tree.write(dst_dir+'/data-file.xml')
