#!/usr/bin/env python

# Reads a QE save directory, extracts useful information, and create a file
# suitable for band structure plotting.
#
# Felipe H. da Jornada (2015)

import numpy as np
import sys
import os
import xml.etree.cElementTree as ET

units = {'hartree':27.21139, 'rydberg':13.605695, 'ev':1.0}
dirname = sys.argv[1]
fname_out = sys.argv[2]
tree = ET.parse(os.path.join(dirname, 'data-file.xml'))

EF = float(tree.find('BAND_STRUCTURE_INFO/FERMI_ENERGY').text)*units['hartree']
nk = int(tree.find('BAND_STRUCTURE_INFO/NUMBER_OF_K-POINTS').text)
nb = int(tree.find('BAND_STRUCTURE_INFO/NUMBER_OF_BANDS').text)
kpts = np.empty((nk,3), dtype='float')
evals = np.empty((nk,nb), dtype='float')
print('EF=%.4f eV'%(EF))
print('nk=%d'%(nk))
print('nb=%d'%(nb))

evals_root = tree.find('EIGENVALUES')
for el in tree.find('BRILLOUIN_ZONE'):
	items = el.tag.split('.')
	if items[0] != 'K-POINT': continue
	ik = int(items[1])-1
	kpts[ik,:] = np.array(map(float, el.attrib['XYZ'].split()))

	evals_fname = evals_root.find('%s/DATAFILE'%(el.tag)).attrib['iotk_link']
	evals_fname = os.path.join(dirname, evals_fname)
	evals_tree = ET.parse(evals_fname)
	evals_k = map(float, evals_tree.find('EIGENVALUES').text.split())
	try:
		unit_text = evals_tree.find('UNITS_FOR_ENERGIES').attrib['UNITS'].lower()
		fact = units[unit_text]
	except:
		fact = 1.0
	del evals_tree
	evals[ik,:] = np.array(evals_k)*fact

data = np.column_stack((kpts, evals))
f_out = open(fname_out, 'w')
f_out.write('#EF={:.4f}\n#nb={}\n#nk={}\n'.format(EF,nb,nk))
np.savetxt(f_out, data, fmt='%.12f')
f_out.close()
