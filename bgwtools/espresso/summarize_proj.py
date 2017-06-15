#!/usr/bin/env python

# This script will accumulate the projections over all atoms of the specified
# type. The output file contains a matrix of dimensions (nk,nb), where each
# value is \sum_atom \sum_p(atom) |<p(atom)|ib,ik>|^2
#
# Make sure you change the parameter nignore and atoms below!
#
# Felipe H. da Jornada (2015)

#Number of lines to ignore. The last line to be ignored should contain the
#last atom in the file.
nignore = 14

#Atomic species over which we accumulate the projection.
atoms = ('Ta', 'S')

import numpy as np

def main(fname, fname_out):
    f = open(fname)

    # Ignore first "nignore" lines in the file
    for i in range(nignore):
        f.readline()
    nproj, nk, nb = map(int, f.readline().split())
    f.readline()

    atoms_ = [atom.strip().lower() for atom in atoms]

    projs = np.empty((nproj,nk,nb))
    should_acc = np.zeros((nproj,), dtype=bool)
    for iproj in range(nproj):
        header = f.readline()
        iproj_ = int(header[0:5])
        ia = int(header[5:10])
        at = header[10:13]
        n = int(header[13:18])
        l = int(header[18:23])
        m = int(header[23:28])
        should_acc[iproj] = at.strip().lower() in atoms_
        print iproj_, ia, at, n, l, m

        data = []
        for il in range(nk*nb):
            data.append(float(f.readline()[16:-1]))
        projs[iproj,:,:] = np.array(data).reshape((nk,nb))

    print
    projs = projs[should_acc,:,:].sum(axis=0)
    np.savetxt(fname_out, projs)


if __name__=="__main__":
    import sys
    fname_in, fname_out = sys.argv[1:3]
    main(fname_in, fname_out)
