#!/usr/bin/env python

import h5py
import numpy as np
from bgwtools.IO.wfn import wfnIO


def main(fname, fname_out):
    f_out = open(fname_out, 'w')

    f = h5py.File(fname)
    nk = f['mf_header/kpoints/nrk'][()]
    nb = f['mf_header/kpoints/mnband'][()]
    ngk = f['mf_header/kpoints/ngk'][()]
    f_out.write('{}\n{}\n'.format(nk,nb))
    print(nk)
    print(nb)
    order = np.arange(nb)

    gvecs = f['wfns/gvecs'][()]
    dims = [nb]+[np.amax(gvecs[:,idim]) - np.amin(gvecs[:,idim]) + 1 for idim in range(3)]
    box_last = np.zeros(dims, dtype=np.complex128)
    box = np.zeros(dims, dtype=np.complex128)

    def get_coeffs(box, offsetk, ng, order):
        gvec = gvecs[offsetk:offsetk+ng]
        cg = f['wfns/coeffs'][:,0,offsetk:offsetk+ng,:].view(np.complex128)[...,0]
        box[...] = 0
        for ig in range(ng):
            box[:,gvec[ig,0],gvec[ig,1],gvec[ig,2]] = cg[order,ig]

    order = np.arange(nb)
    tol = np.sqrt(0.5)
    ng = ngk[0]
    get_coeffs(box_last, 0, ng, order)
    box_last = box_last.conj()
    offsetk = ng
    print('ik = {}'.format(0))
    f_out.write(' '.join((order+1).astype(str)))
    f_out.write('\n')

    for ik in xrange(1,nk):
        print('ik = {}'.format(ik))
        ng = ngk[ik]
        get_coeffs(box, offsetk, ng, order)
        print('  - done reading wfns')

        olap = (box_last * box).sum(-1).sum(-1).sum(-1)
        olap = np.abs(olap)
        cond = olap < tol
        ind = np.where(cond)[0]
        ind_l = []
        ind_r = []

        submat = np.abs(np.dot(box_last.reshape(nb,-1)[cond], box.reshape(nb,-1)[cond].T))
        print('  - submatrix size: {}'.format(submat.shape))
        for i in range(len(submat)):
            coord = np.unravel_index(np.argmax(submat), submat.shape)
            ind_l.append(coord[0])
            ind_r.append(coord[1])
            submat[coord[0],:] = 0
            submat[:,coord[1]] = 0

        ind_l = ind[ind_l]
        ind_r = ind[ind_r]

        order_ = np.arange(nb)
        order_[ind_l]  = np.arange(nb)[ind_r]
        order[ind_l] = order[ind_r]
  
        f_out.write(' '.join((order+1).astype(str)))
        f_out.write('\n')
        box_last[...] = box.conj()[order_]
        offsetk += ng

    f_out.close()

    '''
    wfn = wfnIO(fname)
    gvecs = wfn.get_gvectors_bufer()
    data = wfn.get_data_buffer()

    for ik in xrange(wfn.nk):
        wfn.read_gvectors(gvecs)
        cond = np.all(gvecs[:,:wfn.ngk[ik]]==0, axis=0)
        idx = np.nonzero(cond)[0][0]
        print ik, idx, gvecs[:,idx]
        for ib in xrange(wfn.nbands):
            wfn.read_data(data)
            heads[ik,ib] = data[idx,0]

        break
    '''


if __name__=="__main__":
    import argparse

    desc = 'A utility to reorder the bands from a WFN file.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('wfn', help='a WFN.h5 file')
    parser.add_argument('output', help='output file with bands in the correct order')
    args = parser.parse_args()

    main(args.wfn, args.output)
