#!/usr/bin/env python

# Calculates the G=0 component of the WFN

# Felipe Homrich da Jornada (Aug 2013)

import numpy as np
from bgwtools.IO.wfn import wfnIO
from bgwtools.common.common import get_numpy_flavor


class GaugeCalculator:

    def __init__(self, fname):
        self.wfn = wfnIO(fname)

    def get_gauges(self, bands, verbose=True):
        wfn = self.wfn

        nbands = len(bands)
        ng_max = np.amax(wfn.ngk)
        gvec = np.empty((3, ng_max), dtype=int, order='F')
        data = np.empty((ng_max, wfn.ns), dtype=get_numpy_flavor(wfn.flavor), order='F')
        gauges = np.empty((nbands,wfn.nk), dtype=get_numpy_flavor(wfn.flavor))
        k_done = 0
        for ik in range(wfn.nk):
            k_done += 1
            if verbose:
                print k_done,'/',wfn.nk
            #wfn.read_gvectors(gvec)
            #Assume G-vectors are order by kinetic energy
            wfn.read_gvectors()
            ib_ = 0
            for ib in range(wfn.nbands):
                if ib in bands:
                    wfn.read_data(data)
                    gauges[ib_,ik] = data[0]
                    ib_ += 1
                else:
                    wfn.read_data()
        return gauges

if __name__ == '__main__':
    import sys

    fname = sys.argv[1]
    gauge_calc = GaugeCalculator(fname)
    gauges = gauge_calc.get_gauges([3,4])
    print gauges
    np.save('gauges', gauges)

