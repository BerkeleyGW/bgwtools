#!/usr/bin/env python

# This script opens an arbitrary number of epsmat files and fits epsinv to a
# spline. The parameters are:
#  Gz_max: maximum value of |Gz|
#  avgcut_xy: maximum value of |q_xy|^2
#
# Felipe H. da Jornada

import numpy as np
from scipy.spatial import cKDTree
from scipy.interpolate import splrep, splev

# Not exactly the BZ, but good enough
def _min_range_abs(x):
    '''Moves vector x to the [-0.5, 0.5) range, and apply abs().'''
    return np.fabs(x - np.floor(x + 0.5))

class EpsmatModeler:
    def __init__(self, wfn, Gz_max, avgcut_xy):
        self.bdot = wfn.bdot.copy()
        self.Lz = wfn.alat * wfn.avec[2,2]
        self.M = np.linalg.cholesky(self.bdot).T
        self.Gz_max = Gz_max
        # We are storing eps(Gz,Gz), where Gz is stored in this order:
        self.Gzs = np.arange(-self.Gz_max, self.Gz_max+1)
        self.nG = len(self.Gzs)
        # These are the gvecs we want from the epsmat files, in the appropriate order:
        self.gvecs_want = np.zeros((3,self.nG), dtype='int')
        self.gvecs_want[2,:] = self.Gzs
        self.avgcut_xy = avgcut_xy
        self.qs = []
        self.qlens = []
        self.eps = []
        self.nq = 0
        self.eps00 = []

    def add_epsmat(self, epsmat):
        print 'Dealing with epsmat file "%s"'%(epsmat.fname)
        #construct a KDTree to quickly select the G-vectors we care about
        tree = cKDTree(epsmat.gvec_k.T)
        dists, igs_glob = tree.query(self.gvecs_want.T)
        assert(np.all(dists==0))

        for iq in range(epsmat.nq):
            qq = epsmat.qpt[:,iq]
            qq = _min_range_abs(qq)
            qlen = np.sqrt(np.sum(np.dot(self.M, qq)**2))
            print '(%d/%d) |q| = %.08f Bohr^-1'%(iq+1, epsmat.nq, qlen),
            if qlen>self.avgcut_xy:
                print '[ignored]'
                epsmat.read_qpt(in_place=False, ignore=True)
                continue

            self.qs.append(qq)
            self.qlens.append(qlen)

            buf = epsmat.read_qpt(in_place=False)
            igs_local = epsmat.isort_i[iq][igs_glob] - 1
            eps = buf[igs_local, igs_local]
            assert(np.all(np.fabs(np.imag(eps))<1e-15))
            eps = np.real(eps)

            self.eps.append(eps)
            print '[OK]'
        print 'Done'
        print

    def commit_data(self):
        self.qs = np.array(self.qs).T
        self.qlens = np.array(self.qlens)
        self.eps = np.array(self.eps).T
        order = np.argsort(self.qlens)
        self.qs = self.qs[:,order]
        self.qlens = self.qlens[order]
        self.eps = self.eps[:,order]

    def get_vcoul(self, qlen, Gz, trunc=True):
        if isinstance(Gz, np.ndarray):
            Gz = Gz[:, np.newaxis]
            qlen = self.qlens[np.newaxis, :]
        G2 = Gz**2 * self.bdot[2,2]
        cos_term = 1 - 2*(Gz%2)
        zc = self.Lz/2.0

        # (nq, nG) array with truncated Coulomb interaction
        # otherwise, (nq) array
        if trunc:
            return 8.*np.pi/(qlen**2 + G2) * ( 1. - np.exp(-qlen*zc)*cos_term )
        else:
            return 8.*np.pi/(qlen**2 + G2)
        

    def model(self, model=0):
        import matplotlib.pyplot as plt

        if model==0:
            ys = self.eps
        else:
            vT = self.get_vcoul(self.qlens, self.Gzs)
            chi = (self.eps - 1.0)/vT
            #ys = np.array(chi, dtype='float')
            ys = chi

        #for ig, Gz in zip([self.Gz_max], [0]):
        for ig, Gz in zip(np.arange(self.nG), self.Gzs):
            #print ig, Gz
            x = self.qlens
            y = ys[ig]

            if (Gz==0):
                # Let's force epsinv(0) = 1 => chi(0) = 0
                x = np.append(0, self.qlens)
                y = np.append(0, y)
                #x = np.append(0, np.append(self.qlens[0], np.append(self.qlens[0]*2, self.qlens[1:])))
                #y = np.append(0, np.append(y[0], np.append(y[0]*4, y[1:])))
                #x = np.append(0, np.append(self.qlens[0], np.append(self.qlens[0]*2, self.qlens[1:])))
                #y = np.append(0, np.append(y[0], np.append(y[0]*4, y[1:])))

            #For Gz>=2, we use a smoothing spline interpolation (i.e., s>0)
            #This would not be necessary if Ms. Qiu calculated all epsmat files
            #with the same cutoff ;)
            if np.fabs(Gz)<1:
                s=0
            else:
                s=len(x)
            tck = splrep(x, y, s=s)
            x_intp = np.linspace(0, np.amax(x), 10000) + 1e-12
            y_intp = splev(x_intp, tck)

            if model>0:
                vT_intp = self.get_vcoul(x_intp, Gz)
                y_intp = 1.0 + vT_intp*y_intp

            if Gz<0: continue # only need to plot the second part

            #lines = plt.plot(self.qlens, self.eps[ig], 'o')
            #lines = plt.plot(x, y, 'o')
            #color = lines[0].get_color()
            #plt.plot(x_intp, y_intp, '-', lw=2, color=color, label='Gz=%d'%(Gz))
            plt.plot(x_intp, y_intp, '-', lw=2, label='Gz=%d'%(Gz))
        plt.legend(prop={'size':14})
        plt.xlabel('$|q|$', size=18)
        plt.ylabel(r'$\varepsilon^{-1}$', size=18)
        plt.show()
    
    def get_bgw_params(self):
        print self.Gzs

if __name__=="__main__":
    from bgwtools.IO.wfn import wfnIO
    from bgwtools.IO.epsmat import epsmatIO
    import cPickle

    from optparse import OptionParser
    usage = ("\n" +
             "(1) %prog [options] wfn epsmat1 [epsmat2 ...]\n" + 
             "(2) %prog model.pck")
    parser = OptionParser(usage)
    parser.add_option("--Gz_max", default=0, type="int",
                      help="keep G vectors up to Gz<=|Gz_max|")
    parser.add_option("--avgcut_xy", default=0.0, type="float", 
                      help="""keep q-points up |q|^2<=|avgcut_xy|.
                      If 0.0, keep all q-points.""")
    parser.add_option("-k", default=3, type="int", 
                      help="order of the spline polynomial.")
    parser.add_option("-d", "--dump", default='model.pck', metavar="FILE", type="str", 
                      help="dumps model to FILE.")
    (options, args) = parser.parse_args()

    if len(args)<1:
        parser.error('Insuficient number of arguments.')
    elif len(args)==1:
        print "Reading model from file '%s'"%(args[0])
        f = open(args[0], 'rb')
        epsmat_modeler = cPickle.load(f)
        f.close()
    else:
        print "Parsing wfn file '%s'"%(args[0])
        wfn = wfnIO(args[0])
        epsmat_modeler = EpsmatModeler(wfn, options.Gz_max, options.avgcut_xy)
        for arg in args[1:]:
            print "Parsing epsmat file '%s'"%(arg)
            epsmat = epsmatIO(arg, read_all=False)
            epsmat_modeler.add_epsmat(epsmat)
        epsmat_modeler.commit_data()

    epsmat_modeler.model(model=1)
    if len(args)>1:
        print "Dumping model to file '%s'"%(options.dump)
        f = open(options.dump, 'wb')
        cPickle.dump(epsmat_modeler, f)
        f.close()
