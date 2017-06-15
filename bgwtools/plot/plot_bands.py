#!/usr/bin/env python
#coding:utf-8

# Script to plot bandstructure with color proportional to the character of a
# particular state. The plot will be saved into the file "bs.png."
# Don`t forget to adjust the energy range!
#
# Felipe H. da Jornada (2015)

import numpy as np

try:
    import mpl_fonts
    import seaborn as sns
    sns.set(font='Helvetica World', style='ticks', context='talk')
    plt = sns.plt
    mpl = sns.mpl
    cm = sns.cubehelix_palette(8, start=2.8, rot=0, dark=0.2, light=1, as_cmap=True)
except:
    sys.stdout.write('WARNING: could not load seaborn. Falling back to vanilla matplotlib\n')
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    cm = plt.cm.Blues

from matplotlib.artist import allow_rasterization
from matplotlib import rcParams

# We`ll hijack the "draw" method from matplotlib.collections.LineCollection
# to allow us to specify linecaps and linejoins, otherwise the plot looks ugly.

@allow_rasterization
def my_draw(self, renderer):
    if not self.get_visible():
        return
    renderer.open_group(self.__class__.__name__, self.get_gid())

    self.update_scalarmappable()

    transform, transOffset, offsets, paths = self._prepare_points()

    gc = renderer.new_gc()
    self._set_gc_clip(gc)
    gc.set_snap(self.get_snap())

    if self._hatch:
        gc.set_hatch(self._hatch)

    if self.get_sketch_params() is not None:
        gc.set_sketch_params(*self.get_sketch_params())

    if self.get_path_effects():
        from matplotlib.patheffects import PathEffectRenderer
        renderer = PathEffectRenderer(self.get_path_effects(), renderer)

    # If the collection is made up of a single shape/color/stroke,
    # it can be rendered once and blitted multiple times, using
    # `draw_markers` rather than `draw_path_collection`.  This is
    # *much* faster for Agg, and results in smaller file sizes in
    # PDF/SVG/PS.

    trans = self.get_transforms()
    facecolors = self.get_facecolor()
    edgecolors = self.get_edgecolor()
    do_single_path_optimization = False
    if (len(paths) == 1 and len(trans) <= 1 and
        len(facecolors) == 1 and len(edgecolors) == 1 and
        len(self._linewidths) == 1 and
        self._linestyles == [(None, None)] and
        len(self._antialiaseds) == 1 and len(self._urls) == 1 and
        self.get_hatch() is None):
        if len(trans):
            combined_transform = (transforms.Affine2D(trans[0]) +
                                  transform)
        else:
            combined_transform = transform
        extents = paths[0].get_extents(combined_transform)
        width, height = renderer.get_canvas_width_height()
        if (extents.width < width and
            extents.height < height):
            do_single_path_optimization = True

    #HACK
    join = rcParams['lines.dash_joinstyle']
    cap = rcParams['lines.solid_capstyle']
    gc.set_joinstyle(join)
    gc.set_capstyle(cap)
    #/HACK

    if do_single_path_optimization:
        gc.set_foreground(tuple(edgecolors[0]))
        gc.set_linewidth(self._linewidths[0])
        gc.set_linestyle(self._linestyles[0])
        gc.set_antialiased(self._antialiaseds[0])
        gc.set_url(self._urls[0])
        renderer.draw_markers(
            gc, paths[0], combined_transform.frozen(),
            mpath.Path(offsets), transOffset, tuple(facecolors[0]))
    else:
        renderer.draw_path_collection(
            gc, transform.frozen(), paths,
            self.get_transforms(), offsets, transOffset,
            self.get_facecolor(), self.get_edgecolor(),
            self._linewidths, self._linestyles,
            self._antialiaseds, self._urls,
            self._offset_position)

    gc.restore()
    renderer.close_group(self.__class__.__name__)

#Hijack matplotlib, inject new "draw" method
mpl.collections.LineCollection.draw = my_draw

LineCollection = mpl.collections.LineCollection

def plot(fname, fname_order, fname_proj, Emin, Emax):
    with open(fname) as f:
        tmp = f.readline()
        # Assume the file will contain the line #EF=value
        if 'EF' in tmp:
            E0 = float(tmp.split('=')[-1])
        else:
            E0 = 0.0
            print 'WARNING: could not read Fermi energy from file.'
        data = np.loadtxt(f)

    print 'Setting Fermi energy to {:.3f}. Make sure this is correct!'.format(E0)

    # x axis (kk) is actually the length kk(i) = \sum_{j<i} |delta k|_j
    k_cart = data[:,:3]
    dk = np.zeros(data.shape[0])
    dk[0] = 0
    dk[1:,] = np.sqrt(np.sum((k_cart[1:]-k_cart[0:-1])**2, axis=1))
    kk = np.cumsum(dk)

    bands = data[:,3:] - E0
    nk, nb = bands.shape

    # Open file with projection?
    if fname_proj is not None:
        data_proj = np.loadtxt(fname_proj)

    # Open file with band ordering?
    if fname_order is not None:
        with open(fname_order) as f:
            f.readline()
            f.readline()
            order = np.loadtxt(f).astype(int) - 1
    else:
        order = np.empty((nk,nb), dtype=int)
        order[:,:] = np.arange(nb)[None,:]

    # Reorder bands
    for ik in range(nk):
        bands[ik] = bands[ik,order[ik]]

    # Reorder projections, is applicable
    if fname_proj is not None:
        for ik in range(nk):
            data_proj[ik] = data_proj[ik,order[ik]]

    # Plot! We use a fancy scheme based on the hijacked LineCollection
    # to draw very efficiently the BS with colors.
    if fname_proj is not None:
        norm = mpl.colors.NoNorm()
        proj = (data_proj[1:] + data_proj[:-1])*0.5
        nbins = 32
        bins = np.linspace(0, 1, nbins+1)
        ind = np.digitize(proj.ravel(), bins).reshape(proj.shape) - 1
        for ibin in range(nbins):
            cond = ind==ibin
            npts = np.sum(cond)
            print ibin, npts
            if npts==0:
                continue
            ind_k, ind_b = np.where(cond)
            x0 = kk[ind_k]
            x1 = kk[ind_k+1]
            y0 = bands[ind_k,ind_b]
            y1 = bands[ind_k+1,ind_b]
            segments = []
            projs2 = []
            for ipt in range(npts):
                segments.append( ((x0[ipt],y0[ipt]), (x1[ipt],y1[ipt])) )

            lc = LineCollection(segments, cmap=cm, norm=norm)
            proj2 = proj[ind_k,ind_b]
            lc.set_array(proj2*0.9 + 0.1)
            zorder = (proj2.mean()*100 - 100).astype(int)
            lc.set_zorder(zorder)
            lc.set_lw(2.5)
            #lc.set_solid_capstyle('round')
            # Fancy scheme if the user has the projections
            plt.gca().add_collection(lc)

    else:
        for ib in range(nb):
            if fname_proj is not None:
                #We`ll never execute this branch, it`s here for debugging purposes
                if 0:
                    proj = data_proj[:,ib]
                    proj = (proj[1:] + proj[:-1])*0.5
                    color = cm(proj*0.9 + 0.1)
                    zorder = (proj*100 - 100).astype(int)
                    alpha = proj*0.5 + 0.5
                    for ik in range(nk-1):
                        plt.plot(kk[ik:ik+2], bands[ik:ik+2, ib], color=color[ik],
                            zorder=zorder[ik], alpha=alpha[ik])
                else:
                    proj = data_proj[:,ib].mean()
                    color = cm(proj*0.9 + 0.1)
                    zorder = (proj*100 - 100).astype(int)
                    alpha = proj*0.5 + 0.5
                    plt.plot(kk, bands[:, ib], color=color,
                        zorder=zorder, alpha=alpha)
            else:
                # Simple plot if the user doesn`t have the projections
                plt.plot(kk, bands[:, ib], '-', lw=2.5, color=cm(1.))

    opts = dict(color='#aaaaaa', ls='solid', lw=1, zorder=5)
    plt.axvline(kk[50], **opts)
    plt.axvline(kk[100], **opts)
    plt.xticks([kk[0], kk[50], kk[100], kk[150]], [u'Γ', u'Μ', u'Κ', u'Γ'])
    plt.xlim(kk[0], kk[-1])
    plt.ylim(Emin, Emax)
    plt.ylabel('Energy (eV)')
    plt.axhline(ls='-', color='#aaaaaa', lw=1)


if __name__ == "__main__":
    from argparse import ArgumentParser

    desc = 'Plots band structure, optionally with band reorder and colored by projection.'
    parser = ArgumentParser(description=desc)
    parser.add_argument('fname', help=(
        'Band structure file gernerated by espresso2bs.py (bs.dat)'))
    parser.add_argument('--fname_proj',
        help='Projection file generated by summarize_proj.py (bs_proj.dat)')
    parser.add_argument('--fname_order',
        help='Band ordering file generated by path_reorder.py')
    parser.add_argument('--Emin', default=-7, type=float,
        help='Minimum energy window')
    parser.add_argument('--Emax', default=5, type=float,
        help='Maximum energy window')
    args = parser.parse_args()

    plt.figure(figsize=(5,5))
    plot(args.fname, args.fname_order, args.fname_proj, args.Emin, args.Emax)

    plt.tight_layout()
    plt.savefig('bs.png')
    plt.show()
