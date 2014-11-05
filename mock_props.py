#!/usr/bin/env python

#Duncan Campbell
#April 8, 2014
#Yale University
#Plot some properties from the mock

import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu
import sys

def main():
    
    if len(sys.argv)>1: catalogue = sys.argv[1]
    else: catalogue = 'Mr19_age_distribution_matching_mock'
    
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    plotpath = cu.get_plot_path()        + 'analysis/hearin_mocks/'

    print 'opening mock catalogue:', catalogue+'.hdf5'
    #open catalogue
    f1 = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f1.get(catalogue)
    mock = np.array(mock)
    print 'length:', len(mock)
    for name in mock.dtype.names: print '     ', name
    
    #calculate the satellite fraction
    N_sat = float(np.sum((mock['ID_host']!=-1)))
    N_cen = float(np.sum((mock['ID_host']==-1)))
    fsat = N_sat/(N_sat+N_cen)
    print 'satellite fraction:', fsat
    
    #subsample the galaxies to make plotting easier
    N=10000
    selection = np.random.permutation(len(mock))[0:N]
    mock = mock[selection]
    
    
    host =  np.where(mock['ID_host']==-1)[0]
    sub = np.where(mock['ID_host']!=-1)[0]
    
    color    = mock['g-r']
    LHS      = 0.7 - 0.032*(mock['M_r,0.1']+16.5) #Weinmann 2006
    blue     = np.where((color<LHS) & (mock['M_r,0.1']!=-99))[0]
    red      = np.where((color>LHS) & (mock['M_r,0.1']!=-99))[0]
    
    fig = plt.figure(figsize=(3.3,3.3))
    ax = fig.add_subplot(1,1,1)
    fig.subplots_adjust(left=0.2, right=0.9, bottom=0.2,top=0.9)
    p1, = ax.plot(10.0**mock['M_host'][host],mock['V_peak'][host],'o',ms=1, alpha=1, color='orange', mew=0, rasterized=True)
    p2, = ax.plot(10.0**mock['M_host'][sub],mock['V_peak'][sub],'o',ms=1, alpha=1, color='green', mew=0, rasterized=True)
    ax.plot([10**10,10**16],[100,100],'--',color='black')
    ax.plot([2*10**11,2*10**11],[1,800],'--',color='black')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim([50,2000])
    ax.set_xlim([10**10,10**15])
    ax.set_xlabel(r'$M_{\rm host}$ $[M_{\odot}/h]$')
    ax.set_ylabel(r'$V_{\rm peak}$')
    ax.legend((p1,p2),('centrals','satellites'),loc=2, numpoints=1, fontsize=10, markerscale=2, frameon=False)
    plt.show()
    filename1 = 'mock_props_censata_MhVpeak.pdf'
    fig.savefig(plotpath+filename1, dpi=400)
    
    fig = plt.figure(figsize=(3.3,3.3))
    ax = fig.add_subplot(1,1,1)
    fig.subplots_adjust(left=0.2, right=0.9, bottom=0.2,top=0.9)
    p1, = ax.plot(10.0**mock['M_host'][red],mock['V_peak'][red],'o',ms=1, alpha=1, color='red', mew=0, rasterized=True)
    p2, = ax.plot(10.0**mock['M_host'][blue],mock['V_peak'][blue],'o',ms=1, alpha=1, color='blue', mew=0, rasterized=True)
    ax.plot([10**10,10**16],[100,100],'--',color='black')
    ax.plot([2*10**11,2*10**11],[1,800],'--',color='black')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim([50,2000])
    ax.set_xlim([10**10,10**15])
    ax.set_xlabel(r'$M_{\rm host}$ $[M_{\odot}/h]$')
    ax.set_ylabel(r'$V_{\rm peak}$')
    ax.legend((p1,p2),('red subsample','blue subsample'),loc=2, numpoints=1, fontsize=10, markerscale=2, frameon=False)
    plt.show()
    filename2 = 'mock_props_color_MhVpeak.pdf'
    fig.savefig(plotpath+filename2, dpi=400)


import numpy as np
from matplotlib import pyplot as plt


def scatter_contour(x, y,
                    levels=10,
                    threshold=100,
                    log_counts=False,
                    histogram2d_args=None,
                    plot_args=None,
                    contour_args=None,
                    filled_contour=True,
                    ax=None):
    """Scatter plot with contour over dense regions

    Parameters
    ----------
    x, y : arrays
        x and y data for the contour plot
    levels : integer or array (optional, default=10)
        number of contour levels, or array of contour levels
    threshold : float (default=100)
        number of points per 2D bin at which to begin drawing contours
    log_counts :boolean (optional)
        if True, contour levels are the base-10 logarithm of bin counts.
    histogram2d_args : dict
        keyword arguments passed to numpy.histogram2d
        see doc string of numpy.histogram2d for more information
    plot_args : dict
        keyword arguments passed to plt.plot.  By default it will use
        dict(marker='.', linestyle='none').
        see doc string of pylab.plot for more information
    contour_args : dict
        keyword arguments passed to plt.contourf or plt.contour
        see doc string of pylab.contourf for more information
    filled_contour : bool
        If True (default) use filled contours. Otherwise, use contour outlines.
    ax : pylab.Axes instance
        the axes on which to plot.  If not specified, the current
        axes will be used

    Returns
    -------
    points, contours :
       points is the return value of ax.plot()
       contours is the return value of ax.contour or ax.contourf
    """
    x = np.asarray(x)
    y = np.asarray(y)

    default_plot_args = dict(marker='.', linestyle='none')

    if plot_args is not None:
        default_plot_args.update(plot_args)
    plot_args = default_plot_args

    if histogram2d_args is None:
        histogram2d_args = {}

    if contour_args is None:
        contour_args = {}

    if ax is None:
        ax = plt.gca()

    H, xbins, ybins = np.histogram2d(x, y, **histogram2d_args)

    Nx = len(xbins)
    Ny = len(ybins)

    if log_counts:
        H = np.log10(1 + H)
        threshold = np.log10(1 + threshold)

    levels = np.asarray(levels)

    if levels.size == 1:
        levels = np.linspace(threshold, H.max(), levels)

    extent = [xbins[0], xbins[-1], ybins[0], ybins[-1]]

    i_min = np.argmin(levels)

    # draw a zero-width line: this gives us the outer polygon to
    # reduce the number of points we draw
    # somewhat hackish... we could probably get the same info from
    # the full contour plot below.
    outline = ax.contour(H.T, levels[i_min:i_min + 1],
                         linewidths=0, extent=extent)

    if filled_contour:
        contours = ax.contourf(H.T, levels, extent=extent, **contour_args)
    else:
        contours = ax.contour(H.T, levels, extent=extent, **contour_args)

    X = np.hstack([x[:, None], y[:, None]])

    if len(outline.allsegs[0]) > 0:
        outer_poly = outline.allsegs[0][0]
        try:
            # this works in newer matplotlib versions
            from matplotlib.path import Path
            points_inside = Path(outer_poly).contains_points(X)
        except:
            # this works in older matplotlib versions
            import matplotlib.nxutils as nx
            points_inside = nx.points_inside_poly(X, outer_poly)

        Xplot = X[~points_inside]
    else:
        Xplot = X

    points = ax.plot(Xplot[:, 0], Xplot[:, 1], zorder=1, **plot_args)

    return points, contours

    
if __name__ == '__main__':
    main() 
