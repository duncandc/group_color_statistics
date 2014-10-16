#Duncan Campbell
#April 7, 2014
#Calculate the halo transiton probability function.

from __future__ import division
import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
import custom_utilities as cu
import sys

def main():

    if len(sys.argv)>1: group_cat = sys.argv[1]
    else: group_cat = 'tinker_mocks'
    if len(sys.argv)>2: catalogue = sys.argv[2]
    else: catalogue = 'Mr19_age_distribution_matching_mock'

    filename=catalogue+'_HTP'

    if group_cat == 'tinker_mocks':
        filepath = cu.get_output_path() + 'processed_data/tinker_groupcat/mock_runs/4th_run/custom_catalogues/'
        savepath = cu.get_output_path() + 'analysis/tinker_groupcat/'
        plotpath = cu.get_plot_path() + 'analysis/tinker_groupcat/'
        catalogue = catalogue + '_clf_groups_M19'
    if group_cat == 'berlind_mocks':
        filepath = cu.get_output_path() + 'processed_data/berlind_groupcat/mock_runs/4th_run/custom_catalogues/'
        savepath = cu.get_output_path() + 'analysis/berlind_groupcat/'
        plotpath = cu.get_plot_path() + 'analysis/berlind_groupcat/'
        catalogue = catalogue + '_groups'
    if group_cat == 'yang_mocks':
        filepath = cu.get_output_path() + 'processed_data/yang_groupcat/mock_runs/4th_run/version_2/custom_catalogues/'
        savepath = cu.get_output_path() + 'analysis/yang_groupcat/'
        plotpath = cu.get_plot_path() + 'analysis/yang_groupcat/'
        catalogue = catalogue + '_groups'

    print 'opening group catalogue:', catalogue
    #open catalogue
    f  =  h5py.File(filepath+catalogue+'.hdf5', 'r') #open catalogue file
    GC = f.get(catalogue)

    #determine central/satellite designation
    group_centrals   = np.where(GC['RANK']==1)[0]
    group_satellites = np.where(GC['RANK']!=1)[0]
    halo_centrals    = np.where(GC['HALO_RANK']==1)[0]
    halo_satellites  = np.where(GC['HALO_RANK']!=1)[0]

    #define group/halo mass bins
    mass_bins = np.arange(11,15.2,0.2)
    bin_width = mass_bins[1]-mass_bins[0]
    bin_centers = (mass_bins[:-1]-mass_bins[1:])/2.0

    #calculate the HTP
    H_cc, H_sc, H_ss, H_cs = htp(GC['MGROUP'], GC['HALO_M'], group_centrals, group_satellites, halo_centrals, halo_satellites, mass_bins) 

    #plot the results
    fig = plot_htp(H_cc, H_sc, H_ss, H_cs, mass_bins)
    plt.show(block=True)
    #print plotpath+filename+'.eps'
    #fig.savefig(plotpath+filename+'.eps') 

def htp(inf_mass, halo_mass, inf_cen, inf_sat, halo_cen, halo_sat, mass_bins):
    #returns the halo transition probability given:
    #    inf_mass: inferred group mass containing the galaxy
    #    halo_mass: halo mass containing the galaxy
    #    inf_cen: indices of the inffered centrals of the groups
    #    inf_sat: indices of the inffered satellites of the groups
    #    halo_cen: indices of the halo centrals
    #    halo_sat: indices of the halo satellites
    import numpy as np

    #check to make sure the inputs of of the right size
    if len(inf_mass)!=len(halo_mass): return 'error'
    
    counts, inds = cu.histogram2d(inf_mass, halo_mass, mass_bins, mass_bins)
    #sums = np.sum(H,axis=1) #sum horizontal
    #sums = np.sum(H,axis=0) #sum vertical

    H_cc = np.zeros((len(mass_bins),len(mass_bins)))
    H_ss = np.zeros((len(mass_bins),len(mass_bins)))
    H_sc = np.zeros((len(mass_bins),len(mass_bins)))
    H_cs = np.zeros((len(mass_bins),len(mass_bins)))
    
    for i in range(0,len(mass_bins)):
        for j in range(0,len(mass_bins)):
            selection = np.where((inds[:,0]==i) & (inds[:,1]==j))[0]
            Nsample = len(selection)
            selection = np.in1d(selection,halo_cen)
            selection = halo_cen[selection]
            Nsub = np.sum(np.in1d(selection,inf_cen))
            p = Nsub/Nsample
            H_cc[i,j] = p
    
    for i in range(0,len(mass_bins)):
        for j in range(0,len(mass_bins)):
            selection = np.where((inds[:,0]==i) & (inds[:,1]==j))[0]
            Nsample = len(selection)
            selection = np.in1d(selection,halo_sat)
            selection = halo_sat[selection]
            Nsub = np.sum(np.in1d(selection,inf_sat))
            p = Nsub/Nsample
            H_ss[i,j] = p
    
    for i in range(0,len(mass_bins)):
        for j in range(0,len(mass_bins)):
            selection = np.where((inds[:,0]==i) & (inds[:,1]==j))[0]
            Nsample = len(selection)
            selection = np.in1d(selection,halo_cen)
            selection = halo_cen[selection]
            Nsub = np.sum(np.in1d(selection,inf_sat))
            p = Nsub/Nsample
            H_sc[i,j] = p
    
    for i in range(0,len(mass_bins)):
        for j in range(0,len(mass_bins)):
            selection = np.where((inds[:,0]==i) & (inds[:,1]==j))[0]
            Nsample = len(selection)
            selection = np.in1d(selection,halo_sat)
            selection = halo_sat[selection]
            Nsub = np.sum(np.in1d(selection,inf_cen))
            p = Nsub/Nsample
            H_cs[i,j] = p
                    
    return H_cc, H_sc, H_ss, H_cs

def plot_htp(H_cc, H_sc, H_ss, H_cs, mass_bins):
    
    xedges=mass_bins
    yedges=mass_bins

    fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0,wspace=0.05)
    fig.subplots_adjust(right=0.85, bottom=0.2)
    axes = axes.flatten()

    ax = axes[0]
    X, Y = np.meshgrid(xedges, yedges)
    ax.pcolormesh(X, Y, H_cc.T, vmin=0.001, vmax=1,cmap='Oranges', norm=matplotlib.colors.LogNorm())
    ax.plot([min(mass_bins),max(mass_bins)],[min(mass_bins),max(mass_bins)],color='white')
    ax.set_ylabel(r'$log(M_{halo}) [M_{\odot}h^{-1}]$')
    ax.set_title('group centrals')
    ax.text(13.25,11.5,'halo centrals', color='black')

    ax = axes[1]
    X, Y = np.meshgrid(xedges, yedges)
    ax.pcolormesh(X, Y, H_ss.T, vmin=0.001, vmax=1, cmap='Oranges', norm=matplotlib.colors.LogNorm())
    ax.plot([min(mass_bins),max(mass_bins)],[min(mass_bins),max(mass_bins)],color='white')
    ax.set_title('group satellites')
    ax.text(11.5,14.5,'halo satellites', color='black')
    
    ax = axes[2]
    X, Y = np.meshgrid(xedges, yedges)
    ax.pcolormesh(X, Y, H_cs.T, vmin=0.001, vmax=1, cmap='Oranges', norm=matplotlib.colors.LogNorm())
    ax.plot([min(mass_bins),max(mass_bins)],[min(mass_bins),max(mass_bins)],color='white')
    ax.set_ylabel(r'$log(M_{halo}) [M_{\odot}h^{-1}]$')
    ax.set_yticks([11.5,12,12.5,13,13.5,14,14.5])
    ax.xaxis.set_ticklabels([])
    ax.text(13.25,11.5,'halo satellites', color='black')
    ax.set_xlabel(r'$log(M_{group}) [M_{\odot}h^{-1}]$')
    ax.set_xticks([11.5,12,12.5,13,13.5,14,14.5])

    ax = axes[3]
    X, Y = np.meshgrid(xedges, yedges)
    p = ax.pcolormesh(X, Y, H_sc.T, vmin=0.001, vmax=1,cmap='Oranges', norm=matplotlib.colors.LogNorm())
    ax.plot([min(mass_bins),max(mass_bins)],[min(mass_bins),max(mass_bins)],color='white')
    ax.xaxis.set_ticklabels([])
    ax.set_yticks([11.5,12,12.5,13,13.5,14,14.5])
    ax.text(11.5,14.5,'halo centrals', color='black')
    ax.set_xlabel(r'$log(M_{group}) [M_{\odot}h^{-1}]$')
    ax.set_xticks([11.5,12,12.5,13,13.5,14,14.5])

    cbar_ax = fig.add_axes([0.875, 0.2, 0.025, 0.7])
    fig.colorbar(mappable=p,cax=cbar_ax)
    cbar_ax.set_ylabel('frequency')
    cbar_ax.set_yticklabels(['','0.01','0.1','1'])

    return fig 

if __name__ == '__main__':
    main()
