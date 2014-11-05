#Duncan Campbell
#April 7, 2014
#Calculate the halo transiton probability function.

def main():
    import numpy as np
    import h5py
    import matplotlib.pyplot as plt
    import custom_utilities as cu
    import sys

    if len(sys.argv)>1: group_cat = sys.argv[1]
    else: group_cat = 'tinker_mocks'
    if len(sys.argv)>2: catalogue = sys.argv[2]
    else: catalogue = 'Mr19_age_distribution_matching_mock'

    filename=catalogue+'_HTP_mass'

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
    mass_bins = np.arange(11,15.2,0.1)
    bin_width = mass_bins[1]-mass_bins[0]
    bin_centers = (mass_bins[:-1]-mass_bins[1:])/2.0

    #calculate the HTP
    H = htp(GC['MGROUP'], GC['HALO_M'], group_centrals, group_satellites, halo_centrals, halo_satellites, mass_bins) 

    #plot the results
    fig = plot_htp(H, mass_bins)
    plt.show(block=True)
    print plotpath+filename+'.eps'
    fig.savefig(plotpath+filename+'.eps') 

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
    
    H, xedges, yedges = np.histogram2d(inf_mass, halo_mass, bins=mass_bins)
    sums = np.sum(H,axis=1)

    for i in range(0,len(H)):
        H[i]=H[i]/(sums[i])
    remove = np.isnan(H)
    H[remove]=0 

    return H

def plot_htp(H, mass_bins):
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    
    xedges=mass_bins
    yedges=mass_bins

    sums = np.sum(H,axis=1)

    fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(3.3,3.3))
    fig.subplots_adjust(hspace=0,wspace=0.05)
    fig.subplots_adjust(left=0.2, right=0.8, bottom=0.2,top=0.9)

    ax = axes
    X, Y = np.meshgrid(xedges, yedges)
    p = ax.pcolormesh(X, Y, H.T, vmin=0.001, vmax=1,cmap='Oranges', norm=matplotlib.colors.LogNorm())
    ax.plot([min(mass_bins),max(mass_bins)],[min(mass_bins),max(mass_bins)],color='white')
    ax.set_ylabel(r'$log(M_{halo}) [M_{\odot}h^{-1}]$')
    ax.set_yticks([11.5,12,12.5,13,13.5,14,14.5])
    ax.set_xlabel(r'$log(M_{group}) [M_{\odot}h^{-1}]$')
    ax.set_xticks([11.5,12,12.5,13,13.5,14,14.5])
    
    cbar_ax = fig.add_axes([0.875, 0.2, 0.025, 0.7])
    fig.colorbar(mappable=p,cax=cbar_ax)
    cbar_ax.set_ylabel('frequency')
    cbar_ax.set_yticklabels(['','0.01','0.1','1'])
   
    return fig 

if __name__ == '__main__':
    main()
