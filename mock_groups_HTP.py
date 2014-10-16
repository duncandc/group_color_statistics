#Duncan Campbell
#April 7, 2014
#Calculate the halo transiton probability function.

def main():
    import numpy as np
    import h5py
    import matplotlib.pyplot as plt
    import custom_utilities as cu
    import sys

    group_cat = sys.argv[1]  # tinker_mocks, berlind_mocks, yang_mocks
    catalogue = sys.argv[2]  # Mr19_age_distribution_matching_mock, Mr19_age_distribution_matching_mock_sys_empty_shuffle_satrel_shuffle

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
    
    #Pcc
    selection = np.in1d(inf_cen,halo_cen)
    selection = inf_cen[selection]
    H_cc, xedges, yedges = np.histogram2d(inf_mass[selection], halo_mass[selection], bins=mass_bins)
    sums_cc = np.sum(H_cc,axis=1)
    
    #Pss
    selection = np.in1d(inf_sat, halo_sat)
    selection = inf_sat[selection]
    H_ss, xedges, yedges = np.histogram2d(inf_mass[selection], halo_mass[selection], bins=mass_bins)
    sums_ss = np.sum(H_ss,axis=1)

    #Psc
    selection = np.in1d(inf_cen, halo_sat)
    selection = inf_cen[selection]
    H_sc, xedges, yedges = np.histogram2d(inf_mass[selection], halo_mass[selection], bins=mass_bins)
    sums_sc = np.sum(H_sc,axis=1)

    #Pcs
    selection = np.in1d(inf_sat, halo_cen)
    selection = inf_sat[selection]
    H_cs, xedges, yedges = np.histogram2d(inf_mass[selection], halo_mass[selection], bins=mass_bins)
    sums_cs = np.sum(H_cs,axis=1)

    #Pcc
    for i in range(0,len(H_cc)):
        H_cc[i]=H_cc[i]/(sums_cc[i]+sums_sc[i])
    remove = np.isnan(H_cc)
    H_cc[remove]=0 

    #Pss
    for i in range(0,len(H_ss)):
        H_ss[i]=H_ss[i]/(sums_ss[i]+sums_cs[i])
    remove = np.isnan(H_ss)
    H_ss[remove]=0 

    #Pcs
    for i in range(0,len(H_cs)):
        H_cs[i]=H_cs[i]/(sums_cs[i]+sums_ss[i])
    remove = np.isnan(H_cs)
    H_cs[remove]=0

    #Psc
    for i in range(0,len(H_cs)):
        H_sc[i]=H_sc[i]/(sums_sc[i]+sums_cc[i])
    remove = np.isnan(H_sc)
    H_sc[remove]=0

    return H_cc, H_sc, H_ss, H_cs

def plot_htp(H_cc, H_sc, H_ss, H_cs, mass_bins):
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    
    xedges=mass_bins
    yedges=mass_bins

    sums_cc = np.sum(H_cc,axis=1)
    sums_cs = np.sum(H_cs,axis=1)
    sums_ss = np.sum(H_ss,axis=1)
    sums_sc = np.sum(H_sc,axis=1)

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
    ax.pcolormesh(X, Y, H_sc.T, vmin=0.001, vmax=1, cmap='Oranges', norm=matplotlib.colors.LogNorm())
    ax.plot([min(mass_bins),max(mass_bins)],[min(mass_bins),max(mass_bins)],color='white')
    ax.set_ylabel(r'$log(M_{halo}) [M_{\odot}h^{-1}]$')
    ax.set_yticks([11.5,12,12.5,13,13.5,14,14.5])
    ax.xaxis.set_ticklabels([])
    ax.text(13.25,11.5,'halo satellites', color='black')

    ax = axes[3]
    X, Y = np.meshgrid(xedges, yedges)
    p = ax.pcolormesh(X, Y, H_cs.T, vmin=0.001, vmax=1,cmap='Oranges', norm=matplotlib.colors.LogNorm())
    ax.plot([min(mass_bins),max(mass_bins)],[min(mass_bins),max(mass_bins)],color='white')
    ax.xaxis.set_ticklabels([])
    ax.set_yticks([11.5,12,12.5,13,13.5,14,14.5])
    ax.text(11.5,14.5,'halo centrals', color='black')

    cbar_ax = fig.add_axes([0.875, 0.2, 0.025, 0.7])
    fig.colorbar(mappable=p,cax=cbar_ax)
    cbar_ax.set_ylabel('frequency')
    cbar_ax.set_yticklabels(['','0.01','0.1','1'])

    hist_1_ax = fig.add_axes([0.125, 0.1, 0.354, 0.1])
    hist_2_ax = fig.add_axes([0.4965, 0.1, 0.354, 0.1])

    N = np.array(sums_cc, dtype=np.float)
    N_1 = N
    N = np.array(sums_sc, dtype=np.float)
    N_2 = N
    hist_1_ax.plot(mass_bins[:-1], N_1/(N_1+N_2), color='black')
    hist_1_ax.set_ylim([0,1])
    hist_1_ax.set_xlim([11,15])
    hist_1_ax.set_ylabel('C')
    hist_1_ax.set_xlabel(r'$log(M_{group}) [M_{\odot}h^{-1}]$')
    hist_1_ax.set_yticks([0,0.25,0.5,0.75])
    hist_1_ax.set_xticks([11,11.5,12,12.5,13,13.5,14,14.5])

    N = np.array(sums_ss, dtype=np.float)
    N_1 = N
    N = np.array(sums_cs, dtype=np.float)
    N_2 = N
    hist_2_ax.plot(mass_bins[:-1], N_1/(N_1+N_2), color='black')
    hist_2_ax.set_ylim([0,1])
    hist_2_ax.set_xlim([11,15])
    hist_2_ax.set_xlabel(r'$log(M_{group}) [M_{\odot}h^{-1}]$')
    hist_2_ax.set_yticks([0,0.25,0.5,0.75])
    hist_2_ax.yaxis.set_ticklabels([])
    hist_2_ax.set_xticks([11,11.5,12,12.5,13,13.5,14,14.5,15.0])

    return fig 

if __name__ == '__main__':
    main()
