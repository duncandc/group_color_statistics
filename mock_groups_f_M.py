import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
import custom_utilities as cu
from astropy import cosmology
import sys
from f_prop import f_prop

def main():
    catalogue = sys.argv[1]
    plotpath = cu.get_plot_path() + 'analysis/groupcats/'

    filename = catalogue+'_mock_group_f_red_sat_M.eps'

    fig1, axes = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, figsize=(6.95, 3.3))
    fig1.subplots_adjust(hspace=0, wspace=0.05)
    fig1.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.9)
    axes = axes.flatten()
    ax=axes[0]
    ax.set_xlim([11,15])
    ax.set_ylim([0,1])
    ax.set_ylabel(r'$f_{red}$')
    ax.set_xlabel(r'$log(M/[M_{\odot}h^{-1}])$')
    ax.set_xticks([9.6,9.8,10.0,10.2,10.4,10.6])
    ax.set_title(r'Berlind FoF groups')
    ax=axes[1]
    ax.set_xlim([11,15])
    ax.set_ylim([0,1])
    ax.set_xlabel(r'$log(M/[M_{\odot}h^{-1}])$')
    ax.set_xticks([12.5,13,13.5,14,14.5])
    ax.set_xlim([12,15])
    ax.set_title(r'Tinker SO groups')
    ax=axes[2]
    ax.set_xlim([11,15])
    ax.set_ylim([0,1])
    ax.set_xlabel(r'$log(M/[M_{\odot}h^{-1}])$')
    ax.set_xticks([12.5,13,13.5,14,14.5])
    ax.set_xlim([12,15])
    ax.set_title(r'Yang SO groups')

    N_boots = 50

    bins = np.arange(11,15,0.2)
    bin_centers = (bins[:-1]+bins[1:])/2.0

    #open mock catalogue
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    print 'opening mock catalogue:', catalogue+'.hdf5'
    f1 = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f1.get(catalogue)

    centrals   = np.array(mock['ID_host']==-1)
    satellites = np.array(mock['ID_host']!=-1)
    
    #galaxy color
    color = mock['g-r']
    LHS   = 0.21-0.03*mock['M_r,0.1']
    blue  = np.where(color<LHS)[0] #indices of blue galaxies
    red   = np.where(color>LHS)[0] #indicies of red galaxies

    f_red_cen = f_prop(mock['M_host'],bins,red,blue,centrals)
    f_red_sat = f_prop(mock['M_host'],bins,red,blue,satellites)
    
    ax=axes[0]
    p1a, = ax.plot(bin_centers,f_red_cen,color='orange')
    p2a, = ax.plot(bin_centers,f_red_sat,color='green')
    ax=axes[1]
    p1a, = ax.plot(bin_centers,f_red_cen,color='orange')
    p2a, = ax.plot(bin_centers,f_red_sat,color='green')
    ax=axes[2]
    p1a, = ax.plot(bin_centers,f_red_cen,color='orange')
    p2a, = ax.plot(bin_centers,f_red_sat,color='green')


    ###################################################
    # Ideal Groups
    ###################################################
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/ideal_groups/'
    print 'opening mock catalogue:', catalogue+'_groups.hdf5'
    f1 = h5py.File(filepath_mock+catalogue+'_groups.hdf5', 'r') #open catalogue file
    GC = f1.get(catalogue+'_groups')

    f_red_cen = np.zeros((len(bin_centers)))
    f_red_sat = np.zeros((len(bin_centers)))
    f_sat_red = np.zeros((len(bin_centers)))
    f_sat_blue = np.zeros((len(bin_centers)))

    #with mass estimate but not cen/sat
    centrals_ind   = np.where(GC['HALO_RANK']==0)[0]
    satellites_ind = np.where(GC['HALO_RANK']!=0)[0]
    centrals_bool   = (GC['HALO_RANK']==0)
    satellites_bool = (GC['HALO_RANK']!=0)

    centrals_mock_ind   = np.where(GC['HALO_RANK']==1)[0]
    satellites_mock_ind = np.where(GC['HALO_RANK']!=1)[0]
    N_sat_mock = float(len(satellites_mock_ind))
    N_cen_mock = float(len(centrals_mock_ind))
    
    #galaxy color
    color = GC['M_g,0.1']-GC['M_r,0.1']
    LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
    blue_ind  = np.where(color<LHS)[0] #indices of blue galaxies
    red_ind   = np.where(color>LHS)[0] #indicies of red galaxies
    blue_bool = (color<LHS) #indices of blue galaxies
    red_bool  = (color>LHS) #indicies of red galaxies

    f_red_cen  = f_prop(GC['MGROUP'],bins,red_ind,blue_ind,centrals_bool)
    f_red_sat  = f_prop(GC['MGROUP'],bins,red_ind,blue_ind,satellites_bool)
    f_sat_red  = f_prop(GC['MGROUP'],bins,satellites_ind,centrals_ind,red_bool)
    f_sat_blue = f_prop(GC['MGROUP'],bins,satellites_ind,centrals_ind,blue_bool)

    ax=axes[0]
    p1a, = ax.plot(bin_centers,f_red_cen,color='orange', alpha=0.5)
    p2a, = ax.plot(bin_centers,f_red_sat,color='green', alpha=0.5)
    ax=axes[1]
    p1a, = ax.plot(bin_centers,f_red_cen,color='orange', alpha=0.5)
    p2a, = ax.plot(bin_centers,f_red_sat,color='green', alpha=0.5)
    ax=axes[2]
    p1a, = ax.plot(bin_centers,f_red_cen,color='orange', alpha=0.5)
    p2a, = ax.plot(bin_centers,f_red_sat,color='green', alpha=0.5)
    
    #with cen/sat but not mass estimate
    centrals_ind   = np.where(GC['RANK']==0)[0]
    satellites_ind = np.where(GC['RANK']!=0)[0]
    centrals_bool   = (GC['RANK']==0)
    satellites_bool = (GC['RANK']!=0)

    centrals_mock_ind   = np.where(GC['HALO_RANK']==1)[0]
    satellites_mock_ind = np.where(GC['HALO_RANK']!=1)[0]
    N_sat_mock = float(len(satellites_mock_ind))
    N_cen_mock = float(len(centrals_mock_ind))
    
    #galaxy color
    color = GC['M_g,0.1']-GC['M_r,0.1']
    LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
    blue_ind  = np.where(color<LHS)[0] #indices of blue galaxies
    red_ind   = np.where(color>LHS)[0] #indicies of red galaxies
    blue_bool = (color<LHS) #indices of blue galaxies
    red_bool  = (color>LHS) #indicies of red galaxies

    f_red_cen  = f_prop(GC['HALO_M'],bins,red_ind,blue_ind,centrals_bool)
    f_red_sat  = f_prop(GC['HALO_M'],bins,red_ind,blue_ind,satellites_bool)
    f_sat_red  = f_prop(GC['HALO_M'],bins,satellites_ind,centrals_ind,red_bool)
    f_sat_blue = f_prop(GC['HALO_M'],bins,satellites_ind,centrals_ind,blue_bool)

    ax=axes[0]
    p1a, = ax.plot(bin_centers,f_red_cen,'--',color='orange', alpha=0.5)
    p2a, = ax.plot(bin_centers,f_red_sat,'--',color='green', alpha=0.5)
    ax=axes[1]
    p1a, = ax.plot(bin_centers,f_red_cen,'--',color='orange', alpha=0.5)
    p2a, = ax.plot(bin_centers,f_red_sat,'--',color='green', alpha=0.5)
    ax=axes[2]
    p1a, = ax.plot(bin_centers,f_red_cen,'--',color='orange', alpha=0.5)
    p2a, = ax.plot(bin_centers,f_red_sat,'--',color='green', alpha=0.5)


    ###################################################
    groupcat='berlind'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    f_red_cen = np.zeros((N_boots,len(bin_centers)))
    f_red_sat = np.zeros((N_boots,len(bin_centers)))
    f_sat_red = np.zeros((N_boots,len(bin_centers)))
    f_sat_blue = np.zeros((N_boots,len(bin_centers)))
    for boot in range(0,N_boots):
        #open catalogue
        catalogue_1 = group_catalogue+'_'+str(boot)
        f =  h5py.File(filepath_cat+'bootstraps/'+catalogue_1+'.hdf5', 'r') #open catalogue file
        GC = f.get(catalogue_1)

        centrals_ind   = np.where(GC['RANK']==1)[0]
        satellites_ind = np.where(GC['RANK']!=1)[0]
        centrals_bool   = (GC['RANK']==1)
        satellites_bool = (GC['RANK']!=1)

        centrals_mock_ind   = np.where(GC['HALO_RANK']==1)[0]
        satellites_mock_ind = np.where(GC['HALO_RANK']!=1)[0]
        N_sat_mock = float(len(satellites_mock_ind))
        N_cen_mock = float(len(centrals_mock_ind))
    
        #galaxy color
        color = GC['M_g,0.1']-GC['M_r,0.1']
        LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
        blue_ind  = np.where(color<LHS)[0] #indices of blue galaxies
        red_ind   = np.where(color>LHS)[0] #indicies of red galaxies
        blue_bool = (color<LHS) #indices of blue galaxies
        red_bool  = (color>LHS) #indicies of red galaxies

        f_red_cen[boot,:]  = f_prop(GC['MGROUP'],bins,red_ind,blue_ind,centrals_bool)
        f_red_sat[boot,:]  = f_prop(GC['MGROUP'],bins,red_ind,blue_ind,satellites_bool)
        f_sat_red[boot,:]  = f_prop(GC['MGROUP'],bins,satellites_ind,centrals_ind,red_bool)
        f_sat_blue[boot,:] = f_prop(GC['MGROUP'],bins,satellites_ind,centrals_ind,blue_bool)

    err_cen = np.nanstd(f_red_cen, axis=0)
    err_sat = np.nanstd(f_red_sat, axis=0)
    f_red_cen = np.nanmean(f_red_cen, axis=0)
    f_red_sat = np.nanmean(f_red_sat, axis=0)

    err_sat_red = np.nanstd(f_sat_red, axis=0)
    err_sat_blue = np.nanstd(f_sat_blue, axis=0)
    f_sat_red = np.nanmean(f_sat_red, axis=0)
    f_sat_blue = np.nanmean(f_sat_blue, axis=0)

    ax=axes[0]
    p3a=ax.errorbar(bin_centers, f_red_cen, yerr=err_cen,fmt='o',color='orange', mec='none', ms=3)
    p4a=ax.errorbar(bin_centers, f_red_sat, yerr=err_sat,fmt='o',color='green', mec='none', ms=3)
    ax.legend((p1a,p2a),('halo cen','halo sat'), loc='lower right', fontsize=10, numpoints=1, frameon=False)


    ###################################################
    groupcat='tinker'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    f_red_cen = np.zeros((N_boots,len(bin_centers)))
    f_red_sat = np.zeros((N_boots,len(bin_centers)))
    f_sat_red = np.zeros((N_boots,len(bin_centers)))
    f_sat_blue = np.zeros((N_boots,len(bin_centers)))
    for boot in range(0,N_boots):
        #open catalogue
        catalogue_1 = group_catalogue+'_'+str(boot)
        f =  h5py.File(filepath_cat+'bootstraps/'+catalogue_1+'.hdf5', 'r') #open catalogue file
        GC = f.get(catalogue_1)

        centrals_ind   = np.where(GC['RANK']==1)[0]
        satellites_ind = np.where(GC['RANK']!=1)[0]
        centrals_bool   = (GC['RANK']==1)
        satellites_bool = (GC['RANK']!=1)

        centrals_mock_ind   = np.where(GC['HALO_RANK']==1)[0]
        satellites_mock_ind = np.where(GC['HALO_RANK']!=1)[0]
        N_sat_mock = float(len(satellites_mock_ind))
        N_cen_mock = float(len(centrals_mock_ind))
    
        #galaxy color
        color = GC['M_g,0.1']-GC['M_r,0.1']
        LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
        blue_ind  = np.where(color<LHS)[0] #indices of blue galaxies
        red_ind   = np.where(color>LHS)[0] #indicies of red galaxies
        blue_bool = (color<LHS) #indices of blue galaxies
        red_bool  = (color>LHS) #indicies of red galaxies

        f_red_cen[boot,:]  = f_prop(GC['MGROUP'],bins,red_ind,blue_ind,centrals_bool)
        f_red_sat[boot,:]  = f_prop(GC['MGROUP'],bins,red_ind,blue_ind,satellites_bool)
        f_sat_red[boot,:]  = f_prop(GC['MGROUP'],bins,satellites_ind,centrals_ind,red_bool)
        f_sat_blue[boot,:] = f_prop(GC['MGROUP'],bins,satellites_ind,centrals_ind,blue_bool)

    err_cen = np.nanstd(f_red_cen, axis=0)
    err_sat = np.nanstd(f_red_sat, axis=0)
    f_red_cen = np.nanmean(f_red_cen, axis=0)
    f_red_sat = np.nanmean(f_red_sat, axis=0)

    err_sat_red = np.nanstd(f_sat_red, axis=0)
    err_sat_blue = np.nanstd(f_sat_blue, axis=0)
    f_sat_red = np.nanmean(f_sat_red, axis=0)
    f_sat_blue = np.nanmean(f_sat_blue, axis=0)

    ax=axes[1]
    p3a=ax.errorbar(bin_centers, f_red_cen, yerr=err_cen,fmt='o',color='orange', mec='none', ms=3)
    p4a=ax.errorbar(bin_centers, f_red_sat, yerr=err_sat,fmt='o',color='green', mec='none', ms=3)
    ax.legend((p3a,p4a),('groups cen','groups sat'), loc='lower right', fontsize=10, numpoints=1, frameon=False)

  
    ###################################################
    groupcat='yang'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    f_red_cen = np.zeros((N_boots,len(bin_centers)))
    f_red_sat = np.zeros((N_boots,len(bin_centers)))
    f_sat_red = np.zeros((N_boots,len(bin_centers)))
    f_sat_blue = np.zeros((N_boots,len(bin_centers)))
    for boot in range(0,N_boots):
        #open catalogue
        catalogue_1 = group_catalogue+'_'+str(boot)
        f =  h5py.File(filepath_cat+'bootstraps/'+catalogue_1+'.hdf5', 'r') #open catalogue file
        GC = f.get(catalogue_1)

        centrals_ind   = np.where(GC['RANK']==1)[0]
        satellites_ind = np.where(GC['RANK']!=1)[0]
        centrals_bool   = (GC['RANK']==1)
        satellites_bool = (GC['RANK']!=1)

        centrals_mock_ind   = np.where(GC['HALO_RANK']==1)[0]
        satellites_mock_ind = np.where(GC['HALO_RANK']!=1)[0]
        N_sat_mock = float(len(satellites_mock_ind))
        N_cen_mock = float(len(centrals_mock_ind))
    
        #galaxy color
        color = GC['M_g,0.1']-GC['M_r,0.1']
        LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
        blue_ind  = np.where(color<LHS)[0] #indices of blue galaxies
        red_ind   = np.where(color>LHS)[0] #indicies of red galaxies
        blue_bool = (color<LHS) #indices of blue galaxies
        red_bool  = (color>LHS) #indicies of red galaxies

        f_red_cen[boot,:]  = f_prop(GC['MGROUP'],bins,red_ind,blue_ind,centrals_bool)
        f_red_sat[boot,:]  = f_prop(GC['MGROUP'],bins,red_ind,blue_ind,satellites_bool)
        f_sat_red[boot,:]  = f_prop(GC['MGROUP'],bins,satellites_ind,centrals_ind,red_bool)
        f_sat_blue[boot,:] = f_prop(GC['MGROUP'],bins,satellites_ind,centrals_ind,blue_bool)

    err_cen = np.nanstd(f_red_cen, axis=0)
    err_sat = np.nanstd(f_red_sat, axis=0)
    f_red_cen = np.nanmean(f_red_cen, axis=0)
    f_red_sat = np.nanmean(f_red_sat, axis=0)

    err_sat_red = np.nanstd(f_sat_red, axis=0)
    err_sat_blue = np.nanstd(f_sat_blue, axis=0)
    f_sat_red = np.nanmean(f_sat_red, axis=0)
    f_sat_blue = np.nanmean(f_sat_blue, axis=0)

    ax=axes[2]
    p3a=ax.errorbar(bin_centers, f_red_cen, yerr=err_cen,fmt='o',color='orange', mec='none', ms=3)
    p4a=ax.errorbar(bin_centers, f_red_sat, yerr=err_sat,fmt='o',color='green', mec='none', ms=3)

   
    plt.show()

    fig1.savefig(plotpath+filename)
  
def get_gc_path(groupcat):

    if groupcat=='berlind':
        s = cu.get_output_path() + 'processed_data/berlind_groupcat/mock_runs/4th_run/custom_catalogues/'
    if groupcat=='tinker':
        s = cu.get_output_path() + 'processed_data/tinker_groupcat/mock_runs/4th_run/custom_catalogues/'
    if groupcat=='yang':
        s = cu.get_output_path() + 'processed_data/yang_groupcat/mock_runs/4th_run/version_2/custom_catalogues/'
    return s

def get_gc_name(groupcat,catalogue):

    if groupcat=='tinker': s = catalogue+'_clf_groups_M19'
    if groupcat=='berlind': s = catalogue+'_groups'
    if groupcat=='yang': s = catalogue+'_groups'
    return s

if __name__ == '__main__':
    main() 
