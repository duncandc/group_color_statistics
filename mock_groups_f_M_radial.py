import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
import custom_utilities as cu
from astropy import cosmology
import sys
from f_prop import f_prop

def main():
    catalogue = 'Mr19_age_distribution_matching_mock'
    if len(sys.argv)>1: catalogue = sys.argv[1]
    plotpath = cu.get_plot_path() + 'analysis/groupcats/'

    filename = catalogue+'_mock_group_f_red_sat_M.eps'
    fig1, axes = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, figsize=(6.95, 3.3))
    fig1.subplots_adjust(hspace=0, wspace=0.05)
    fig1.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.9)
    axes = axes.flatten()
    ax=axes[0]
    ax.set_xlim([0,1.5])
    ax.set_ylim([0,1])
    ax.set_ylabel(r'$f_{red}$')
    ax.set_xlabel(r'$R_{\rm proj} [{\rm Mpc}]$')
    ax.set_title(r'$13<log(M)<13.5$')
    ax=axes[1]
    ax.set_xlim([0,1.5])
    ax.set_ylim([0,1])
    ax.set_xlabel(r'$R_{\rm proj} [{\rm Mpc}]$')
    ax.set_title(r'$13.5<log(M)<14$')
    ax=axes[2]
    ax.set_xlim([0,1.5])
    ax.set_ylim([0,1])
    ax.set_xlabel(r'$R_{\rm proj} [{\rm Mpc}]$')
    ax.set_title(r'$14<log(M)<14.5$')
    
    N_boots = 50

    bins = np.arange(0,1.5,0.1)
    bin_centers = (bins[:-1]+bins[1:])/2.0

    #open mock catalogue
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    print 'opening mock catalogue:', catalogue+'.hdf5'
    f1 = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f1.get(catalogue)
    print mock.dtype.names

    centrals   = np.array(mock['ID_host']==-1)
    satellites = np.array(mock['ID_host']!=-1)
    
    #galaxy color
    color = mock['g-r']
    LHS   = 0.21-0.03*mock['M_r,0.1']
    blue  = np.where(color<LHS)[0] #indices of blue galaxies
    red   = np.where(color>LHS)[0] #indicies of red galaxies

    selection = ((mock['M_host']>13.0) & (mock['M_host']<13.5))
    selection = satellites & selection
    f_red_1 = f_prop(mock['R_proj'][selection],bins,red,blue,selection)

    selection = ((mock['M_host']>13.5) & (mock['M_host']<14.0))
    selection = satellites & selection
    f_red_2 = f_prop(mock['R_proj'][selection],bins,red,blue,selection)

    selection = ((mock['M_host']>14.0) & (mock['M_host']<14.5))
    selection = satellites & selection
    f_red_3 = f_prop(mock['R_proj'][selection],bins,red,blue,selection)
    
    ax=axes[0]
    p1, = ax.plot(bin_centers,f_red_1,color='orange')
    ax=axes[1]
    p2, = ax.plot(bin_centers,f_red_2,color='orange')
    ax=axes[2]
    p3, = ax.plot(bin_centers,f_red_3,color='orange')

    
    ###################################################
    groupcat='berlind'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    f_red_1 = np.zeros((N_boots,len(bin_centers)))
    f_red_2 = np.zeros((N_boots,len(bin_centers)))
    f_red_3 = np.zeros((N_boots,len(bin_centers)))
    for boot in range(0,N_boots):
        #open catalogue
        catalogue_1 = group_catalogue+'_'+str(boot)
        f =  h5py.File(filepath_cat+'bootstraps/'+catalogue_1+'.hdf5', 'r') #open catalogue file
        GC = f.get(catalogue_1)

        centrals_ind   = np.where(GC['RANK']==1)[0]
        satellites_ind = np.where(GC['RANK']!=1)[0]
        centrals_bool   = (GC['RANK']==1)
        satellites_bool = (GC['RANK']!=1)
    
        #galaxy color
        color = GC['M_g,0.1']-GC['M_r,0.1']
        LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
        blue_ind  = np.where(color<LHS)[0] #indices of blue galaxies
        red_ind   = np.where(color>LHS)[0] #indicies of red galaxies
        blue_bool = (color<LHS) #indices of blue galaxies
        red_bool  = (color>LHS) #indicies of red galaxies

        selection = ((GC['MGROUP']>13.0) & (GC['MGROUP']<13.5))
        selection = satellites_bool & selection
        f_red_1[boot,:]  = f_prop(GC['RPROJ']/1000.0,bins,red_ind,blue_ind,selection)
        selection = ((GC['MGROUP']>13.5) & (GC['MGROUP']<14.0))
        selection = satellites_bool & selection
        f_red_2[boot,:]  = f_prop(GC['RPROJ']/1000.0,bins,red_ind,blue_ind,selection)
        selection = ((GC['MGROUP']>14.0) & (GC['MGROUP']<14.5))
        selection = satellites_bool & selection
        f_red_3[boot,:]  = f_prop(GC['RPROJ']/1000.0,bins,red_ind,blue_ind,selection)
        

    err_1 = np.nanstd(f_red_1, axis=0)
    f_red_1 = np.nanmean(f_red_1, axis=0)
    err_2 = np.nanstd(f_red_2, axis=0)
    f_red_2 = np.nanmean(f_red_2, axis=0)
    err_3 = np.nanstd(f_red_3, axis=0)
    f_red_3 = np.nanmean(f_red_3, axis=0)

    ax=axes[0]
    p0a=ax.errorbar(bin_centers, f_red_1, yerr=err_1, fmt='o',color='red', mec='none', ms=3)
    ax=axes[1]
    p1a=ax.errorbar(bin_centers, f_red_2, yerr=err_2, fmt='o',color='red', mec='none', ms=3)
    ax=axes[2]
    p2a=ax.errorbar(bin_centers, f_red_3, yerr=err_3, fmt='o',color='red', mec='none', ms=3)

    ###################################################
    groupcat='tinker'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    f_red_1 = np.zeros((N_boots,len(bin_centers)))
    f_red_2 = np.zeros((N_boots,len(bin_centers)))
    f_red_3 = np.zeros((N_boots,len(bin_centers)))
    for boot in range(0,N_boots):
        #open catalogue
        catalogue_1 = group_catalogue+'_'+str(boot)
        f =  h5py.File(filepath_cat+'bootstraps/'+catalogue_1+'.hdf5', 'r') #open catalogue file
        GC = f.get(catalogue_1)

        centrals_ind   = np.where(GC['RANK']==1)[0]
        satellites_ind = np.where(GC['RANK']!=1)[0]
        centrals_bool   = (GC['RANK']==1)
        satellites_bool = (GC['RANK']!=1)
    
        #galaxy color
        color = GC['M_g,0.1']-GC['M_r,0.1']
        LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
        blue_ind  = np.where(color<LHS)[0] #indices of blue galaxies
        red_ind   = np.where(color>LHS)[0] #indicies of red galaxies
        blue_bool = (color<LHS) #indices of blue galaxies
        red_bool  = (color>LHS) #indicies of red galaxies

        selection = ((GC['MGROUP']>13.0) & (GC['MGROUP']<13.5))
        selection = satellites_bool & selection
        f_red_1[boot,:]  = f_prop(GC['RPROJ']/1000.0,bins,red_ind,blue_ind,selection)
        selection = ((GC['MGROUP']>13.5) & (GC['MGROUP']<14.0))
        selection = satellites_bool & selection
        f_red_2[boot,:]  = f_prop(GC['RPROJ']/1000.0,bins,red_ind,blue_ind,selection)
        selection = ((GC['MGROUP']>14.0) & (GC['MGROUP']<14.5))
        selection = satellites_bool & selection
        f_red_3[boot,:]  = f_prop(GC['RPROJ']/1000.0,bins,red_ind,blue_ind,selection)
        

    err_1 = np.nanstd(f_red_1, axis=0)
    f_red_1 = np.nanmean(f_red_1, axis=0)
    err_2 = np.nanstd(f_red_2, axis=0)
    f_red_2 = np.nanmean(f_red_2, axis=0)
    err_3 = np.nanstd(f_red_3, axis=0)
    f_red_3 = np.nanmean(f_red_3, axis=0)

    ax=axes[0]
    p0b=ax.errorbar(bin_centers, f_red_1, yerr=err_1, fmt='o',color='green', mec='none', ms=3)
    ax=axes[1]
    p1b=ax.errorbar(bin_centers, f_red_2, yerr=err_2, fmt='o',color='green', mec='none', ms=3)
    ax=axes[2]
    p2b=ax.errorbar(bin_centers, f_red_3, yerr=err_3, fmt='o',color='green', mec='none', ms=3)

    ###################################################
    groupcat='yang'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    f_red_1 = np.zeros((N_boots,len(bin_centers)))
    f_red_2 = np.zeros((N_boots,len(bin_centers)))
    f_red_3 = np.zeros((N_boots,len(bin_centers)))
    for boot in range(0,N_boots):
        #open catalogue
        catalogue_1 = group_catalogue+'_'+str(boot)
        f =  h5py.File(filepath_cat+'bootstraps/'+catalogue_1+'.hdf5', 'r') #open catalogue file
        GC = f.get(catalogue_1)

        centrals_ind   = np.where(GC['RANK']==1)[0]
        satellites_ind = np.where(GC['RANK']!=1)[0]
        centrals_bool   = (GC['RANK']==1)
        satellites_bool = (GC['RANK']!=1)
    
        #galaxy color
        color = GC['M_g,0.1']-GC['M_r,0.1']
        LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
        blue_ind  = np.where(color<LHS)[0] #indices of blue galaxies
        red_ind   = np.where(color>LHS)[0] #indicies of red galaxies
        blue_bool = (color<LHS) #indices of blue galaxies
        red_bool  = (color>LHS) #indicies of red galaxies

        selection = ((GC['MGROUP']>13.0) & (GC['MGROUP']<13.5))
        selection = satellites_bool & selection
        f_red_1[boot,:]  = f_prop(GC['RPROJ']/1000.0,bins,red_ind,blue_ind,selection)
        selection = ((GC['MGROUP']>13.5) & (GC['MGROUP']<14.0))
        selection = satellites_bool & selection
        f_red_2[boot,:]  = f_prop(GC['RPROJ']/1000.0,bins,red_ind,blue_ind,selection)
        selection = ((GC['MGROUP']>14.0) & (GC['MGROUP']<14.5))
        selection = satellites_bool & selection
        f_red_3[boot,:]  = f_prop(GC['RPROJ']/1000.0,bins,red_ind,blue_ind,selection)
        

    err_1 = np.nanstd(f_red_1, axis=0)
    f_red_1 = np.nanmean(f_red_1, axis=0)
    err_2 = np.nanstd(f_red_2, axis=0)
    f_red_2 = np.nanmean(f_red_2, axis=0)
    err_3 = np.nanstd(f_red_3, axis=0)
    f_red_3 = np.nanmean(f_red_3, axis=0)

    ax=axes[0]
    p0c=ax.errorbar(bin_centers, f_red_1, yerr=err_1, fmt='o',color='blue', mec='none', ms=3)
    ax=axes[1]
    p1c=ax.errorbar(bin_centers, f_red_2, yerr=err_2, fmt='o',color='blue', mec='none', ms=3)
    ax=axes[2]
    p2c=ax.errorbar(bin_centers, f_red_3, yerr=err_3, fmt='o',color='blue', mec='none', ms=3)

    plt.show()

    
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
