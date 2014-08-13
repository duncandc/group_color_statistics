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

    bins = np.arange(11,15,0.2)
    bin_centers = (bins[:-1]+bins[1:])/2.0

    plt.figure()
    ###################################################
    groupcat='berlind'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    catalogue_1 = group_catalogue
    f =  h5py.File(filepath_cat+catalogue_1+'.hdf5', 'r') #open catalogue file
    GC = f.get(catalogue_1)

    centrals_ind   = np.where(GC['RANK']==1)[0]
    satellites_ind = np.where(GC['RANK']!=1)[0]

    result = np.histogram(GC['MGROUP'][satellites_ind],bins=bins)[0]
    p1,=plt.plot(bin_centers,np.log10(result),'-o', color='orange')
    ###################################################
    groupcat='tinker'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    catalogue_1 = group_catalogue
    f =  h5py.File(filepath_cat+catalogue_1+'.hdf5', 'r') #open catalogue file
    GC = f.get(catalogue_1)

    centrals_ind   = np.where(GC['RANK']==1)[0]
    satellites_ind = np.where(GC['RANK']!=1)[0]

    result = np.histogram(GC['MGROUP'][satellites_ind],bins=bins)[0]
    p2,=plt.plot(bin_centers,np.log10(result),'-o', color='green')
    ###################################################
    groupcat='yang'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    catalogue_1 = group_catalogue
    f =  h5py.File(filepath_cat+catalogue_1+'.hdf5', 'r') #open catalogue file
    GC = f.get(catalogue_1)

    centrals_ind   = np.where(GC['RANK']==1)[0]
    satellites_ind = np.where(GC['RANK']!=1)[0]

    result = np.histogram(GC['MGROUP'][satellites_ind],bins=bins)[0]
    p3,=plt.plot(bin_centers+np.log10(0.7),np.log10(result),'-o',color='cyan')

    
    plt.xlabel(r'$log(M)$')
    plt.ylabel(r'$log(N)$')
    plt.legend((p1,p2,p3),('berlind','tinker','yang'),loc='lower right')

    plt.show()

    bins = np.arange(11,15,0.2)
    bin_centers = (bins[:-1]+bins[1:])/2.0

    plt.figure()
    ###################################################
    groupcat='berlind'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    catalogue_1 = group_catalogue
    f =  h5py.File(filepath_cat+catalogue_1+'.hdf5', 'r') #open catalogue file
    GC = f.get(catalogue_1)

    centrals_ind   = np.where(GC['RANK']==1)[0]
    satellites_ind = np.where(GC['RANK']!=1)[0]

    result = np.histogram(GC['MGROUP'][centrals_ind],bins=bins)[0]
    p1,=plt.plot(bin_centers,np.log10(result),'-o', color='orange')
    ###################################################
    groupcat='tinker'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    catalogue_1 = group_catalogue
    f =  h5py.File(filepath_cat+catalogue_1+'.hdf5', 'r') #open catalogue file
    GC = f.get(catalogue_1)

    centrals_ind   = np.where(GC['RANK']==1)[0]
    satellites_ind = np.where(GC['RANK']!=1)[0]

    result = np.histogram(GC['MGROUP'][centrals_ind],bins=bins)[0]
    p2,=plt.plot(bin_centers,np.log10(result),'o', color='green')
    ###################################################
    groupcat='yang'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    catalogue_1 = group_catalogue
    f =  h5py.File(filepath_cat+catalogue_1+'.hdf5', 'r') #open catalogue file
    GC = f.get(catalogue_1)

    print GC.dtype.names

    centrals_ind   = np.where(GC['RANK']==1)[0]
    satellites_ind = np.where(GC['RANK']!=1)[0]

    halo_centrals_ind = np.where(GC['HALO_RANK']==1)[0]

    result = np.histogram(GC['MGROUP'][centrals_ind],bins=bins)[0]
    p3,=plt.plot(bin_centers,np.log10(result),'o',color='cyan')
    result = np.histogram(GC['HALO_M'][halo_centrals_ind],bins=bins)[0]
    p4,=plt.plot(bin_centers,np.log10(result),'-',color='red')

    
    plt.xlabel(r'$log(M)$')
    plt.ylabel(r'$log(N)$')
    plt.legend((p1,p2,p3,p4),('berlind FoF groups','tinker SO groups','yang SO groups','Bolshoi mock haloes'),loc='upper right')

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
