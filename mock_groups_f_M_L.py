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
    LHS   = 0.7 - 0.032*(mock['M_r,0.1']+16.5) #Weinmann 2006
    blue  = np.where(color<LHS)[0] #indices of blue galaxies
    red   = np.where(color>LHS)[0] #indicies of red galaxies
    blue_bool  = np.array(color<LHS) #indices of blue galaxies
    red_bool   = np.array(color>LHS) #indicies of red galaxies

    lumbins = np.arange(9.6,11.0,0.2)
    lumbin_centers = (lumbins[:-1]+lumbins[1:])/2.0
    print len(lumbin_centers)
    print lumbin_centers
    inds = np.digitize(solar_lum(mock['M_r,0.1'],4.64),bins=lumbins)
    
    bininds=np.array(inds==0)
    result1 = np.histogram(mock['M_host'][(centrals & bininds)],bins=bins)[0]
    central_mean_mass_0 = np.mean(mock['M_host'][(centrals & bininds)])
    f_red_cen_0 = float(len(mock['M_host'][red_bool&(centrals & bininds)]))/float(len(mock['M_host'][(centrals & bininds)]))
    f_red_sat_0 = f_prop(mock['M_host'],bins,red,blue,(satellites & bininds))

    bininds=np.array(inds==2)
    result1 = np.histogram(mock['M_host'][(centrals & bininds)],bins=bins)[0]
    central_mean_mass_2 = np.mean(mock['M_host'][(centrals & bininds)])
    f_red_cen_2 = float(len(mock['M_host'][red_bool&(centrals & bininds)]))/float(len(mock['M_host'][(centrals & bininds)]))
    f_red_sat_2 = f_prop(mock['M_host'],bins,red,blue,(satellites & bininds))

    bininds=np.array(inds==4)
    result1 = np.histogram(mock['M_host'][(centrals & bininds)],bins=bins)[0]
    central_mean_mass_4 = np.mean(mock['M_host'][(centrals & bininds)])
    f_red_cen_4 = float(len(mock['M_host'][red_bool&(centrals & bininds)]))/float(len(mock['M_host'][(centrals & bininds)]))
    f_red_sat_4 = f_prop(mock['M_host'],bins,red,blue,(satellites & bininds))

    bininds=np.array(inds==6)
    result1 = np.histogram(mock['M_host'][(centrals & bininds)],bins=bins)[0]
    central_mean_mass_6 = np.mean(mock['M_host'][(centrals & bininds)])
    f_red_cen_6 = float(len(mock['M_host'][red_bool&(centrals & bininds)]))/float(len(mock['M_host'][(centrals & bininds)]))
    f_red_sat_6 = f_prop(mock['M_host'],bins,red,blue,(satellites & bininds))
    
    plt.figure()
    plt.plot(central_mean_mass_0,f_red_cen_0,'o',color='blue')
    p0,=plt.plot(bin_centers[bin_centers>central_mean_mass_0],f_red_sat_0[bin_centers>central_mean_mass_0],color='blue')

    plt.plot(central_mean_mass_2,f_red_cen_2,'o',color='green')
    p2,=plt.plot(bin_centers[bin_centers>central_mean_mass_2],f_red_sat_2[bin_centers>central_mean_mass_2],color='green')

    plt.plot(central_mean_mass_4,f_red_cen_4,'o',color='red')
    p4,=plt.plot(bin_centers[bin_centers>central_mean_mass_4],f_red_sat_4[bin_centers>central_mean_mass_4],color='red')

    plt.plot(central_mean_mass_6,f_red_cen_6,'o',color='black')
    p6,=plt.plot(bin_centers[bin_centers>central_mean_mass_6],f_red_sat_6[bin_centers>central_mean_mass_6],color='black')

    plt.legend((p0,p2,p4,p6),(lumbin_centers[0],lumbin_centers[2],lumbin_centers[4],lumbin_centers[6]),loc='lower right')
    plt.xlabel(r'$\log(M)$')
    plt.ylabel(r'$f_{red}$')

    ###################################################
    N_boots=50
    groupcat='yang'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    f_red_cen_0 = np.zeros((N_boots,))
    f_red_sat_0 = np.zeros((N_boots,len(bin_centers)))
    f_red_cen_2 = np.zeros((N_boots,))
    f_red_sat_2 = np.zeros((N_boots,len(bin_centers)))
    f_red_cen_4 = np.zeros((N_boots,))
    f_red_sat_4 = np.zeros((N_boots,len(bin_centers)))
    f_red_cen_6 = np.zeros((N_boots,))
    f_red_sat_6 = np.zeros((N_boots,len(bin_centers)))
    central_mean_mass_0 = np.zeros((N_boots,))
    central_mean_mass_2 = np.zeros((N_boots,))
    central_mean_mass_4 = np.zeros((N_boots,))
    central_mean_mass_6 = np.zeros((N_boots,))
    for boot in range(0,N_boots):
        #open catalogue
        catalogue_1 = group_catalogue+'_'+str(boot)
        f =  h5py.File(filepath_cat+'bootstraps/'+catalogue_1+'.hdf5', 'r') #open catalogue file
        GC = f.get(catalogue_1)

        centrals_ind   = np.where(GC['RANK']==1)[0]
        satellites_ind = np.where(GC['RANK']!=1)[0]
        centrals   = np.array(GC['RANK']==1)
        satellites = np.array(GC['RANK']!=1)
    
        #galaxy color
        color = GC['M_g,0.1']-GC['M_r,0.1']
        LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
        blue  = np.where(color<LHS)[0] #indices of blue galaxies
        red   = np.where(color>LHS)[0] #indicies of red galaxies
        blue_bool = (color<LHS) #indices of blue galaxies
        red_bool  = (color>LHS) #indicies of red galaxies
    
        inds = np.digitize(solar_lum(GC['M_r,0.1'],4.64),bins=lumbins)
    
        bininds=np.array(inds==0)
        result1 = np.histogram(GC['MGROUP'][(centrals & bininds)],bins=bins)[0]
        central_mean_mass_0[boot] = np.mean(GC['MGROUP'][(centrals & bininds)])
        f_red_cen_0[boot] = float(len(GC['MGROUP'][red_bool&(centrals & bininds)]))/float(len(GC['MGROUP'][(centrals & bininds)]))
        f_red_sat_0[boot,:] = f_prop(GC['MGROUP'],bins,red,blue,(satellites & bininds))

        bininds=np.array(inds==2)
        result1 = np.histogram(GC['MGROUP'][(centrals & bininds)],bins=bins)[0]
        central_mean_mass_2[boot] = np.mean(GC['MGROUP'][(centrals & bininds)])
        f_red_cen_2[boot] = float(len(GC['MGROUP'][red_bool&(centrals & bininds)]))/float(len(GC['MGROUP'][(centrals & bininds)]))
        f_red_sat_2[boot:,] = f_prop(GC['MGROUP'],bins,red,blue,(satellites & bininds))

        bininds=np.array(inds==4)
        result1 = np.histogram(GC['MGROUP'][(centrals & bininds)],bins=bins)[0]
        central_mean_mass_4[boot] = np.mean(GC['MGROUP'][(centrals & bininds)])
        f_red_cen_4[boot] = float(len(GC['MGROUP'][red_bool&(centrals & bininds)]))/float(len(GC['MGROUP'][(centrals & bininds)]))
        f_red_sat_4[boot,:] = f_prop(GC['MGROUP'],bins,red,blue,(satellites & bininds))

        bininds=np.array(inds==6)
        result1 = np.histogram(GC['MGROUP'][(centrals & bininds)],bins=bins)[0]
        central_mean_mass_6[boot] = np.mean(GC['MGROUP'][(centrals & bininds)])
        f_red_cen_6[boot] = float(len(GC['MGROUP'][red_bool&(centrals & bininds)]))/float(len(GC['MGROUP'][(centrals & bininds)]))
        f_red_sat_6[boot,:] = f_prop(GC['MGROUP'],bins,red,blue,(satellites & bininds))
    
    err_cen_0 = np.nanstd(f_red_cen_0, axis=0)
    err_sat_0 = np.nanstd(f_red_sat_0, axis=0)
    f_red_cen_0 = np.nanmean(f_red_cen_0, axis=0)
    f_red_sat_0 = np.nanmean(f_red_sat_0, axis=0)
    central_mean_mass_err_0 = np.nanstd(central_mean_mass_0, axis=0)
    central_mean_mass_0 = np.nanmean(central_mean_mass_0, axis=0)
    err_cen_2 = np.nanstd(f_red_cen_2, axis=0)
    err_sat_2 = np.nanstd(f_red_sat_2, axis=0)
    f_red_cen_2 = np.nanmean(f_red_cen_2, axis=0)
    f_red_sat_2 = np.nanmean(f_red_sat_2, axis=0)
    central_mean_mass_err_2 = np.nanstd(central_mean_mass_2, axis=0)
    central_mean_mass_2 = np.nanmean(central_mean_mass_2, axis=0)
    err_cen_4 = np.nanstd(f_red_cen_4, axis=0)
    err_sat_4 = np.nanstd(f_red_sat_4, axis=0)
    f_red_cen_4 = np.nanmean(f_red_cen_4, axis=0)
    f_red_sat_4 = np.nanmean(f_red_sat_4, axis=0)
    central_mean_mass_err_4 = np.nanstd(central_mean_mass_4, axis=0)
    central_mean_mass_4 = np.nanmean(central_mean_mass_4, axis=0)
    err_cen_6 = np.nanstd(f_red_cen_6, axis=0)
    err_sat_6 = np.nanstd(f_red_sat_6, axis=0)
    f_red_cen_6 = np.nanmean(f_red_cen_6, axis=0)
    f_red_sat_6 = np.nanmean(f_red_sat_6, axis=0)
    central_mean_mass_err_6 = np.nanstd(central_mean_mass_6, axis=0)
    central_mean_mass_6 = np.nanmean(central_mean_mass_6, axis=0)

    plt.errorbar(central_mean_mass_0,f_red_cen_0,yerr=err_cen_0,fmt='o',color='blue')
    p0=plt.errorbar(bin_centers[bin_centers>central_mean_mass_0],f_red_sat_0[bin_centers>central_mean_mass_0],yerr=err_sat_0[bin_centers>central_mean_mass_0],\
                    fmt='o',color='blue', alpha=0.5)

    plt.errorbar(central_mean_mass_2,f_red_cen_2,yerr=err_cen_2,fmt='o',color='green')
    p2=plt.errorbar(bin_centers[bin_centers>central_mean_mass_2],f_red_sat_2[bin_centers>central_mean_mass_2],yerr=err_sat_2[bin_centers>central_mean_mass_2],\
                    fmt='o',color='green',alpha=0.5)

    plt.errorbar(central_mean_mass_4,f_red_cen_4,yerr=err_cen_4,fmt='o',color='red')
    p4=plt.errorbar(bin_centers[bin_centers>central_mean_mass_4],f_red_sat_4[bin_centers>central_mean_mass_4],yerr=err_sat_4[bin_centers>central_mean_mass_4],\
                    fmt='o',color='red',alpha=0.5)

    plt.errorbar(central_mean_mass_6,f_red_cen_6,yerr=err_cen_6,fmt='o',color='black')
    p6=plt.errorbar(bin_centers[bin_centers>central_mean_mass_6],f_red_sat_6[bin_centers>central_mean_mass_6],yerr=err_sat_6[bin_centers>central_mean_mass_6],\
                    fmt='o',color='black',alpha=0.5)

    plt.ylim([0.2,1])
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

def solar_lum(M,Msol):
    L = ((Msol-M)/2.5)
    return L

if __name__ == '__main__':
    main() 
