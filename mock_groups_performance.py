#!/usr/bin/env python

#Duncan Campbell
#October 15, 2014
#Yale University
#Plot classic measures of group finder performance.

#load packages
import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu
from astropy import cosmology
import sys

def main():

    catalogue = sys.argv[1]
    
    plotpath = cu.get_plot_path() + 'analysis/groupcats/'


    ######################################################################################
    #Satellite and central purity and completeness
    ######################################################################################
    #define outputs
    filename1 = catalogue+'_performance_M'
    filename2 = catalogue+'_performance_L'

    #First make plots as a function of galaxy Mgroup or Mhalo
    #set up plots
    fig1, axes = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, figsize=(6.95, 3.3))
    fig1.subplots_adjust(hspace=0, wspace=0.05)
    fig1.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.9)
    axes = axes.flatten()

    axes[0].set_xlim(12,15)
    axes[0].set_ylim(0.2,1)
    axes[0].set_xticks([12.5,13,13.5,14,14.5])
    axes[0].set_xlabel(r'$log(M/[M_{\odot}h^{-1}])$')
    axes[0].set_ylabel(r'$C,P$')
    axes[0].set_title('Berlind FoF groups')
    axes[1].set_title('Tinker SO groups')
    axes[1].set_xlabel(r'$log(M/M_{\odot}h)$')
    axes[2].set_title('Yang SO groups')
    axes[2].set_xlabel(r'$log(M/M_{\odot}h)$')
    
    #open the Berlind group catalogue
    groupcat='berlind'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    #open catalogue
    f =  h5py.File(filepath_cat+group_catalogue+'.hdf5', 'r') #open catalogue file
    GC = f.get(group_catalogue)

    #define mass bins
    bins = np.arange(12,15,0.2)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    result1 = sp_in_mass_bins(GC,bins)
    result2 = sc_in_mass_bins(GC,bins)
    result3 = cp_in_mass_bins(GC,bins)
    result4 = cc_in_mass_bins(GC,bins)

    ax=axes[0]
    p1,=ax.plot(bin_centers,result1, color='green',ls='-')
    p2,=ax.plot(bin_centers,result2, color='green',ls='--')
    p3,=ax.plot(bin_centers,result3, color='orange',ls='-')
    p4,=ax.plot(bin_centers,result4, color='orange',ls='--')
    ax.legend((p1,p2),('satellite purity','satellite \ncompleteness'),loc='lower right', fontsize=10, numpoints=1, frameon=False)
    
    #open the Tinker group catalogue
    groupcat='tinker'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    #open catalogue
    f =  h5py.File(filepath_cat+group_catalogue+'.hdf5', 'r') #open catalogue file
    GC = f.get(group_catalogue)

    result1 = sp_in_mass_bins(GC,bins)
    result2 = sc_in_mass_bins(GC,bins)
    result3 = cp_in_mass_bins(GC,bins)
    result4 = cc_in_mass_bins(GC,bins)

    ax=axes[1]
    p1,=ax.plot(bin_centers,result1, color='green', ls='-')
    p2,=ax.plot(bin_centers,result2, color='green',ls='--')
    p3,=ax.plot(bin_centers,result3, color='orange',ls='-')
    p4,=ax.plot(bin_centers,result4, color='orange',ls='--')
    ax.legend((p3,p4),('central purity','central \ncompleteness'),loc='lower right', fontsize=10, numpoints=1, frameon=False)

    #open the Yang group catalogue
    groupcat='yang'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    #open catalogue
    f =  h5py.File(filepath_cat+group_catalogue+'.hdf5', 'r') #open catalogue file
    GC = f.get(group_catalogue)

    result1 = sp_in_mass_bins(GC,bins)
    result2 = sc_in_mass_bins(GC,bins)
    result3 = cp_in_mass_bins(GC,bins)
    result4 = cc_in_mass_bins(GC,bins)

    ax=axes[2]
    ax.plot(bin_centers,result1, color='green', ls='-')
    ax.plot(bin_centers,result2, color='green',ls='--')
    ax.plot(bin_centers,result3, color='orange',ls='-')
    ax.plot(bin_centers,result4, color='orange',ls='--')
    plt.show(block=False)

    #save plot
    print plotpath+filename1+'.pdf'
    fig1.savefig(plotpath+filename1+'.pdf')

    #Second make plots as a function of galaxy Luminosity
    #define luminosity bins
    bins = np.arange(9.4,10.8,0.1)
    bin_centers = (bins[:-1]+bins[1:])/2.0

    #set up plots
    fig2, axes = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, figsize=(6.95, 3.3))
    fig2.subplots_adjust(hspace=0, wspace=0.05)
    fig2.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.9)
    axes = axes.flatten()

    axes[0].set_xlim(9,11)
    axes[0].set_ylim(0.2,1)
    axes[0].set_xticks([9.5,10,10.5])
    axes[0].set_xlabel(r'$log(L/L_{\odot})$')
    axes[0].set_ylabel(r'$C,P$')
    axes[0].set_title('Berlind FoF groups')
    axes[1].set_title('Tinker SO groups')
    axes[1].set_xlabel(r'$log(L/L_{\odot})$')
    axes[2].set_title('Yang SO groups')
    axes[2].set_xlabel(r'$log(L/L_{\odot})$')
    
    #open the Berlind group catalogue
    groupcat='berlind'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    #open catalogue
    f =  h5py.File(filepath_cat+group_catalogue+'.hdf5', 'r') #open catalogue file
    GC = f.get(group_catalogue)

    result1 = sp_in_L_bins(GC,bins)
    result2 = sc_in_L_bins(GC,bins)
    result3 = cp_in_L_bins(GC,bins)
    result4 = cc_in_L_bins(GC,bins)

    ax=axes[0]
    p1,=ax.plot(bin_centers,result1, color='green',ls='-')
    p2,=ax.plot(bin_centers,result2, color='green',ls='--')
    p3,=ax.plot(bin_centers,result3, color='orange',ls='-')
    p4,=ax.plot(bin_centers,result4, color='orange',ls='--')
    ax.legend((p1,p2),('satellite purity','satellite \ncompleteness'),loc='lower left', fontsize=10, numpoints=1, frameon=False)
    
    #open the Tinker group catalogue
    groupcat='tinker'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    #open catalogue
    f =  h5py.File(filepath_cat+group_catalogue+'.hdf5', 'r') #open catalogue file
    GC = f.get(group_catalogue)

    result1 = sp_in_L_bins(GC,bins)
    result2 = sc_in_L_bins(GC,bins)
    result3 = cp_in_L_bins(GC,bins)
    result4 = cc_in_L_bins(GC,bins)

    ax=axes[1]
    p1,=ax.plot(bin_centers,result1, color='green', ls='-')
    p2,=ax.plot(bin_centers,result2, color='green',ls='--')
    p3,=ax.plot(bin_centers,result3, color='orange',ls='-')
    p4,=ax.plot(bin_centers,result4, color='orange',ls='--')
    ax.legend((p3,p4),('central purity','central \ncompleteness'),loc='lower left', fontsize=10, numpoints=1, frameon=False)

    #open the Yang group catalogue
    groupcat='yang'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    #open catalogue
    f =  h5py.File(filepath_cat+group_catalogue+'.hdf5', 'r') #open catalogue file
    GC = f.get(group_catalogue)

    result1 = sp_in_L_bins(GC,bins)
    result2 = sc_in_L_bins(GC,bins)
    result3 = cp_in_L_bins(GC,bins)
    result4 = cc_in_L_bins(GC,bins)

    ax=axes[2]
    ax.plot(bin_centers,result1, color='green', ls='-')
    ax.plot(bin_centers,result2, color='green',ls='--')
    ax.plot(bin_centers,result3, color='orange',ls='-')
    ax.plot(bin_centers,result4, color='orange',ls='--')
    plt.show(block=False)

    #save plots
    print plotpath+filename2+'.pdf'
    fig2.savefig(plotpath+filename2+'.pdf')
    
    ######################################################################################
    ##### look at fraction of satellites that are more luminous than the central
    ######################################################################################
    
    #open ideal group catalogue
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/ideal_groups/'
    print 'opening mock catalogue:', catalogue+'_groups.hdf5'
    #open catalogue
    f1 = h5py.File(filepath_mock+catalogue+'_groups.hdf5', 'r') #open catalogue file
    GC = f1.get(catalogue+'_groups')
    
    #define mass bins
    bins = np.arange(10,15,0.2)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    
    #which are true centrals and satellites?
    true_centrals   = np.where(GC['HALO_RANK']==0)[0]
    true_satellites = np.where(GC['HALO_RANK']==1)[0]
    true_centrals_bool = (GC['HALO_RANK']==0)
    
    #which are identified centrals and satellites?
    identified_centrals   = np.where(GC['RANK']==0)[0]
    identified_satellites = np.where(GC['RANK']!=0)[0]
    identified_centrals_bool = (GC['RANK']==0)
    
    #which identified centrals are not true centrals?
    #fu = np.setdiff1d(true_centrals,identified_centrals)
    fu = np.in1d(true_centrals,identified_centrals)
    fu = true_centrals[fu]
    
    #which centrals are identified correctly?
    good = np.intersect1d(true_centrals,identified_centrals)
    
    f_fu = f_prop(GC['HALO_M'],bins,fu,true_centrals,true_centrals_bool)
    
    plt.figure()
    plt.plot(bin_centers, f_fu)
    plt.ylim([0,1])
    plt.show(block=True)
    

def f_prop(prop,prop_bins,group_1,group_2,mask):
    #returns the fraction of galaxies in group_1 as a function
    #of a property. 
    import numpy as np
    f = np.zeros(len(prop_bins)-1)

    group_1_mask = mask[group_1]
    group_2_mask = mask[group_2]

    group_1 = group_1[group_1_mask].copy()
    group_2 = group_2[group_2_mask].copy()

    result = np.digitize(prop,bins=prop_bins)
    for i in range(0,len(prop_bins)-1):
        ind = np.where(result==i+1)[0]
        group_1_gal = ind[np.in1d(ind,group_1)]
        group_2_gal = ind[np.in1d(ind,group_2)]
        N=len(group_1_gal)+len(group_2_gal)
        if N>0:
            f[i] = float(len(group_1_gal))/(float(N))
        else: f[i]=0.0

    return f

def sp_in_mass_bins(GC,bins):
    x=np.zeros(len(bins)-1)
    result = np.digitize(GC['HALO_M'],bins=bins)
    for i in range(0,len(x)):
        ind = (result==i+1)
        x[i]=sat_purity(GC[ind])
    return x

def sc_in_mass_bins(GC,bins):
    x=np.zeros(len(bins)-1)
    result = np.digitize(GC['HALO_M'],bins=bins)
    for i in range(0,len(x)):
        ind = (result==i+1)
        x[i]=sat_comp(GC[ind])
    return x

def cp_in_mass_bins(GC,bins):
    x=np.zeros(len(bins)-1)
    result = np.digitize(GC['HALO_M'],bins=bins)
    for i in range(0,len(x)):
        ind = (result==i+1)
        x[i]=cen_purity(GC[ind])
    return x

def cc_in_mass_bins(GC,bins):
    x=np.zeros(len(bins)-1)
    result = np.digitize(GC['HALO_M'],bins=bins)
    for i in range(0,len(x)):
        ind = (result==i+1)
        x[i]=cen_comp(GC[ind])
    return x

def sp_in_L_bins(GC,bins):
    x=np.zeros(len(bins)-1)
    result = np.digitize(solar_lum(GC['M_r,0.1'],4.64),bins=bins)
    for i in range(0,len(x)):
        ind = (result==i+1)
        x[i]=sat_purity(GC[ind])
    return x

def sc_in_L_bins(GC,bins):
    x=np.zeros(len(bins)-1)
    result = np.digitize(solar_lum(GC['M_r,0.1'],4.64),bins=bins)
    for i in range(0,len(x)):
        ind = (result==i+1)
        x[i]=sat_comp(GC[ind])
    return x

def cp_in_L_bins(GC,bins):
    x=np.zeros(len(bins)-1)
    result = np.digitize(solar_lum(GC['M_r,0.1'],4.64),bins=bins)
    for i in range(0,len(x)):
        ind = (result==i+1)
        x[i]=cen_purity(GC[ind])
    return x

def cc_in_L_bins(GC,bins):
    x=np.zeros(len(bins)-1)
    result = np.digitize(solar_lum(GC['M_r,0.1'],4.64),bins=bins)
    for i in range(0,len(x)):
        ind = (result==i+1)
        x[i]=cen_comp(GC[ind])
    return x

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

def sat_purity(GC):
    Nss=float(len(np.where((GC['RANK']!=1) & (GC['HALO_RANK']!=1))[0]))
    Nsc=float(len(np.where((GC['RANK']!=1) & (GC['HALO_RANK']==1))[0]))
    x=Nss/(Nss+Nsc)
    return x

def cen_purity(GC):
    Ncc=float(len(np.where((GC['RANK']==1) & (GC['HALO_RANK']==1))[0]))
    Ncs=float(len(np.where((GC['RANK']==1) & (GC['HALO_RANK']!=1))[0]))
    y=Ncc/(Ncc+Ncs)
    return y

def sat_comp(GC):
    Nss=float(len(np.where((GC['RANK']!=1) & (GC['HALO_RANK']!=1))[0]))
    Ncs=float(len(np.where((GC['RANK']==1) & (GC['HALO_RANK']!=1))[0]))
    x=Nss/(Nss+Ncs)
    return x

def cen_comp(GC):
    Ncc=float(len(np.where((GC['RANK']==1) & (GC['HALO_RANK']==1))[0]))
    Nsc=float(len(np.where((GC['RANK']!=1) & (GC['HALO_RANK']==1))[0]))
    y=Ncc/(Ncc+Nsc)
    return y

def sat_bias(GC):
    Nsc=float(len(np.where((GC['RANK']!=1) & (GC['HALO_RANK']==1))[0]))
    Ncs=float(len(np.where((GC['RANK']==1) & (GC['HALO_RANK']!=1))[0]))
    b=Nsc/(Ncs)
    return b

def solar_lum(M,Msol):
    L = ((Msol-M)/2.5)
    return L

if __name__ == '__main__':
    main() 
