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

    if len(sys.argv)>1: catalogue = sys.argv[1]
    else: catalogue = 'Mr19_age_distribution_matching_mock'
    
    plotpath = cu.get_plot_path() + 'analysis/groupcats/'

    ######################################################################################
    #Satellite and central purity and completeness
    ######################################################################################
    #define outputs
    filename1 = catalogue+'_performance_M'
    filename2 = catalogue+'_performance_L'
    filename3 = catalogue+'_performance_Tsat'

    #First make plots as a function of galaxy Mgroup and Mhalo
    #set up plots
    fig1, axes = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, figsize=(6.95, 3.3))
    fig1.subplots_adjust(hspace=0, wspace=0.05)
    fig1.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.9)
    axes = axes.flatten()
    axes[0].set_xlim(10**12,10**15)
    axes[0].set_ylim(0.0,1.0)
    axes[0].set_xscale('log')
    axes[0].set_xlabel(r'$M$ $[M_{\odot}h^{-1}]$')
    axes[0].set_ylabel(r'$C,P$')
    axes[0].set_title('Berlind et al. groups')
    axes[1].set_title('Tinker et al. groups')
    axes[1].set_xlabel(r'$M$ $[M_{\odot}h^{-1}]$')
    axes[1].set_xscale('log')
    axes[2].set_title('Yang et al. groups')
    axes[2].set_xlabel(r'$M$ $[M_{\odot}h^{-1}]$')
    axes[2].set_xscale('log')
    axes[2].set_xticks([10**12,10**13,10**14,10**15])
    axes[2].set_xticklabels([r'$10^{12}$',r'$10^{13}$',r'$10^{14}$',r' '])
    
    #define mass bins
    bins = np.arange(12,15,0.2)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    
    #open the Berlind group catalogue
    groupcat='berlind'
    filepath_cat=get_gc_path(groupcat)
    group_catalogue=get_gc_name(groupcat,catalogue)
    #open catalogue
    f =  h5py.File(filepath_cat+group_catalogue+'.hdf5', 'r') #open catalogue file
    GC = f.get(group_catalogue)

    result1 = sp_in_mass_bins(GC,bins)
    result2 = sc_in_mass_bins(GC,bins)
    result3 = cp_in_mass_bins(GC,bins)
    result4 = cc_in_mass_bins(GC,bins)
    result5_berlind = Tsat_in_mass_bins(GC,bins)
    result5_berlind_global = sat_bias(GC)

    ax=axes[0]
    p1,=ax.plot(10**bin_centers,result1, color='green',ls='-')
    p2,=ax.plot(10**bin_centers,result2, color='green',ls='--')
    p3,=ax.plot(10**bin_centers,result3, color='orange',ls='-')
    p4,=ax.plot(10**bin_centers,result4, color='orange',ls='--')
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
    result5_tinker = Tsat_in_mass_bins(GC,bins)
    result5_tinker_global = sat_bias(GC)

    ax=axes[1]
    p1,=ax.plot(10**bin_centers,result1, color='green', ls='-')
    p2,=ax.plot(10**bin_centers,result2, color='green',ls='--')
    p3,=ax.plot(10**bin_centers,result3, color='orange',ls='-')
    p4,=ax.plot(10**bin_centers,result4, color='orange',ls='--')
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
    result5_yang = Tsat_in_mass_bins(GC,bins)
    result5_yang_global = sat_bias(GC)

    ax=axes[2]
    ax.plot(10**bin_centers,result1, color='green', ls='-')
    ax.plot(10**bin_centers,result2, color='green',ls='--')
    ax.plot(10**bin_centers,result3, color='orange',ls='-')
    ax.plot(10**bin_centers,result4, color='orange',ls='--')
    
    '''
    #open ideal group catalogue
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/ideal_groups/'
    print 'opening mock catalogue:', catalogue+'_groups.hdf5'
    #open catalogue
    f1 = h5py.File(filepath_mock+catalogue+'_groups.hdf5', 'r') #open catalogue file
    GC = f1.get(catalogue+'_groups')
    GC = np.array(GC)
    
    #look at purity and completeness for perfect groups
    satellites = (GC['RANK']>0)
    GC['RANK'][satellites]=2
    centrals = (GC['RANK']==0)
    GC['RANK'][centrals]=1
    
    satellites = (GC['HALO_RANK']!=0)
    GC['HALO_RANK'][satellites]=2
    centrals = (GC['HALO_RANK']==0)
    GC['HALO_RANK'][centrals]=1
    
    #define mass bins
    bins = np.arange(12,15,0.2)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    
    result1 = sp_in_mass_bins(GC,bins)
    result2 = sc_in_mass_bins(GC,bins)
    result3 = cp_in_mass_bins(GC,bins)
    result4 = cc_in_mass_bins(GC,bins)
    result5 = Tsat_in_mass_bins(GC,bins)
    
    ax=axes[0]
    ax.plot(10**bin_centers,result2, color='green', ls='-')
    ax.plot(10**bin_centers,result4, color='orange',ls='-')
    '''
    
    plt.show()

    #save plot
    print plotpath+filename1+'.pdf'
    fig1.savefig(plotpath+filename1+'.pdf')

    #purity and completeness as a function of galaxy Luminosity
    #set up plots
    fig2, axes = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, figsize=(6.95, 3.3))
    fig2.subplots_adjust(hspace=0, wspace=0.05)
    fig2.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.9)
    axes = axes.flatten()
    axes[0].set_xlim(9,11)
    axes[0].set_ylim(0.0,1)
    axes[0].set_xticks([9.5,10,10.5])
    axes[0].set_xlabel(r'$log(L/L_{\odot})$')
    axes[0].set_ylabel(r'$C,P$')
    axes[0].set_title('Berlind et al. groups')
    axes[1].set_title('Tinker et al. groups')
    axes[1].set_xlabel(r'$log(L/L_{\odot})$')
    axes[2].set_title('Yang et al. groups')
    axes[2].set_xlabel(r'$log(L/L_{\odot})$')
    
    #define luminosity bins
    bins = np.arange(9.4,10.8,0.1)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    
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
    
    
    #plot the satellite transition  factor as a function of Mgroup and Mhalo
    
    #define mass bins
    bins = np.arange(12,15,0.2)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    fig = plt.figure(figsize=(3.3,3.3))
    ax = fig.add_subplot(1,1,1)
    fig.subplots_adjust(left=0.2, right=0.9, bottom=0.2,top=0.9)
    p1, = ax.plot(10.0**bin_centers, result5_berlind, '--', color='cyan')
    ax.plot([1.2*10.0**bin_centers[-1]], [result5_berlind_global], 'o', color='cyan', mec='none')
    p2, = ax.plot(10.0**bin_centers, result5_tinker, '-', color='grey')
    ax.plot([1.2*10.0**bin_centers[-1]], [result5_tinker_global], 'o', color='grey', mec='none')
    p3, = ax.plot(10.0**bin_centers, result5_yang, ':', color='purple')
    ax.plot([1.2*10.0**bin_centers[-1]], [result5_yang_global], 'o', color='purple', mec='none')
    ax.legend((p1,p2,p3),('Berlind et al. groups','Tinker et al. groups','Yang et al. groups'), frameon=False, fontsize=8)
    ax.plot([10.0**11,10**16],[1,1], '--', color='black')
    ax.set_xlabel(r'$M$ $[M_{\odot}h^{-1}]$')
    ax.set_ylabel(r'$T_{\rm sat}$')
    ax.set_xscale('log')
    ax.set_ylim([0.5,4])
    ax.set_xlim([10**12,10**15])
    plt.show()
    fig.savefig(plotpath+filename3+'.pdf')


    ######################################################################################
    ##### look at fraction of satellites that are more luminous than the central
    ######################################################################################
    
    filename1 = catalogue+'_performance_fu_censat_Mh'
    filename2 = catalogue+'_performance_fu_censat_Mg'
    filename3 = catalogue+'_performance_fuMg'
    filename4 = catalogue+'_performance_fuMh'
    filename5 = catalogue+'_performance_fuMhMg'
    
    #open ideal group catalogue
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/ideal_groups/'
    print 'opening mock catalogue:', catalogue+'_groups.hdf5'
    #open catalogue
    f1 = h5py.File(filepath_mock+catalogue+'_groups.hdf5', 'r') #open catalogue file
    GC = f1.get(catalogue+'_groups')
    GC = np.array(GC)

    #define mass bins
    bins = np.arange(10,15,0.2)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    
    #which are true centrals and satellites?
    true_centrals   = np.where(GC['HALO_RANK']==0)[0]
    true_satellites = np.where(GC['HALO_RANK']==1)[0]
    true_centrals_bool = (GC['HALO_RANK']==0)
    true_satellites_bool = (GC['HALO_RANK']==1)
    
    #which are identified centrals and satellites?
    identified_centrals   = np.where(GC['RANK']==0)[0]
    identified_satellites = np.where(GC['RANK']!=0)[0]
    identified_centrals_bool = (GC['RANK']==0)
    identified_satellites_bool = (GC['RANK']!=0)
    
    #which identified centrals are not true centrals?
    fu_cen = (np.in1d(true_centrals,identified_centrals)==False)
    nfu_cen = (np.in1d(true_centrals,identified_centrals))
    fu_cen = true_centrals[fu_cen]
    
    #which identified centrals are not true centrals?
    fu_icen = (np.in1d(identified_centrals,true_satellites))
    nfu_icen = (np.in1d(identified_centrals,true_centrals))
    fu_icen = identified_centrals[fu_icen]
    
    #which identified centrals are actually satellites?
    fu_sat = (np.in1d(true_satellites,identified_satellites)==False)
    nfu_sat = (np.in1d(true_satellites,identified_satellites))
    fu_sat = true_satellites[fu_sat]
    
    #which identified centrals are actually satellites?
    fu_isat = (np.in1d(identified_satellites,true_centrals))
    nfu_isat = (np.in1d(identified_satellites,true_satellites))
    fu_isat = identified_satellites[fu_isat]
    
    f_fu_icen = f_prop(GC['HALO_M'],bins,fu_icen,identified_centrals[nfu_icen],identified_centrals_bool)
    f_fu_isat = f_prop(GC['HALO_M'],bins,fu_isat,identified_satellites[nfu_isat],identified_satellites_bool)
    f_fu_cen = f_prop(GC['HALO_M'],bins,fu_cen,true_centrals[nfu_cen],true_centrals_bool)
    f_fu_sat = f_prop(GC['HALO_M'],bins,fu_sat,true_satellites[nfu_sat],true_satellites_bool)
    
    f_fu_icen_im = f_prop(GC['MGROUP'],bins,fu_icen,identified_centrals[nfu_icen],identified_centrals_bool)
    f_fu_isat_im = f_prop(GC['MGROUP'],bins,fu_isat,identified_satellites[nfu_isat],identified_satellites_bool)
    f_fu_cen_im = f_prop(GC['MGROUP'],bins,fu_cen,true_centrals[nfu_cen],true_centrals_bool)
    f_fu_sat_im = f_prop(GC['MGROUP'],bins,fu_sat,true_satellites[nfu_sat],true_satellites_bool)
    
    fig1 = plt.figure(figsize=(3.3,3.3))
    l1, = plt.plot(bin_centers, f_fu_icen, '--', color='black')
    l2, = plt.plot(bin_centers, f_fu_cen, '-', color='black')
    l3, = plt.plot(bin_centers, f_fu_isat, '--', color='green')
    l4, = plt.plot(bin_centers, f_fu_sat, '-', color='green')
    plt.ylim([0,1])
    plt.xlim([11,15])
    plt.xlabel(r'$M_{\rm halo}$')
    plt.ylabel(r'$f_{fu}$')
    plt.legend([l1,l2,l3,l4],['inferred centrals','true centrals','inferred satellites','true satellites'])
    plt.show(block=False)
    fig1.savefig(plotpath+filename1+'.pdf')
    
    fig2 = plt.figure(figsize=(3.3,3.3))
    l5, = plt.plot(bin_centers, f_fu_icen_im, '--', color='black')
    l6, = plt.plot(bin_centers, f_fu_cen_im, '-', color='black',alpha=0.5)
    l7, = plt.plot(bin_centers, f_fu_isat_im, '--', color='green')
    l8, = plt.plot(bin_centers, f_fu_sat_im, '-', color='green',alpha=0.5)
    plt.ylim([0,1])
    plt.xlim([11,15])
    plt.xlabel(r'$M_{\rm group}$')
    plt.ylabel(r'$f_{fu}$')
    plt.legend([l1,l2,l3,l4],['inferred centrals','true centrals','inferred satellites','true satellites'])
    plt.show(block=False)
    fig2.savefig(plotpath+filename2+'.pdf')
    
    inds = np.digitize(GC['MGROUP'][identified_centrals],bins = bins)-1
    
    f_fu = np.zeros((len(bins)-1,))
    for i in range(0,len(bins)-1):
        selection = (inds==i)
        selection = identified_centrals[selection]
        N_bad = np.not_equal(GC['RANK'][selection],GC['HALO_RANK'][selection])
        N_bad = np.sum(N_bad)
        N_good = np.equal(GC['RANK'][selection],GC['HALO_RANK'][selection])
        N_good = np.sum(N_good)
        N_tot = len(selection)
        print N_tot, N_good+N_bad
        f_fu[i] = N_bad/float(N_good+N_bad)
    
    fig3 = plt.figure(figsize=(3.3,3.3))
    ax = fig3.add_subplot(1,1,1)
    fig3.subplots_adjust(left=0.2, right=0.9, bottom=0.2,top=0.9)
    p1, = ax.plot(10.0**bin_centers,f_fu)
    ax.set_ylim([0,1])
    ax.set_xlim([10.0**11,10.0**15])
    ax.set_xscale('log')
    ax.set_xlabel(r'$M_{\rm group}$ $[M_{\odot}/h]$')
    ax.set_ylabel(r'$f_{fu}$')
    plt.show(block=False)
    fig3.savefig(plotpath+filename3+'.pdf')
    
    inds = np.digitize(GC['HALO_M'][identified_centrals],bins = bins)-1
    f_fu = np.zeros((len(bins)-1,))
    for i in range(0,len(bins)-1):
        selection = (inds==i)
        selection = identified_centrals[selection]
        bad = np.not_equal(GC['RANK'][selection],GC['HALO_RANK'][selection])
        N_bad = np.sum(bad)
        good = np.equal(GC['RANK'][selection],GC['HALO_RANK'][selection])
        N_good = np.sum(good)
        N_tot = len(selection)
        print N_tot, N_good+N_bad
        f_fu[i] = N_bad/float(N_good+N_bad)
    
    fig4 = plt.figure(figsize=(3.3,3.3))
    ax = fig4.add_subplot(1,1,1)
    fig4.subplots_adjust(left=0.2, right=0.9, bottom=0.2,top=0.9)
    p2, = ax.plot(10.0**bin_centers,f_fu)
    ax.set_ylim([0,1])
    ax.set_xlim([10.0**11,10.0**15])
    ax.set_xscale('log')
    ax.set_xlabel(r'$M_{\rm halo}$ $[M_{\odot}/h]$')
    ax.set_ylabel(r'$f_{inversion}$')
    plt.show(block=False)
    fig4.savefig(plotpath+filename4+'.pdf')
    
    #look at purity and completeness for perfect groups
    satellites = (GC['RANK']>0)
    GC['RANK'][satellites]=2
    centrals = (GC['RANK']==0)
    GC['RANK'][centrals]=1
    
    satellites = (GC['HALO_RANK']!=0)
    GC['HALO_RANK'][satellites]=2
    centrals = (GC['HALO_RANK']==0)
    GC['HALO_RANK'][centrals]=1
    
    #define mass bins
    bins = np.arange(12,15,0.2)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    
    result1 = sp_in_mass_bins(GC,bins)
    result2 = sc_in_mass_bins(GC,bins)
    result3 = cp_in_mass_bins(GC,bins)
    result4 = cc_in_mass_bins(GC,bins)
    result5 = Tsat_in_mass_bins(GC,bins)

    fig = plt.figure(figsize=(3.3,3.3))
    ax = fig.add_subplot(1,1,1)
    p1,=ax.plot(10**bin_centers,result1, color='green', ls='-')
    p2,=ax.plot(10**bin_centers,result2, color='green',ls='--')
    p3,=ax.plot(10**bin_centers,result3, color='orange',ls='-')
    p4,=ax.plot(10**bin_centers,result4, color='orange',ls='--')
    ax.set_ylim([0,1])
    ax.set_xscale('log')
    ax.legend((p1,p2, p3, p4),('satellite purity','satellite \ncompleteness','central purity','central \ncompleteness'),loc='lower right', fontsize=10, numpoints=1, frameon=False)
    plt.show()
    
    fig = plt.figure(figsize=(3.3,3.3))
    ax = fig.add_subplot(1,1,1)
    fig.subplots_adjust(left=0.2, right=0.9, bottom=0.2,top=0.9)
    p1, = ax.plot(10**bin_centers, result5, '--', color='black')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([0.1,100])
    plt.show()

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

def Tsat_in_mass_bins(GC,bins):
    x=np.zeros(len(bins)-1)
    result_2 = np.digitize(GC['HALO_M'],bins=bins)
    result_1 = np.digitize(GC['MGROUP'],bins=bins)
    #result_1 = result_2
    #result_2 = result_1
    for i in range(0,len(x)):
        ind_1 = (result_1==i+1)
        ind_2 = (result_2==i+1)
        GC_1 = GC[ind_1]
        GC_2 = GC[ind_2]
        Nsc=float(len(np.where((GC_1['RANK']!=1) & (GC_1['HALO_RANK']==1))[0]))
        Ncs=float(len(np.where((GC_2['RANK']==1) & (GC_2['HALO_RANK']!=1))[0]))
        Nss_a=float(len(np.where((GC_1['RANK']!=1) & (GC_1['HALO_RANK']!=1))[0]))
        Nss_b=float(len(np.where((GC_2['RANK']==1) & (GC_2['HALO_RANK']==1))[0]))
        x[i]=((Nsc+Nss_a)/(Ncs+Nss_a))
    return x

def sat_bias(C):
    Nsc=float(len(np.where((C['RANK']!=1) & (C['HALO_RANK']==1))[0]))
    Ncs=float(len(np.where((C['RANK']==1) & (C['HALO_RANK']!=1))[0]))
    Nss=float(len(np.where((C['RANK']!=1) & (C['HALO_RANK']!=1))[0]))
    b=((Nsc+Nss)/(Ncs+Nss))
    return b

def sp_in_mass_bins(GC,bins):
    x=np.zeros(len(bins)-1)
    #result = np.digitize(GC['HALO_M'],bins=bins)
    result = np.digitize(GC['MGROUP'],bins=bins)
    for i in range(0,len(x)):
        ind = (result==i+1)
        x[i]=sat_purity(GC[ind])
    return x

def sc_in_mass_bins(GC,bins):
    x=np.zeros(len(bins)-1)
    result = np.digitize(GC['HALO_M'],bins=bins)
    #result = np.digitize(GC['MGROUP'],bins=bins)
    for i in range(0,len(x)):
        ind = (result==i+1)
        x[i]=sat_comp(GC[ind])
    return x

def cp_in_mass_bins(GC,bins):
    x=np.zeros(len(bins)-1)
    #result = np.digitize(GC['HALO_M'],bins=bins)
    result = np.digitize(GC['MGROUP'],bins=bins)
    for i in range(0,len(x)):
        ind = (result==i+1)
        x[i]=cen_purity(GC[ind])
    return x

def cc_in_mass_bins(GC,bins):
    x=np.zeros(len(bins)-1)
    result = np.digitize(GC['HALO_M'],bins=bins)
    #result = np.digitize(GC['MGROUP'],bins=bins)
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
    print Ncc, Ncs
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

def solar_lum(M,Msol):
    L = ((Msol-M)/2.5)
    return L

if __name__ == '__main__':
    main() 
