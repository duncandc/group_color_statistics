#!/usr/bin/env python

#Duncan Campbell
#April 8, 2014
#Yale University
#Plot the CLF of the group finder runs on a mock and the intrinsic mock results. Includes 
#bootstrap error bars.

#load packages
import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
import custom_utilities as cu
import sys
from CLF import clf

def main():
    catalogue = sys.argv[1]

    savepath = cu.get_output_path() + 'analysis/groupcats/'
    plotpath = cu.get_plot_path()   + 'analysis/groupcats/'
    filename = catalogue + '_CLF'

    N_boots=50

    #set up ploting parameters
    fig,axes = plt.subplots(nrows=3, ncols=3, sharex=True, sharey=True, figsize=(6.95, 7.92))
    fig.subplots_adjust(hspace=0.05,wspace=0)
    fig.subplots_adjust(right=0.95, bottom=0.083, left=0.1, top=0.958)
    axes = axes.flatten()
    #plot 1
    axes[0].set_ylabel(r'$\phi(L)dlogL/group$')
    axes[0].set_title(r'$12.5<\log(M/[M_{\odot}h^{-1}])<13.0$')
    axes[0].set_yscale('log', nonposy='clip')
    axes[0].set_ylim(0.01,50)
    axes[0].set_xlim([9,11.5])
    #plot 2
    axes[1].set_title(r'$13.0<\log(M/[M_{\odot}h^{-1}])<13.5$')
    axes[1].set_yscale('log', nonposy='clip')
    axes[1].set_ylim(0.01,50)
    axes[1].set_xlim([9,11.5])
    #plot 3
    axes[2].set_title(r'$13.5<\log(M/[M_{\odot}h^{-1}])<14.0$')
    axes[2].set_yscale('log', nonposy='clip')
    axes[2].set_ylim(0.01,50)
    axes[2].set_xlim([9,11.5])
    axes[2].text(11.05,8,'Berlind FoF groups', rotation=90)
    #plot 4
    axes[3].set_ylabel(r'$\phi(L)dlogL/group$')
    axes[3].set_yscale('log', nonposy='clip')
    axes[3].set_ylim(0.01,50)
    axes[3].set_xlim([9,11.5])
    #plot 5
    axes[4].set_yscale('log', nonposy='clip')
    axes[4].set_ylim(0.01,50)
    axes[4].set_xlim([9,11.5])
    #plot 6
    axes[5].set_yscale('log', nonposy='clip')
    axes[5].set_ylim(0.01,50)
    axes[5].set_xlim([9,11.5])
    axes[5].text(11.05,8,'Tinker SO groups', rotation=90)
    #plot 7
    axes[6].set_ylabel(r'$\phi(L)dlogL/group$')
    axes[6].set_xlabel(r'$log(L/[L_{\odot}h^{-2}])$')
    axes[6].set_yscale('log', nonposy='clip')
    axes[6].set_ylim(0.01,50)
    axes[6].set_xlim([9,11.5])
    #plot 8
    axes[7].set_xlabel(r'$log(L/[L_{\odot}h^{-2}])$')
    axes[7].set_yscale('log', nonposy='clip')
    axes[7].set_ylim(0.01,50)
    axes[7].set_xlim([9,11.5])
    #plot 9
    axes[8].set_xlabel(r'$log(L/[L_{\odot}h^{-2}])$')
    axes[8].set_yscale('log', nonposy='clip')
    axes[8].set_ylim(0.01,50)
    axes[8].set_xlim([9.5,11.0])
    axes[8].text(11.05,6,'Yang SO groups', rotation=90)
    axes[8].set_xticklabels(["","9.6", "9.8","10.0", "10.2","10.4","10.6","10.8",""])
    #axes[8].set_yticklabels(["", "$0.1$", "$1$","$10$"])

    bins = np.arange(9.45,12.0,0.15)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    S_r = 4.64

    ##########################
    #calculate for mock HOD
    ########################## 
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    print 'opening mock catalogue:', catalogue+'.hdf5'
    #open catalogue
    f1 = h5py.File(filepath_mock+catalogue+'_extended.hdf5', 'r') #open catalogue file
    mock = f1.get(catalogue+'_extended')
    mock = np.array(mock)
    print 'length:', len(mock)
    for name in mock.dtype.names: print '     ', name

    occupied = np.where(mock['M_r,0.1']!=-99)
    empty    = np.where(mock['M_r,0.1']==-99)
    host     = np.where(mock['ID_host']==-1)[0]
    color    = mock['g-r']
    LHS      = 0.7 - 0.032*(mock['M_r,0.1']+16.5) #Weinmann 2006
    blue     = np.array((color<LHS) & (mock['M_r,0.1']!=-99))
    red      = np.array((color>LHS) & (mock['M_r,0.1']!=-99))
    centrals   = np.array(mock['ID_host']==-1)
    satellites = np.array(mock['ID_host']!=-1)

    host_ID = mock['ID_host']
    host_ID[host]=mock['ID_halo'][host]
    L = solar_lum(mock['M_r,0.1'],S_r)
    
    mass_bin = [12.5,13.0]
    phi = clf(host_ID,L,bins,mock['M_host'],mass_bin)
    mask = red & centrals
    phi_red_cen = clf(host_ID,L,bins,mock['M_host'],mass_bin,mask)
    mask = blue & centrals
    phi_blue_cen = clf(host_ID,L,bins,mock['M_host'],mass_bin,mask)
    mask = red & satellites
    phi_red_sat = clf(host_ID,L,bins,mock['M_host'],mass_bin,mask)
    mask = blue & satellites
    phi_blue_sat = clf(host_ID,L,bins,mock['M_host'],mass_bin,mask)
    print phi_blue_sat
    #plot the results
    ax=axes[0]
    p1a,=ax.plot(bin_centers, phi, color='black')
    ax=axes[0]
    p2a,=ax.plot(bin_centers, phi_red_cen, '--', color='red')
    ax=axes[0]
    p3a,=ax.plot(bin_centers, phi_blue_cen, '--', color='blue')
    ax=axes[3]
    p1a,=ax.plot(bin_centers, phi, color='black')
    ax=axes[3]
    p2a,=ax.plot(bin_centers, phi_red_cen, '--', color='red')
    ax=axes[3]
    p3a,=ax.plot(bin_centers, phi_blue_cen, '--', color='blue')
    ax=axes[6]
    p1a,=ax.plot(bin_centers, phi, color='black')
    ax=axes[6]
    p2a,=ax.plot(bin_centers, phi_red_cen, '--', color='red')
    ax=axes[6]
    p3a,=ax.plot(bin_centers, phi_blue_cen, '--', color='blue')
    ax=axes[0]
    p2a,=ax.plot(bin_centers, phi_red_sat, color='red')
    ax=axes[0]
    p3a,=ax.plot(bin_centers, phi_blue_sat, color='blue')
    ax=axes[3]
    p2a,=ax.plot(bin_centers, phi_red_sat, color='red')
    ax=axes[3]
    p3a,=ax.plot(bin_centers, phi_blue_sat, color='blue')
    ax=axes[6]
    p2a,=ax.plot(bin_centers, phi_red_sat, color='red')
    ax=axes[6]
    p3a,=ax.plot(bin_centers, phi_blue_sat, color='blue')

    mass_bin = [13.0,13.5]
    phi = clf(host_ID,L,bins,mock['M_host'],mass_bin)
    mask = red & centrals
    phi_red_cen = clf(host_ID,L,bins,mock['M_host'],mass_bin,mask)
    mask = blue & centrals
    phi_blue_cen = clf(host_ID,L,bins,mock['M_host'],mass_bin,mask)
    mask = red & satellites
    phi_red_sat = clf(host_ID,L,bins,mock['M_host'],mass_bin,mask)
    mask = blue & satellites
    phi_blue_sat = clf(host_ID,L,bins,mock['M_host'],mass_bin,mask)
    print phi_blue_sat
    #plot the results
    ax=axes[1]
    p1a,=ax.plot(bin_centers, phi, color='black')
    ax=axes[1]
    p2a,=ax.plot(bin_centers, phi_red_cen, '--', color='red')
    ax=axes[1]
    p3a,=ax.plot(bin_centers, phi_blue_cen, '--', color='blue')
    ax=axes[4]
    p1a,=ax.plot(bin_centers, phi, color='black')
    ax=axes[4]
    p2a,=ax.plot(bin_centers, phi_red_cen, '--', color='red')
    ax=axes[4]
    p3a,=ax.plot(bin_centers, phi_blue_cen, '--', color='blue')
    ax=axes[7]
    p1a,=ax.plot(bin_centers, phi, color='black')
    ax=axes[7]
    p2a,=ax.plot(bin_centers, phi_red_cen, '--', color='red')
    ax=axes[7]
    p3a,=ax.plot(bin_centers, phi_blue_cen, '--', color='blue')
    ax=axes[1]
    p2a,=ax.plot(bin_centers, phi_red_sat, color='red')
    ax=axes[1]
    p3a,=ax.plot(bin_centers, phi_blue_sat, color='blue')
    ax=axes[4]
    p2a,=ax.plot(bin_centers, phi_red_sat, color='red')
    ax=axes[4]
    p3a,=ax.plot(bin_centers, phi_blue_sat, color='blue')
    ax=axes[7]
    p2a,=ax.plot(bin_centers, phi_red_sat, color='red')
    ax=axes[7]
    p3a,=ax.plot(bin_centers, phi_blue_sat, color='blue')

    mass_bin = [13.5,14.0]
    phi = clf(host_ID,L,bins,mock['M_host'],mass_bin)
    mask = red & centrals
    phi_red_cen = clf(host_ID,L,bins,mock['M_host'],mass_bin,mask)
    mask = blue & centrals
    phi_blue_cen = clf(host_ID,L,bins,mock['M_host'],mass_bin,mask)
    mask = red & satellites
    phi_red_sat = clf(host_ID,L,bins,mock['M_host'],mass_bin,mask)
    mask = blue & satellites
    phi_blue_sat = clf(host_ID,L,bins,mock['M_host'],mass_bin,mask)
    print phi_blue_sat
    #plot the results
    ax=axes[2]
    p1a,=ax.plot(bin_centers, phi, color='black')
    ax=axes[2]
    p2a,=ax.plot(bin_centers, phi_red_cen, '--', color='red')
    ax=axes[2]
    p3a,=ax.plot(bin_centers, phi_blue_cen, '--', color='blue')
    ax=axes[5]
    p1a,=ax.plot(bin_centers, phi, color='black')
    ax=axes[5]
    p2a,=ax.plot(bin_centers, phi_red_cen, '--', color='red')
    ax=axes[5]
    p3a,=ax.plot(bin_centers, phi_blue_cen, '--', color='blue')
    ax=axes[8]
    p1a,=ax.plot(bin_centers, phi, color='black')
    ax=axes[8]
    p2a,=ax.plot(bin_centers, phi_red_cen, '--', color='red')
    ax=axes[8]
    p3a,=ax.plot(bin_centers, phi_blue_cen, '--', color='blue')
    ax=axes[2]
    p2a,=ax.plot(bin_centers, phi_red_sat, color='red')
    ax=axes[2]
    p3a,=ax.plot(bin_centers, phi_blue_sat, color='blue')
    ax=axes[5]
    p2a,=ax.plot(bin_centers, phi_red_sat, color='red')
    ax=axes[5]
    p3a,=ax.plot(bin_centers, phi_blue_sat, color='blue')
    ax=axes[8]
    p2a,=ax.plot(bin_centers, phi_red_sat, color='red')
    ax=axes[8]
    p3a,=ax.plot(bin_centers, phi_blue_sat, color='blue')


    ##########################
    #calculate for ideal groups
    ########################## 
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/ideal_groups/'
    print 'opening mock catalogue:', catalogue+'_groups.hdf5'
    f1 = h5py.File(filepath_mock+catalogue+'_groups.hdf5', 'r') #open catalogue file
    GC = f1.get(catalogue+'_groups')

    #centrals   = np.array(GC['HALO_RANK']==0)
    #satellites = np.array(GC['HALO_RANK']>0)
    centrals   = np.array(GC['RANK']==0)
    satellites = np.array(GC['RANK']>0)
    color = GC['M_g,0.1']-GC['M_r,0.1']
    LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
    blue  = np.array(color<LHS) 
    red   = np.array(color>LHS)

    host_ID = GC['GROUP_ID']
    L = solar_lum(GC['M_r,0.1'],S_r)
    
    #halo_mass = GC['HALO_M']
    halo_mass = GC['MGROUP']
    
    mass_bin = [12.5,13.0]
    phi_1 = clf(host_ID,L,bins,halo_mass,mass_bin)
    mask=red&centrals
    phi_red_cen_1 = clf(host_ID,L,bins,halo_mass,mass_bin,mask)
    mask=blue&centrals
    phi_blue_cen_1 = clf(host_ID,L,bins,halo_mass,mass_bin,mask)
    mask=red&satellites
    phi_red_sat_1 = clf(host_ID,L,bins,halo_mass,mass_bin,mask)
    mask=blue&satellites
    phi_blue_sat_1 = clf(host_ID,L,bins,halo_mass,mass_bin,mask)

    mass_bin = [13.0,13.5]
    phi_2 = clf(host_ID,L,bins,halo_mass,mass_bin)
    mask=red&centrals
    phi_red_cen_2 = clf(host_ID,L,bins,halo_mass,mass_bin,mask)
    mask=blue&centrals
    phi_blue_cen_2 = clf(host_ID,L,bins,halo_mass,mass_bin,mask)
    mask=red&satellites
    phi_red_sat_2 = clf(host_ID,L,bins,halo_mass,mass_bin,mask)
    mask=blue&satellites
    phi_blue_sat_2= clf(host_ID,L,bins,halo_mass,mass_bin,mask)

    mass_bin = [13.5,14.0]
    phi_3 = clf(host_ID,L,bins,halo_mass,mass_bin)
    mask=red&centrals
    phi_red_cen_3 = clf(host_ID,L,bins,halo_mass,mass_bin,mask)
    mask=blue&centrals
    phi_blue_cen_3 = clf(host_ID,L,bins,halo_mass,mass_bin,mask)
    mask=red&satellites
    phi_red_sat_3 = clf(host_ID,L,bins,halo_mass,mass_bin,mask)
    mask=blue&satellites
    phi_blue_sat_3 = clf(host_ID,L,bins,halo_mass,mass_bin,mask)
    
    #plot the results
    ax=axes[0]
    p1a,=ax.plot(bin_centers, phi_1, color='black', alpha=0.5)
    ax=axes[0]
    p2a,=ax.plot(bin_centers, phi_red_cen_1, '--', color='red', alpha=0.5)
    ax=axes[0]
    p3a,=ax.plot(bin_centers, phi_blue_cen_1, '--', color='blue', alpha=0.5)
    ax=axes[3]
    p1a,=ax.plot(bin_centers, phi_1, color='black', alpha=0.5)
    ax=axes[3]
    p2a,=ax.plot(bin_centers, phi_red_cen_1, '--', color='red', alpha=0.5)
    ax=axes[3]
    p3a,=ax.plot(bin_centers, phi_blue_cen_1, '--', color='blue', alpha=0.5)
    ax=axes[6]
    p1a,=ax.plot(bin_centers, phi_1, color='black', alpha=0.5)
    ax=axes[6]
    p2a,=ax.plot(bin_centers, phi_red_cen_1, '--', color='red', alpha=0.5)
    ax=axes[6]
    p3a,=ax.plot(bin_centers, phi_blue_cen_1, '--', color='blue', alpha=0.5)
    ax=axes[0]
    p2a,=ax.plot(bin_centers, phi_red_sat_1, color='red', alpha=0.5)
    ax=axes[0]
    p3a,=ax.plot(bin_centers, phi_blue_sat_1, color='blue', alpha=0.5)
    ax=axes[3]
    p2a,=ax.plot(bin_centers, phi_red_sat_1, color='red', alpha=0.5)
    ax=axes[3]
    p3a,=ax.plot(bin_centers, phi_blue_sat_1, color='blue', alpha=0.5)
    ax=axes[6]
    p2a,=ax.plot(bin_centers, phi_red_sat_1, color='red', alpha=0.5)
    ax=axes[6]
    p3a,=ax.plot(bin_centers, phi_blue_sat_1, color='blue', alpha=0.5)

    ax=axes[1]
    p1a,=ax.plot(bin_centers, phi_2, color='black', alpha=0.5)
    ax=axes[1]
    p2a,=ax.plot(bin_centers, phi_red_cen_2, '--', color='red', alpha=0.5)
    ax=axes[1]
    p3a,=ax.plot(bin_centers, phi_blue_cen_2, '--', color='blue', alpha=0.5)
    ax=axes[4]
    p1a,=ax.plot(bin_centers, phi_2, color='black', alpha=0.5)
    ax=axes[4]
    p2a,=ax.plot(bin_centers, phi_red_cen_2, '--', color='red', alpha=0.5)
    ax=axes[4]
    p3a,=ax.plot(bin_centers, phi_blue_cen_2, '--', color='blue', alpha=0.5)
    ax=axes[7]
    p1a,=ax.plot(bin_centers, phi_2, color='black', alpha=0.5)
    ax=axes[7]
    p2a,=ax.plot(bin_centers, phi_red_cen_2, '--', color='red', alpha=0.5)
    ax=axes[7]
    p3a,=ax.plot(bin_centers, phi_blue_cen_2, '--', color='blue', alpha=0.5)
    ax=axes[1]
    p2a,=ax.plot(bin_centers, phi_red_sat_2, color='red', alpha=0.5)
    ax=axes[1]
    p3a,=ax.plot(bin_centers, phi_blue_sat_2, color='blue', alpha=0.5)
    ax=axes[4]
    p2a,=ax.plot(bin_centers, phi_red_sat_2, color='red', alpha=0.5)
    ax=axes[4]
    p3a,=ax.plot(bin_centers, phi_blue_sat_2, color='blue', alpha=0.5)
    ax=axes[7]
    p2a,=ax.plot(bin_centers, phi_red_sat_2, color='red', alpha=0.5)
    ax=axes[7]
    p3a,=ax.plot(bin_centers, phi_blue_sat_2, color='blue', alpha=0.5)

    ax=axes[2]
    p1a,=ax.plot(bin_centers, phi_3, color='black', alpha=0.5)
    ax=axes[2]
    p2a,=ax.plot(bin_centers, phi_red_cen_3, '--', color='red', alpha=0.5)
    ax=axes[2]
    p3a,=ax.plot(bin_centers, phi_blue_cen_3, '--', color='blue', alpha=0.5)
    ax=axes[5]
    p1a,=ax.plot(bin_centers, phi_3, color='black', alpha=0.5)
    ax=axes[5]
    p2a,=ax.plot(bin_centers, phi_red_cen_3, '--', color='red', alpha=0.5)
    ax=axes[5]
    p3a,=ax.plot(bin_centers, phi_blue_cen_3, '--', color='blue', alpha=0.5)
    ax=axes[8]
    p1a,=ax.plot(bin_centers, phi_3, color='black', alpha=0.5)
    ax=axes[8]
    p2a,=ax.plot(bin_centers, phi_red_cen_3, '--', color='red', alpha=0.5)
    ax=axes[8]
    p3a,=ax.plot(bin_centers, phi_blue_cen_3, '--', color='blue', alpha=0.5)
    ax=axes[2]
    p2a,=ax.plot(bin_centers, phi_red_sat_3, color='red', alpha=0.5)
    ax=axes[2]
    p3a,=ax.plot(bin_centers, phi_blue_sat_3, color='blue', alpha=0.5)
    ax=axes[5]
    p2a,=ax.plot(bin_centers, phi_red_sat_3, color='red', alpha=0.5)
    ax=axes[5]
    p3a,=ax.plot(bin_centers, phi_blue_sat_3, color='blue', alpha=0.5)
    ax=axes[8]
    p2a,=ax.plot(bin_centers, phi_red_sat_3, color='red', alpha=0.5)
    ax=axes[8]
    p3a,=ax.plot(bin_centers, phi_blue_sat_3, color='blue', alpha=0.5)


    ##########################
    #calculate mock groups CLF
    ##########################
    group_cat = 'berlind_groupcat'
    filepath_GC   = cu.get_output_path() + 'processed_data/'+group_cat+'/mock_runs/4th_run/custom_catalogues/'
    if group_cat == 'tinker_groupcat': group_catalogue = catalogue + '_clf_groups_M19'
    if group_cat == 'berlind_groupcat': group_catalogue = catalogue + '_groups'
    if group_cat == 'yang_groupcat': group_catalogue = catalogue + '_groups'
    #set up arrays to store results
    phi_1 = np.zeros((N_boots,len(bins)-1))
    phi_red_cen_1 = np.zeros((N_boots,len(bins)-1))
    phi_blue_cen_1 = np.zeros((N_boots,len(bins)-1))
    phi_red_sat_1 = np.zeros((N_boots,len(bins)-1))
    phi_blue_sat_1 = np.zeros((N_boots,len(bins)-1))
    phi_2 = np.zeros((N_boots,len(bins)-1))
    phi_red_cen_2 = np.zeros((N_boots,len(bins)-1))
    phi_blue_cen_2 = np.zeros((N_boots,len(bins)-1))
    phi_red_sat_2 = np.zeros((N_boots,len(bins)-1))
    phi_blue_sat_2 = np.zeros((N_boots,len(bins)-1))
    phi_3 = np.zeros((N_boots,len(bins)-1))
    phi_red_cen_3 = np.zeros((N_boots,len(bins)-1))
    phi_blue_cen_3 = np.zeros((N_boots,len(bins)-1))
    phi_red_sat_3 = np.zeros((N_boots,len(bins)-1))
    phi_blue_sat_3 = np.zeros((N_boots,len(bins)-1))
    for boot in range(0,N_boots):
        print boot
        f2 = h5py.File(filepath_GC+'bootstraps/'+group_catalogue+'_'+str(boot)+'.hdf5', 'r')
        GC = f2.get(group_catalogue+'_'+str(boot))
        GC = np.array(GC)

        centrals   = np.array(GC['RANK']==1)
        satellites = np.array(GC['RANK']==0)
        color = GC['M_g,0.1']-GC['M_r,0.1']
        LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
        blue  = np.array(color<LHS) 
        red   = np.array(color>LHS)

        host_ID = GC['GROUP_ID']
        L = solar_lum(GC['M_r,0.1'],S_r)
    
        mass_bin = [12.5,13.0]
        phi_1[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin)
        mask=red&centrals
        phi_red_cen_1[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=blue&centrals
        phi_blue_cen_1[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=red&satellites
        phi_red_sat_1[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=blue&satellites
        phi_blue_sat_1[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)

        mass_bin = [13.0,13.5]
        phi_2[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin)
        mask=red&centrals
        phi_red_cen_2[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=blue&centrals
        phi_blue_cen_2[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=red&satellites
        phi_red_sat_2[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=blue&satellites
        phi_blue_sat_2[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)

        mass_bin = [13.5,14.0]
        phi_3[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin)
        mask=red&centrals
        phi_red_cen_3[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=blue&centrals
        phi_blue_cen_3[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=red&satellites
        phi_red_sat_3[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=blue&satellites
        phi_blue_sat_3[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)

    phi_err_1 = np.std(phi_1,axis=0)
    phi_1 = np.mean(phi_1,axis=0)
    phi_red_cen_err_1 = np.std(phi_red_cen_1,axis=0)
    phi_red_cen_1 = np.mean(phi_red_cen_1,axis=0)
    phi_blue_cen_err_1 = np.std(phi_blue_cen_1,axis=0)
    phi_blue_cen_1 = np.mean(phi_blue_cen_1,axis=0)
    phi_red_sat_err_1 = np.std(phi_red_sat_1,axis=0)
    phi_red_sat_1 = np.mean(phi_red_sat_1,axis=0)
    phi_blue_sat_err_1 = np.std(phi_blue_sat_1,axis=0)
    phi_blue_sat_1 = np.mean(phi_blue_sat_1,axis=0)
    phi_err_2 = np.std(phi_2,axis=0)
    phi_2 = np.mean(phi_2,axis=0)
    phi_red_cen_err_2 = np.std(phi_red_cen_2,axis=0)
    phi_red_cen_2 = np.mean(phi_red_cen_2,axis=0)
    phi_blue_cen_err_2 = np.std(phi_blue_cen_2,axis=0)
    phi_blue_cen_2 = np.mean(phi_blue_cen_2,axis=0)
    phi_red_sat_err_2 = np.std(phi_red_sat_2,axis=0)
    phi_red_sat_2 = np.mean(phi_red_sat_2,axis=0)
    phi_blue_sat_err_2 = np.std(phi_blue_sat_2,axis=0)
    phi_blue_sat_2 = np.mean(phi_blue_sat_2,axis=0)
    phi_err_3 = np.std(phi_3,axis=0)
    phi_3 = np.mean(phi_3,axis=0)
    phi_red_cen_err_3 = np.std(phi_red_cen_3,axis=0)
    phi_red_cen_3 = np.mean(phi_red_cen_3,axis=0)
    phi_blue_cen_err_3 = np.std(phi_blue_cen_3,axis=0)
    phi_blue_cen_3 = np.mean(phi_blue_cen_3,axis=0)
    phi_red_sat_err_3 = np.std(phi_red_sat_3,axis=0)
    phi_red_sat_3 = np.mean(phi_red_sat_3,axis=0)
    phi_blue_sat_err_3 = np.std(phi_blue_sat_3,axis=0)
    phi_blue_sat_3 = np.mean(phi_blue_sat_3,axis=0)

    ax=axes[0]
    ax.errorbar(bin_centers,phi_1,yerr=phi_err_1, fmt='o', color='black', ms=4, mec='none')
    ax.plot(bin_centers,phi_red_cen_1, 'o', mec='red',mew=1.1, ms=4, mfc='none')
    ax.plot(bin_centers,phi_blue_cen_1, 'o', mec='blue',mew=1.1, ms=4, mfc='none')
    ax.plot(bin_centers,phi_red_sat_1, 'o', color='red', ms=4, mec='none')
    ax.plot(bin_centers,phi_blue_sat_1, 'o', color='blue', ms=4, mec='none')
    ax=axes[1]
    ax.errorbar(bin_centers,phi_2,yerr=phi_err_2, fmt='o', color='black', ms=4, mec='none')
    ax.plot(bin_centers,phi_red_cen_2, 'o', mec='red',mew=1.1, ms=4, mfc='none')
    ax.plot(bin_centers,phi_blue_cen_2, 'o', mec='blue',mew=1.1, ms=4, mfc='none')
    ax.plot(bin_centers,phi_red_sat_2, 'o', color='red', ms=4, mec='none')
    ax.plot(bin_centers,phi_blue_sat_2, 'o', color='blue', ms=4, mec='none')
    ax=axes[2]
    ax.errorbar(bin_centers,phi_3,yerr=phi_err_3, fmt='o', color='black', ms=4, mec='none')
    ax.plot(bin_centers,phi_red_cen_3, 'o', mec='red',mew=1.1, ms=4, mfc='none')
    ax.plot(bin_centers,phi_blue_cen_3, 'o', mec='blue',mew=1.1, ms=4, mfc='none')
    ax.plot(bin_centers,phi_red_sat_3, 'o', color='red', ms=4, mec='none')
    ax.plot(bin_centers,phi_blue_sat_3, 'o', color='blue', ms=4, mec='none')


    ##########################
    #calculate mock groups HOD
    ##########################
    group_cat = 'tinker_groupcat'
    filepath_GC   = cu.get_output_path() + 'processed_data/'+group_cat+'/mock_runs/4th_run/custom_catalogues/'
    if group_cat == 'tinker_groupcat': group_catalogue = catalogue + '_clf_groups_M19'
    if group_cat == 'berlind_groupcat': group_catalogue = catalogue + '_groups'
    if group_cat == 'yang_groupcat': group_catalogue = catalogue + '_groups'
    #set up arrays to store results
    phi_1 = np.zeros((N_boots,len(bins)-1))
    phi_red_cen_1 = np.zeros((N_boots,len(bins)-1))
    phi_blue_cen_1 = np.zeros((N_boots,len(bins)-1))
    phi_red_sat_1 = np.zeros((N_boots,len(bins)-1))
    phi_blue_sat_1 = np.zeros((N_boots,len(bins)-1))
    phi_2 = np.zeros((N_boots,len(bins)-1))
    phi_red_cen_2 = np.zeros((N_boots,len(bins)-1))
    phi_blue_cen_2 = np.zeros((N_boots,len(bins)-1))
    phi_red_sat_2 = np.zeros((N_boots,len(bins)-1))
    phi_blue_sat_2 = np.zeros((N_boots,len(bins)-1))
    phi_3 = np.zeros((N_boots,len(bins)-1))
    phi_red_cen_3 = np.zeros((N_boots,len(bins)-1))
    phi_blue_cen_3 = np.zeros((N_boots,len(bins)-1))
    phi_red_sat_3 = np.zeros((N_boots,len(bins)-1))
    phi_blue_sat_3 = np.zeros((N_boots,len(bins)-1))
    for boot in range(0,N_boots):
        print boot
        f2 = h5py.File(filepath_GC+'bootstraps/'+group_catalogue+'_'+str(boot)+'.hdf5', 'r')
        GC = f2.get(group_catalogue+'_'+str(boot))
        GC = np.array(GC)

        centrals   = np.array(GC['RANK']==1)
        satellites = np.array(GC['RANK']==0)
        color = GC['M_g,0.1']-GC['M_r,0.1']
        LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
        blue  = np.array(color<LHS) 
        red   = np.array(color>LHS)

        host_ID = GC['GROUP_ID']
        L = solar_lum(GC['M_r,0.1'],S_r)
    
        mass_bin = [12.5,13.0]
        phi_1[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin)
        mask=red&centrals
        phi_red_cen_1[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=blue&centrals
        phi_blue_cen_1[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=red&satellites
        phi_red_sat_1[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=blue&satellites
        phi_blue_sat_1[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)

        mass_bin = [13.0,13.5]
        phi_2[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin)
        mask=red&centrals
        phi_red_cen_2[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=blue&centrals
        phi_blue_cen_2[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=red&satellites
        phi_red_sat_2[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=blue&satellites
        phi_blue_sat_2[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)

        mass_bin = [13.5,14.0]
        phi_3[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin)
        mask=red&centrals
        phi_red_cen_3[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=blue&centrals
        phi_blue_cen_3[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=red&satellites
        phi_red_sat_3[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=blue&satellites
        phi_blue_sat_3[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)

    phi_err_1 = np.std(phi_1,axis=0)
    phi_1 = np.mean(phi_1,axis=0)
    phi_red_cen_err_1 = np.std(phi_red_cen_1,axis=0)
    phi_red_cen_1 = np.mean(phi_red_cen_1,axis=0)
    phi_blue_cen_err_1 = np.std(phi_blue_cen_1,axis=0)
    phi_blue_cen_1 = np.mean(phi_blue_cen_1,axis=0)
    phi_red_sat_err_1 = np.std(phi_red_sat_1,axis=0)
    phi_red_sat_1 = np.mean(phi_red_sat_1,axis=0)
    phi_blue_sat_err_1 = np.std(phi_blue_sat_1,axis=0)
    phi_blue_sat_1 = np.mean(phi_blue_sat_1,axis=0)
    phi_err_2 = np.std(phi_2,axis=0)
    phi_2 = np.mean(phi_2,axis=0)
    phi_red_cen_err_2 = np.std(phi_red_cen_2,axis=0)
    phi_red_cen_2 = np.mean(phi_red_cen_2,axis=0)
    phi_blue_cen_err_2 = np.std(phi_blue_cen_2,axis=0)
    phi_blue_cen_2 = np.mean(phi_blue_cen_2,axis=0)
    phi_red_sat_err_2 = np.std(phi_red_sat_2,axis=0)
    phi_red_sat_2 = np.mean(phi_red_sat_2,axis=0)
    phi_blue_sat_err_2 = np.std(phi_blue_sat_2,axis=0)
    phi_blue_sat_2 = np.mean(phi_blue_sat_2,axis=0)
    phi_err_3 = np.std(phi_3,axis=0)
    phi_3 = np.mean(phi_3,axis=0)
    phi_red_cen_err_3 = np.std(phi_red_cen_3,axis=0)
    phi_red_cen_3 = np.mean(phi_red_cen_3,axis=0)
    phi_blue_cen_err_3 = np.std(phi_blue_cen_3,axis=0)
    phi_blue_cen_3 = np.mean(phi_blue_cen_3,axis=0)
    phi_red_sat_err_3 = np.std(phi_red_sat_3,axis=0)
    phi_red_sat_3 = np.mean(phi_red_sat_3,axis=0)
    phi_blue_sat_err_3 = np.std(phi_blue_sat_3,axis=0)
    phi_blue_sat_3 = np.mean(phi_blue_sat_3,axis=0)

    ax=axes[3]
    ax.errorbar(bin_centers,phi_1,yerr=phi_err_1, fmt='o', color='black', ms=4, mec='none')
    ax.plot(bin_centers,phi_red_cen_1, 'o', mec='red',mew=1.1, ms=4, mfc='none')
    ax.plot(bin_centers,phi_blue_cen_1, 'o', mec='blue',mew=1.1, ms=4, mfc='none')
    ax.plot(bin_centers,phi_red_sat_1, 'o', color='red', ms=4, mec='none')
    ax.plot(bin_centers,phi_blue_sat_1, 'o', color='blue', ms=4, mec='none')
    ax=axes[4]
    ax.errorbar(bin_centers,phi_2,yerr=phi_err_2, fmt='o', color='black', ms=4, mec='none')
    ax.plot(bin_centers,phi_red_cen_2, 'o', mec='red',mew=1.1, ms=4, mfc='none')
    ax.plot(bin_centers,phi_blue_cen_2, 'o', mec='blue',mew=1.1, ms=4, mfc='none')
    ax.plot(bin_centers,phi_red_sat_2, 'o', color='red', ms=4, mec='none')
    ax.plot(bin_centers,phi_blue_sat_2, 'o', color='blue', ms=4, mec='none')
    ax=axes[5]
    ax.errorbar(bin_centers,phi_3,yerr=phi_err_3, fmt='o', color='black', ms=4, mec='none')
    ax.plot(bin_centers,phi_red_cen_3, 'o', mec='red',mew=1.1, ms=4, mfc='none')
    ax.plot(bin_centers,phi_blue_cen_3, 'o', mec='blue',mew=1.1, ms=4, mfc='none')
    ax.plot(bin_centers,phi_red_sat_3, 'o', color='red', ms=4, mec='none')
    ax.plot(bin_centers,phi_blue_sat_3, 'o', color='blue', ms=4, mec='none')


    ##########################
    #calculate mock groups CLF
    ##########################
    group_cat = 'yang_groupcat'
    filepath_GC   = cu.get_output_path() + 'processed_data/'+group_cat+'/mock_runs/4th_run/version_2/custom_catalogues/'
    if group_cat == 'tinker_groupcat': group_catalogue = catalogue + '_clf_groups_M19'
    if group_cat == 'berlind_groupcat': group_catalogue = catalogue + '_groups'
    if group_cat == 'yang_groupcat': group_catalogue = catalogue + '_groups'
    #set up arrays to store results
    phi_1 = np.zeros((N_boots,len(bins)-1))
    phi_red_cen_1 = np.zeros((N_boots,len(bins)-1))
    phi_blue_cen_1 = np.zeros((N_boots,len(bins)-1))
    phi_red_sat_1 = np.zeros((N_boots,len(bins)-1))
    phi_blue_sat_1 = np.zeros((N_boots,len(bins)-1))
    phi_2 = np.zeros((N_boots,len(bins)-1))
    phi_red_cen_2 = np.zeros((N_boots,len(bins)-1))
    phi_blue_cen_2 = np.zeros((N_boots,len(bins)-1))
    phi_red_sat_2 = np.zeros((N_boots,len(bins)-1))
    phi_blue_sat_2 = np.zeros((N_boots,len(bins)-1))
    phi_3 = np.zeros((N_boots,len(bins)-1))
    phi_red_cen_3 = np.zeros((N_boots,len(bins)-1))
    phi_blue_cen_3 = np.zeros((N_boots,len(bins)-1))
    phi_red_sat_3 = np.zeros((N_boots,len(bins)-1))
    phi_blue_sat_3 = np.zeros((N_boots,len(bins)-1))
    for boot in range(0,N_boots):
        print boot
        f2 = h5py.File(filepath_GC+'bootstraps/'+group_catalogue+'_'+str(boot)+'.hdf5', 'r')
        GC = f2.get(group_catalogue+'_'+str(boot))
        GC = np.array(GC)

        centrals   = np.array(GC['RANK']==1)
        satellites = np.array(GC['RANK']==0)
        color = GC['M_g,0.1']-GC['M_r,0.1']
        LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
        blue  = np.array(color<LHS) 
        red   = np.array(color>LHS)

        host_ID = GC['GROUP_ID']
        L = solar_lum(GC['M_r,0.1'],S_r)
    
        mass_bin = [12.5,13.0]
        phi_1[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin)
        mask=red&centrals
        phi_red_cen_1[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=blue&centrals
        phi_blue_cen_1[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=red&satellites
        phi_red_sat_1[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=blue&satellites
        phi_blue_sat_1[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)

        mass_bin = [13.0,13.5]
        phi_2[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin)
        mask=red&centrals
        phi_red_cen_2[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=blue&centrals
        phi_blue_cen_2[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=red&satellites
        phi_red_sat_2[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=blue&satellites
        phi_blue_sat_2[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)

        mass_bin = [13.5,14.0]
        phi_3[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin)
        mask=red&centrals
        phi_red_cen_3[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=blue&centrals
        phi_blue_cen_3[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=red&satellites
        phi_red_sat_3[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)
        mask=blue&satellites
        phi_blue_sat_3[boot,:] = clf(host_ID,L,bins,GC['MGROUP'],mass_bin,mask)

    phi_err_1 = np.std(phi_1,axis=0)
    phi_1 = np.mean(phi_1,axis=0)
    phi_red_cen_err_1 = np.std(phi_red_cen_1,axis=0)
    phi_red_cen_1 = np.mean(phi_red_cen_1,axis=0)
    phi_blue_cen_err_1 = np.std(phi_blue_cen_1,axis=0)
    phi_blue_cen_1 = np.mean(phi_blue_cen_1,axis=0)
    phi_red_sat_err_1 = np.std(phi_red_sat_1,axis=0)
    phi_red_sat_1 = np.mean(phi_red_sat_1,axis=0)
    phi_blue_sat_err_1 = np.std(phi_blue_sat_1,axis=0)
    phi_blue_sat_1 = np.mean(phi_blue_sat_1,axis=0)
    phi_err_2 = np.std(phi_2,axis=0)
    phi_2 = np.mean(phi_2,axis=0)
    phi_red_cen_err_2 = np.std(phi_red_cen_2,axis=0)
    phi_red_cen_2 = np.mean(phi_red_cen_2,axis=0)
    phi_blue_cen_err_2 = np.std(phi_blue_cen_2,axis=0)
    phi_blue_cen_2 = np.mean(phi_blue_cen_2,axis=0)
    phi_red_sat_err_2 = np.std(phi_red_sat_2,axis=0)
    phi_red_sat_2 = np.mean(phi_red_sat_2,axis=0)
    phi_blue_sat_err_2 = np.std(phi_blue_sat_2,axis=0)
    phi_blue_sat_2 = np.mean(phi_blue_sat_2,axis=0)
    phi_err_3 = np.std(phi_3,axis=0)
    phi_3 = np.mean(phi_3,axis=0)
    phi_red_cen_err_3 = np.std(phi_red_cen_3,axis=0)
    phi_red_cen_3 = np.mean(phi_red_cen_3,axis=0)
    phi_blue_cen_err_3 = np.std(phi_blue_cen_3,axis=0)
    phi_blue_cen_3 = np.mean(phi_blue_cen_3,axis=0)
    phi_red_sat_err_3 = np.std(phi_red_sat_3,axis=0)
    phi_red_sat_3 = np.mean(phi_red_sat_3,axis=0)
    phi_blue_sat_err_3 = np.std(phi_blue_sat_3,axis=0)
    phi_blue_sat_3 = np.mean(phi_blue_sat_3,axis=0)

    ax=axes[6]
    ax.errorbar(bin_centers,phi_1,yerr=phi_err_1, fmt='o', color='black', ms=4, mec='none')
    ax.plot(bin_centers,phi_red_cen_1, 'o', mec='red',mew=1.1, ms=4, mfc='none')
    ax.plot(bin_centers,phi_blue_cen_1, 'o', mec='blue',mew=1.1, ms=4, mfc='none')
    ax.plot(bin_centers,phi_red_sat_1, 'o', color='red', ms=4, mec='none')
    ax.plot(bin_centers,phi_blue_sat_1, 'o', color='blue', ms=4, mec='none')
    ax=axes[7]
    ax.errorbar(bin_centers,phi_2,yerr=phi_err_2, fmt='o', color='black', ms=4, mec='none')
    ax.plot(bin_centers,phi_red_cen_2, 'o', mec='red',mew=1.1, ms=4, mfc='none')
    ax.plot(bin_centers,phi_blue_cen_2, 'o', mec='blue',mew=1.1, ms=4, mfc='none')
    ax.plot(bin_centers,phi_red_sat_2, 'o', color='red', ms=4, mec='none')
    ax.plot(bin_centers,phi_blue_sat_2, 'o', color='blue', ms=4, mec='none')
    ax=axes[8]
    ax.errorbar(bin_centers,phi_3,yerr=phi_err_3, fmt='o', color='black', ms=4, mec='none')
    ax.plot(bin_centers,phi_red_cen_3, 'o', mec='red', mew=1.1, ms=4, mfc='none')
    ax.plot(bin_centers,phi_blue_cen_3, 'o', mec='blue', mew=1.1, ms=4, mfc='none')
    ax.plot(bin_centers,phi_red_sat_3, 'o', color='red', ms=4, mec='none')
    ax.plot(bin_centers,phi_blue_sat_3, 'o', color='blue', ms=4, mec='none')

    plt.show(block=False)
    print plotpath+filename+'.pdf'
    fig.savefig(plotpath+filename+'.pdf')

def solar_lum(M,Msol):
    L = ((Msol-M)/2.5)
    return L    

if __name__ == '__main__':
    main()
