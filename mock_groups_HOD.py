#!/usr/bin/env python

#Duncan Campbell
#April 8, 2014
#Yale University
#Plot the HOD of the group finder runs on a mock and the intrinsic mock results. Includes 
#bootstrap error bars.

#load packages
import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
import custom_utilities as cu
import sys
from HOD import hod

def main():
    
    if len(sys.argv)>1: catalogue = sys.argv[1]
    else: catalogue = 'Mr19_age_distribution_matching_mock'

    savepath = cu.get_output_path() + 'analysis/groupcats/'
    plotpath = cu.get_plot_path()   + 'analysis/groupcats/'
    filename = catalogue + '_HOD'

    N_boots=50

    #set up ploting parameters
    fig,axes = plt.subplots(nrows=3, ncols=3, sharex=True, sharey=True, figsize=(6.95, 7.92))
    fig.subplots_adjust(hspace=0.05,wspace=0)
    fig.subplots_adjust(right=0.95, bottom=0.083, left=0.1, top=0.958)
    axes = axes.flatten()
    #plot 1
    axes[0].set_ylabel(r'$\langle N(M)\rangle$')
    axes[0].text(10**11.5,50,'all galaxies\n' r'$M_r<-19+5\log(h)$')
    axes[0].set_yscale('log', nonposy='clip')
    axes[0].set_xscale('log', nonposy='clip')
    axes[0].set_ylim([0.1,200])
    axes[0].set_xlim([10**11,10**15])
    #plot 2
    axes[1].text(10**11.5,50,'red galaxies\n' r'$M_r<-19+5\log(h)$')
    axes[1].set_yscale('log', nonposy='clip')
    axes[1].set_xscale('log', nonposy='clip')
    axes[1].set_ylim([0.1,200])
    axes[1].set_xlim([10**11,10**15])
    #plot 3
    axes[2].text(10**11.5,50,'blue galaxies\n' r'$M_r<-19+5\log(h)$')
    axes[2].text(10**15.15,25,'Berlind FoF groups', rotation=90)
    axes[2].set_yscale('log', nonposy='clip')
    axes[2].set_xscale('log', nonposy='clip')
    axes[2].set_ylim([0.1,200])
    axes[2].set_xlim([10**11,10**15])
    #plot 4
    axes[3].set_ylabel(r'$\langle N(M)\rangle$')
    axes[3].set_yscale('log', nonposy='clip')
    axes[3].set_xscale('log', nonposy='clip')
    axes[3].set_ylim([0.1,200])
    axes[3].set_xlim([10**11,10**15])
    #plot 5
    axes[4].set_yscale('log', nonposy='clip')
    axes[4].set_xscale('log', nonposy='clip')
    axes[4].set_ylim([0.1,200])
    axes[4].set_xlim([10**11,10**15])
    #plot 6
    axes[5].text(10**15.15,25,'Tinker SO groups', rotation=90)
    axes[5].set_yscale('log', nonposy='clip')
    axes[5].set_xscale('log', nonposy='clip')
    axes[5].set_ylim([0.1,200])
    axes[5].set_xlim([10**11,10**15])
    #plot 7
    axes[6].set_ylabel(r'$\langle N(M)\rangle$')
    axes[6].set_xlabel(r'$M$ $[M_{\odot}/h]$')
    axes[6].set_yscale('log', nonposy='clip')
    axes[6].set_xscale('log', nonposy='clip')
    axes[6].set_ylim([0.1,200])
    axes[6].set_xlim([10**11,10**15])
    #plot 8
    axes[7].set_xlabel(r'$M$ $[M_{\odot}/h]$')
    axes[7].set_yscale('log', nonposy='clip')
    axes[7].set_xscale('log', nonposy='clip')
    axes[7].set_ylim([0.1,200])
    axes[7].set_xlim([10**11,10**15])
    #plot 9
    axes[8].set_xlabel(r'$M$ $[M_{\odot}/h]$')
    axes[8].text(10**15.15,25,'Yang SO groups', rotation=90)
    axes[8].set_yscale('log', nonposy='clip')
    axes[8].set_xscale('log', nonposy='clip')
    axes[8].set_ylim([0.1,200])
    axes[8].set_xlim([10**11,10**15])
    axes[8].set_xticks([10**11,10**12,10**13,10**14,10**15])
    axes[8].set_xticklabels([r' ', r'$10^{12}$', r'$10^{13}$', r'$10^{14}$',r' '])
    axes[0].set_yticks([0.1,1.0,10.0,100.0])
    axes[0].set_yticklabels([r' ', r'$1$', r'$10$', r'$100$'])

    bins = np.arange(10,15,0.2)
    bin_centers = (bins[:-1]+bins[1:])/2.0

    ##########################
    #calculate for mock HOD
    ########################## 
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    print 'opening mock catalogue:', catalogue+'.hdf5'
    #open catalogue
    f1 = h5py.File(filepath_mock+catalogue+'_extended.hdf5', 'r') #open catalogue file
    mock = f1.get(catalogue+'_extended')
    mock = np.array(mock)
    #print 'length:', len(mock)
    #for name in mock.dtype.names: print '     ', name

    bins = np.arange(10,15,0.2)
    bin_centers = (bins[:-1]+bins[1:])/2.0

    occupied = np.where(mock['M_r,0.1']!=-99)
    empty    = np.where(mock['M_r,0.1']==-99)
    host     =  np.where(mock['ID_host']==-1)[0]
    color    = mock['g-r']
    LHS      = 0.7 - 0.032*(mock['M_r,0.1']+16.5) #Weinmann 2006
    blue     = np.where((color<LHS) & (mock['M_r,0.1']!=-99))[0]
    red      = np.where((color>LHS) & (mock['M_r,0.1']!=-99))[0]

    Ngal = np.zeros(len(mock))
    Ngal[occupied]=1
    host_ID = mock['ID_host']
    host_ID[host]=mock['ID_halo'][host]
    
    #All galaxies
    mask = np.zeros(len(mock))
    mask[:] = 1
    avg_N = hod(host_ID,mock['M_host'],Ngal,bins,mask)
    #Red Galaxies
    mask = np.zeros(len(mock))
    mask[red] = 1
    avg_N_red = hod(host_ID,mock['M_host'],Ngal,bins,mask)
    #Blue Galaxies
    mask = np.zeros(len(mock))
    mask[blue] = 1
    avg_N_blue = hod(host_ID,mock['M_host'],Ngal,bins,mask)

    #plot the results
    ax=axes[0]
    p1a,=ax.plot(10**bin_centers, avg_N, color='black')
    ax=axes[1]
    p2a,=ax.plot(10**bin_centers, avg_N_red, color='red')
    ax=axes[2]
    p3a,=ax.plot(10**bin_centers, avg_N_blue, color='blue')
    ax=axes[3]
    ax.plot(10**bin_centers, avg_N, color='black')
    ax=axes[4]
    ax.plot(10**bin_centers, avg_N_red, color='red')
    ax=axes[5]
    ax.plot(10**bin_centers, avg_N_blue, color='blue')
    ax=axes[6]
    ax.plot(10**bin_centers, avg_N, color='black')
    ax=axes[7]
    ax.plot(10**bin_centers, avg_N_red, color='red')
    ax=axes[8]
    ax.plot(10**bin_centers, avg_N_blue, color='blue')
    
    ##########################
    #calculate ideal groups HOD
    ########################## 
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/ideal_groups/'
    print 'opening mock catalogue:', catalogue+'_groups.hdf5'
    #open catalogue
    f1 = h5py.File(filepath_mock+catalogue+'_groups.hdf5', 'r') #open catalogue file
    GC = f1.get(catalogue+'_groups')
    GC = np.array(GC)
    #print 'length:', len(GC)
    #for name in GC.dtype.names: print '     ', name

    bins = np.arange(10,15,0.2)
    bin_centers = (bins[:-1]+bins[1:])/2.0

    centrals   = np.where(GC['HALO_RANK']==1)[0]
    satellites = np.where(GC['HALO_RANK']==0)[0]
    color = GC['M_g,0.1']-GC['M_r,0.1']
    LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
    blue  = np.where(color<LHS)[0] 
    red   = np.where(color>LHS)[0]

    Ngal = np.zeros(len(GC))
    Ngal[:] = 1
    host_ID = GC['GROUP_ID']
    
    #mfactor=-1.0*np.log10(1.0/0.7)
    mfactor = 0.0
    
    #All galaxies
    mask = np.zeros(len(GC))
    mask[:] = 1
    avg_N_all = hod(host_ID,GC['MGROUP']+mfactor,Ngal,bins,mask)
    #Red Galaxies
    mask = np.zeros(len(GC))
    mask[red] = 1
    avg_N_red = hod(host_ID,GC['MGROUP']+mfactor,Ngal,bins,mask)
    #Blue Galaxies
    mask = np.zeros(len(GC))
    mask[blue] = 1
    avg_N_blue = hod(host_ID,GC['MGROUP']+mfactor,Ngal,bins,mask)

    #plot the results
    alpha=1
    ax=axes[0]
    pppp, =ax.plot(10**bin_centers, avg_N_all, ':', color='black', alpha=alpha, fillstyle='none')
    ax=axes[1]
    p2c,=ax.plot(10**bin_centers, avg_N_red, ':', color='red', alpha=alpha, fillstyle='none')
    ax=axes[2]
    p3a,=ax.plot(10**bin_centers, avg_N_blue, ':', color='blue', alpha=alpha, fillstyle='none')
    ax=axes[3]
    ax.plot(10**bin_centers, avg_N_all, ':', color='black', alpha=alpha, fillstyle='none')
    ax=axes[4]
    ax.plot(10**bin_centers, avg_N_red, ':', color='red', alpha=alpha, fillstyle='none')
    ax=axes[5]
    ax.plot(10**bin_centers, avg_N_blue, ':', color='blue', alpha=alpha, fillstyle='none')
    ax=axes[6]
    ax.plot(10**bin_centers, avg_N_all, ':', color='black', alpha=alpha, fillstyle='none')
    ax=axes[7]
    ax.plot(10**bin_centers, avg_N_red, ':', color='red', alpha=alpha, fillstyle='none')
    ax=axes[8]
    ax.plot(10**bin_centers, avg_N_blue, ':', color='blue', alpha=alpha, fillstyle='none')

    ##########################
    #calculate mock groups HOD
    ##########################
    group_cat = 'berlind_groupcat'
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    filepath_GC   = cu.get_output_path() + 'processed_data/'+group_cat+'/mock_runs/4th_run/custom_catalogues/'
    if group_cat == 'tinker_groupcat': group_catalogue = catalogue + '_clf_groups_M19'
    if group_cat == 'berlind_groupcat': group_catalogue = catalogue + '_groups'
    if group_cat == 'yang_groupcat': group_catalogue = catalogue + '_groups'
    #set up arrays to store results
    avg_N_all  = np.zeros((N_boots,len(bins)-1))
    avg_N_red  = np.zeros((N_boots,len(bins)-1))
    avg_N_blue = np.zeros((N_boots,len(bins)-1))
    for boot in range(0,N_boots):
        #print 'opening mock group catalogue:', filepath_GC+'bootstraps/'+group_catalogue+'_'+str(boot)+'.hdf5'
        f2 = h5py.File(filepath_GC+'bootstraps/'+group_catalogue+'_'+str(boot)+'.hdf5', 'r') #open catalogue file
        GC = f2.get(group_catalogue+'_'+str(boot))
        GC = np.array(GC)
        #print len(GC)
        #for name in GC.dtype.names: print '     ', name

        centrals   = np.where(GC['RANK']==1)[0]
        satellites = np.where(GC['RANK']==0)[0]
        color = GC['M_g,0.1']-GC['M_r,0.1']
        LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
        blue  = np.where(color<LHS)[0] 
        red   = np.where(color>LHS)[0]

        Ngal = np.zeros(len(GC))
        Ngal[:] = 1
        host_ID = GC['GROUP_ID']
    
        #All galaxies
        mask = np.zeros(len(GC))
        mask[:] = 1
        avg_N_all[boot,:] = hod(host_ID,GC['MGROUP'],Ngal,bins,mask)
        #Red Galaxies
        mask = np.zeros(len(GC))
        mask[red] = 1
        avg_N_red[boot,:] = hod(host_ID,GC['MGROUP'],Ngal,bins,mask)
        #Blue Galaxies
        mask = np.zeros(len(GC))
        mask[blue] = 1
        avg_N_blue[boot,:] = hod(host_ID,GC['MGROUP'],Ngal,bins,mask)
    
    #calculate the error bars
    avg_N_all_err  = np.std(avg_N_all,axis=0)
    avg_N_red_err  = np.std(avg_N_red,axis=0)
    avg_N_blue_err = np.std(avg_N_blue,axis=0)
    avg_N_all  = np.mean(avg_N_all,axis=0)
    avg_N_red  = np.mean(avg_N_red,axis=0)
    avg_N_blue = np.mean(avg_N_blue,axis=0)

    #plot the results with error bars ontop of the mock result
    ax = axes[0]
    p1b = ax.errorbar(10**bin_centers, avg_N_all, yerr=avg_N_all_err, fmt='o', color='black', ms=3, mec='none')
    ax.legend((p1a,p1b),('mock','group catalogue'), loc='lower right', fontsize=10, numpoints=1, frameon=False )
    ax = axes[1]
    p2b = ax.errorbar(10**bin_centers, avg_N_red, yerr=avg_N_red_err, fmt='o', color='red', ms=3, mec='none')
    ax = axes[2]
    p3b = ax.errorbar(10**bin_centers, avg_N_blue, yerr=avg_N_blue_err, fmt='o', color='blue', ms=3, mec='none')

    ##########################
    #calculate mock groups HOD
    ##########################
    group_cat = 'tinker_groupcat'
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    filepath_GC   = cu.get_output_path() + 'processed_data/'+group_cat+'/mock_runs/4th_run/custom_catalogues/'
    if group_cat == 'tinker_groupcat': group_catalogue = catalogue + '_clf_groups_M19'
    if group_cat == 'berlind_groupcat': group_catalogue = catalogue + '_groups'
    if group_cat == 'yang_groupcat': group_catalogue = catalogue + '_groups'
    #set up arrays to store results
    avg_N_all  = np.zeros((N_boots,len(bins)-1))
    avg_N_red  = np.zeros((N_boots,len(bins)-1))
    avg_N_blue = np.zeros((N_boots,len(bins)-1))
    for boot in range(0,N_boots):
        #print 'opening mock group catalogue:', filepath_GC+'bootstraps/'+group_catalogue+'_'+str(boot)+'.hdf5'
        f2 = h5py.File(filepath_GC+'bootstraps/'+group_catalogue+'_'+str(boot)+'.hdf5', 'r') #open catalogue file
        GC = f2.get(group_catalogue+'_'+str(boot))
        GC = np.array(GC)
        #print len(GC)
        #for name in GC.dtype.names: print '     ', name

        centrals   = np.where(GC['RANK']==1)[0]
        satellites = np.where(GC['RANK']==0)[0]
        color = GC['M_g,0.1']-GC['M_r,0.1']
        LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
        blue  = np.where(color<LHS)[0] 
        red   = np.where(color>LHS)[0]

        Ngal = np.zeros(len(GC))
        Ngal[:] = 1
        host_ID = GC['GROUP_ID']
    
        #All galaxies
        mask = np.zeros(len(GC))
        mask[:] = 1
        avg_N_all[boot,:] = hod(host_ID,GC['MGROUP'],Ngal,bins,mask)
        #Red Galaxies
        mask = np.zeros(len(GC))
        mask[red] = 1
        avg_N_red[boot,:] = hod(host_ID,GC['MGROUP'],Ngal,bins,mask)
        #Blue Galaxies
        mask = np.zeros(len(GC))
        mask[blue] = 1
        avg_N_blue[boot,:] = hod(host_ID,GC['MGROUP'],Ngal,bins,mask)
    
    #calculate the error bars
    avg_N_all_err  = np.std(avg_N_all,axis=0)
    avg_N_red_err  = np.std(avg_N_red,axis=0)
    avg_N_blue_err = np.std(avg_N_blue,axis=0)
    avg_N_all  = np.mean(avg_N_all,axis=0)
    avg_N_red  = np.mean(avg_N_red,axis=0)
    avg_N_blue = np.mean(avg_N_blue,axis=0)

    #plot the results with error bars ontop of the mock result
    ax = axes[3]
    p1b = ax.errorbar(10**bin_centers, avg_N_all, yerr=avg_N_all_err, fmt='o', color='black', ms=3, mec='none')
    ax.legend([pppp],['limit of ``perfect" \n group membership'], loc='lower right', fontsize=10, numpoints=1, frameon=False )
    ax = axes[4]
    p2b = ax.errorbar(10**bin_centers, avg_N_red, yerr=avg_N_red_err, fmt='o', color='red', ms=3, mec='none')
    ax = axes[5]
    p3b = ax.errorbar(10**bin_centers, avg_N_blue, yerr=avg_N_blue_err, fmt='o', color='blue', ms=3, mec='none')

    ##########################
    #calculate mock groups HOD
    ##########################
    group_cat = 'yang_groupcat'
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    filepath_GC   = cu.get_output_path() + 'processed_data/'+group_cat+'/mock_runs/4th_run/version_2/custom_catalogues/'
    if group_cat == 'tinker_groupcat': group_catalogue = catalogue + '_clf_groups_M19'
    if group_cat == 'berlind_groupcat': group_catalogue = catalogue + '_groups'
    if group_cat == 'yang_groupcat': group_catalogue = catalogue + '_groups'
    #set up arrays to store results
    avg_N_all  = np.zeros((N_boots,len(bins)-1))
    avg_N_red  = np.zeros((N_boots,len(bins)-1))
    avg_N_blue = np.zeros((N_boots,len(bins)-1))
    for boot in range(0,N_boots):
        #print 'opening mock group catalogue:', filepath_GC+'bootstraps/'+group_catalogue+'_'+str(boot)+'.hdf5'
        f2 = h5py.File(filepath_GC+'bootstraps/'+group_catalogue+'_'+str(boot)+'.hdf5', 'r') #open catalogue file
        GC = f2.get(group_catalogue+'_'+str(boot))
        GC = np.array(GC)
        #print len(GC)
        #for name in GC.dtype.names: print '     ', name

        centrals   = np.where(GC['RANK']==1)[0]
        satellites = np.where(GC['RANK']==0)[0]
        color = GC['M_g,0.1']-GC['M_r,0.1']
        LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
        blue  = np.where(color<LHS)[0] 
        red   = np.where(color>LHS)[0]

        Ngal = np.zeros(len(GC))
        Ngal[:] = 1
        host_ID = GC['GROUP_ID']
    
        #All galaxies
        mask = np.zeros(len(GC))
        mask[:] = 1
        avg_N_all[boot,:] = hod(host_ID,GC['MGROUP'],Ngal,bins,mask)
        #Red Galaxies
        mask = np.zeros(len(GC))
        mask[red] = 1
        avg_N_red[boot,:] = hod(host_ID,GC['MGROUP'],Ngal,bins,mask)
        #Blue Galaxies
        mask = np.zeros(len(GC))
        mask[blue] = 1
        avg_N_blue[boot,:] = hod(host_ID,GC['MGROUP'],Ngal,bins,mask)
    
    #calculate the error bars
    avg_N_all_err  = np.std(avg_N_all,axis=0)
    avg_N_red_err  = np.std(avg_N_red,axis=0)
    avg_N_blue_err = np.std(avg_N_blue,axis=0)
    avg_N_all  = np.mean(avg_N_all,axis=0)
    avg_N_red  = np.mean(avg_N_red,axis=0)
    avg_N_blue = np.mean(avg_N_blue,axis=0)

    #plot the results with error bars ontop of the mock result
    ax = axes[6]
    p1b = ax.errorbar(10**(bin_centers), avg_N_all, yerr=avg_N_all_err, fmt='o', color='black', ms=3, mec='none')
    ax = axes[7]
    p2b = ax.errorbar(10**(bin_centers), avg_N_red, yerr=avg_N_red_err, fmt='o', color='red', ms=3, mec='none')
    ax = axes[8]
    p3b = ax.errorbar(10**(bin_centers), avg_N_blue, yerr=avg_N_blue_err, fmt='o', color='blue', ms=3, mec='none')

    plt.show(block=True)
    print plotpath+filename+'.pdf'
    fig.savefig(plotpath+filename+'.pdf')

if __name__ == '__main__':
    main() 
