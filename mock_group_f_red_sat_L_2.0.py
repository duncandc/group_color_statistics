#!/usr/bin/env python

#Duncan Campbell
#September 2, 2012
#Yale University
#Plot the red fraction of satellites and centrals as a function of galaxy luminosity.
#Includes group finder runs on mock results, intrinsic mock results, and bootstrapped 
#error bars.

import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu
import sys

def main():

    #process user input
    if len(sys.argv)==2:
        catalogue = sys.argv[1]
    else:
        catalogue = 'Mr19_age_distribution_matching_mock_sys_empty_shuffle_satrel_shuffle'
    print "running for:", catalogue

    #setup figure
    fig1, axes = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True, 
                              figsize=(6.95, 6.6-0.99))
    fig1.subplots_adjust(hspace=0, wspace=0.05)
    fig1.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
    axes = axes.flatten()

    #get figure saving information
    plotpath = cu.get_plot_path() + 'analysis/groupcats/'
    filename = catalogue+'_f_red_sat_L'

    #number of bootstraps per group catalogue.
    N_boots = 3

#run code for intrinsic mock results
##########################################################################################

    filepath_mock = cu.get_output_path() + \
        'processed_data/hearin_mocks/custom_catalogues/'

    print 'opening mock catalogue:', catalogue+'.hdf5'
    #open catalogue
    f1 = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f1.get(catalogue)
    
    centrals   = np.where(mock['ID_host']==-1)
    satellites = np.where(mock['ID_host']!=-1)
    
    #galaxy color
    color = mock['g-r']
    LHS   = 0.21-0.03*mock['M_r,0.1']
    blue  = np.where(color<LHS)[0] #indices of blue galaxies
    red   = np.where(color>LHS)[0] #indicies of red galaxies
    
    S_r = 4.64
    L = solar_lum(mock['M_r,0.1'],S_r)

    bins = np.arange(9.5,10.8,0.1)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    result = np.digitize(L,bins=bins)
     
    f_red_cen = np.zeros(len(bin_centers))
    f_red_sat = np.zeros(len(bin_centers))
    f_sat_red = np.zeros(len(bin_centers))
    f_sat_blue = np.zeros(len(bin_centers))
    for i in range(0,len(bins)-1):
        print i
        ind = np.where(result==i+1)[0]
        centrals_in_bin = np.in1d(ind,centrals)
        centrals_in_bin = ind[centrals_in_bin]
        satellites_in_bin = np.in1d(ind,satellites)
        satellites_in_bin = ind[satellites_in_bin]

        red_centrals = np.in1d(centrals_in_bin,red)
        red_centrals = centrals_in_bin[red_centrals]
        red_satellites = np.in1d(satellites_in_bin,red)
        red_satellites = satellites_in_bin[red_satellites]

        blue_centrals = np.in1d(centrals_in_bin,blue)
        blue_centrals = centrals_in_bin[blue_centrals]
        blue_satellites = np.in1d(satellites_in_bin,blue)
        blue_satellites = satellites_in_bin[blue_satellites]

        if (len(red_centrals)+len(blue_centrals)) > 0:
            f_red_cen[i] = float(len(red_centrals))/(len(red_centrals)+len(blue_centrals))
        if (len(red_satellites)+len(blue_satellites)) > 0:
            f_red_sat[i] = float(len(red_satellites))/(len(red_satellites)+len(blue_satellites))

        #f_sat_red[i] = float(len(red_satellites))/(len(centrals_in_bin)+len(satellites_in_bin))
        #f_sat_blue[i] = float(len(blue_satellites))/(len(centrals_in_bin)+len(satellites_in_bin))
        f_sat_red[i] = float(len(red_satellites))/(len(red_centrals)+len(red_satellites))
        f_sat_blue[i] = float(len(blue_satellites))/(len(blue_centrals)+len(blue_satellites))

    ax=axes[0]
    ax.set_xlim([9.5,10.8])
    ax.set_ylim([0,1])
    ax.set_ylabel(r'$f_{red}$')
    p1a, = ax.plot(bin_centers,f_red_cen,color='orange')
    p2a, = ax.plot(bin_centers,f_red_sat,color='green')

    ax=axes[3]
    p1b, = ax.plot(bin_centers,f_sat_red,color='red')
    p2b, = ax.plot(bin_centers,f_sat_blue,color='blue')
    
    ax=axes[1]
    ax.set_xlim([9.5,10.8])
    ax.set_ylim([0,1])
    p1a, = ax.plot(bin_centers,f_red_cen,color='orange')
    p2a, = ax.plot(bin_centers,f_red_sat,color='green')
    
    ax=axes[4]
    p1b, = ax.plot(bin_centers,f_sat_red,color='red')
    p2b, = ax.plot(bin_centers,f_sat_blue,color='blue')
    
    ax=axes[2]
    ax.set_xlim([9.5,10.8])
    ax.set_ylim([0,1])
    p1a, = ax.plot(bin_centers,f_red_cen,color='orange')
    p2a, = ax.plot(bin_centers,f_red_sat,color='green')
    
    ax=axes[5]
    p1b, = ax.plot(bin_centers,f_sat_red,color='red')
    p2b, = ax.plot(bin_centers,f_sat_blue,color='blue')

#run code for Berlind group catalogue
##########################################################################################

    group_cat = 'berlind'
    filepath_cat  = cu.get_output_path() + \
        'processed_data/'+group_cat+'_groupcat/mock_runs/4th_run/custom_catalogues/'
    group_catalogue = catalogue+'_groups'

    f_red_cen = np.zeros((N_boots,len(bin_centers)))
    f_red_sat = np.zeros((N_boots,len(bin_centers)))
    f_sat_red = np.zeros((N_boots,len(bin_centers)))
    f_sat_blue = np.zeros((N_boots,len(bin_centers)))
    for boot in range(0,N_boots):

        #open catalogue
        catalogue_1 = group_catalogue+'_'+str(boot)
        print 'opening group catalogue:', group_catalogue
        f =  h5py.File(filepath_cat+'bootstraps/'+catalogue_1+'.hdf5', 'r') #open catalogue file
        GC = f.get(catalogue_1)

        centrals   = np.where(GC['RPROJ']==0)
        satellites = np.where(GC['RPROJ']>0)
    
        #galaxy color
        color = GC['M_g,0.1']-GC['M_r,0.1']
        LHS   = 0.21-0.03*GC['M_r,0.1']
        blue  = np.where(color<LHS)[0] #indices of blue galaxies
        red   = np.where(color>LHS)[0] #indicies of red galaxies
    
        S_r = 4.64
        L = solar_lum(GC['M_r,0.1'],S_r)
    
        result = np.digitize(L,bins=bins)

        for i in range(0,len(bins)-1):
            print i
            ind = np.where(result==i+1)[0]
            centrals_in_bin = np.in1d(ind,centrals)
            centrals_in_bin = ind[centrals_in_bin]
            satellites_in_bin = np.in1d(ind,satellites)
            satellites_in_bin = ind[satellites_in_bin]

            red_centrals = np.in1d(centrals_in_bin,red)
            red_centrals = centrals_in_bin[red_centrals]
            red_satellites = np.in1d(satellites_in_bin,red)
            red_satellites = satellites_in_bin[red_satellites]

            blue_centrals = np.in1d(centrals_in_bin,blue)
            blue_centrals = centrals_in_bin[blue_centrals]
            blue_satellites = np.in1d(satellites_in_bin,blue)
            blue_satellites = satellites_in_bin[blue_satellites]

            if (len(red_centrals)+len(blue_centrals)) > 0:
                f_red_cen[boot,i] = float(len(red_centrals))/(len(red_centrals)+len(blue_centrals))
            if (len(red_satellites)+len(blue_satellites)) > 0:
                f_red_sat[boot,i] = float(len(red_satellites))/(len(red_satellites)+len(blue_satellites))

            #f_sat_red[boot,i] = float(len(red_satellites))/(len(centrals_in_bin)+len(satellites_in_bin))
            #f_sat_blue[boot,i] = float(len(blue_satellites))/(len(centrals_in_bin)+len(satellites_in_bin))
            f_sat_red[boot,i] = float(len(red_satellites))/(len(red_centrals)+len(red_satellites))
            f_sat_blue[boot,i] = float(len(blue_satellites))/(len(blue_centrals)+len(blue_satellites))


    f_sat_red_111 = f_sat_red
    f_sat_blue_111 = f_sat_blue

    f_cen = np.nanmean(f_red_cen, axis=0)
    f_sat = np.nanmean(f_red_sat, axis=0)
    err_cen = np.nanstd(f_red_cen, axis=0)
    err_sat = np.nanstd(f_red_sat, axis=0)

    f_sat_red = np.nanmean(f_sat_red, axis=0)
    f_sat_blue = np.nanmean(f_sat_blue, axis=0)
    err_sat_red = np.nanstd(f_sat_red_111, axis=0)
    err_sat_blue = np.nanstd(f_sat_blue_111, axis=0)

    ax=axes[0]
    p3a=ax.errorbar(bin_centers, f_cen, yerr=err_cen,fmt='o',color='orange', mec='none', ms=3)
    p4a=ax.errorbar(bin_centers, f_sat, yerr=err_sat,fmt='o',color='green', mec='none', ms=3)
    ax.legend((p1a,p2a),('halo central','halo satellite'), loc='lower right', fontsize=10,  numpoints=1, frameon=False)
    ax.set_title(r'Berlind FoF groups')

    ax=axes[3]
    p3b=ax.errorbar(bin_centers, f_sat_red, yerr=err_sat_red, fmt='o', color='red', mec='none', ms=3)
    p4b=ax.errorbar(bin_centers+0.01, f_sat_blue, yerr=err_sat_blue, fmt='o', color='blue', mec='none', ms=3)
    ax.set_xlabel(r'$\log(L)$ $[L_{\odot}]$')
    ax.set_xticks([9.6,9.8,10.0,10.2,10.4,10.6])
    ax.set_xlim([9.5,10.7])
    ax.set_ylabel(r'$f_{sat}$')
    ax.legend((p1b,p2b),('halo red cen/sat','halo blue cen/sat'), loc='upper right', fontsize=10, numpoints=1, frameon=False)

#run code for Tinker group catalogue
##########################################################################################

    group_cat = 'tinker'
    filepath_cat  = cu.get_output_path() + \
        'processed_data/'+group_cat+'_groupcat/mock_runs/4th_run/custom_catalogues/'
     group_catalogue = catalogue+'_clf_groups_M19'

    f_red_cen = np.zeros((N_boots,len(bin_centers)))
    f_red_sat = np.zeros((N_boots,len(bin_centers)))
    f_sat_red = np.zeros((N_boots,len(bin_centers)))
    f_sat_blue = np.zeros((N_boots,len(bin_centers)))
    for boot in range(0,N_boots):

        #open catalogue
        catalogue_1 = group_catalogue+'_'+str(boot)
        print 'opening group catalogue:', group_catalogue
        f =  h5py.File(filepath_cat+'bootstraps/'+catalogue_1+'.hdf5', 'r') #open catalogue file
        GC = f.get(catalogue_1)

        centrals   = np.where(GC['RPROJ']==0)
        satellites = np.where(GC['RPROJ']>0)
    
        #galaxy color
        color = GC['M_g,0.1']-GC['M_r,0.1']
        LHS   = 0.21-0.03*GC['M_r,0.1']
        blue  = np.where(color<LHS)[0] #indices of blue galaxies
        red   = np.where(color>LHS)[0] #indicies of red galaxies
    
        S_r = 4.64
        L = solar_lum(GC['M_r,0.1'],S_r)
    
        result = np.digitize(L,bins=bins)

        for i in range(0,len(bins)-1):
            print i
            ind = np.where(result==i+1)[0]
            centrals_in_bin = np.in1d(ind,centrals)
            centrals_in_bin = ind[centrals_in_bin]
            satellites_in_bin = np.in1d(ind,satellites)
            satellites_in_bin = ind[satellites_in_bin]

            red_centrals = np.in1d(centrals_in_bin,red)
            red_centrals = centrals_in_bin[red_centrals]
            red_satellites = np.in1d(satellites_in_bin,red)
            red_satellites = satellites_in_bin[red_satellites]

            blue_centrals = np.in1d(centrals_in_bin,blue)
            blue_centrals = centrals_in_bin[blue_centrals]
            blue_satellites = np.in1d(satellites_in_bin,blue)
            blue_satellites = satellites_in_bin[blue_satellites]

            if (len(red_centrals)+len(blue_centrals)) > 0:
                f_red_cen[boot,i] = float(len(red_centrals))/(len(red_centrals)+len(blue_centrals))
            if (len(red_satellites)+len(blue_satellites)) > 0:
                f_red_sat[boot,i] = float(len(red_satellites))/(len(red_satellites)+len(blue_satellites))

            #f_sat_red[boot,i] = float(len(red_satellites))/(len(centrals_in_bin)+len(satellites_in_bin))
            #f_sat_blue[boot,i] = float(len(blue_satellites))/(len(centrals_in_bin)+len(satellites_in_bin))
            f_sat_red[boot,i] = float(len(red_satellites))/(len(red_centrals)+len(red_satellites))
            f_sat_blue[boot,i] = float(len(blue_satellites))/(len(blue_centrals)+len(blue_satellites))


    f_sat_red_111 = f_sat_red
    f_sat_blue_111 = f_sat_blue

    f_cen = np.nanmean(f_red_cen, axis=0)
    f_sat = np.nanmean(f_red_sat, axis=0)
    err_cen = np.nanstd(f_red_cen, axis=0)
    err_sat = np.nanstd(f_red_sat, axis=0)

    f_sat_red = np.nanmean(f_sat_red, axis=0)
    f_sat_blue = np.nanmean(f_sat_blue, axis=0)
    err_sat_red = np.nanstd(f_sat_red_111, axis=0)
    err_sat_blue = np.nanstd(f_sat_blue_111, axis=0)

    ax=axes[1]
    p3a=ax.errorbar(bin_centers, f_cen, yerr=err_cen,fmt='o',color='orange', mec='none', ms=3)
    p4a=ax.errorbar(bin_centers, f_sat, yerr=err_sat,fmt='o',color='green', mec='none', ms=3)
    ax.legend((p3a,p4a),('group central','group satellite'), loc='lower right', fontsize=10, numpoints=1, frameon=False)
    ax.set_title(r'Tinker SO groups')

    ax=axes[4]
    p3b=ax.errorbar(bin_centers, f_sat_red, yerr=err_sat_red, fmt='o', color='red', mec='none', ms=3)
    p4b=ax.errorbar(bin_centers+0.01, f_sat_blue, yerr=err_sat_blue, fmt='o', color='blue', mec='none', ms=3)
    ax.set_xlabel(r'$\log(L)$ $[L_{\odot}]$')
    ax.set_xticks([9.6,9.8,10.0,10.2,10.4,10.6])
    ax.set_xlim([9.5,10.7])
    ax.legend((p3b,p4b),('group red cen/sat','group blue cen/sat'), loc='upper right', fontsize=10, numpoints=1, frameon=False)

#run code for Yang group catalogue
##########################################################################################

    group_cat = 'yang'
    filepath_cat  = cu.get_output_path() + \
        'processed_data/'+group_cat+'_groupcat/mock_runs/4th_run/custom_catalogues/'
    group_catalogue = catalogue+'_groups'
    
    f_red_cen = np.zeros((N_boots,len(bin_centers)))
    f_red_sat = np.zeros((N_boots,len(bin_centers)))
    f_sat_red = np.zeros((N_boots,len(bin_centers)))
    f_sat_blue = np.zeros((N_boots,len(bin_centers)))
    for boot in range(0,N_boots):

        #open catalogue
        catalogue_1 = group_catalogue+'_'+str(boot)
        print 'opening group catalogue:', group_catalogue
        f =  h5py.File(filepath_cat+'bootstraps/'+catalogue_1+'.hdf5', 'r') #open catalogue file
        GC = f.get(catalogue_1)

        centrals   = np.where(GC['RPROJ']==0)
        satellites = np.where(GC['RPROJ']>0)
    
        #galaxy color
        color = GC['M_g,0.1']-GC['M_r,0.1']
        LHS   = 0.21-0.03*GC['M_r,0.1']
        blue  = np.where(color<LHS)[0] #indices of blue galaxies
        red   = np.where(color>LHS)[0] #indicies of red galaxies
    
        S_r = 4.64
        L = solar_lum(GC['M_r,0.1'],S_r)
    
        result = np.digitize(L,bins=bins)

        for i in range(0,len(bins)-1):
            print i
            ind = np.where(result==i+1)[0]
            centrals_in_bin = np.in1d(ind,centrals)
            centrals_in_bin = ind[centrals_in_bin]
            satellites_in_bin = np.in1d(ind,satellites)
            satellites_in_bin = ind[satellites_in_bin]

            red_centrals = np.in1d(centrals_in_bin,red)
            red_centrals = centrals_in_bin[red_centrals]
            red_satellites = np.in1d(satellites_in_bin,red)
            red_satellites = satellites_in_bin[red_satellites]

            blue_centrals = np.in1d(centrals_in_bin,blue)
            blue_centrals = centrals_in_bin[blue_centrals]
            blue_satellites = np.in1d(satellites_in_bin,blue)
            blue_satellites = satellites_in_bin[blue_satellites]

            if (len(red_centrals)+len(blue_centrals)) > 0:
                f_red_cen[boot,i] = float(len(red_centrals))/(len(red_centrals)+len(blue_centrals))
            if (len(red_satellites)+len(blue_satellites)) > 0:
                f_red_sat[boot,i] = float(len(red_satellites))/(len(red_satellites)+len(blue_satellites))

            #f_sat_red[boot,i] = float(len(red_satellites))/(len(centrals_in_bin)+len(satellites_in_bin))
            #f_sat_blue[boot,i] = float(len(blue_satellites))/(len(centrals_in_bin)+len(satellites_in_bin))
            f_sat_red[boot,i] = float(len(red_satellites))/(len(red_centrals)+len(red_satellites))
            f_sat_blue[boot,i] = float(len(blue_satellites))/(len(blue_centrals)+len(blue_satellites))


    f_sat_red_111 = f_sat_red
    f_sat_blue_111 = f_sat_blue

    f_cen = np.nanmean(f_red_cen, axis=0)
    f_sat = np.nanmean(f_red_sat, axis=0)
    err_cen = np.nanstd(f_red_cen, axis=0)
    err_sat = np.nanstd(f_red_sat, axis=0)

    f_sat_red = np.nanmean(f_sat_red, axis=0)
    f_sat_blue = np.nanmean(f_sat_blue, axis=0)
    err_sat_red = np.nanstd(f_sat_red_111, axis=0)
    err_sat_blue = np.nanstd(f_sat_blue_111, axis=0)

    ax=axes[2]
    p3a=ax.errorbar(bin_centers, f_cen, yerr=err_cen,fmt='o',color='orange', mec='none', ms=3)
    p4a=ax.errorbar(bin_centers, f_sat, yerr=err_sat,fmt='o',color='green', mec='none', ms=3)
    #ax.legend((p1a,p2a,p3a,p4a),('halo cen','halo sat','groups cen','groups sat'), loc='lower right', fontsize=10)
    ax.set_title(r'Yang SO groups')
    ax.set_yticks([0,0.2,0.4,0.6,0.8])

    ax=axes[5]
    p3b=ax.errorbar(bin_centers, f_sat_red, yerr=err_sat_red, fmt='o', color='red', mec='none', ms=3)
    p4b=ax.errorbar(bin_centers+0.01, f_sat_blue, yerr=err_sat_blue, fmt='o', color='blue', mec='none', ms=3)
    ax.set_xlabel(r'$\log(L)$ $[L_{\odot}]$')
    ax.set_xticks([9.6,9.8,10.0,10.2,10.4,10.6])
    ax.set_xlim([9.5,10.7])
    ax.set_yticks([0,0.2,0.4,0.6,0.8])

#plot results, show, and save.
##########################################################################################

    plt.show()
    fig1.savefig(plotpath+filename+'.pdf', dpi=400, bbox_inches='tight')


def solar_lum(M,Msol):
    L = ((Msol-M)/2.5)
    return L  

if __name__ == '__main__':
    main() 
