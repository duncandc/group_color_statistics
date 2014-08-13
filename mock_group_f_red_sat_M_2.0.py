import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
import custom_utilities as cu
from astropy import cosmology
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable

def main():
    catalogue = sys.argv[1]
    plotpath = cu.get_plot_path() + 'analysis/groupcats/'

    fig1, axes = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, figsize=(6.95, 3.3))
    fig1.subplots_adjust(hspace=0, wspace=0.05)
    fig1.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.9)
    axes = axes.flatten()
    filename = catalogue+'_mock_group_f_red_sat_M.eps'

    N_boots = 50

    ###################################################################################################################################

    group_cat = 'berlind'
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    filepath_cat  = cu.get_output_path() + 'processed_data/'+group_cat+'_groupcat/mock_runs/4th_run/custom_catalogues/'

    if group_cat=='tinker': group_catalogue = catalogue+'_clf_groups_M19'
    if group_cat=='berlind': group_catalogue = catalogue+'_groups'

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

    #bins = np.arange(9.5,10.8,0.1)
    bins = np.arange(11,15,0.2)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    result = np.digitize(mock['M_host'],bins=bins)
     
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
        if (len(centrals_in_bin))>0:
            #f_sat_red[i] = float(len(red_satellites))/(len(centrals_in_bin)+len(satellites_in_bin))
            #f_sat_blue[i] = float(len(blue_satellites))/(len(centrals_in_bin)+len(satellites_in_bin))
            f_sat_red[i] = float(len(red_satellites))/(len(red_centrals)+len(red_satellites))
            f_sat_blue[i] = float(len(blue_satellites))/(len(blue_centrals)+len(blue_satellites))

    ax=axes[0]
    ax.set_xlim([11,15])
    ax.set_ylim([0,1])
    ax.set_ylabel(r'$f_{red}$')
    ax.set_xlabel(r'$log(M/[M_{\odot}h^{-1}])$')
    p1a, = ax.plot(bin_centers,f_red_cen,color='orange')
    p2a, = ax.plot(bin_centers,f_red_sat,color='green')

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
    
        #S_r = 4.64
        #L = solar_lum(GC['M_r,0.1'],S_r)
    
        result = np.digitize(GC['MGROUP'],bins=bins)

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
            if (len(centrals_in_bin))>0:
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
    ax.legend((p1a,p2a),('halo cen','halo sat'), loc='lower right', fontsize=10, numpoints=1, frameon=False)
    ax.set_xticks([9.6,9.8,10.0,10.2,10.4,10.6])
    ax.set_title(r'Berlind FoF groups')

#######################################################################################################################
    group_cat = 'tinker'
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    filepath_cat  = cu.get_output_path() + 'processed_data/'+group_cat+'_groupcat/mock_runs/4th_run/custom_catalogues/'

    if group_cat=='tinker': group_catalogue = catalogue+'_clf_groups_M19'
    if group_cat=='berlind': group_catalogue = catalogue+'_groups'

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

    #bins = np.arange(9.5,10.8,0.1)
    bins = np.arange(11,15,0.2)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    result = np.digitize(mock['M_host'],bins=bins)
     
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
        if (len(centrals_in_bin))>0:
            #f_sat_red[i] = float(len(red_satellites))/(len(centrals_in_bin)+len(satellites_in_bin))
            #f_sat_blue[i] = float(len(blue_satellites))/(len(centrals_in_bin)+len(satellites_in_bin))
            f_sat_red[i] = float(len(red_satellites))/(len(red_centrals)+len(red_satellites))
            f_sat_blue[i] = float(len(blue_satellites))/(len(blue_centrals)+len(blue_satellites))

    ax=axes[1]
    ax.set_xlim([11,15])
    ax.set_ylim([0,1])
    #ax.set_xticks([9.6,9.8,10.0,10.2,10.4,10.6])
    #ax.set_ylabel(r'$f_{red}$')
    ax.set_xlabel(r'$log(M/[M_{\odot}h^{-1}])$')
    p1a, = ax.plot(bin_centers,f_red_cen,color='orange')
    p2a, = ax.plot(bin_centers,f_red_sat,color='green')

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
    
        #S_r = 4.64
        #L = solar_lum(GC['M_r,0.1'],S_r)
    
        result = np.digitize(GC['MGROUP'],bins=bins)

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
            if (len(centrals_in_bin))>0:
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
    ax.set_xticks([12.5,13,13.5,14,14.5])
    ax.set_xlim([12,15])
    ax.set_title(r'Tinker SO groups')
    ax.legend((p3a,p4a),('groups cen','groups sat'), loc='lower right', fontsize=10, numpoints=1, frameon=False)

    '''
    #calculate with correction
    x=0.8
    y=0.93
    f_sat_1 = (y*f_sat+(x-1)*f_cen)/(x*y+(x-1)*(1-y))
    f_cen_1 = (f_cen-(1-y)*f_sat_1)/y
    p4a_1, = ax.plot(bin_centers,f_sat_1,'o',color='green', mfc='none', mec='green')
    p4a_1, = ax.plot(bin_centers,f_cen_1,'o',color='orange', mfc='none', mec='orange')
    '''
 

#######################################################################################################################
    group_cat = 'yang'
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    filepath_cat  = cu.get_output_path() + 'processed_data/'+group_cat+'_groupcat/mock_runs/4th_run/custom_catalogues/'

    if group_cat=='tinker': group_catalogue = catalogue+'_clf_groups_M19'
    if group_cat=='berlind': group_catalogue = catalogue+'_groups'
    if group_cat=='yang': group_catalogue = catalogue+'_groups'

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

    #bins = np.arange(9.5,10.8,0.1)
    bins = np.arange(11,15,0.2)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    result = np.digitize(mock['M_host'],bins=bins)
     
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
        if (len(centrals_in_bin))>0:
            #f_sat_red[i] = float(len(red_satellites))/(len(centrals_in_bin)+len(satellites_in_bin))
            #f_sat_blue[i] = float(len(blue_satellites))/(len(centrals_in_bin)+len(satellites_in_bin))
            f_sat_red[i] = float(len(red_satellites))/(len(red_centrals)+len(red_satellites))
            f_sat_blue[i] = float(len(blue_satellites))/(len(blue_centrals)+len(blue_satellites))

    ax=axes[2]
    ax.set_xlim([11,15])
    ax.set_ylim([0,1])
    #ax.set_xticks([9.6,9.8,10.0,10.2,10.4,10.6])
    #ax.set_ylabel(r'$f_{red}$')
    ax.set_xlabel(r'$log(M/[M_{\odot}h^{-1}])$')
    p1a, = ax.plot(bin_centers,f_red_cen,color='orange')
    p2a, = ax.plot(bin_centers,f_red_sat,color='green')

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
    
        #S_r = 4.64
        #L = solar_lum(GC['M_r,0.1'],S_r)
    
        result = np.digitize(GC['MGROUP'],bins=bins)

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
            if (len(centrals_in_bin))>0:
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
    ax.set_xticks([12.5,13,13.5,14,14.5])
    ax.set_xlim([12,15])
    ax.set_title(r'Yang SO groups')
   
    plt.show()

    fig1.savefig(plotpath+filename)
  
          
def solar_lum(M,Msol):
    L = ((Msol-M)/2.5)
    return L  

if __name__ == '__main__':
    main() 
