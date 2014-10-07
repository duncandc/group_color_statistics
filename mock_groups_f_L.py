import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
import custom_utilities as cu
from astropy import cosmology
import sys

def main():

    catalogue = sys.argv[1]

    #define outputs
    plotpath = cu.get_plot_path() + 'analysis/groupcats/'
    filename = catalogue+'_f_L'

    #set up plot
    fig1, axes = plt.subplots(nrows=2,ncols=3,sharex=True,sharey=True,figsize=(6.95,6.6-1.5))
    fig1.subplots_adjust(hspace=0, wspace=0.05)
    fig1.subplots_adjust(left=0.1, right=0.95, bottom=0.125, top=0.95)
    axes = axes.flatten()
    ax=axes[0]
    ax.set_xlim([9.5,10.8])
    ax.set_ylim([0,1])
    ax.set_ylabel(r'$f_{red}$')
    ax.set_title(r'Berlind FoF groups')
    ax.set_yticks([0,0.2,0.4,0.6,0.8])
    ax=axes[3]
    ax.set_xlabel(r'$\log(L/[L_{\odot}h^{-2}])$')
    ax.set_xticks([9.6,9.8,10.0,10.2,10.4,10.6])
    ax.set_xlim([9.5,10.7])
    ax.set_ylabel(r'$f_{sat}$')
    ax.set_yticks([0,0.2,0.4,0.6,0.8])
    ax=axes[1]
    ax.set_xlim([9.5,10.8])
    ax.set_ylim([0,1])
    ax.set_title(r'Tinker SO groups')
    ax=axes[4]
    ax.set_xlabel(r'$\log(L/[L_{\odot}h^{-2}])$')
    ax.set_xticks([9.6,9.8,10.0,10.2,10.4,10.6])
    ax.set_xlim([9.5,10.7])
    ax.set_yticks([0,0.2,0.4,0.6,0.8])
    ax=axes[2]
    ax.set_xlim([9.5,10.8])
    ax.set_ylim([0,1])
    ax.set_title(r'Yang SO groups')
    ax=axes[5]
    ax.set_xlabel(r'$\log(L/[L_{\odot}h^{-2}])$')
    ax.set_xticks([9.6,9.8,10.0,10.2,10.4,10.6])
    ax.set_xlim([9.5,10.7])

    #define luminosity bins
    bins = np.arange(9.5,10.8,0.1)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    #define constants
    S_r = 4.64
    N_boots = 50 #number of group cat bootstraps

    #open mock catalogue
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'

    print 'opening mock catalogue:', catalogue+'.hdf5'
    f1 = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f1.get(catalogue)
    
    #central/satellites
    centrals_ind   = np.where(mock['ID_host']==-1)[0]
    satellites_ind = np.where(mock['ID_host']!=-1)[0]
    centrals_bool   = (mock['ID_host']==-1)
    satellites_bool = (mock['ID_host']!=-1)

    f_sat_mock = float(len(satellites_ind))/(len(satellites_ind)+len(centrals_ind))
    print f_sat_mock

    #galaxy color
    color = mock['g-r']
    LHS   = 0.7 - 0.032*(mock['M_r,0.1']+16.5) #Weinmann 2006
    blue_ind  = np.where(color<LHS)[0] #indices of blue galaxies
    red_ind   = np.where(color>LHS)[0] #indicies of red galaxies
    blue_bool = (color<LHS) #indices of blue galaxies
    red_bool  = (color>LHS) #indicies of red galaxies

    L = solar_lum(mock['M_r,0.1'],S_r)

    f_red_cen = f_prop(L,bins,red_ind,blue_ind,centrals_bool)
    f_red_sat = f_prop(L,bins,red_ind,blue_ind,satellites_bool)
    f_sat_red = f_prop(L,bins,satellites_ind,centrals_ind,red_bool)
    f_sat_blue = f_prop(L,bins,satellites_ind,centrals_ind,blue_bool)

    ax=axes[0]
    p1a, = ax.plot(bin_centers,f_red_cen,color='orange')
    p2a, = ax.plot(bin_centers,f_red_sat,color='green')
    ax.legend((p1a,p2a),('halo central','halo satellite'), loc='lower right', fontsize=10,  numpoints=1, frameon=False)
    ax=axes[1]
    p1a, = ax.plot(bin_centers,f_red_cen,color='orange')
    p2a, = ax.plot(bin_centers,f_red_sat,color='green')
    ax=axes[2]
    p1a, = ax.plot(bin_centers,f_red_cen,color='orange')
    p2a, = ax.plot(bin_centers,f_red_sat,color='green')
    ax=axes[3]
    p1b, = ax.plot(bin_centers,f_sat_red,color='red')
    p2b, = ax.plot(bin_centers,f_sat_blue,color='blue')
    ax.legend((p1b,p2b),('halo red sat/cen','halo blue sat/cen'), loc='upper right', fontsize=10, numpoints=1, frameon=False)
    ax=axes[4]
    p1b, = ax.plot(bin_centers,f_sat_red,color='red')
    p2b, = ax.plot(bin_centers,f_sat_blue,color='blue')
    ax=axes[5]
    p1b, = ax.plot(bin_centers,f_sat_red,color='red')
    p2b, = ax.plot(bin_centers,f_sat_blue,color='blue')

    f_red_cen_mock = f_red_cen 
    f_red_sat_mock = f_red_sat 
    f_sat_red_mock = f_sat_red 
    f_sat_blue_mock = f_sat_blue
    
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

    centrals_ind   = np.where(GC['RANK']==0)[0]
    satellites_ind = np.where(GC['RANK']!=0)[0]
    centrals_bool   = (GC['RANK']==0)
    satellites_bool = (GC['RANK']!=0)
    Nsat = float(len(satellites_ind))
    Ncen = float(len(centrals_ind))
    N_orphan = len(np.where(GC['CEN_IND']==-99)[0])

    centrals_mock_ind   = np.where(GC['HALO_RANK']==0)[0]
    satellites_mock_ind = np.where(GC['HALO_RANK']==1)[0]
    N_sat_mock = float(len(satellites_mock_ind))
    N_cen_mock = float(len(centrals_mock_ind))
    
    #calculate global satellite fractions
    f_sat_mock = N_sat_mock/(N_sat_mock+N_cen_mock)
    f_sat_ideal_groups = N_sat_ideal_groups/(N_sat_ideal_groups+N_cen_ideal_groups)
    
    f_sat_err = len(np.where((GC['RANK']>0).astype(int)!=GC['HALO_RANK'])[0])/(float(N_sat_ideal_groups))
    print f_sat_mock, f_sat_ideal_groups, f_sat_err
    
    
    #galaxy color
    color = GC['M_g,0.1']-GC['M_r,0.1']
    LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
    blue_ind  = np.where(color<LHS)[0] #indices of blue galaxies
    red_ind   = np.where(color>LHS)[0] #indicies of red galaxies
    blue_bool = (color<LHS) #indices of blue galaxies
    red_bool  = (color>LHS) #indicies of red galaxies
    
    L = solar_lum(GC['M_r,0.1'],S_r)
    
    result = np.digitize(L,bins=bins)
    f_red_cen  = f_prop(L,bins,red_ind,blue_ind,centrals_bool)
    f_red_sat  = f_prop(L,bins,red_ind,blue_ind,satellites_bool)
    f_sat_red  = f_prop(L,bins,satellites_ind,centrals_ind,red_bool)
    f_sat_blue = f_prop(L,bins,satellites_ind,centrals_ind,blue_bool)
    
    ax=axes[0]
    p1a, = ax.plot(bin_centers,f_red_cen,color='orange', alpha=0.5)
    p2a, = ax.plot(bin_centers,f_red_sat,color='green', alpha=0.5)
    ax=axes[1]
    p1a, = ax.plot(bin_centers,f_red_cen,color='orange', alpha=0.5)
    p2a, = ax.plot(bin_centers,f_red_sat,color='green', alpha=0.5)
    ax=axes[2]
    p1a, = ax.plot(bin_centers,f_red_cen,color='orange', alpha=0.5)
    p2a, = ax.plot(bin_centers,f_red_sat,color='green', alpha=0.5)
    ax=axes[3]
    p1b, = ax.plot(bin_centers,f_sat_red,color='red', alpha=0.5)
    p2b, = ax.plot(bin_centers,f_sat_blue,color='blue', alpha=0.5)
    ax=axes[4]
    p1b, = ax.plot(bin_centers,f_sat_red,color='red', alpha=0.5)
    p2b, = ax.plot(bin_centers,f_sat_blue,color='blue', alpha=0.5)
    ax=axes[5]
    p1b, = ax.plot(bin_centers,f_sat_red,color='red', alpha=0.5)
    p2b, = ax.plot(bin_centers,f_sat_blue,color='blue', alpha=0.5)

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
    
        L = solar_lum(GC['M_r,0.1'],S_r)
    
        result = np.digitize(L,bins=bins)
        f_red_cen[boot,:]  = f_prop(L,bins,red_ind,blue_ind,centrals_bool)
        f_red_sat[boot,:]  = f_prop(L,bins,red_ind,blue_ind,satellites_bool)
        f_sat_red[boot,:]  = f_prop(L,bins,satellites_ind,centrals_ind,red_bool)
        f_sat_blue[boot,:] = f_prop(L,bins,satellites_ind,centrals_ind,blue_bool)
    
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
    ax=axes[3]
    p3b=ax.errorbar(bin_centers, f_sat_red, yerr=err_sat_red, fmt='o', color='red', mec='none', ms=3)
    p4b=ax.errorbar(bin_centers+0.01, f_sat_blue, yerr=err_sat_blue, fmt='o', color='blue', mec='none', ms=3)

    '''
    #calculate with correction
    catalogue_1 = group_catalogue
    f =  h5py.File(filepath_cat+catalogue_1+'.hdf5', 'r') #open catalogue file
    GC = f.get(catalogue_1)
    centrals_ind   = np.where(GC['RANK']==1)[0]
    satellites_ind = np.where(GC['RANK']!=1)[0]
    centrals_bool   = (GC['RANK']==1)
    satellites_bool = (GC['RANK']!=1)
    centrals_mock_ind   = np.where(GC['HALO_RANK']==1)[0]
    satellites_mock_ind = np.where(GC['HALO_RANK']!=1)[0]
    N_sat_mock = float(len(satellites_mock_ind))
    N_cen_mock = float(len(centrals_mock_ind))
    N_sat = float(len(satellites_ind))
    N_cen = float(len(centrals_ind))
    b_sat=sat_bias(GC)
    Ps=sat_purity(GC)
    Pc=cen_purity(GC)
    Cs=sat_comp(GC)
    Cc=cen_comp(GC)
    print Ps, Cs, Pc, Cc, b_sat
    print 'test:',(1-Ps)*N_sat_mock*Cs*Pc==(1-Pc)*N_cen_mock*b_sat*Cc*Ps
    print 'test:',N_sat*Ps/Cs==N_sat_mock
    print 'test:',N_cen*Pc/Cc==N_cen_mock
    print 'test:',N_cen+N_sat==N_cen_mock+N_sat_mock
    f_red_sat_1 = (f_red_sat-f_red_cen*(((1-Ps)*(1-Cc))/((1-Pc)*Cc*b_sat)))/(1-(((1-Cs)*(1-Cc))/(Ps*Pc)))/Ps
    f_red_cen_1 = (f_red_cen-f_red_sat_1*(((1-Cs)*(1-Pc)*Ps*b_sat)/((1-Ps)*Cs)))/Pc
    ax=axes[0]
    p4a_1, = ax.plot(bin_centers,f_red_sat_1,'o',color='green', mfc='none', mec='green', ms=3)
    p4a_1, = ax.plot(bin_centers,f_red_cen_1,'o',color='orange', mfc='none', mec='orange', ms=3)
    f_sat_red_1 =f_sat_red*b_sat
    f_sat_blue_1 =f_sat_blue*b_sat
    ax=axes[3]
    N_red_sat = f_red_sat_1 * N_sat*Ps/Cs
    N_red_cen = f_red_cen_1 * N_cen*Pc/Cc
    N_blue_sat = (1-f_red_sat_1) * N_sat*Ps/Cs
    N_blue_cen = (1-f_red_cen_1) * N_cen*Pc/Cc
    f_sat_red_1 = N_red_sat/(N_red_cen+N_red_sat)
    f_sat_blue_1 = N_blue_sat/(N_blue_cen+N_blue_sat)
    p4b_1, = ax.plot(bin_centers,f_sat_red_1,'o',color='red', mfc='none', mec='red', ms=3)
    p4b_1, = ax.plot(bin_centers,f_sat_blue_1,'o',color='blue', mfc='none', mec='blue', ms=3)
    '''

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
    
        #galaxy color
        color = GC['M_g,0.1']-GC['M_r,0.1']
        LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
        blue_ind  = np.where(color<LHS)[0] #indices of blue galaxies
        red_ind   = np.where(color>LHS)[0] #indicies of red galaxies
        blue_bool = (color<LHS) #indices of blue galaxies
        red_bool  = (color>LHS) #indicies of red galaxies
    
        L = solar_lum(GC['M_r,0.1'],S_r)
    
        result = np.digitize(L,bins=bins)
        f_red_cen[boot,:]  = f_prop(L,bins,red_ind,blue_ind,centrals_bool)
        f_red_sat[boot,:]  = f_prop(L,bins,red_ind,blue_ind,satellites_bool)
        f_sat_red[boot,:]  = f_prop(L,bins,satellites_ind,centrals_ind,red_bool)
        f_sat_blue[boot,:] = f_prop(L,bins,satellites_ind,centrals_ind,blue_bool)
    
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
    ax.legend((p3a,p4a),('group central','group satellite'), loc='lower right', fontsize=10, numpoints=1, frameon=False)
    ax=axes[4]
    p3b=ax.errorbar(bin_centers, f_sat_red, yerr=err_sat_red, fmt='o', color='red', mec='none', ms=3)
    p4b=ax.errorbar(bin_centers+0.01, f_sat_blue, yerr=err_sat_blue, fmt='o', color='blue', mec='none', ms=3)
    ax.legend((p3b,p4b),('group red cen/sat','group blue cen/sat'), loc='upper right', fontsize=10, numpoints=1, frameon=False)

    '''
    #calculate with correction
    catalogue_1 = group_catalogue
    f =  h5py.File(filepath_cat+catalogue_1+'.hdf5', 'r') #open catalogue file
    GC = f.get(catalogue_1)
    centrals_ind   = np.where(GC['RANK']==1)[0]
    satellites_ind = np.where(GC['RANK']!=1)[0]
    centrals_bool   = (GC['RANK']==1)
    satellites_bool = (GC['RANK']!=1)
    centrals_mock_ind   = np.where(GC['HALO_RANK']==1)[0]
    satellites_mock_ind = np.where(GC['HALO_RANK']!=1)[0]
    N_sat_mock = float(len(satellites_mock_ind))
    N_cen_mock = float(len(centrals_mock_ind))
    N_sat = float(len(satellites_ind))
    N_cen = float(len(centrals_ind))
    x=0.8
    y=0.93
    f_sat = (1-y)/((1-x)+(1-y))
    f_red_sat_1 = (y*f_red_sat+(x-1)*f_red_cen)/(x*y+(x-1)*(1-y))
    f_red_cen_1 = (f_red_cen-(1-y)*f_red_sat_1)/y
    ax=axes[1]
    p4a_1, = ax.plot(bin_centers,f_red_sat_1,'x',color='green', mfc='none', mec='green', ms=3)
    p4a_1, = ax.plot(bin_centers,f_red_cen_1,'x',color='orange', mfc='none', mec='orange', ms=3)
    N_red_sat = f_red_sat_1 * N_sat
    N_red_cen = f_red_cen_1 * N_cen
    N_blue_sat = (1-f_red_sat_1) * N_sat
    N_blue_cen = (1-f_red_cen_1) * N_cen
    f_sat_red_1 = N_red_sat/(N_red_cen+N_red_sat)
    f_sat_blue_1 = N_blue_sat/(N_blue_cen+N_blue_sat)
    ax=axes[4]
    p4b_1, = ax.plot(bin_centers,f_sat_red_1,'x',color='red', mfc='none', mec='red', ms=3)
    p4b_1, = ax.plot(bin_centers,f_sat_blue_1,'x',color='blue', mfc='none', mec='blue', ms=3)

    #calculate with correction
    b_sat=sat_bias(GC)
    Ps=sat_purity(GC)
    Pc=cen_purity(GC)
    Cs=sat_comp(GC)
    Cc=cen_comp(GC)
    print Ps, Cs, Pc, Cc, b_sat
    print 'test:',(1-Ps)*N_sat_mock*Cs*Pc==(1-Pc)*N_cen_mock*b_sat*Cc*Ps
    print 'test:',N_sat*Ps/Cs==N_sat_mock
    print 'test:',N_cen*Pc/Cc==N_cen_mock
    print 'test:',N_cen+N_sat==N_cen_mock+N_sat_mock
    f_red_sat_1 = (f_red_sat-f_red_cen*(((1-Ps)*(1-Cc))/((1-Pc)*Cc*b_sat)))/(1-(((1-Cs)*(1-Cc))/(Ps*Pc)))/Ps
    f_red_cen_1 = (f_red_cen-f_red_sat_1*(((1-Cs)*(1-Pc)*Ps*b_sat)/((1-Ps)*Cs)))/Pc
    ax=axes[1]
    p4a_1, = ax.plot(bin_centers,f_red_sat_1,'o',color='green', mfc='none', mec='green', ms=3)
    p4a_1, = ax.plot(bin_centers,f_red_cen_1,'o',color='orange', mfc='none', mec='orange', ms=3)
    f_sat_red_1 =f_sat_red*b_sat
    f_sat_blue_1 =f_sat_blue*b_sat
    ax=axes[4]
    N_red_sat = f_red_sat_1 * N_sat*Ps/Cs
    N_red_cen = f_red_cen_1 * N_cen*Pc/Cc
    N_blue_sat = (1-f_red_sat_1) * N_sat*Ps/Cs
    N_blue_cen = (1-f_red_cen_1) * N_cen*Pc/Cc
    f_sat_red_1 = N_red_sat/(N_red_cen+N_red_sat)
    f_sat_blue_1 = N_blue_sat/(N_blue_cen+N_blue_sat)
    p4b_1, = ax.plot(bin_centers,f_sat_red_1,'o',color='red', mfc='none', mec='red', ms=3)
    p4b_1, = ax.plot(bin_centers,f_sat_blue_1,'o',color='blue', mfc='none', mec='blue', ms=3)
    '''

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
    
        #galaxy color
        color = GC['M_g,0.1']-GC['M_r,0.1']
        LHS   = 0.7 - 0.032*(GC['M_r,0.1']+16.5) #Weinmann 2006
        blue_ind  = np.where(color<LHS)[0] #indices of blue galaxies
        red_ind   = np.where(color>LHS)[0] #indicies of red galaxies
        blue_bool = (color<LHS) #indices of blue galaxies
        red_bool  = (color>LHS) #indicies of red galaxies
    
        L = solar_lum(GC['M_r,0.1'],S_r)
    
        result = np.digitize(L,bins=bins)
        f_red_cen[boot,:]  = f_prop(L,bins,red_ind,blue_ind,centrals_bool)
        f_red_sat[boot,:]  = f_prop(L,bins,red_ind,blue_ind,satellites_bool)
        f_sat_red[boot,:]  = f_prop(L,bins,satellites_ind,centrals_ind,red_bool)
        f_sat_blue[boot,:] = f_prop(L,bins,satellites_ind,centrals_ind,blue_bool)
    
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
    ax=axes[5]
    p3b=ax.errorbar(bin_centers, f_sat_red, yerr=err_sat_red, fmt='o', color='red', mec='none', ms=3)
    p4b=ax.errorbar(bin_centers+0.01, f_sat_blue, yerr=err_sat_blue, fmt='o', color='blue', mec='none', ms=3)
   
    '''
    #calculate with correction
    catalogue_1 = group_catalogue
    f =  h5py.File(filepath_cat+catalogue_1+'.hdf5', 'r') #open catalogue file
    GC = f.get(catalogue_1)
    centrals_ind   = np.where(GC['RANK']==1)[0]
    satellites_ind = np.where(GC['RANK']!=1)[0]
    centrals_bool   = (GC['RANK']==1)
    satellites_bool = (GC['RANK']!=1)
    centrals_mock_ind   = np.where(GC['HALO_RANK']==1)[0]
    satellites_mock_ind = np.where(GC['HALO_RANK']!=1)[0]
    N_sat_mock = float(len(satellites_mock_ind))
    N_cen_mock = float(len(centrals_mock_ind))
    N_sat = float(len(satellites_ind))
    N_cen = float(len(centrals_ind))
    b_sat=sat_bias(GC)
    Ps=sat_purity(GC)
    Pc=cen_purity(GC)
    Cs=sat_comp(GC)
    Cc=cen_comp(GC)
    print Ps, Cs, Pc, Cc, b_sat
    print 'test:',(1-Ps)*N_sat_mock*Cs*Pc==(1-Pc)*N_cen_mock*b_sat*Cc*Ps
    print 'test:',N_sat*Ps/Cs==N_sat_mock
    print 'test:',N_cen*Pc/Cc==N_cen_mock
    print 'test:',N_cen+N_sat==N_cen_mock+N_sat_mock
    f_red_sat_1 = (f_red_sat-f_red_cen*(((1-Ps)*(1-Cc))/((1-Pc)*Cc*b_sat)))/(1-(((1-Cs)*(1-Cc))/(Ps*Pc)))/Ps
    f_red_cen_1 = (f_red_cen-f_red_sat_1*(((1-Cs)*(1-Pc)*Ps*b_sat)/((1-Ps)*Cs)))/Pc
    ax=axes[2]
    p4a_1, = ax.plot(bin_centers,f_red_sat_1,'o',color='green', mfc='none', mec='green', ms=3)
    p4a_1, = ax.plot(bin_centers,f_red_cen_1,'o',color='orange', mfc='none', mec='orange', ms=3)
    f_sat_red_1 =f_sat_red*b_sat
    f_sat_blue_1 =f_sat_blue*b_sat
    ax=axes[5]
    N_red_sat = f_red_sat_1 * N_sat*Ps/Cs
    N_red_cen = f_red_cen_1 * N_cen*Pc/Cc
    N_blue_sat = (1-f_red_sat_1) * N_sat*Ps/Cs
    N_blue_cen = (1-f_red_cen_1) * N_cen*Pc/Cc
    f_sat_red_1 = N_red_sat/(N_red_cen+N_red_sat)
    f_sat_blue_1 = N_blue_sat/(N_blue_cen+N_blue_sat)
    p4b_1, = ax.plot(bin_centers,f_sat_red_1,'o',color='red', mfc='none', mec='red', ms=3)
    p4b_1, = ax.plot(bin_centers,f_sat_blue_1,'o',color='blue', mfc='none', mec='blue', ms=3)
    '''

    plt.show(block=False)
    print plotpath+filename+'.eps'
    fig1.savefig(plotpath+filename+'.eps')
    

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
    
def solar_lum(M,Msol):
    L = ((Msol-M)/2.5)
    return L

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

    

if __name__ == '__main__':
    main() 
