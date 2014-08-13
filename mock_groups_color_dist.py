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
    blue  = np.array(color<LHS) #indices of blue galaxies
    red   = np.array(color>LHS) #indicies of red galaxies

    massbins = np.arange(11,15,0.2)
    massbin_centers = (massbins[:-1]+massbins[1:])/2.0
    mass_inds = np.digitize(mock['M_host'],bins=massbins)

    lumbins = np.arange(9.6,11.8,0.2)
    lumbin_centers = (lumbins[:-1]+lumbins[1:])/2.0
    lum_inds = np.digitize(solar_lum(mock['M_r,0.1'],4.64),bins=lumbins)

    colorbins = np.arange(0,1.5,0.05)
    colorbin_centers = (colorbins[:-1]+colorbins[1:])/2.0

    lum_selection = np.array(lum_inds==5)
    selection = lum_selection & centrals
    result_cen = np.histogram(mock['M_host'][selection],bins=massbins)[0]
    result_cen = result_cen/np.float(np.sum(result_cen))
    selection = lum_selection & satellites
    result_sat = np.histogram(mock['M_host'][selection],bins=massbins)[0]
    result_sat = result_sat/np.float(np.sum(result_sat))

    lower_cen = max(np.where(np.cumsum(result_cen)<0.25)[0])
    upper_cen = min(np.where(np.cumsum(result_cen)>.75)[0])
    print lower_cen, upper_cen
    
    '''
    plt.figure()
    plt.step(massbins[:-1],result_cen,where='pre')
    plt.step(massbins[:-1],result_sat,where='pre')
    plt.show()
    '''

    plt.figure()
    for i in range(lower_cen,upper_cen):
        mass_selection = np.array(mass_inds==i)
        selection = mass_selection & lum_selection
        selection = selection & centrals
        colors = np.sort(mock['g-r'][selection])
        N = np.cumsum(np.zeros(len(colors))+1.0)
        N = N/np.float(len(N))
        plt.plot(colors,N,color='black')

    c = [(0.1,1,0,0.1),(0.1,1,0,0.2),(0.1,1,0,0.3),(0.1,1,0,0.4),(0.1,1,0,0.5),\
         (0.1,1,0,0.6),(0.1,1,0,0.7),(0.1,1,0,0.8),(0.1,1,0,0.9),(0.1,1,0,1),(0.1,1,0,1),(0.1,1,0,1)]
    print colors
    for i in range(upper_cen,len(massbins)-1):
        print i
        mass_selection = np.array(mass_inds==i)
        selection = mass_selection & lum_selection
        selection = selection & satellites
        colors = np.sort(mock['g-r'][selection])
        N = np.cumsum(np.zeros(len(colors))+1.0)
        N = N/np.float(len(N))
        print c[i-upper_cen]
        plt.plot(colors,N,color=c[i-upper_cen])

    plt.xlim([0,1.5])
    plt.ylim([0,1])
    plt.show()


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

if __name__ == '__main__':
    main()
