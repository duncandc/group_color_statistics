import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
import custom_utilities as cu
from astropy import cosmology
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def main():
    catalogue = 'Mr19_age_distribution_matching_mock'
    if len(sys.argv)>1: catalogue = sys.argv[1]

    #open mock catalogue
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    print 'opening mock catalogue:', catalogue+'.hdf5'
    f1 = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f1.get(catalogue)
    
    print mock.dtype.names

    centrals   = np.where(mock['ID_host']==-1)[0] #central galaxies
    satellites = np.where(mock['ID_host']!=-1)[0] #satellite galaxies
    LHS   = 0.7 - 0.032*(mock['M_r,0.1']+16.5) #Weinmann 2006 color cut
    blue  = np.array(mock['g-r']<LHS)[0] #blue galaxies
    red   = np.array(mock['g-r']>LHS)[0] #red galaxies
    L = solar_lum(mock['M_r,0.1'],4.64)

    lumbins = np.arange(9.6,11.0,0.15)
    lumbin_centers = (lumbins[:-1]+lumbins[1:])/2.0
    luminds = np.digitize(solar_lum(mock['M_r,0.1'],4.64),bins=lumbins)

    massbins = np.arange(11,15,0.2)
    massbin_centers = (massbins[:-1]+massbins[1:])/2.0
    massinds = np.digitize(mock['M_host'],bins=massbins)

    colorbins = np.arange(0,1.5,0.1)
    colorbin_centers = (colorbins[:-1]+colorbins[1:])/2.0
    colorinds = np.digitize(mock['g-r'],bins=colorbins)

    print np.all(np.in1d(mock['ID_host'][satellites],mock['ID_halo'][centrals]))
    result = np.in1d(mock['ID_host'][satellites],mock['ID_halo'][centrals])
    print len(np.where(result==True)[0])
    print len(np.where(result==False)[0])
    print len(np.where(result==False)[0]) + len(np.where(result==True)[0])

    sys.exit()

    f = np.zeros(len(massbin_centers))
    for i in range(0,len(massbin_centers)):
        print i, massbin_centers[i]
        selection = np.where(massinds==i)[0]
        central_selection = np.in1d(centrals,selection)
        central_selection = centrals[central_selection]
        satellite_selection = np.in1d(satellites,selection)
        satellite_selection = satellites[satellite_selection]
        print 'test:', np.all(np.in1d(mock['ID_host'][satellite_selection],mock['ID_halo'][central_selection]))
        print np.in1d(mock['ID_host'][satellite_selection],mock['ID_halo'][central_selection])
        if len(satellite_selection>0):
            index    = np.argsort(mock['ID_halo'][central_selection])
            sorted_x = mock['ID_halo'][central_selection[index]]
            ind      = np.searchsorted(sorted_x,mock['ID_host'][satellite_selection])
            print 'test:', len(ind)==len(satellite_selection)
            ind      = index[ind]
            dM =  mock['M_r,0.1'][central_selection[ind]] - mock['M_r,0.1'][satellite_selection]
            print 'test:', np.all(mock['ID_halo'][central_selection[ind]]==mock['ID_host'][satellite_selection])
    
            mag_bins = np.arange(-5,2,0.1)
            mag_bin_centers = (mag_bins[:-1]+mag_bins[1:])/2.0
            result = np.histogram(dM, bins=mag_bins)[0]
            #plt.figure()
            #plt.plot(mag_bin_centers,result)
            #plt.show()

            keep = np.where(mag_bin_centers>=0)[0]
            print np.sum(result[keep]), np.sum(result), np.sum(result[keep])/np.float(np.sum(result))
            f[i] = np.sum(result[keep])/np.float(np.sum(result))
    plt.figure()
    plt.plot(massbin_centers,f)
    plt.show()

def solar_lum(M,Msol):
    L = ((Msol-M)/2.5)
    return L    

if __name__ == '__main__':
    main() 
