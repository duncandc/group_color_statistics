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
    
    if len(sys.argv)>1: catalogue = sys.argv[1]
    else: catalogue = 'Mr19_age_distribution_matching_mock'

    #open mock catalogue
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    print 'opening mock catalogue:', catalogue+'.hdf5'
    f1 = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
    mock1 = f1.get(catalogue)
    
    #open mock catalogue
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/ideal_mocks/
    catalogue = catalogue+'_groups'
    print 'opening mock catalogue:', catalogue+'.hdf5'
    f2 = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
    mock2 = f2.get(catalogue)
    
    print mock1.dtype.names
    print mock2.dtype.names

    ind = np.where(mock2['MGROUP']==-99.9)[0]
    print len(ind)

    #sys.exit()

    plt.figure()
    plt.plot(mock2['MGROUP'], mock2['HALO_M'],'.', ms=1, alpha=0.25)
    plt.xlim([10,16])
    plt.show()

    print min(mock2['RANK']), max(mock2['RANK'])
    print min(mock2['HALO_RANK']), max(mock2['HALO_RANK'])

    b_sats = np.where(mock2['HALO_RANK']==1)[0]

    print len(np.where(mock2['HALO_RANK']==-99)[0])
    print len(np.where(mock2['HALO_RANK']==0)[0])
    print len(np.where(mock2['HALO_RANK']==1)[0])

    dM = mock2['M_r,0.1'][b_sats]-mock2['M_r,0.1'][mock2['CEN_IND'][b_sats]]
    print dM
    
    mag_bins = np.arange(-5,2,0.1)
    mag_bin_centers = (mag_bins[:-1]+mag_bins[1:])/2.0
    print mag_bin_centers
    result = np.histogram(dM, bins=mag_bins)[0]
    plt.figure()
    plt.plot(mag_bin_centers,result)
    plt.show()

if __name__ == '__main__':
    main() 
