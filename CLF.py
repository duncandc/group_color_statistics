#!/usr/bin/env python

#Duncan Campbell
#April 7, 2014
#Yale University
#Calculate the conditonal luminosity function

def main():
    import numpy as np
    import h5py
    import matplotlib.pyplot as plt
    import custom_utilities as cu
    import sys

    catalogue = sys.argv[1]

    bins = np.arange(7.5,11.5,0.1)
    bin_centers = (bins[:-1]+bins[1:])/2.0
    S_r = 4.64

    print 'opening mock catalogue:', catalogue+'.hdf5'
    #open catalogue
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    f1 = h5py.File(filepath_mock+catalogue+'_extended.hdf5', 'r') #open catalogue file
    mock = f1.get(catalogue+'_extended')
    mock = np.array(mock)
    print 'length:', len(mock)
    for name in mock.dtype.names: print '     ', name
    
    occupied = np.where(mock['M_r,0.1']!=-99)
    empty = np.where(mock['M_r,0.1']==-99)
    host =  np.where(mock['ID_host']==-1)[0]
    L = solar_lum(mock['M_r,0.1'],S_r)
    color = mock['g-r']
    LHS = 0.7 - 0.032*(mock['M_r,0.1']+16.5) #Weinmann 2006
    blue  = np.array((color<LHS) & (mock['M_r,0.1']!=-99))
    red   = np.array((color>LHS) & (mock['M_r,0.1']!=-99))

    host_ID = mock['ID_host']
    host_ID[host]=mock['ID_halo'][host]
    
    mass_bin = [12.9,13.2]
    phi = clf(host_ID,L,bins,mock['M_host'],mass_bin)
    mask = red
    phi_red = clf(host_ID,L,bins,mock['M_host'],mass_bin,mask)
    mask = blue
    phi_blue = clf(host_ID,L,bins,mock['M_host'],mass_bin,mask)

    plt.figure()
    plt.plot(bin_centers,phi, color='black')
    plt.plot(bin_centers,phi_red, color='red')
    plt.plot(bin_centers,phi_blue, color='blue')
    plt.yscale('log')
    plt.ylim([0.1,100])
    plt.ylabel(r'$\phi(L)dL/{\rm group}$')
    plt.xlabel('L')
    plt.show()

def clf(halo_ID, L, bins, mass, mass_bin, mask='none'):
    import numpy as np
    
    selection = np.array((mass>mass_bin[0]) & (mass<mass_bin[1]))
    IDs, inds = np.unique(halo_ID[selection], return_index=True)
    N_halo = float(len(IDs))
    print N_halo
    
    
    if mask!='none':
        L = L[np.array(selection&mask)].copy()
    else: L = L[selection].copy()

    dL = np.fabs(bins[1]-bins[0])
    avg_N  = np.zeros((len(bins)-1)) #store result

    result = np.digitize(L,bins)
    for i in range(0,len(bins)-1):
        inds = np.where(result==i+1)[0]
        N_gal = float(len(inds))
        avg_N[i] = N_gal/N_halo

    phi = avg_N/dL
    return phi

def solar_lum(M,Msol):
    L = ((Msol-M)/2.5)
    return L    

if __name__ == '__main__':
    main()
