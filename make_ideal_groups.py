#!/usr/bin/env python

#Duncan Campbell
#October 3, 2014
#Yale University
#use the halo memberships from the mock to create an ideal group catalogue. Group masses 
#are estimated using abundance matching. Central/satellite designation is 
#brightest/non-brightest.

#load packages
import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
import custom_utilities as cu
from astropy import cosmology
from astropy.io import ascii
from scipy import interpolate

def main():
    
    if len(sys.argv)>1: catalogue = sys.argv[1]
    else: catalogue = 'Mr19_age_distribution_matching_mock'
    
    savepath = cu.get_output_path() + 'processed_data/hearin_mocks/ideal_groups/'

    #open mock catalogue
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    print 'opening mock catalogue:', catalogue+'.hdf5'
    f1 = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f1.get(catalogue)
    mock = np.array(mock)

    #create a new catalogue to store results
    dtype=[('ID','>i8'),('k_1','>i8'),('k_2','>i8'),\
           ('RA','>f8'),('DEC','>f8'),('Z','>f8'),('red','>i8'),\
           ('M_u,0.1','>f8'),('M_g,0.1','>f8'),('M_r,0.1','>f8'),('M_i,0.1','>f8'),('M_z,0.1','>f8'),\
           ('MSTAR','>f8'),\
           ('GROUP_ID','>i8'),('MGROUP','>f8'),('ZGROUP','>f8'),('R200','>f8'),\
           ('CEN_IND','>i8'),('RANK','>i8'),('RPROJ','>f8'),\
           ('N_sat','>i8'),('N_sat_red','>i8'),('N_sat_blue','>i8'),\
           ('HALO_M','>f8'),('HALO_RANK','>i8')]
    dtype = np.dtype(dtype)
    data = np.recarray((len(mock),), dtype=dtype)
    data.fill(-99.9) #empty value indicator

    centrals   = np.where(mock['ID_host']==-1)[0] #central galaxies
    satellites = np.where(mock['ID_host']!=-1)[0] #satellite galaxies
    result     = np.in1d(mock['ID_host'][satellites],mock['ID_halo'][centrals])
    special    = np.where(result==False)[0] #satellites without a host
    special    = satellites[special]
    satellites = satellites[result] #satellites with a host

    index    = np.argsort(mock['ID_halo'][centrals]) #indices which sort halo IDs of centrals
    sorted_x = mock['ID_halo'][centrals[index]] #sorted halo IDs of centrals
    ind      = np.searchsorted(sorted_x,mock['ID_host'][satellites])
    ind      = index[ind] #index of centrals for each satellite
 
    data['ID'] = mock['ID_halo']
 
    data['GROUP_ID']      = mock['ID_host'] #satellites
    #data['GROUP_ID'][ind] = mock['ID_halo'][ind] #centrals with satellites
    
    #result = (data['GROUP_ID']==-1)
    result = (mock['ID_host']==-1)
    data['GROUP_ID'][result] = mock['ID_halo'][result] #centrals with no satellites

    data['M_g,0.1'] = mock['M_r,0.1']+mock['g-r']
    data['M_r,0.1'] = mock['M_r,0.1']
    data['HALO_M']  = mock['M_host']

    #determine cen/sat designation in xyz mock
    result = np.where(mock['ID_host']==-1)[0]
    data['HALO_RANK'][result] = 0 #central
    result = np.where(mock['ID_host']!=-1)[0]
    data['HALO_RANK'][result] = 1 #satellite

    #calculate galaxy colors
    color = data['M_g,0.1']-data['M_r,0.1']
    LHS = 0.7 - 0.032*(data['M_r,0.1']+16.5) #Weinmann 2006
    blue = np.where(color<LHS)[0] #indices of blue galaxies
    red = np.where(color>LHS)[0]  #indicies of red galaxies
    
    #record color designation
    data['red'][red]  = 1
    data['red'][blue] = 0

    for group_id in np.unique(data['GROUP_ID']):
        members  = np.where(data['GROUP_ID']==group_id)[0]
        central = np.where(data['HALO_RANK'][members]==0)[0]
        central  = members[central]
        satellites = np.where(data['HALO_RANK'][members]!=0)[0]
        satellites = members[satellites]
        rank = np.argsort(data['M_r,0.1'][members])
        data['RANK'][members[rank]]    = np.arange(0,len(rank))
        #record number of satellites in the group
        data['N_sat'][members]   = len(satellites)
        sat_red  = np.where(np.in1d(satellites,red)==True)[0]
        data['N_sat_red'][members]  = len(sat_red)
        sat_blue = np.where(np.in1d(satellites,blue)==True)[0]
        data['N_sat_blue'][members] = len(sat_blue)
        #record other group information
        data['CEN_IND'][members] = central

    centrals = np.where(data['RANK']==0)[0]

    #read in mass halo function
    filepath = '/scratch/dac29/fortran_code/mass_functions/'
    filename = 'Bolshoi_Massfunc.dat'
    names = ['dM','dn','nsum']
    dndM = ascii.read(filepath+filename, delimiter='\s', names=names, data_start=0)
    dndM = np.array(dndM)

    #calculate group total r-band luminosities
    S_r = 4.64
    group_L = np.zeros((len(data),),dtype=np.float)

    for group_id in np.unique(data['GROUP_ID']):
        members  = np.where(data['GROUP_ID']==group_id)[0]
        group_L[members] = np.log10(np.sum(10.0**(solar_lum(data['M_r,0.1'][members], S_r))))
    tot_lum = group_L[centrals]

    #calculate abundance matched masses for groups
    mock_volume = 250.0**3.0

    #caclulate the group luminosity function
    N_gal = np.cumsum(np.zeros(len(centrals))+1) #cumulative number of groups
    n_gal = N_gal/mock_volume #number density
    L_gal = np.sort(tot_lum)[::-1] #group luminosity
    ind   = np.argsort(tot_lum)[::-1]

    #integrate halo mass function
    n_halo  = dndM['nsum'][::-1] #cumulative number density
    M_halo  = dndM['dM'][::-1] #halo mass

    #interpolate the halo mass function
    #x = np.log10(n_halo)
    x = np.log10(n_halo)
    y = M_halo
    f = interpolate.interp1d(x, y, kind='linear', bounds_error='False', fill_value=0.0)

    data['MGROUP'][centrals[ind]] = f(np.log10(n_gal))

    for group_id in np.unique(data['GROUP_ID']):
        members  = np.where(data['GROUP_ID']==group_id)[0]
        central = np.where(data['RANK'][members]==0)[0]
        central  = members[central]
        data['MGROUP'][members] = data['MGROUP'][central]

    print 'saving hdf5 version of the catalogue...'
    filename = catalogue+'_groups'
    f = h5py.File(savepath+filename+'.hdf5', 'w')
    dset = f.create_dataset(filename, data=data)
    f.close()

def solar_lum(M,Msol):
    L = ((Msol-M)/2.5)
    return L 


if __name__ == '__main__':
    main() 
