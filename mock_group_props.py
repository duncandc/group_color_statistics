#!/usr/bin/env python

import h5py
import astropy
from astropy import table
from astropy.io import ascii
import numpy as np
import custom_utilities as cu
import matplotlib.pyplot as plt
import os.path
import sys

def main():

    mocks = ['Mr19_age_distribution_matching_mock',\
             'Mr19_age_distribution_matching_mock_cen_shuffle',\
             'Mr19_age_distribution_matching_mock_satsys_shuffle',\
             'Mr19_age_distribution_matching_mock_sys_empty_shuffle',\
             'Mr19_age_distribution_matching_mock_sys_empty_shuffle_satsys_shuffle',\
             'Mr19_age_distribution_matching_mock_sys_empty_shuffle_cen_shuffle',\
             'Mr19_age_distribution_matching_mock_sys_empty_shuffle_satrel_shuffle']
    short_names = ['Mr19_age_distribution_matching_mock',\
                   'cen_shuffle',\
                   'satsys_shuffle',\
                   'sys_empty_shuffle',\
                   'sys_empty_shuffle_satsys_shuffle',\
                   'sys_empty_shuffle_cen_shuffle',\
                   'sys_empty_shuffle_satrel_shuffle']
    
    dtype = [('name','a50'),('group_finder','a25'),('Ngal','>i8'),('Ngroups','>i8'),('fsat','>f2'),\
             ('z_min','>f2'),('z_max','>f2')]
    dtype = np.dtype(dtype)
    data = np.recarray((3*len(mocks),), dtype=dtype)
    data.fill('-')

    N_mocks = len(mocks)

    for i in range(0,len(mocks)):
        mock = mocks[i]
        print mock
        catalogue = mock
        filepath = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
        f = h5py.File(filepath+catalogue+'.hdf5', 'r') #open catalogue file
        GC = f.get(catalogue)

        catalogue = mock + '_clf_groups_M19'
        filepath = cu.get_output_path() + 'processed_data/tinker_groupcat/mock_runs/4th_run/custom_catalogues/'
        f  =  h5py.File(filepath+catalogue+'.hdf5', 'r') #open catalogue file
        tinker_GC = f.get(catalogue)

        catalogue = mock + '_groups'
        filepath = cu.get_output_path() + 'processed_data/berlind_groupcat/mock_runs/4th_run/custom_catalogues/'
        f  =  h5py.File(filepath+catalogue+'.hdf5', 'r') #open catalogue file
        berlind_GC = f.get(catalogue)

        catalogue = mock + '_groups'
        filepath = cu.get_output_path() + 'processed_data/yang_groupcat/mock_runs/4th_run/version_2/custom_catalogues/'
        if os.path.isfile(filepath+catalogue+'.hdf5'):
            f  =  h5py.File(filepath+catalogue+'.hdf5', 'r') #open catalogue file
            yang_GC = f.get(catalogue)

        name = short_names[i]
        group_finder = 'tinker'
        Ngal = len(tinker_GC)
        Ngroups = len(np.where(tinker_GC['RPROJ']==0)[0])
        Nsat = len(np.where(tinker_GC['RPROJ']!=0)[0])
        fsat = float(Nsat)/(float(Ngal))
        zmin = min(tinker_GC['Z'])
        zmax = max(tinker_GC['Z'])
        data['name'][i] = name
        data['group_finder'][i] = group_finder
        data['Ngal'][i] = Ngal
        data['Ngroups'][i] = Ngroups
        data['fsat'][i] = fsat
        data['z_min'][i] = zmin
        data['z_max'][i] = zmax

        name = short_names[i]
        group_finder = 'berlind'
        Ngal = len(berlind_GC)
        Ngroups = len(np.where(berlind_GC['RPROJ']==0)[0])
        Nsat = len(np.where(berlind_GC['RPROJ']!=0)[0])
        fsat = float(Nsat)/(float(Ngal))
        zmin = min(berlind_GC['Z'])
        zmax = max(berlind_GC['Z'])
        data['name'][N_mocks+i] = name
        data['group_finder'][N_mocks+i] = group_finder
        data['Ngal'][N_mocks+i] = Ngal
        data['Ngroups'][N_mocks+i] = Ngroups
        data['fsat'][N_mocks+i] = fsat
        data['z_min'][N_mocks+i] = zmin
        data['z_max'][N_mocks+i] = zmax


        if os.path.isfile(filepath+catalogue+'.hdf5'):
            name = short_names[i]
            group_finder = 'yang'
            Ngal = len(yang_GC)
            Ngroups = len(np.where(yang_GC['RPROJ']==0)[0])
            Nsat = len(np.where(yang_GC['RPROJ']!=0)[0])
            fsat = float(Nsat)/(float(Ngal))
            zmin = min(yang_GC['Z'])
            zmax = max(yang_GC['Z'])
            data['name'][2*N_mocks+i] = name
            data['group_finder'][2*N_mocks+i] = group_finder
            data['Ngal'][2*N_mocks+i] = Ngal
            data['Ngroups'][2*N_mocks+i] = Ngroups
            data['fsat'][2*N_mocks+i] = fsat
            data['z_min'][2*N_mocks+i] = zmin
            data['z_max'][2*N_mocks+i] = zmax

    print 'saving ascii version of the catalogue...'
    filename = 'mock_group_properties'
    savepath = './'
    data_table = table.table.Table(data=data)
    ascii.write(data_table, savepath+filename+'.dat')
    print data_table

    catalogue = sys.argv[1]

    ##########################
    #open mock
    ########################## 
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    print 'opening mock catalogue:', catalogue+'.hdf5'
    #open catalogue
    f1 = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f1.get(catalogue)
    mock = np.array(mock)
    #print 'length:', len(mock)
    #for name in mock.dtype.names: print '     ', name

    Ngal = len(mock)
    Ngroups = len(np.where(mock['ID_host']==-1)[0])
    Nsat = len(np.where(mock['ID_host']!=-1)[0])
    fsat = float(Nsat)/(float(Ngal))
    print fsat


if __name__ == '__main__':
    main()
