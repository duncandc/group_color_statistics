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

    centrals   = np.array(mock['ID_host']==-1) #central galaxies
    satellites = np.array(mock['ID_host']!=-1) #satellite galaxies
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

    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(L[satellites],mock['g-r'][satellites],mock['M_host'][satellites],marker='.',alpha=0.1)
    #ax.set_ylim([0,1.5])
    #plt.show()

    fig1, axes = plt.subplots(nrows=3, ncols=3, sharex=True, sharey=True, figsize=(10,10))
    fig1.subplots_adjust(hspace=0, wspace=0.0)
    fig1.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.9)
    axes = axes.flatten()

    #satellites
    ax = axes[0]
    selection = np.array(luminds==0)
    selection = selection & satellites
    ax.plot(mock['M_host'][selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='green')
    ax.text(11.5,1.1,lumbin_centers[0])
    ax = axes[1]
    selection = np.array(luminds==1)
    selection = selection & satellites
    ax.plot(mock['M_host'][selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='green')
    ax.text(11.5,1.1,lumbin_centers[1])
    ax = axes[2]
    selection = np.array(luminds==2)
    selection = selection & satellites
    ax.plot(mock['M_host'][selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='green')
    ax.text(11.5,1.1,lumbin_centers[2])
    ax = axes[3]
    selection = np.array(luminds==3)
    selection = selection & satellites
    ax.plot(mock['M_host'][selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='green')
    ax.text(11.5,1.1,lumbin_centers[3])
    ax = axes[4]
    selection = np.array(luminds==4)
    selection = selection & satellites
    ax.plot(mock['M_host'][selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='green')
    ax.text(11.5,1.1,lumbin_centers[4])
    ax = axes[5]
    selection = np.array(luminds==5)
    selection = selection & satellites
    ax.plot(mock['M_host'][selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='green')
    ax.text(11.5,1.1,lumbin_centers[5])
    ax = axes[6]
    selection = np.array(luminds==6)
    selection = selection & satellites
    ax.plot(mock['M_host'][selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='green')
    ax.text(11.5,1.1,lumbin_centers[6])
    ax = axes[7]
    selection = np.array(luminds==7)
    selection = selection & satellites
    ax.plot(mock['M_host'][selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='green')
    ax.text(11.5,1.1,lumbin_centers[7])
    ax = axes[8]
    selection = np.array(luminds==8)
    selection = selection & satellites
    ax.plot(mock['M_host'][selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='green')
    ax.text(11.5,1.1,lumbin_centers[8])
    
    #centrals
    ax = axes[0]
    selection = np.array(luminds==0)
    selection = selection & centrals
    ax.plot(mock['M_host'][selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    ax = axes[1]
    selection = np.array(luminds==1)
    selection = selection & centrals
    ax.plot(mock['M_host'][selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    ax = axes[2]
    selection = np.array(luminds==2)
    selection = selection & centrals
    ax.plot(mock['M_host'][selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    ax = axes[3]
    selection = np.array(luminds==3)
    selection = selection & centrals
    ax.plot(mock['M_host'][selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    ax = axes[4]
    selection = np.array(luminds==4)
    selection = selection & centrals
    ax.plot(mock['M_host'][selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    ax = axes[5]
    selection = np.array(luminds==5)
    selection = selection & centrals
    ax.plot(mock['M_host'][selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    ax = axes[6]
    selection = np.array(luminds==6)
    selection = selection & centrals
    ax.plot(mock['M_host'][selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    ax = axes[7]
    selection = np.array(luminds==7)
    selection = selection & centrals
    ax.plot(mock['M_host'][selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    ax = axes[8]
    selection = np.array(luminds==8)
    selection = selection & centrals
    ax.plot(mock['M_host'][selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')

    ax.set_xlim([11,15])
    ax.set_ylim([0,1.25])
    plt.show()

    fig2, axes = plt.subplots(nrows=3, ncols=3, sharex=True, sharey=True, figsize=(10,10))
    fig2.subplots_adjust(hspace=0, wspace=0.0)
    fig2.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.9)
    axes = axes.flatten()

    #satellites
    ax = axes[0]
    selection = np.array(massinds==5)
    selection = selection & satellites
    ax.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='green')
    ax.text(9.5,1.1,massbin_centers[5])
    ax = axes[1]
    selection = np.array(massinds==6)
    selection = selection & satellites
    ax.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='green')
    ax.text(9.5,1.1,massbin_centers[6])
    ax = axes[2]
    selection = np.array(massinds==7)
    selection = selection & satellites
    ax.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='green')
    ax.text(9.5,1.1,massbin_centers[7])
    ax = axes[3]
    selection = np.array(massinds==8)
    selection = selection & satellites
    ax.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='green')
    ax.text(9.5,1.1,massbin_centers[8])
    ax = axes[4]
    selection = np.array(massinds==9)
    selection = selection & satellites
    ax.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='green')
    ax.text(9.5,1.1,massbin_centers[9])
    ax = axes[5]
    selection = np.array(massinds==10)
    selection = selection & satellites
    ax.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='green')
    ax.text(9.5,1.1,massbin_centers[10])
    ax = axes[6]
    selection = np.array(massinds==11)
    selection = selection & satellites
    ax.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='green')
    ax.text(9.5,1.1,massbin_centers[11])
    ax = axes[7]
    selection = np.array(massinds==12)
    selection = selection & satellites
    ax.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='green')
    ax.text(9.5,1.1,massbin_centers[12])
    ax = axes[8]
    selection = np.array(massinds==13)
    selection = selection & satellites
    ax.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='green')
    ax.text(9.5,1.1,massbin_centers[13])
    
    #centrals
    ax = axes[0]
    selection = np.array(massinds==5)
    selection = selection & centrals
    ax.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    ax = axes[1]
    selection = np.array(massinds==6)
    selection = selection & centrals
    ax.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    ax = axes[2]
    selection = np.array(massinds==7)
    selection = selection & centrals
    ax.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    ax = axes[3]
    selection = np.array(massinds==8)
    selection = selection & centrals
    ax.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    ax = axes[4]
    selection = np.array(massinds==9)
    selection = selection & centrals
    ax.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    ax = axes[5]
    selection = np.array(massinds==10)
    selection = selection & centrals
    ax.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    ax = axes[6]
    selection = np.array(massinds==11)
    selection = selection & centrals
    ax.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    ax = axes[7]
    selection = np.array(massinds==12)
    selection = selection & centrals
    ax.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    ax = axes[8]
    selection = np.array(massinds==13)
    selection = selection & centrals
    ax.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')

    ax.set_xlim([9.5,11])
    ax.set_ylim([0,1.25])
    plt.show()

    plt.figure()
    selection = np.array(massinds<=8) & np.array(massinds>6)
    selection = selection & centrals
    plt.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    plt.xlim([9.5,11])
    plt.ylim([0,1.25])
    plt.figure()
    selection = np.array(massinds>8)
    selection = selection & centrals
    plt.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    plt.xlim([9.5,11])
    plt.ylim([0,1.25])

    plt.figure()
    selection = np.array(massinds<5)
    selection = selection & centrals
    plt.plot(L[selection],mock['g-r'][selection],'.',ms=1, alpha=0.2, color='orange')
    plt.xlim([9.5,11])
    plt.ylim([0,1.25])
    plt.show()
    
    


def solar_lum(M,Msol):
    L = ((Msol-M)/2.5)
    return L    

if __name__ == '__main__':
    main() 
