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
    filepath_mock = cu.get_output_path() + 'processed_data/hearin_mocks/custom_catalogues/'
    plotpath = cu.get_plot_path()        + 'analysis/hearin_mocks/'

 
    print 'opening mock catalogue:', catalogue+'.hdf5'
    #open catalogue
    f1 = h5py.File(filepath_mock+catalogue+'.hdf5', 'r') #open catalogue file
    mock = f1.get(catalogue)
    mock = np.array(mock)
    print 'length:', len(mock)
    for name in mock.dtype.names: print '     ', name

    color = mock['g-r']
    LHS = 0.7 - 0.032*(mock['M_r,0.1']+16.5) #Weinmann 2006
    blue = np.where(color<LHS)[0] #indices of blue galaxies
    red = np.where(color>LHS)[0] #indicies of red galaxies

    #create subsample for plotting purposes
    N=10000.0
    N_red = float(len(red))
    N_blue = float(len(blue))
    N_red_1= int(N_red/(N_red+N_blue)*N)
    N_blue_1 = int(N_blue/(N_red+N_blue)*N)
    red_sub = np.random.random_integers(0,len(red)-1,N_red_1)
    blue_sub = np.random.random_integers(0,len(blue)-1,N_blue_1)

    x = np.arange(-25,-15,0.1)
    y = 0.7 - 0.032*(x+16.5) #Weinmann 2006

    #make density plot
    H, xedges, yedges = np.histogram2d(mock['M_r,0.1'],mock['g-r'],bins=50)
    xbins =  (xedges[:-1] + xedges[1:])/2.0
    ybins =  (yedges[:-1] + yedges[1:])/2.0

    fig = plt.figure(figsize=(3.3,3.3))
    ax = fig.add_subplot(1,1,1)
    fig.subplots_adjust(left=0.2, right=0.9, bottom=0.2,top=0.9)
    ax.scatter(mock['M_r,0.1'][red[red_sub]], mock['g-r'][red[red_sub]], s=1, marker='o', color='red', edgecolor='none')#, rasterized=True)
    ax.scatter(mock['M_r,0.1'][blue[blue_sub]], mock['g-r'][blue[blue_sub]], s=1, marker='o', color='blue', edgecolor='none')#, rasterized=True)
    #ax.contour(xbins,ybins,H.T, colors='0.75', linewidths=1)
    ax.plot(x,y, '-', color='black', linewidth=1)
    ax.set_xlabel(r'$M_r-5log(h)$')
    ax.set_ylabel(r'$(g-r)$')
    ax.set_xlim([-18.5,-22.5])
    ax.set_ylim([0,1.5])
    ax.set_yticks([0,0.5,1,1.5])
    ax.set_yticks(np.arange(0,1.5,0.1), minor=True)
    ax.set_xticks(np.arange(-22,-18,1))
    ax.set_xticks(np.arange(-22.5,-18.5,0.5), minor=True)
    plt.show(block=False)
    filename = catalogue+'_Mr_colors'
    fig.savefig(plotpath+filename+'.eps', dpi=400)

    fig = plt.figure(figsize=(3.3,3.3))
    ax = fig.add_subplot(1,1,1)
    fig.subplots_adjust(left=0.2, right=0.9, bottom=0.2,top=0.9)
    ax.plot(mock['M_host'][red], mock['V_peak'][red],'.', color='red', ms=0.2)
    ax.plot(mock['M_host'][blue], mock['V_peak'][blue],'.', color='blue', ms=0.2)
    ax.set_yscale('log')
    ax.set_ylim([90,1500])
    ax.set_xlim([10,15])
    ax.set_xlabel(r'$log(M_{\rm halo})$ $[M_{\odot}/h]$')
    ax.set_ylabel(r'$V_{\rm peak}$ $[{\rm km/s}]$')
    plt.show(block=False)
    filename = catalogue+'_Mhost_Vmax_color'
    #fig.savefig(plotpath+filename+'.eps')



    
if __name__ == '__main__':
    main() 
