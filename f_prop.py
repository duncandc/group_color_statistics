#Duncan Campbell
#April 7, 2014
#Calculate the fraction of galaxies with a specified prop as a function of a specified variable.

def main():
    print 'main'

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

if __name__ == '__main__':
    main()
