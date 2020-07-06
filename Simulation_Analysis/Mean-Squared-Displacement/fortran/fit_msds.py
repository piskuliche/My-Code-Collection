#!/usr/bin/env python
'''
Python program for reading in a msd and block data and outputting the block undertainty of the trajectory.
The program also fits the mean squared displacements and prints the diffusion coefficient to a file.
'''

# Define the Linear Fitting Function
def linear(x, m, b):
    return m*x + b

def D_from_MSD(time, msd, blmsd, startskip, endskip, nblocks,molname,label,t_val):
    from scipy.optimize import curve_fit
    import numpy as np
    """ 
    Function that fits an msd vector and an nblock msd vector over a range startskip:endskip, and reports 
    the Diffusion coefficient and the uncertainty in the diffusion coefficient according to a confidence interval
    defined by the input t_val
    """
    # Fit the msd
    popt, pcov = curve_fit(linear, time[startskip:endskip], msd[startskip:endskip])
    
    # output fitted msd to file
    fittedout = open('fitmsd_'+molname+'.log','w')
    for i in time:
        fittedout.write("%s %s\n" % (i, linear(i, popt[0], popt[1])))
    fittedout.close()
    
    # Fit block msds
    popt_bl = []
    for i in range(nblocks):
        tmp_popt, tmp_pcov = curve_fit(linear, time[startskip:endskip], blmsd[i][startskip:endskip])
        popt_bl.append(tmp_popt[0])
    # Calculate diffusion coeffs
    d_tot = popt[0]/6*10**(-4)
    d_std = np.std(popt_bl)*t_val/(6*np.sqrt(nblocks))*10**(-4)
    f = open('bl_diff_'+str(molname)+'.dat','w')
    for i in range(nblocks):
        dbl = popt_bl[i]/6*10**(-4)
        print("block %d diffusion: %2.4e" %(i,dbl))
        f.write('%d %s' % (i,dbl))
    f.close()

    print("%s %s" % (d_tot, d_std))

    diffout = open('cm_'+molname+'_diffusion.log','w')
    diffout.write('%s %2.4e %2.4e\n' % (label, d_tot, d_std))
    diffout.close()

    return d_tot, d_std




if __name__ == "__main__":
    import numpy as np
    import sys
    import math
    import scipy.stats as stats

    # Read in parameters
    if len(sys.argv) != 8:
        print("Usage: fit_msds.py file nblocks molname startskip endskip prepend startindex") 
        exit(0)
    filename = str(sys.argv[1])
    nblocks  = int(sys.argv[2])
    molname  = str(sys.argv[3])
    startskip = int(sys.argv[4])
    endskip   = int(sys.argv[5])
    prepend   = str(sys.argv[6])
    startindex = int(sys.argv[7])
    t_val=stats.t.ppf(0.975,nblocks-1)

    # Define the Linear Fitting Function
    def linear(x, m, b):
        return m*x + b

    # Read in the data
    time, msd = np.genfromtxt(filename, usecols = (0,startindex+1), unpack=True)
    blockmsd = []
    for i in range(nblocks):
        j = 2+startindex + i
        blockmsd.append(np.genfromtxt(filename, usecols = j, unpack=True))

    # Calculate block uncertainty, output to a file.
    bl_std = np.std(blockmsd,axis=0)
    bl_err = t_val*bl_std
    np.savetxt(molname+'_'+'msd.dat', np.c_[time, msd, bl_err])

    # Fit the correlation functions
    endskip = len(msd) - endskip
    if endskip < 0:
        print("Error: Endskip is too big")
        exit(1)
    label = molname
    
    d_tot, d_std = D_from_MSD(time, msd, blockmsd, startskip, endskip, nblocks,molname, t_val)

    output = open('msd_'+molname+'.log', 'a')
    output.write("%s %s %s\n" % (prepend, d_tot, d_std))
    print("Diffusion coefficient and error are %s %s" % (d_tot, d_std))
    output.close()
