#!/usr/bin/env python
"""
This is a code that is meant to plot, and manipulate the gofr's output by the fortran code.
"""

def read_input(fname, selec1, selec2):
    with open(fname, 'r') as f:
        lines=f.readlines()
        for line  in lines:
            if "L" in line:
                L=float(line[2:-2])
            if "selec1" in line:
                tmp=int(line[7:-2])
            if "selec2" in line:
                tmp2=int(line[7:-2])
            if "nblocks" in line:
                nblocks=int(line[8:-2])
    # Doesn't override command line
    if int(selec1) == 0:
        selec1 = str(tmp)
    if int(selec2) == 0:
        selec2 = str(tmp2)

    if L == 0.0:
        npt = 1
        L = np.genfromtxt("L.dat", usecols=0)
    else:
        npt=0
    return L, selec1, selec2, npt


def read_gofr(t_val,selec1, selec2):
    """
    This function reads in the gofr data (two column distance gofr formatted)
    from a file (i.e. pairdist_1_1.dat)
    and then the block gofrs from the block flies (two column dist gofr formt)
    from the files (i.e. bl_1_pairdist_1_1.dat)
    It then calculates the error in the RDF from the blocks, according to a 95% confidence
    interval associated with Student's t-distribution
    Writes to a file (i.e RDF_1_1_298.dat)
    """
    print("Reading in the RDF")
    # Reads in RDF
    fstring="pairdist_"+selec1+"_"+selec2+".dat"
    dist, gofr = np.genfromtxt(fstring, usecols=(0,1), unpack=True)

    # Reads in Blcok RDF
    bl_gofr = []
    for i in range(1,nblocks+1):
        bstring=("bl_%d_"%i)
        bl_gofr.append(np.genfromtxt(bstring+fstring, usecols=1, unpack=True))
    # Calculates the confidence interval (note t_val includes a factor of 1/sqrt(nblocks) implicitly)
    err_gofr = np.std(bl_gofr,axis=0)*t_val
    np.savetxt('RDF_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dist,gofr,err_gofr])
    return dist, gofr, bl_gofr, err_gofr

def read_egofr(t_val,selec1, selec2,item):
    """
    This function reads in the gofr data (two column distance gofr formatted)
    from a file (i.e. epairdist_1_1.dat)
    and then the block gofrs from the block flies (two column dist gofr formt)
    from the files (i.e. bl_1_epairdist_1_1.dat)
    It then calculates the error in the RDF from the blocks, according to a 95% confidence
    interval associated with Student's t-distribution
    Writes to a file (i.e eRDF_1_1_298.dat)
    """
    print("Reading in the RDF weighted by e")
    # Read in eRDF
    fstring=str(item)+"pairdist_"+selec1+"_"+selec2+".dat"
    dist, egofr,e2gofr = np.genfromtxt(fstring, usecols=(0,1,2), unpack=True)

    # Reads in block eRDF
    bl_egofr = []
    for i in range(1,nblocks+1):
        bstring=("bl_%d_"%i)
        bl_egofr.append(np.genfromtxt(bstring+fstring, usecols=1, unpack=True))

    err_egofr = np.std(bl_egofr,axis=0)*t_val
    np.savetxt(str(item)+'RDF_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dist,egofr,err_egofr])
    return egofr, bl_egofr, err_egofr,e2gofr

def read_volgofr(t_val,selec1, selec2):
    """
    This function reads the gofr data (two column distance gofr formatted)
    from a file (i.e. volpairdist_1_1.dat)
    and then the block gofrs from the block files (two column dist gofr format)
    from the files (i.e. bl_1_volpairdist_1_1.dat)
    It then calculates the error in the RDF from the blocks, according to a 95% confidence interval assoc
    with the Student's t-distribution
    Writes to a file (i.e. volRDF_1_1_298.dat)
    """
    print("Reading in the RDF weighted by vol")
    # Read in volRDF
    # Read in eRDF
    fstring="volpairdist_"+selec1+"_"+selec2+".dat"
    dist, volgofr, vol2gofr,vol3gofr = np.genfromtxt(fstring, usecols=(0,1,2,3), unpack=True)
    # Need to add in the extra factor of -beta
    kb_JK=1.380648E-23
    AperM=1E10
    PaperBar=1E5
    InvSi_to_invBar = (1/(kb_JK*T))*(1/AperM)**3.*PaperBar 
    volgofr=volgofr*InvSi_to_invBar # Units = bar^-1
    vol2gofr=vol2gofr*InvSi_to_invBar**2.0 # Units = bar^-2
    vol3gofr=vol3gofr*InvSi_to_invBar**3.0 # Unites = bar^-3
    # Reads in block volRDF
    bl_volgofr = []
    for i in range(1,nblocks+1):
        bstring=("bl_%d_"%i)
        bl_volgofr.append(np.genfromtxt(bstring+fstring, usecols=1, unpack=True)*InvSi_to_invBar)

    err_volgofr = np.std(bl_volgofr,axis=0)*t_val
    np.savetxt('volRDF_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dist,volgofr,err_volgofr])
    np.savetxt('vol2RDF_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dist,vol2gofr])
    np.savetxt('vol3RDF_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dist,vol3gofr])
    return volgofr, bl_volgofr, err_volgofr,vol2gofr,vol3gofr

def pred_gofr_Tdep(gofr, egofr,e2gofr, bl_gofr, bl_egofr, T, Tp):
    """
    This function calculates the first order taylor series prediction based on the gofr and the derivative
    Note: egofr = <dH delta(r-r')> = - d(RDF)/d(beta)
    Thus the TS expansion is:
    RDF(beta) = RDF(beta_o) + d(RDF)/d(beta)|beta_o (beta-beta_o)
              = RDF(beta_o) - <dH delta(r-r')>|beta_o (beta-beta_o)
    This prediction is saved to a file (i.e. gofr_pred_1_1_235.0_from_298.dat)
    """
    print("Predicting the RDF")
    pred = []
    pred2 = []
    #pred3 = []
    taylorfactor=1/(kb*Tp) - 1/(kb*T) # (beta-beta_o)
    # Makes the prediction for each point in the gofr
    for i in range(len(gofr)):
        pred.append(gofr[i]-egofr[i]*taylorfactor)
        pred2.append(pred[i]+0.5*e2gofr[i]*taylorfactor**2.)
        #pred3.append(pred2[i]+(1/6)*e3gofr[i]*taylorfactor**2.)
    # Block averaging
    bl_pred = []
    for b in range(nblocks):
        tmp_pred = []
        for i in range(len(bl_gofr[b])):
            tmp_pred.append(bl_gofr[b][i]-bl_egofr[b][i]*taylorfactor)
        bl_pred.append(tmp_pred)
    err_pred = np.std(bl_pred, axis=0)*t_val
    # Saves the result to a file (i.e. gofr_pred_1_1_235.0_from_298.dat)
    if graphmovie == 0:
        np.savetxt('gofr_pred_T_'+str(selec1)+"_"+str(selec2)+"_"+str(Tp)+'_from_'+str(int(T))+'.dat', np.c_[dist,pred,err_pred,pred2])
        for b in range(nblocks):
            np.savetxt('bl_'+str(b)+'_gofr_pred_T_'+str(selec1)+'_'+str(selec2)+'_'+str(Tp)+'_from_'+str(int(T))+'.dat',np.c_[dist, bl_pred[b]])
    else:
        np.savetxt('gofr_frame_'+str(int(Tp))+'.dat', np.c_[dist,pred,err_pred])
    print(len(gofr), len(pred))
    pred_sconfig_Tdep(dist,pred,bl_pred, Tp)
    return pred

def pred_gofr_Pdep(gofr, volgofr,vol2gofr,vol3gofr, bl_gofr, bl_volgofr, P, Pp):
    """
    This function calculates the first order taylor series prediction based on the gofr and the derivative
    Note: egofr = <dH delta(r-r')> = - d(RDF)/d(beta)
    Thus the TS expansion is:
    RDF(beta) = RDF(beta_o) + d(RDF)/d(p)|p_0 (P-P_o)
              = RDF(beta_o) - <dH delta(r-r')>|P_o (P-P_o)
    This prediction is saved to a file (i.e. gofr_pred_P_1_1_100.0_from_1.0.dat)
    """
    print("Predicting the RDF")
    pred = []
    pred2 = []
    pred3 = []
    taylorfactor=(Pp-P) # (beta-beta_o)
    # Makes the prediction for each point in the gofr
    for i in range(len(gofr)):
        pred.append(gofr[i]-volgofr[i]*taylorfactor)
        pred2.append(pred[i]+0.5*vol2gofr[i]*taylorfactor**2.)
        pred3.append(pred2[i]+(1/6)*vol3gofr[i]*taylorfactor**3.)
    # Block averaging
    bl_pred = []
    for b in range(nblocks):
        tmp_pred = []
        for i in range(len(bl_gofr[b])):
            tmp_pred.append(bl_gofr[b][i]-bl_volgofr[b][i]*taylorfactor)
        bl_pred.append(tmp_pred)
    err_pred = np.std(bl_pred, axis=0)*t_val
    # Saves the result to a file (i.e. gofr_pred_1_1_235.0_from_298.dat)
    if graphmovie == 0:
        np.savetxt('gofr_pred_P_'+str(selec1)+"_"+str(selec2)+"_"+str(Pp)+'_from_'+str(int(P))+'.dat', np.c_[dist,pred,err_pred,pred2,pred3])
    else:
        np.savetxt('gofr_pframe_'+str(int(Pp))+'.dat', np.c_[dist,pred,err_pred])
    return pred


def calc_nofr(gofr, dist, bl_gofr, L, N):
    """
    This piece of code calculates the coordination number.
    This is achieved via a trapezoid rule integration
    This also calculates the block values of the coordination number
    and the uncertainty according to a 95% confidence interval
    as determined by Student's t-distribution
    """
    print("Calculating the coordination number")
    if npt == 0:
        nofr= integrate.cumtrapz(gofr*dist**2.*N/L**3.*4*np.pi, dist, initial=0)
    else:
        nofr= integrate.cumtrapz(gofr*dist**2.*N/np.average(L)**3.*4*np.pi, dist, initial=0)

    bl_nofr=[] 
    for i in range(1,nblocks+1):
        if npt == 0:
            bl_nofr.append(integrate.cumtrapz(bl_gofr[i-1]*dist**2.*N/L**3.*4*np.pi, dist, initial=0))
        else:
            bstart=int((i-1)*(len(L)-1)/nblocks)
            bend=int(i*(len(L)-1)/nblocks)
            bl_nofr.append(integrate.cumtrapz(bl_gofr[i-1]*dist**2.*N/np.average(L[bstart:bend])**3.*4*np.pi, dist, initial=0))

    err_nofr = np.std(bl_nofr,axis=0)*t_val # Note: t_val includes the 1/sqrt(nblocks) factor implicitly
    np.savetxt('nofr_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dist, nofr, err_nofr])
    return nofr, bl_nofr

def calc_dnofr_db(egofr, dist, bl_egofr, L, N):
    """
    This piece of code calculates the derivative of the coordination number.
    This is achieved via a trapezoid rule integration
    This also calculates the block values of the coordination number
    and the uncertainty according to a 95% confidence interval
    as determined by Student's t-distribution
    """
    print("Calculating the derivative of nofr")
    if npt == 0:
        dnofr= integrate.cumtrapz(-egofr*dist**2.*N/L**3.*4*np.pi, dist, initial=0)
    else:
        dnofr= integrate.cumtrapz(-egofr*dist**2.*N/np.average(L)**3.*4*np.pi, dist, initial=0)

    bl_dnofr=[]
    for i in range(1,nblocks+1):
        if npt == 0:
            bl_dnofr.append(integrate.cumtrapz(-bl_egofr[i-1]*dist**2.*N/L**3.*4*np.pi, dist, initial=0))
        else:
            bstart=int((i-1)*(len(L)-1)/nblocks)
            bend=int(i*(len(L)-1)/nblocks)
            bl_dnofr.append(integrate.cumtrapz(-bl_egofr[i-1]*dist**2.*N/np.average(L[bstart:bend])**3.*4*np.pi, dist, initial=0))

    err_dnofr = np.std(bl_dnofr,axis=0)*t_val # Note: t_val includes the 1/sqrt(nblocks) factor implicitly
    np.savetxt('dnofr_db_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dist, dnofr, err_dnofr])
    return dnofr, bl_dnofr

def pred_nofr_Tdep(nofr, dnofr, bl_nofr, bl_dnofr, T, Tp):
    """
    This function calculates the first order taylor series prediction based on the nofr and the derivative
    This prediction is saved to a file (i.e. nofr_pred_1_1_235.0_from_298.dat)
    """
    print("Predicting NofR")
    pred = []
    taylorfactor=1/(kb*Tp) - 1/(kb*T) # (beta-beta_o)
    # Makes the prediction for each point in the gofr
    for i in range(len(gofr)):
        pred.append(nofr[i]+dnofr[i]*taylorfactor)
    # Block averaging
    bl_pred = []
    for b in range(nblocks):
        tmp_pred = []
        for i in range(len(bl_nofr[b])):
            tmp_pred.append(bl_nofr[b][i]+bl_dnofr[b][i]*taylorfactor)
        bl_pred.append(tmp_pred)
    err_pred = np.std(bl_pred, axis=0)*t_val
    # Saves the result to a file (i.e. gofr_pred_1_1_235.0_from_298.dat)
    np.savetxt('nofr_pred_'+str(selec1)+"_"+str(selec2)+"_"+str(Tp)+'_from_'+str(int(T))+'.dat', np.c_[dist,pred,err_pred])
    return pred

def calc_dnofr_dp(volgofr, dist, bl_volgofr, L, N):
    """
    This piece of code calculates the derivative of the coordination number.
    This is achieved via a trapezoid rule integration
    This also calculates the block values of the coordination number
    and the uncertainty according to a 95% confidence interval
    as determined by Student's t-distribution
    """
    print("Calculating the derivative of nofr")
    if npt == 0:
        dpnofr= integrate.cumtrapz(-volgofr*dist**2.*N/L**3.*4*np.pi, dist, initial=0)
    else:
        dpnofr= integrate.cumtrapz(-volgofr*dist**2.*N/np.average(L)**3.*4*np.pi, dist, initial=0)

    bl_dpnofr=[]
    for i in range(1,nblocks+1):
        if npt == 0:
            bl_dpnofr.append(integrate.cumtrapz(-bl_volgofr[i-1]*dist**2.*N/L**3.*4*np.pi, dist, initial=0))
        else:
            bstart=int((i-1)*(len(L)-1)/nblocks)
            bend=int(i*(len(L)-1)/nblocks)
            bl_dpnofr.append(integrate.cumtrapz(-bl_volgofr[i-1]*dist**2.*N/np.average(L[bstart:bend])**3.*4*np.pi, dist, initial=0))

    err_dpnofr = np.std(bl_dpnofr,axis=0)*t_val # Note: t_val includes the 1/sqrt(nblocks) factor implicitly
    np.savetxt('dnofr_dp_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dist, dpnofr, err_dpnofr])
    return dpnofr, bl_dpnofr

def pred_nofr_pdep(nofr, dpnofr, bl_nofr, bl_dpnofr, P, Pp):
    """
    This function calculates the first order taylor series prediction based on the nofr and the derivative
    This prediction is saved to a file (i.e. nofr_pred_P_1_1_1000.0_from_1.dat)
    """
    print("Predicting NofR")
    pred = []
    taylorfactor=Pp-P # (beta-beta_o)
    # Makes the prediction for each point in the gofr
    for i in range(len(gofr)):
        pred.append(nofr[i]+dpnofr[i]*taylorfactor)
    # Block averaging
    bl_pred = []
    for b in range(nblocks):
        tmp_pred = []
        for i in range(len(bl_nofr[b])):
            tmp_pred.append(bl_nofr[b][i]+bl_dpnofr[b][i]*taylorfactor)
        bl_pred.append(tmp_pred)
    err_pred = np.std(bl_pred, axis=0)*t_val
    # Saves the result to a file (i.e. gofr_pred_1_1_235.0_from_298.dat)
    np.savetxt('nofr_pred_P_'+str(selec1)+"_"+str(selec2)+"_"+str(Pp)+'_from_'+str(int(P))+'.dat', np.c_[dist,pred,err_pred])
    return pred


def peak_find(dist, gofr, nofr, thresh):
    """
    This is a cute little script for finding the peaks and minima in the RDF
    It prints the value of the maxima, the minima, and the coordination number at that point. 
    The "thresh" parameter determines the threshold by which a peak must stand out compared to
    its neighbors. 
    
    This code writes these values to both the screen, and to a minimum and maximum file.
    """
    print("Finding peaks in the RDF")
    # Maxima
    mlocs, dict =  find_peaks(gofr, threshold=thresh, height=1)
    maxvals = []
    dmax = []
    nmax = []
    for i in mlocs:
        print('Maximum value at %4.5f %4.5f with coordination number %4.5f' %(dist[i],gofr[i],nofr[i]))
        nmax.append(nofr[i])
        maxvals.append(gofr[i])
        dmax.append(dist[i])
    minvals = []
    dmin = []
    nmin = []
    for i in range(len(mlocs)-1):
        loc = np.where(gofr == min(gofr[mlocs[i]:mlocs[i+1]]))[0][0]
        minvals.append(gofr[loc])
        dmin.append(dist[loc])
        nmin.append(nofr[loc])
        print('Minimum value at %4.5f %4.5f with coordination number %4.5f' %(dist[loc],gofr[loc],nofr[loc]))
    # Write maximums to file (i.e. gofr_maximums_1_1_298.dat)
    f = open('gofr_maximums_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat','w')
    f.write("#distance(Angstroms) gofrmax nofrmax\n")
    for i in range(len(maxvals)):
        f.write("{0:.5f} {1:.5f} {2:.5f}\n".format(dmax[i],maxvals[i],nmax[i]))
    f.close()
    # Write minimums to file (i.e. gofr_minimums_1_1_298.dat
    g = open('gofr_minimums_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat','w')
    g.write("#distance(Angstroms) gofrmin nofrmin\n")
    for i in range(len(minvals)):
        g.write("{0:.5f} {1:.5f} {2:.5f}\n".format(dmin[i],minvals[i],nmin[i]))
    g.close()
    return maxvals, dmax, minvals, dmin

def calc_PMF(dist, gofr, bl_gofr, T):
    """
    This section of the code calculates the potential of mean force
    The PMF is defined as PMF = -kb*T*log[RDF]
    This is written to a file (i.e. pmf_1_1_298.dat)
    """
    print("Calculating the PMF")
    PMF = []
    aofr=[]
    pdist=[]
    c=0
    for i in range(len(gofr)):
        if i > 1 and gofr[i-1] > 0:
            PMF.append(-kb*T*np.log(gofr[i]))
            aofr.append(-kb*T*np.log(gofr[i])-2*kb*T*np.log(dist[i]))
        else:
            c += 1

    bl_pmf = []
    bl_aofr = []
    for b in range(len(bl_gofr)):
        tmp_pmf = []
        tmp_aofr = []
        for i in range(len(bl_gofr[b])):
            if i > 1 and gofr[i-1] > 0:
                tmp_pmf.append(-kb*T*np.log(bl_gofr[b][i]))
                tmp_aofr.append(-kb*T*np.log(bl_gofr[b][i])-2*kb*T*np.log(dist[i]))
        bl_pmf.append(tmp_pmf)
        bl_aofr.append(tmp_aofr)
    err_pmf = np.std(bl_pmf,axis=0)*t_val
    err_aofr = np.std(bl_aofr, axis=0)*t_val
    np.savetxt("pmf_"+str(selec1)+"_"+str(selec2)+'_'+str(int(T))+".dat", np.c_[dist[c:], PMF,err_pmf])
    np.savetxt("A_"+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+".dat", np.c_[dist[c:], aofr, err_aofr])
    return PMF, c, bl_pmf, aofr, bl_aofr

def calc_dPMF_db(PMF, dist, gofr, egofr, bl_gofr, bl_egofr, bl_PMF,T):
    """
    This script calculates the first derivative of the PMF
    The derivative is defined as:
    kb*T*(eRDF/RDF - PMF)
    This is written to a file(i.e. dpmf_1_1_298.dat)
    """
    print("Calculating the PMF T-derivative")
    dPMF = []
    daofr = []
    for i in range(len(PMF)):
        dPMF.append(kb*T*(egofr[i]/gofr[i]-PMF[i]))
        daofr.append(kb*T*(egofr[i]/gofr[i]-PMF[i])+2*kb**2.*T**2.*np.log(dist[i]))
    bl_dPMF = []
    bl_daofr = []
    for b in range(len(bl_PMF)):
        tmp_dPMF = []
        tmp_daofr = []
        for i in range(len(bl_PMF[b])):
            tmp_dPMF.append(kb*T*(bl_egofr[b][i]/bl_gofr[b][i]-bl_PMF[b][i]))
            tmp_daofr.append(kb*T*(bl_egofr[b][i]/bl_gofr[b][i]-bl_PMF[b][i])+2*kb**2.*T**2.*np.log(dist[i]))
        bl_dPMF.append(tmp_dPMF)
        bl_daofr.append(tmp_daofr)
    err_dPMF = np.std(bl_dPMF,axis=0)*t_val
    err_daofr = np.std(bl_daofr, axis=0)*t_val
    np.savetxt("dpmf_"+str(selec1)+"_"+str(selec2)+'_'+str(int(T))+".dat", np.c_[dist, dPMF, err_dPMF])
    np.savetxt("dA_"+str(selec1)+"_"+str(selec2)+'_'+str(int(T))+".dat", np.c_[dist, daofr, err_daofr])
    return dPMF, bl_dPMF, daofr, bl_daofr

def pred_PMF_Tdep(dist, PMF, dPMF, bl_PMF, bl_dPMF, T, Tp):
    """
    This script makes predictions of the PMF at other temperature(s) Tp from the simulation temperature T
    This is done by a standard taylor series expansion out to **one** term
    PMF(beta) = PMF(beta_o) + dPMF/dbeta|beta_o*(beta-beta_o)
    This is written to a file (i.e. pmf_pred_1_1_235.0_from_298.dat)
    """
    print("Predicting the PMF at %s" % Tp)
    pred = []
    taylorfactor=1/(kb*Tp)-1/(kb*T)
    print(taylorfactor)
    for i in range(len(PMF)):
        pred.append(PMF[i]+dPMF[i]*taylorfactor)
    # Block Calculation
    bl_pmfpred = []
    for b in range(len(bl_PMF)):
        tmp_pPMF = []
        for i in range(len(bl_PMF[b])):
            tmp_pPMF.append(bl_PMF[b][i]+bl_dPMF[b][i]*taylorfactor)
        bl_pmfpred.append(tmp_pPMF)
    err_pmfpred = np.std(bl_pmfpred,axis=0)*t_val
    np.savetxt('pmf_pred_'+str(selec1)+'_'+str(selec2)+'_'+str(Tp)+'_from_'+str(int(T))+'.dat', np.c_[dist,pred,err_pmfpred])
    return pred

def pred_AofR_Tdep(dist, aofr, daofr, bl_aofr, bl_daofr, T, Tp):
    """
    This script makes predictions of the aofr at other temperature(s) Tp from the simulation temperature T
    This is done by a standard taylor series expansion out to **one** term
    aofr(beta) = aofr(beta_o) + daofr/dbeta|beta_o*(beta-beta_o)
    This is written to a file (i.e. aofr_pred_1_1_235.0_from_298.dat)
    """
    print("Predicting the aofr at %s" % Tp)
    pred = []
    taylorfactor=1/(kb*Tp)-1/(kb*T)
    print(taylorfactor)
    for i in range(len(PMF)):
        pred.append(aofr[i]+daofr[i]*taylorfactor)
    # Block Calculation
    bl_aofrpred = []
    for b in range(len(bl_PMF)):
        tmp_paofr = []
        for i in range(len(bl_aofr[b])):
            tmp_paofr.append(bl_aofr[b][i]+bl_daofr[b][i]*taylorfactor)
        bl_aofrpred.append(tmp_paofr)
    err_aofrpred = np.std(bl_aofrpred,axis=0)*t_val
    np.savetxt('aofr_pred_'+str(selec1)+'_'+str(selec2)+'_'+str(Tp)+'_from_'+str(int(T))+'.dat', np.c_[dist,pred,err_aofrpred])
    return pred

def calc_dS(PMF, dPMF,bl_PMF,bl_dPMF, dist, T):
    """
    This is a calculation of the entropy from the derivative of the PMF
    dS = 1/(kb*T^2)*(dPMF/dbeta + 2*(kb*T)^2 * log[r])
    Results are written to a file (i.e. dS_1_1.dat)
    Units are (kcal/mol)/K
    """
    print("Calcuating the entropy change")
    dS = []
    for i in range(len(PMF)):
        dS.append(1/(kb*T**2.)*(dPMF[i]+2*(kb*T)**2.*np.log(dist[i])))
    bl_dS = []
    for b in range(len(bl_PMF)):
        tmp_dS = []
        for i in range(len(bl_PMF[b])):
            tmp_dS.append(1/(kb*T**2.)*(bl_dPMF[b][i]+2*(kb*T)**2.*np.log(dist[i])))
        bl_dS.append(tmp_dS)
    err_dS = np.std(bl_dS, axis=0)*t_val
    np.savetxt('dS_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dist,dS,err_dS])
    return dS

def calc_sconfig(dist, gofr, bl_gofr):
    """
    This is a calculation of the configurational entropy from the equation
    s2/N = -4*PI/2*N^2*kb/V*int r^2*dr*(g(r)ln[g(r)]-g(r)+1)
    Results are written to a file (i.e. s2_1_1.dat)
    """
    print("Calculating the configurational entropy")
    val=[]
    dred=[]
    if npt == 0:
        V=L**3.0
    else:
        V=np.average(L)**3.0

    # Total Calculation
    for i in range(len(gofr)):
        if gofr[i] <= 0.0:
            val.append(0.0)
        else:
            val.append(dist[i]**2*(gofr[i]*np.log(gofr[i])-gofr[i]+1))
    s2=-(4*np.pi)/2*N/V*kb*integrate.cumtrapz(val, dist,initial=0)
    np.savetxt('s2_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dist, s2])
    bl_s2 = []
    # Block Calculation
    tb = []
    for b in range(nblocks):
        tmp_val = []
        for i in range(len(bl_gofr[b])):
            if bl_gofr[b][i] <= 0.0: # This makes sure it isn't less than or equal to zero
                tmp_val.append(0)
            else:
                tmp_val.append(dist[i]**2*(bl_gofr[b][i]*np.log(bl_gofr[b][i])-bl_gofr[b][i]+1))
        if npt == 0:
            bl_s2.append(-(4*np.pi)/2*N/V*kb*integrate.cumtrapz(tmp_val, dist,initial=0))
        else:
            bstart=int(b*(len(L)-1)/nblocks)
            bend=int((b+1)*(len(L)-1)/nblocks)
            V=np.average(L[bstart:bend])**3.0
            bl_s2.append(-(4*np.pi)/2*N/V*kb*integrate.cumtrapz(tmp_val, dist,initial=0))

        np.savetxt('bl_'+str(b)+'_s2_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dist, bl_s2[b]])
        tb.append(bl_s2[b][-1])
    err_tb = np.std(tb)*t_val
    print("Configurational Entropy %s %s" % ( np.average(tb), err_tb))
   
    f=open('s2_val_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', 'w')
    f.write('%s %s %s\n' % (T, np.average(tb), err_tb))
    f.close()

    err_s2 = np.std(bl_s2, axis=0)*t_val

    # Writes to a file
    np.savetxt('s2_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dist, s2, err_s2])

    return s2


def calc_dsconfig_db(dist, gofr, bl_gofr, dgofr, bl_dgofr):
    """
    This calculates the beta derivative of the configurational entropy s2.
    \frac{\partial s2}{\partial\beta}=-\frac{4\pi N k_B}{V}\int_0^\infty r^2 dr \frac{\partial g(r)}{\partial \beta}\ln{g(r)}
    """
    print("Calculating the derivative of s2")
    val=[]
    dred=[]
    if npt == 0:
        V=L**3.0
    else:
        V=np.average(L)**3.0

    for i in range(len(gofr)):
        if gofr[i] > 0.0:
            dred.append(dist[i])
            val.append(dist[i]**2*(dgofr[i]*np.log(gofr[i])))
    ds2 = -(4*np.pi)/2*N/V*kb*integrate.cumtrapz(val, dred, initial=0)
    np.savetxt('ds2_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dred,ds2])
    bl_ds2 = []
    dtb = []
    # Block Calculation
    for b in range(nblocks):
        tmp_val = []
        for i in range(len(bl_gofr[b])):
            if gofr[i] > 0.0:
                if bl_gofr[b][i] <= 0.0: # This makes sure it isn't less than or equal to zero
                    tmp_val.append(0)
                else:
                    tmp_val.append(dist[i]**2*bl_dgofr[b][i]*np.log(bl_gofr[b][i]))
        if npt == 0:
            bl_ds2.append(-(4*np.pi)/2*N/V*kb*integrate.cumtrapz(tmp_val, dred,initial=0))
        else:
            bstart=int(b*(len(L)-1)/nblocks)
            bend=int((b+1)*(len(L)-1)/nblocks)
            V=np.average(L[bstart:bend])**3.0
            bl_ds2.append(-(4*np.pi)/2*N/V*kb*integrate.cumtrapz(tmp_val, dred,initial=0))
        dtb.append(bl_ds2[b][-1])

    errdtb = np.std(dtb)*t_val
    print("Derivative of configurational entropy %s %s" % (np.average(dtb), errdtb))
    f=open('ds2_val_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', 'w')
    f.write('%s %s %s\n' % (T, np.average(dtb), errdtb))
    f.close()

    err_ds2 = np.std(bl_ds2, axis=0)*t_val

    # Writes to a file
    np.savetxt('ds2_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[dred, ds2, err_ds2])
    
    return ds2

def pred_sconfig_Tdep(dist,predgofr,bl_predgofr,Tp):
    """This takes the predicted gofr, and then uses it to recalculate the s2"""
    print("Predicting the configurational entropy")
    val=[]
    dred=[]
    if npt == 0:
        V=L**3.0
    else:
        V=np.average(L)**3.0
    # Total Calculation
    for i in range(len(predgofr)):
        if predgofr[i] <= 0.0:
            val.append(0.0)
        else:
            val.append(dist[i]**2*(predgofr[i]*np.log(predgofr[i])-predgofr[i]+1))
    ps2=-(4*np.pi)/2*N/V*kb*integrate.cumtrapz(val, dist,initial=0)[-1]
    blps2=[]
    for b in range(nblocks):
        blval=[]
        for i in range(len(bl_predgofr[b])):
            if bl_predgofr[b][i]<= 0.0:
                blval.append(0)
            else:
                blval.append(dist[i]**2*(bl_predgofr[b][i]*np.log(bl_predgofr[b][i])-bl_predgofr[b][i]+1))
        tmp = -(4*np.pi)/2*N/V*kb*integrate.cumtrapz(blval, dist,initial=0)
        blps2.append(-(4*np.pi)/2*N/V*kb*integrate.cumtrapz(blval, dist,initial=0)[-1])
    ps2err = np.std(blps2)*t_val


    f=open('s2_pred_T_'+str(selec1)+"_"+str(selec2)+"_"+str(Tp)+'_from_'+str(int(T))+'.dat', 'w')
    f.write('%s %s %s\n' % (Tp, ps2,ps2err))
    f.close()
    return

def calc_U(r, g, eg, bl_g, bl_eg,item):
    """
    This calculates the total energy
    """
    U=[]
    for i in range(len(g)):
        if g[i] != 0:
            U.append(eg[i]/g[i])
        else:
            U.append(0)
    bl_U=[]
    for b in range(len(bl_g)):
        tmp_U=[]
        for i in range(len(bl_g[b])):
            if bl_g[b][i] != 0:
                tmp_U.append(bl_eg[b][i]/bl_g[b][i])
            else:
                tmp_U.append(0)
        bl_U.append(tmp_U)
    err_U=np.std(bl_U, axis=0)*t_val
    np.savetxt('U'+str(item)+'_'+str(selec1)+'_'+str(selec2)+'_'+str(int(T))+'.dat', np.c_[r, U, err_U])
    return U, bl_U



def peak_find_simple(r, g, U, A, S,thresh,T):
    """
    """
    print("Finding simple peaks in the RDF")
    # Maxima
    mlocs, dict =  find_peaks(g, threshold=thresh, height=1)
    maxvals = []
    dmax = []
    for i in mlocs:
        maxvals.append(g[i])
        dmax.append(r[i])
    minvals = []
    dmin = []
    for i in range(len(mlocs)-1):
        loc = np.where(g == min(g[mlocs[i]:mlocs[i+1]]))[0][0]
        minvals.append(g[loc])
        dmin.append(r[loc])

    print("There are %d minimums and %d maximums at temperature %s" % (len(minvals), len(maxvals), T))
    outpeak(r, U, A, S, dmin, dmax,T)
    return maxvals, dmax, minvals, dmin

def pred_T(r, g, U, bl_g, bl_U, S, T, Tp):
    """
    This is the optimal way to make predictions of the temperature dependence, assuming that U is temperature independent 
    This is typically a good assumption.
    g(r;\beta)=g(r;\beta_{298})e^{-(\beta-\beta_{298})\Delta U(r)}
    """
    pg      = []
    bl_pg   = []
    err_pg  = []
    pPMF     = []
    bl_pPMF  = []
    err_pPMF = []
    pA       = []
    bl_pA    = []
    err_pA  = []
    r_red   = []
    beta_Tp = 1/(kb*Tp)
    beta_T  = 1/(kb*T)
    c = 0
    for i in range(len(g)):
        pg.append(g[i]*np.exp(-U[i]*(beta_Tp - beta_T)))
        if i > 1 and pg[i-1] > 0:
            r_red.append(r[i])
            pPMF.append(-kb*Tp*np.log(pg[i]))
            pA.append(pPMF[i-c]-2*kb*Tp*np.log(r[i]))
        else:
            c += 1
    
    for b in range(nblocks):
        tmp_pg   = []
        tmp_pPMF = []
        tmp_pA   = []
        for i in range(len(bl_g[b])):
            tmp_pg.append(bl_g[b][i]*np.exp(-bl_U[b][i]*(beta_Tp - beta_T)))
            if i > 1 and pg[i-1] > 0.0:
                tmp_pPMF.append(-kb*Tp*np.log(tmp_pg[i]))
                tmp_pA.append(tmp_pPMF[i-c]-2*kb*Tp*np.log(r[i]))
        bl_pg.append(tmp_pg)
        bl_pPMF.append(tmp_pPMF)
        bl_pA.append(tmp_pA)
    print(len(S),len(pg[c:]),len(U[c:]), len(pA))
    err_pg = np.std(bl_pg, axis=0)*t_val
    err_pPMF = np.std(bl_pPMF, axis=0)*t_val
    err_pA = np.std(bl_pA, axis=0)*t_val
    np.savetxt('pred_gofr_'+str(selec1)+'_'+str(selec2)+'_'+str(Tp)+'_from_'+str(int(T))+'.dat', np.c_[r, pg, err_pg])
    np.savetxt('pred_PMF_'+str(selec1)+'_'+str(selec2)+'_'+str(Tp)+'_from_'+str(int(T))+'.dat', np.c_[r_red, pPMF, err_pPMF])
    np.savetxt('pred_A_'+str(selec1)+'_'+str(selec2)+'_'+str(Tp)+'_from_'+str(int(T))+'.dat', np.c_[r_red, pA, err_pA])
    peak_find_simple(r[c:], pg[c:], U[c:], pA, S, thresh,Tp)

    return pg, bl_pg, err_pg
            

def outpeak(r, U, A, S, minloc,maxloc,T):
    dr = np.where(r == maxloc[0])[0][0]
    db = np.where(r == minloc[0])[0][0]
    dp = np.where(r == maxloc[1])[0][0]
           
    Ar = A[dr]
    Sr = S[dr]
    Ur = U[dr]
    Ab = A[db]
    Sb = S[db]
    Ub = U[db]
    Ap = A[dp]
    Sp = S[dp]
    Up = U[dp]

    Aval = (2*Ab-Ar-Ap)
    TSval = T*(2*Sb-Sr-Sp)
    Uval = (2*Ub-Ur-Up)
    print("***")
    print("Temperature: %s" % T)
    print(r[dr],r[db],r[dp])
    print("Activation Free Energy %s (kcal/mol)" % (2*Ab-Ar-Ap))
    print("Activation Entropy %s (kcal/mol/K)" % (2*Sb-Sr-Sp))
    print("Activation T*S %s (kcal/mol)" % (T*(2*Sb-Sr-Sp)))
    print("Activation U %s (kcal/mol)" % (2*Ub-Ur-Up))
    print("***")
    if T == Tp[0]:
        f = open('activation_params.dat','w')
        f.write('#T A TS U\n')
    else:
        f = open('activation_params.dat', 'a')

    f.write("%s %s %s %s\n" % (T, Aval, TSval, Uval))
    f.close()
    
    return
    

import numpy as np
import argparse, math
import os
#import matplotlib.pyplot as plt
from scipy import stats,integrate
from scipy.signal import find_peaks

kb=0.0019872041 #kcal/mol
convp=1.439E-5 #bar A^3 to kcal/mol


# Command Line Input
parser = argparse.ArgumentParser()
parser.add_argument('-plot', default=0, help='[0] No plot [1] Plots in matplotlib')
parser.add_argument('-f', default="gofr.inp", help="Namelist file for the fortran code, extracts L, selec1, selec2")
parser.add_argument('-selec1',default=0, help='Selection1')
parser.add_argument('-selec2', default=0, help='Selection2')
parser.add_argument('-nblocks',default=5, help='Number of blocks')
parser.add_argument('-N', default=343, help='Number of atoms')
parser.add_argument('-T', default=298.15, help='Temperature (K)')
parser.add_argument('-Tpred', default=[], help='Temperature (K) to predict',action='append')
parser.add_argument('-thresh', default=0.05, help='Threshold for peak finding')
parser.add_argument('-P', default=1.0, help="Pressure (bar)")
parser.add_argument('-Ppred', default=[], help='Pressure (bar) to predict', action='append')
parser.add_argument('-ew', default=0, help='Energy Weighting [0] off [1] on')
parser.add_argument('-vw', default=0, help='Volume Weighting [0] off [1] on')
parser.add_argument('-graphmovie', default=0, help='Make movie [0] off [1] on')
args = parser.parse_args()
selec1  = str(args.selec1)
selec2  = str(args.selec2)
nblocks = int(args.nblocks)
N       = int(args.N)
T       = float(args.T)
Tp      = args.Tpred
P       = float(args.P)
Pp      = args.Ppred
ew      = int(args.ew)
vw      = int(args.vw)
graphmovie = int(args.graphmovie)
pltflag = int(args.plot)

L, selec1, selec2, npt=read_input(str(args.f), selec1, selec2)


for i in range(len(Tp)):
    Tp[i] = float(Tp[i])

#Tp = []

#for i in range(200,360,5):
#    Tp.append(float(i))


if args.thresh is not None:
    thresh=float(args.thresh)
else:
    thresh=0.001

t_val=stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)


# Read GofR
dist, gofr, bl_gofr, err_gofr = read_gofr(t_val,selec1, selec2)

# Calculates the Coordination Number
nofr, bl_nofr = calc_nofr(gofr, dist, bl_gofr, L, N)

# Finds the local minima and maxima
maxvals, dmax, minvals, dmin = peak_find(dist, gofr, nofr,thresh)

# Calculates the PMF
PMF, c, bl_PMF, aofr, bl_aofr= calc_PMF(dist,gofr, bl_gofr, T)

# Calculates the configurational entropy
s2 = calc_sconfig(dist, gofr, bl_gofr)


# Temperature Derivative Calculations
if ew == 1:
    # Reads beta derivative stuff, and makes GofR prediction
    egofr, bl_egofr, err_egofr, e2gofr = read_egofr(t_val, selec1, selec2,"e")
    if os.path.isfile('ljpairdist_'+str(selec1)+'_'+str(selec2)+'.dat'):
        ljgofr, bl_ljgofr, err_ljgofr, lj2gofr = read_egofr(t_val, selec1, selec2,"lj")
        Ulj, bl_Ulj=calc_U(dist, gofr, ljgofr, bl_gofr, bl_ljgofr,"lj")
    if os.path.isfile('kepairdist_'+str(selec1)+'_'+str(selec2)+'.dat'):
        kegofr, bl_kegofr, err_kegofr, ke2gofr = read_egofr(t_val, selec1, selec2,"ke")
        calc_U(dist, gofr, kegofr, bl_gofr, bl_kegofr,"ke")
    if os.path.isfile('volpairdist_'+str(selec1)+'_'+str(selec2)+'.dat'):
        volgofr, bl_volgofr, err_volgofr, vol2gofr = read_egofr(t_val, selec1, selec2,"vol")
        calc_U(dist, gofr, volgofr, bl_gofr, bl_volgofr, "vol")
        # This adds the volume component to the gofr derivative, this makes it so that
        # it actually gets the whole thing correctly!
        egofr = egofr + P*convp*volgofr
        for b in range(nblocks):
            bl_egofr[b] = bl_egofr[b]+P*convp*bl_volgofr[b]

    # Calculation of the internal energy
    U,bl_U=calc_U(dist, gofr, egofr, bl_gofr, bl_egofr,"e")
    pGofR =  []
    if graphmovie == 0:
        for tp in Tp:
            pGofR.append(pred_gofr_Tdep(gofr, egofr, e2gofr, bl_gofr, bl_egofr, T, float(tp)))
    else:
        for tp in range(220,400):
            pGofR.append(pred_gofr_Tdep(gofr, egofr,e2gofr, bl_gofr, bl_egofr, T, float(tp)))
    # Calculates NofR derivative and predicts it!
    dnofr, bl_dnofr = calc_dnofr_db(egofr, dist, bl_egofr, L, N)
    pnofr = []
    for tp in Tp:
        pnofr.append(pred_nofr_Tdep(nofr,dnofr,bl_nofr, bl_dnofr, T, float(tp)))

    bl_cutgofr = [ bl_gofr[b][c:] for b in range(nblocks) ]
    bl_cutegofr = [ bl_egofr[b][c:] for b in range(nblocks) ]
    dPMF,bl_dPMF,daofr,bl_daofr=calc_dPMF_db(PMF, dist[c:], gofr[c:], egofr[c:], bl_cutgofr, bl_cutegofr, bl_PMF, T)
    # PMF Predictions
    pPMF =  []
    paofr = []
    for tp in Tp:
        pPMF.append(pred_PMF_Tdep(dist[c:], PMF, dPMF,bl_PMF,bl_dPMF, T, tp))
        paofr.append(pred_AofR_Tdep(dist[c:], aofr, daofr,bl_aofr,bl_daofr, T, tp))

    # Calculates Entropy Derivative
    dS = calc_dS(PMF, dPMF,bl_PMF,bl_dPMF, dist[c:], T)
    calc_dsconfig_db(dist, gofr, bl_gofr,egofr,bl_egofr)
    outpeak(dist[c:], U[c:], aofr, dS, dmin, dmax, T)
    
    # Better Prediction
    if graphmovie == 0:
        for tp in Tp:
            pred_T(dist, gofr, U, bl_gofr, bl_U, dS, T, tp)
    else:
        for tp in range(150,400, 10):
            pred_T(dist, gofr, U, bl_gofr, bl_U, dS, T, tp)

# Pressure Derivative Calculations!
if vw == 1:
    # Reads in derivative correlation function and Predicts it.
    volgofr, bl_volgofr, err_volgofr, vol2gofr,vol3gofr = read_volgofr(t_val, selec1, selec2)
    pGofR_P = []
    if graphmovie == 0:
        for pp in Pp:
            pGofR_P.append(pred_gofr_Pdep(gofr,volgofr, vol2gofr,vol3gofr, bl_gofr, bl_volgofr, P, float(pp)))
    else:
        for pp in range(1,10000,200):
            pGofr_P.append(pred_gofr_Pdep(gofr, volgofr, vol2gofr,vol3gofr, bl_gofr, bl_volgofr, P, float(pp)))

    # Calculates NofR Derivative and Predicts it
    dpnofr, bl_dpnofr = calc_dnofr_dp(volgofr, dist, bl_volgofr, L, N)
    pnofr_P = []
    for pp in Pp:
        pnofr_P.append(pred_nofr_pdep(nofr,dpnofr,bl_nofr, bl_dpnofr, P, float(pp)))
    # Creates cut arrays for PMF
    bl_cutgofr = [ bl_gofr[b][c:] for b in range(nblocks) ]
    bl_cutvolgofr = [ bl_volgofr[b][c:] for b in range(nblocks) ]
    # Calculates PMF derivative and predicts it
"""
    dPMF,bl_dPMF=calc_dPMF_db(PMF, dist[c:], gofr[c:], egofr[c:], bl_cutgofr, bl_cutvolgofr, bl_PMF, T)
    pPMF_P =  []
    for pp in Pp:
        pPMF_P.append(pred_PMF_Tdep(dist[c:], PMF, dpPMF,bl_PMF,bl_dpPMF, T, tp))
"""
    
if pltflag == 1:
    plt.plot(dmax,maxvals,'bs')
    plt.plot(dmin,minvals,'rs')
    plt.plot(dist,gofr,'-b')
    plt.show()
    if ew == 1:
        plt.plot(dist,egofr)
        plt.show()
    plt.plot(dist,nofr)
    plt.show()
    plt.plot(dist[c:],PMF)
    if ew == 1:
        plt.plot(dist[c:],dPMF)
    plt.show()
    



