import numpy as np
def run_stwham(histograms,minval, maxval,lambdas,eta, Ho,bl):
    """
    This is a variant of ST-WHAM (based on David Stelter's code) that calculates
    the effective temperature and combines the values
    The bl parameter sets which block (or total) is being run.
    """
    kb=0.0019872041
    hist = []
    hist=histograms
    nbins,nwalkers= np.shape(histograms)[0],np.shape(histograms)[1]
    print(nwalkers,nbins)
    checklimit=0.00001
    ### This is a function that calculates STWHAM
    def EffTemp(lambdavalue, H):
        # Evaluates the gREM effective temperature
        Teff = lambdavalue + eta*(H - Ho)
        return Teff

    # Write out the lambdas
    if bl == -1:
        h = np.linspace(minval,maxval,num=200)
        out=[h]
        for i in lambdas:
            tval = EffTemp(i,h)
            out.append(tval)
        out = np.array(out)
        np.savetxt("lambda_funcs.dat",out.T)
            
    # This is a function to run ST-WHAM on histogrammed data for gREM simulations
    pdf      = np.zeros((nbins, nwalkers))
    count    = np.zeros(nwalkers)
    totpdf   = np.zeros(nbins)
    totcount = 0
    binsize  = (maxval-minval)/float(nbins)
    print("ST WHAM Parameters")
    print("binsize = %10.5f" % binsize)
    print("minval  = %10.5f" % minval)
    print("maxval  = %10.5f" % maxval)
    print("nbins   = %10.5f" % nbins)
    # This part builds the pdfs from the histograms
    for rep in range(nwalkers):
        count[rep]=0
        for i in range(nbins):
            pdf[i][rep] = hist[i][rep]  #pdf of each replica
            count[rep] = count[rep] + pdf[i][rep] #count
            totpdf[i] = totpdf[i]+pdf[i][rep] # adds to total pdf
        pdf[:,rep] = pdf[:,rep]/count[rep] # Normalized replica pdf
        totcount = totcount + count[rep] # Total number of data points
    totpdf = totpdf/totcount # Normalized total PDF
    
    ### This next section calculates Ts(H) ###
    # Make sure that numerical derivative is defined
    TH, S         = np.zeros(nbins), np.zeros(nbins)
    betaH , betaW = np.zeros(nbins), np.zeros(nbins)
    hfrac = np.zeros((nbins,nwalkers))
    bstart, bstop = None, None
    for i in range(nbins):
        if (totpdf[i] > checklimit):
            bstart = i + 3
            break
    if (bstart == None): exit("Error: Lower end of logarithm will be undefined")
    for i in range(nbins-1, bstart, -1):
        if (totpdf[i] > checklimit):
            bstop = i - 3
            print("chlim",bl,bstop)
            break
    if (bstop == None): exit("Error: Upper end of logarithm will be undefined")
    # If these checks pass, we should be good to proceed.
    # Calculate Ts(H)
    for i in range(bstart,bstop):
        if (totpdf[i+1] > checklimit and totpdf[i-1]>checklimit): # makes sure it is above checklimit
            betaH[i] = np.log(totpdf[i+1]/totpdf[i-1])/(2*binsize)*kb
        else: # otherwise 0
            betaH[i] = 0
        for rep in range(nwalkers):
            if (totpdf[i] > 0): #positive, not empty
                H      = minval + (i*binsize+binsize/2)
                weight = np.nan_to_num(1/EffTemp(lambdas[rep],H))
                betaW[i] = betaW[i] + ((count[rep]*pdf[i][rep])/(totcount*totpdf[i])*weight)
        TH[i] = 1 / (betaH[i] + betaW[i])

    ### This Section calculates the histogram fraction
    if bl == -1:
        with open("histfrac_STWHAM.out",'w') as fracout:
            for rep in range(nwalkers):
                for i in range(bstart,bstop):
                    hfrac[i][rep] = hist[i][rep]/(totpdf[i]*totcount) # how many counts of total
                    fracout.write("%f %f\n" % (minval + (i*binsize), hfrac[i][rep]))
                fracout.write("\n")
    
    ### This Section calculates the entropy
    def Falpha(i,j,T,binsize):
        # This function interpolates the entropy based on Ts(H)
        f = 0
        for indx in range(i+1,j):
            if (T[indx] == T[indx-1]):
                f = f + binsize/T[indx]
            else:
                f = f + binsize/(T[indx]-T[indx-1])*np.log(T[indx]/T[indx-1])
        return f
    for i in range(bstart,bstop):
        S[i] = Falpha(bstart,i,TH,binsize)
    if bl == -1:
        with open("TandS_STWHAM.out",'w') as tout:
            for i in range(bstart,bstop):
                tout.write("%f %f %f\n" % (minval+(i*binsize), TH[i],S[i]))
    return TH, S, bstart, bstop, binsize
            
