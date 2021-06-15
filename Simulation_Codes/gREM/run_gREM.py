#!/usr/bin/env python
########################################################################################
# This code is a personal version written by Zeke Piskulich
# It takes logfiles from replica exchange simulations and then pulls the data from them
# it takes heavy inspiration from the code written by David Stelter for the gREM method
# in 2015. 
# For questions/concerns reach out to Zeke Piskulich (piskulichz@gmail.com)
# Copyright 2021 - Boston University
########################################################################################

import numpy as np
import argparse

consts={"kb":0.0019872041} # kcal/mol

# Read Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-start', default=0,type=int, help='Starting index for calculating [default=0]')
parser.add_argument('-stop', default=5, type=int, help='Ending index for calculating [default=5]')
parser.add_argument('-workdir', default=".", type=str, help='Working directory location [default=./]')
parser.add_argument('-eta', default=-0.03, type=float, help='Eta value used in the simulations [default=-0.03]')
parser.add_argument('-Ho', default=-30000, type=float, help='H0 value used in the simulations [default=-30000]')
parser.add_argument('-nbins', default=20, type=int, help='Number of bins used for histogramming [default=20]')
args= parser.parse_args()

fstart  = args.start
fend    = args.stop
workdir = args.workdir
eta     = args.eta
Ho      = args.Ho
nbins   = args.nbins

class walkers:
    def __init__(self,workdir,nwalkers,start,stop,ndata,walkloc,eta,Ho):
        self.walkdata,self.repdata={},{}
        self.nwalkers, self.ndata=nwalkers, ndata
        self.start,self.stop=start,stop
        self.walkloc=walkloc
        self.logbase="/log.lammps.w-"
        self.workdir=workdir
        self.hist,self.edges={},{}
        self.minval, self.maxval={},{}
        self.eta, self.Ho = eta, Ho
        print("Run Parameters:")
        print("workdir = %s"%self.workdir)
        print("logbase = %s"%self.logbase)
        print("start, stop = %d, %d" % (self.start, self.stop))
        print("eta, Ho = %10.5f, %10.5f" %(self.eta, self.Ho))
        print("nwalkers, ndata = %d, %d" %(self.nwalkers, self.ndata))
        print("Walker class initiated based on following log format %s" % self.logbase)
        print("Pulling data")
        self.walk_data()
        self.get_replicadata()
        self.write()
        print("Data written.")
    def post_process(self,nbins):
        print("Beginning Post Processing")
        self.calc_Teff()
        for key in self.keys:
            self.hist_data(key,nbins)
        print("Post Processing Completed")
    def walk_data(self):
        # Read the keys and first vals
        self.nruns=self.stop-self.start
        self.keys, firstvals = read_lammps_log(self.workdir+"/"+str(0)+self.logbase, self.start+1,self.ndata)
        perwindow=len(firstvals[self.keys[0]])
        for key in self.keys:
            self.walkdata[key]=np.zeros((self.nwalkers,self.nruns*perwindow))
            self.repdata[key]=np.zeros((self.nwalkers,self.nruns*perwindow))
        for w in range(self.nwalkers):
            for r in range(self.nruns):
                rstart, rstop = r*perwindow, (r+1)*perwindow
                # read in the run
                _, values = read_lammps_log(self.workdir+"/"+str(w)+self.logbase, r+1,self.ndata)
                for key in self.keys:
                    # copy into array
                    if key == "Step":
                        values[key] = (values[key][-1]-values[key][0])*r+values[key]
                    self.walkdata[key][w][rstart:rstop]=values[key]
    def write(self):
        print("Writing walkers and replicas.")
        self.write_array(self.walkdata, "walkers")
        self.write_array(self.repdata, "replica")
        print("Walkers and replicas have been written.")
    def write_array(self,arr,fbase):
        # Used to write the replica, walker files
        for w in range(self.nwalkers):
            # Class function to write out the information from the logfile
            with open(fbase+'-%d-data.dat'%w,'w') as f:
                f.write("#")
                for key in self.keys:
                    f.write("%s " % key)
                f.write("\n")
                for i in range(len(arr[self.keys[0]][w])):
                    f.write("%10.5f " % i)
                    for key in self.keys:
                        f.write("%10.5f " % arr[key][w][i])
                    f.write("\n")
    def get_replicadata(self):
        for r in range(self.nwalkers):
            for key in self.keys:
                self.repdata[key][r]=self.walkdata[key].T[np.where(self.walkloc==r)]
    def calc_Teff(self):
        try:
            self.lambdas=np.genfromtxt("lambdas.dat")
        except:
            exit("Error: need to generate lambdas.dat")
        self.Teff=np.zeros(self.nwalkers)
        self.Havg=np.zeros(self.nwalkers)
        for r in range(self.nwalkers):
            self.Havg[r]=np.average(self.repdata["Enthalpy"][r])
            print("For walker %d the average Enthalpy is %10.5f kcal/mol" % (r,self.Havg[r]))
            self.Teff[r] = self.lambdas[r] + self.eta*(self.Havg[r]-self.Ho)
        np.savetxt("Teff.dat",np.c_[self.Havg,self.Teff, self.lambdas])
    def hist_data(self,key,nbins):
        self.minval[key],self.maxval[key] = np.min(self.repdata[key]),np.max(self.repdata[key])
        histogram,center,edges=[],[],[]
        for r in range(self.nwalkers):
            tmp, bins = np.histogram(self.repdata[key][r],bins=nbins,range=(self.minval[key],self.maxval[key]))
            center = (bins[:-1] + bins[1:])/2
            histogram.append(tmp)
            edges.append(bins)
        self.hist[key]=np.array(histogram).T
        self.edges[key]=np.array(edges).T
        np.savetxt("histogram_"+key+".dat", np.row_stack((center,histogram)).T)

def run_stwham(hist,eta,Ho,nbins,checklimit,lambdas,minval,maxval):
    def EffTemp(lambdavalue, eta, Ho, H):
        # Evaluates the gREM effective temperature
        Teff = lambdavalue + eta*(H - Ho)
        return Teff
    # This is a function to run ST-WHAM on histogrammed data for gREM simulations
    pdf      = np.zeros((nbins, allwalkers.nwalkers))
    count    = np.zeros(allwalkers.nwalkers)
    totpdf   = np.zeros(nbins)
    totcount = 0
    binsize  = (maxval-minval)/float(nbins)
    print("ST WHAM Parameters")
    print("binsize = %10.5f" % binsize)
    print("minval  = %10.5f" % minval)
    print("maxval  = %10.5f" % maxval)
    print("nbins   = %10.5f" % nbins)
    # This part builds the pdfs from the histograms
    for rep in range(allwalkers.nwalkers):
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
    hfrac = np.zeros((nbins,allwalkers.nwalkers))
    bstart, bstop = None, None
    for i in range(nbins):
        if (totpdf[i] > checklimit):
            bstart = i + 3
            break
    if (bstart == None): exit("Error: Lower end of logarithm will be undefined")
    for i in range(nbins-1, bstart, -1):
        if (totpdf[i] > checklimit):
            bstop = i - 3
            break
    if (bstop == None): exit("Error: Upper end of logarithm will be undefined")
    # If these checks pass, we should be good to proceed.
    # Calculate Ts(H)
    for i in range(bstart,bstop):
        if (totpdf[i+1] > checklimit and totpdf[i-1]>checklimit): # makes sure it is above checklimit
            betaH[i] = np.log(totpdf[i+1]/totpdf[i-1])/(2*binsize)*consts["kb"]
        else: # otherwise 0
            betaH[i] = 0
        for rep in range(allwalkers.nwalkers):
            if (totpdf[i] > 0): #positive, not empty
                H      = minval + (i*binsize)
                weight = np.nan_to_num(1/EffTemp(lambdas[rep],eta,Ho,H))
                betaW[i] = betaW[i] + ((count[rep]*pdf[i][rep])/(totcount*totpdf[i])*weight)
        TH[i] = 1 / (betaH[i] + betaW[i])

    ### This Section calculates the histogram fraction
    with open("histfrac_STWHAM.out",'w') as fracout:
        for rep in range(allwalkers.nwalkers):
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

    with open("TandS_STWHAM.out",'w') as tout:
        for i in range(bstart,bstop):
            S[i] = Falpha(bstart,i,TH,binsize)
            tout.write("%f %f %f\n" % (minval+(i*binsize), TH[i],S[i]))
    return
            


def read_lammps_log(logbase,runindex,ndata):
    # This reads a lammps log file and extracts the data from it
    # it saves each of the keys that are written in the "Step ..." line
    # and saves the info as keys and a data array
    keys,logdata=[],{}
    with open(logbase+str(runindex),'r') as lg:
        lines=lg.readlines()
        flag=0
        for line in lines:
            if "Loop" in line: # Stops reading the file
                flag=0
                break
            if flag==1: # reads values
                vals=line.strip().split()
                c=0
                for key in keys:
                    logdata[key].append(float(vals[c]))
                    c+=1
            if "Step" in line and flag==0: # begins read, reads keys
                flag=1
                keys=line.strip().split()
                for key in keys:
                    logdata[key]=[]
    for key in keys:
        logdata[key]=np.array(logdata[key][:-1])
        if ndata != len(logdata[key]):
            rat=int(len(logdata[key])/ndata)
            logdata[key]=logdata[key][::rat]
    return keys,logdata

def read_walker(walklog,start,stop):
    walkdata=np.genfromtxt(walklog+"-"+str(start), skip_header=3, skip_footer=1,dtype=int)
    for run in range(start+1,stop):
        walkdata=np.append(walkdata,np.genfromtxt(walklog+"-"+str(run),skip_header=3,skip_footer=1,dtype=int),axis=0)
    return walkdata

def read_lambdas(lammpsfile):
    vals=[]
    with open(lammpsfile,'r') as f:
        lines=f.readlines()
        for line in lines:
            if "variable lambda world" in line:
                vals = np.array(line.strip().split()[3:],dtype='float')
    np.savetxt("lambdas.dat", np.c_[vals])
    return
if __name__ == "__main__":
    # read in the walker data (which window it is in)
    walkloc=read_walker(workdir+"/log/log.lammps",fstart,fend)[:,1:]
    # read input file
    read_lambdas(workdir+"/start.gREM")
    nreps=walkloc.shape[1]
    nruns = fend-fstart
    ndata = np.shape(walkloc)[0]/nruns
    print("There are %d replicas in the present simulation, with %d runs" % (nreps,nruns))
    allwalkers=walkers(workdir,nreps,fstart,fend,ndata,walkloc,eta,Ho)
    allwalkers.post_process(nbins)
    tol=0.001
    run_stwham(allwalkers.hist["PotEng"],eta,Ho,nbins,tol,allwalkers.lambdas,allwalkers.minval["PotEng"], allwalkers.maxval["PotEng"])


