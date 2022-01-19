#!/usr/bin/env python
import numpy as np
from scipy import stats
import pickle

def read_lammps_log(logbase,runindex,ndata):
    """
    # This reads a lammps log file and extracts the data from it
    # it saves each of the keys that are written in the "Step ..." line
    # and saves the info as keys and a data array
    """
    keys,logdata=[],{}
    with open(logbase+str(runindex),'r') as lg:
        lines=lg.readlines()
        flag=0
        lprev=None
        for line in lines:
            if "Loop" in line: # Stops reading the file
                flag=0
                break
            if "WARNING" in line:
                print("*********************")
                print("Note: There was a warning found")
                print("warning location: %s" % (logbase+str(runindex)))
                print("%s" % line.strip())
                if lprev is not None:
                    print("%s" % lprev.strip())
                print("*********************")
            if flag==1 and "WARNING" not in line: # reads values
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
            lprev = line
    for key in keys:
        logdata[key]=np.array(logdata[key][:-1])
        if ndata != len(logdata[key]):
            rat=int(len(logdata[key])/ndata)
            logdata[key]=logdata[key][::rat]
    return keys,logdata

class walkers:
    def __init__(self,workdir,start,stop,walkloc,eta,Ho):
        self.walkdata,self.repdata={},{}
        self.start,self.stop=start,stop
        self.walkloc=walkloc
        self.nwalkers=np.shape(self.walkloc)[1]
        self.nruns=self.stop-self.start
        self.ndata=np.shape(self.walkloc)[0]/self.nruns
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
    def post_process(self,nbins,nblocks):
        print("Beginning Post Processing")
        self.calc_Teff(nblocks)
        for key in self.keys:
            self.hist_data(key,nbins,nblocks)
        print("Post Processing Completed")
    def walk_data(self):
        # Read the keys and first vals
        self.nruns=self.stop-self.start
        self.keys, firstvals = read_lammps_log(self.workdir+"/"+str(0)+self.logbase, self.start,self.ndata)
        perwindow=len(firstvals[self.keys[0]])
        for key in self.keys:
            self.walkdata[key]=np.zeros((self.nwalkers,self.nruns*perwindow))
            self.repdata[key]=np.zeros((self.nwalkers,self.nruns*perwindow))
        for w in range(self.nwalkers):
            for r in range(self.start,self.stop):
                rstart, rstop = (r-self.start)*perwindow, (r-self.start+1)*perwindow
                # read in the run
                _, values = read_lammps_log(self.workdir+"/"+str(w)+self.logbase, r,self.ndata)
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
        print(np.shape(self.walkdata["PotEng"]))
        for r in range(self.nwalkers):
            for key in self.keys:
                self.repdata[key][r]=self.walkdata[key].T[np.where(self.walkloc==r)]
    def calc_Teff(self,nblocks):
        def block_av(data,nblocks):
            bl_data=np.array_split(data,nblocks)
            bl_av=np.average(bl_data,axis=1)
            return bl_av
        try:
            self.lambdas=np.genfromtxt("lambdas.dat")
        except:
            exit("Error: need to generate lambdas.dat")
        self.Teff=np.zeros(self.nwalkers)
        self.Havg=np.zeros(self.nwalkers)
        bl_Hav = []
        for r in range(self.nwalkers):
            self.Havg[r]=np.average(self.repdata["PotEng"][r])
            bl_Hav.append(block_av(self.repdata["PotEng"][r],nblocks))
            print("For walker %d the average Enthalpy is %10.5f kcal/mol" % (r,self.Havg[r]))
            self.Teff[r] = self.lambdas[r] + self.eta*(self.Havg[r]-self.Ho)
        print(np.shape(self.Havg), np.shape(self.Teff), np.shape(self.lambdas))
        np.savetxt("Teff.dat",np.c_[self.Havg,self.Teff, self.lambdas])
        pickle.dump(bl_Hav,open("bl_Teff.pckl",'wb'),protocol=pickle.HIGHEST_PROTOCOL)
    def hist_data(self,key,nbins,nblocks):
        try:
            self.bl_hist
        except:
            self.bl_hist={}
        # set up blocking
        length = len(self.repdata[key][0])
        nperblock = length/nblocks
        print("nper", nperblock,length)
        nperblock=int(nperblock)
        t_value = stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)
        self.minval[key],self.maxval[key] = np.min(self.repdata[key]),np.max(self.repdata[key])
        histogram,center,edges=[],[],[]
        for r in range(self.nwalkers):
            tmp, bins = np.histogram(self.repdata[key][r],bins=nbins,range=(self.minval[key],self.maxval[key]))
            center = (bins[:-1] + bins[1:])/2
            histogram.append(tmp)
            edges.append(bins)
        self.hist[key]=np.array(histogram).T
        self.edges[key]=np.array(edges).T
        self.bl_hist[key] = []
        for blk in range(nblocks):
            bl_start, bl_end = blk*nperblock, (blk+1)*nperblock
            bl_tmp = []
            for r in range(self.nwalkers):
                tmp,bins = np.histogram(self.repdata[key][r][bl_start:bl_end],bins=nbins, range=(self.minval[key],self.maxval[key]))
                bl_tmp.append(tmp)
            self.bl_hist[key].append(np.array(bl_tmp).T)
        bl_err = np.std(self.bl_hist[key],axis=0)*t_value
        np.savetxt("histogram_"+key+".dat", np.row_stack((center,histogram,bl_err.T)).T)
