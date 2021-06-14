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

# Read Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-start', default=0,type=int, help='Starting index for calculating [default=0]')
parser.add_argument('-stop', default=5, type=int, help='Ending index for calculating [default=5]')
parser.add_argument('-workdir', default=".", type=str, help='Working directory location [default=./]')
args= parser.parse_args()

fstart = args.start
fend   = args.stop
workdir= args.workdir

class walkers:
    def __init__(self,workdir,nwalkers,start,stop,ndata,walkloc):
        self.walkdata={}
        self.repdata={}
        self.nwalkers=nwalkers
        self.ndata=ndata
        self.start,self.stop=start,stop
        self.walkloc=walkloc
        self.logbase="/log.lammps.w-"
        self.workdir=workdir
        print("Walker class initiated based on following log format %s" % self.logbase)
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
        self.write_array(self.walkdata, "walkers")
        self.write_array(self.repdata, "replica")
    def write_array(self,arr,fbase):
        for w in range(self.nwalkers):
            # Class function to write out the information from the logfile
            with open(fbase+'-%d-data.dat'%w,'w') as f:
                f.write("#")
                for key in self.keys:
                    f.write("%s " % key)
                f.write("\n")
                count=0
                for i in range(len(arr[self.keys[0]][w])):
                    f.write("%10.5f " % i)
                    for key in self.keys:
                        f.write("%10.5f " % arr[key][w][i])
                    f.write("\n")
    def get_replicadata(self):
        for r in range(self.nwalkers):
            for key in self.keys:
                self.repdata[key][r]=self.walkdata[key].T[np.where(self.walkloc==r)]


        

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

if __name__ == "__main__":
    # read in the walker data (which window it is in)
    walkloc=read_walker(workdir+"/log/log.lammps",fstart,fend)[:,1:]
    nreps=walkloc.shape[1]
    nruns = fend-fstart
    ndata = np.shape(walkloc)[0]/nruns
    print("There are %d replicas in the present simulation, with %d runs" % (nreps,nruns))
    allwalkers=walkers(workdir,nreps,fstart,fend,ndata,walkloc)
    allwalkers.walk_data()
    allwalkers.get_replicadata()
    allwalkers.write()




