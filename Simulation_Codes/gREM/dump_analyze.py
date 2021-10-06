#!/usr/bin/env python
import argparse,pickle,os
import numpy as np
from scipy import stats
from sortdumps import frame, get_frames,find_leaflet,calc_P2

parser = argparse.ArgumentParser()
parser.add_argument('-dumpid', default=0, type=int, help='ID of dump to analyze (walker-id) [default=0]')
parser.add_argument('-dumpbase', default="None", type=str, help='Base name of the dump file i.e. dump in dump-1.dat. [default=None]')
parser.add_argument('-safe', default=1, type=int, help='Determines whether to identify leaflets every step or not [default=1 : every]')
parser.add_argument('-nb', default=5, type=int, help='Number of blocks [default=5]')
parser.add_argument('-start', default=1, type=int, help='First dump to read [default=1]')
parser.add_argument('-stop', default=11, type=int, help='Last dump to read [default=11]')
parser.add_argument('-rcut', default=12, type=float, help='Cutoff distance in A [default=12]')
parser.add_argument('-hatom',default=4, type=int, help='Atom type for header atom')
parser.add_argument('-workdir',default=".", type=str, help='Workdir [default=.]')
parser.add_argument('-nreps', default=40, type=int, help='Number of replicas [default=40]')
parser.add_argument('-option',default=1, type=int, help='Options: [1] generate properties for single dump [2] join and sort [3] calculate stats')

args = parser.parse_args()

safeleaf    = args.safe
dumpbase    = args.dumpbase
dumpid      = args.dumpid
nblocks     = args.nb
rcut        = args.rcut
hatom       = args.hatom
start       = args.start
stop        = args.stop
workdir     = args.workdir
nreps       = args.nreps
op          = args.option

t_val = stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)

def gen_leaf(frames):
    # This function calls the required function from sort dumps to find the leaflets for each frame
    if safeleaf == 1:
        for fr in range(len(frames)):
            if fr%1000 == 0:
                print("Finding leaflet %d" % fr)
            frames[fr],_=find_leaflet(frames[fr],hatom,rcut)
    else: # bases it off first frame
        frames[0],leaflets=find_leaflet(frames[0],hatom,rcut)
        for fr in range(1,len(frames)):
            frames[fr].add_leaflets(leaflets)

def calc_stats(frames,key,repid):
    # This takes a frame calc_data field (i.e. thickness) and then calculates the average (with uncertainty)
    data=[]
    for fr in frames:
        data.append(fr.calc_data[key])
    data = np.array(data)
    hist,bins = np.histogram(data,bins=100)
    center = (bins[:-1] + bins[1:]) / 2
    np.savetxt("hist/histogram_"+key+"_"+str(repid)+".out",np.c_[center,hist])
    bl_data = np.split(data,nblocks)
    bl_av = np.average(bl_data,axis=1)
    try:
        os.mkdirs("temp")
    except:
        pass
    bl_err = np.std(bl_av,axis=0)*t_val
    f=open("temp/"+key+"_dump_"+str(repid)+".out",'w')
    f.write("%d %10.5f +/- %10.5f\n" %(dumpid,np.average(data),bl_err))
    f.close()
    for block in range(len(bl_av)):
        g=open("temp/bl_"+str(block)+"_"+key+"_dump_"+str(repid)+".out",'w')
        g.write("%d %10.5f\n" %(dumpid,bl_av[block]))
        g.close()
    return

if __name__ == "__main__":
    if dumpbase == "None":
        exit("Error: please specify dumpbase")
    if op == 1:
        print("Reading walker %d" % dumpid)
        walkframes = []
        for fl in range(start,stop):
            print("fl",fl)
            with open(workdir+"/"+str(dumpid)+"/"+dumpbase+"-"+str(fl)+".dat",'r') as f:
                c=0
                leaflets=None
                while True:
                    fr=None
                    try: # read frame, calculate properties, clear xyz data
                        c = c + 1
                        print(c)
                        fr = frame("lmps",f)
                    except:
                        break
                    # Analyze
                    if safeleaf == 1:
                            fr,_ = find_leaflet(fr,hatom,rcut)
                    elif safeleaf == 0 and c == 1:
                        fr, leaflets = find_leaflet(fr,hatom,rcut)
                    elif safeleaf == 0 and c != 1:
                        fr.add_leaflets(leaflets)
                    else:
                        print("Error invalid option for safeleaf")
                    fr.calc_thick([3,4])
                    fr.calc_area()
                    calc_P2(fr, [[5,6],[9,10]])
                    fr.clear_data()
                    #add to tmp
                    walkframes.append(fr)
                walkframes.pop() # deletes final frame which is same as first frame of next dump
        pickle.dump(walkframes,open("walker-"+str(dumpid)+".pckl",'wb'))
    if op == 2:
        wframes = []
        walkloc=pickle.load(open(workdir+'/walkloc.pckl','rb'),encoding='latin1')
        allwalkers=pickle.load(open(workdir+'/allwalkers.pckl','rb'),encoding='latin1')
        for walker in range(nreps):
            wframes.append(pickle.load(open("walker-"+str(walker)+".pckl",'rb')))
        wframes = np.array(wframes)
        print(np.shape(walkloc))
        print(walkloc[0])
        rat = int(np.shape(walkloc)[0]/len(wframes[0]))
        for r in range(nreps):
            tmp=wframes.T[np.where(walkloc[::rat]==r)[1]].T[r]
            pickle.dump(tmp,open("replica-"+dumpbase+"_"+str(r)+".pckl",'wb'))
    if op == 3:
        for rep in range(nreps):
            # Reads frames
            frames=pickle.load(open("replica-"+dumpbase+"_"+str(rep)+".pckl",'rb'),encoding='latin1')
            # Finds Leaflets
            calc_stats(frames,"thickness",rep)
            calc_stats(frames, "area",rep)
            calc_stats(frames, "P2",rep)


            
        


    

