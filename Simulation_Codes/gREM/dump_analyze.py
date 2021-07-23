#!/usr/bin/env python
import argparse,pickle,os
import numpy as np
from scipy import stats
from sortdumps import frame, get_frames,find_leaflet,calc_P2

parser = argparse.ArgumentParser()
parser.add_argument('-dumpid', default=0, type=int, help='ID of dump to analyze [default=0]')
parser.add_argument('-dumpbase', default="None", type=str, help='Base name of the dump file i.e. dump in dump-1.dat. [default=None]')
parser.add_argument('-safe', default=1, type=int, help='Determines whether to identify leaflets every step or not [default=1 : every]')
parser.add_argument('-nb', default=5, type=int, help='Number of blocks [default=5]')
parser.add_argument('-rcut', default=12, type=float, help='Cutoff distance in A [default=12]')
parser.add_argument('-hatom',default=4, type=int, help='Atom type for header atom')
args = parser.parse_args()

safeleaf    = args.safe
dumpbase    = args.dumpbase
dumpid      = args.dumpid
nblocks     = args.nb
rcut        = args.rcut
hatom       = args.hatom

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

def calc_stats(frames,key):
    # This takes a frame calc_data field (i.e. thickness) and then calculates the average (with uncertainty)
    data=[]
    for fr in frames:
        data.append(fr.calc_data[key])
    data = np.array(data)
    hist,bins = np.histogram(data,bins=100)
    center = (bins[:-1] + bins[1:]) / 2
    np.savetxt("hist/histogram_"+key+"_"+str(dumpid)+".out",np.c_[center,hist])
    print(np.average(data),np.std(data))
    bl_data = np.split(data,nblocks)
    bl_av = np.average(bl_data,axis=1)
    try:
        os.mkdirs("temp")
    except:
        pass
    bl_err = np.std(bl_av,axis=0)*t_val
    f=open("temp/"+key+"_dump_"+str(dumpid)+".out",'w')
    f.write("%d %10.5f +/- %10.5f\n" %(dumpid,np.average(data),bl_err))
    f.close()
    for block in range(len(bl_av)):
        g=open("temp/bl_"+str(block)+"_"+key+"_dump_"+str(dumpid)+".out",'w')
        g.write("%d %10.5f\n" %(dumpid,bl_av[block]))
        g.close()
    return

if __name__ == "__main__":
    # Basic example analysis
    if dumpbase != "None":
        # Reads frames
        frames=pickle.load(open(dumpbase+"_"+str(dumpid)+".pckl",'rb'),encoding='latin1')
        print("There are %d frames" % len(frames))
        # Finds Leaflets
        gen_leaf(frames)
        # Calculates thickness and area
        c2values=[]
        for fr in range(len(frames)):
            frames[fr].calc_thick([3,4])
            frames[fr].calc_area()
            calc_P2(frames[fr],[[5,6],[9,10]])
        calc_stats(frames,"thickness")
        calc_stats(frames, "area")
        calc_stats(frames, "P2")

