#!/usr/bin/env python
import argparse,pickle
import numpy as np
from scipy import stats
from sortdumps import frame, get_frames,find_leaflet

parser = argparse.ArgumentParser()
parser.add_argument('-dumpid', default=0, type=int, help='ID of dump to analyze [default=0]')
parser.add_argument('-dumpbase', default="None", type=str, help='Base name of the dump file i.e. dump in dump-1.dat. [default=None]')
parser.add_argument('-safe', default=1, type=int, help='Determines whether to identify leaflets every step or not [default=1 : every]')
parser.add_argument('-nb', default=5, type=int, help='Number of blocks [default=5]')
args = parser.parse_args()

safeleaf    = args.safe
dumpbase    = args.dumpbase
dumpid      = args.dumpid
nblocks     = args.nb

t_val = stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)

def gen_leaf(frames):
    if safeleaf == 1:
        for fr in range(len(frames)):
            if fr%1000 == 0:
                print("Finding leaflet %d" % fr)
            frames[fr],_=find_leaflet(frames[fr],4)
    else: # bases it off first frame
        frames[0],leaflets=find_leaflet(frames[0],4)
        for fr in range(1,len(frames)):
            frames[fr].add_leaflets(leaflets)

def calc_stats(frames,key):
    data=[]
    for fr in frames:
        data.append(fr.calc_data[key])
    data = np.array(data)
    bl_data = np.split(data,nblocks)
    bl_av = np.average(bl_data,axis=1)
    bl_err = np.std(bl_av,axis=0)*t_val
    f=open(key+"_dump_"+str(dumpid)+".out",'w')
    f.write("%d %10.5f +/- %10.5f\n" %(dumpid,np.average(data),bl_err))
    f.close()

if __name__ == "__main__":
    if dumpbase != "None":
        frames=pickle.load(open(dumpbase+"_"+str(dumpid)+".pckl",'rb'),encoding='latin1')
        print("There are %d frames" % len(frames))
        gen_leaf(frames)
        for fr in range(len(frames)):
            frames[fr].calc_thick([3,4])
        calc_stats(frames,"thickness")

