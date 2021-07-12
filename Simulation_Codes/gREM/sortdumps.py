#!/usr/bin/env python
from walker import walkers
import numpy as np
import pickle
from collections import deque


class frame:
    def __init__(self,dumptype,f):
        if dumptype=="lmps":
            self.read_lmpsdump(f)
    def read_lmpsdump(self,f):
        lines=[]
        for i in range(9): lines.append(f.readline())
        self.T = int(lines[1].strip())
        self.N = int(lines[3].strip())
        self.xhi, self.xlo = float(lines[5].strip().split()[1]),float(lines[5].strip().split()[0])
        self.yhi, self.ylo = float(lines[6].strip().split()[1]),float(lines[6].strip().split()[0])
        self.zhi, self.zlo = float(lines[7].strip().split()[1]),float(lines[7].strip().split()[0])
        lx,ly,lz = self.xhi-self.xlo,self.yhi-self.ylo,self.zhi-self.zlo
        self.L=[lx,ly,lz]
        self.fields = lines[8].strip().split()[2:]
        self.data={}
        intfields=["id","mol","type","ix","iy","iz"]
        fltfields=["x","xu","y","yu","z","zu"]
        for field in self.fields:
            self.data[field]=[]
        for i in range(self.N):
            line=f.readline().strip().split()
            for fno in range(len(self.fields)):
                field=self.fields[fno]
                if field in intfields:
                    self.data[field].append(int(line[fno]))
                if field in fltfields:
                    self.data[field].append(float(line[fno]))
        # Sorts by id
        for field in self.fields[1:]:
            _,self.data[field] = zip(*sorted(zip(self.data[self.fields[0]],self.data[field])))
            self.data[field]=np.array(self.data[field])
        self.data[self.fields[0]]=np.array(sorted(self.data[self.fields[0]]))
        return
    def write_lmpsdump(self,f):
        f.write("ITEM: TIMESTEP\n")
        f.write("%d\n"%self.T)
        f.write("ITEM: NUMBER OF ATOMS\n")
        f.write("%d\n"%self.N)
        f.write("ITEM: BOX BOUNDS pp pp pp\n")
        f.write("% 12.8f % 12.8f\n" % (self.xlo,self.xhi))
        f.write("% 12.8f % 12.8f\n" % (self.ylo,self.yhi))
        f.write("% 12.8f % 12.8f\n" % (self.zlo,self.zhi))
        f.write("ITEM: ATOMS ")
        for field in self.fields:
            f.write("%s "%field)
        f.write("\n")
        intfields=["id","mol","type","ix","iy","iz"]
        fltfields=["x","xu","y","yu","z","zu"]
        for atom in range(self.N):
            for field in self.fields:
                if field in intfields: f.write("%5d " % self.data[field][atom])
                if field in fltfields: f.write("% 12.8f " % self.data[field][atom])
            f.write("\n")
    def add_leaflets(self,leaflets):
        self.leaflets=leaflets
        assert np.shape(self.leaflets)[0] == 2, "There are %d leaflets" % np.shape(self.leaflets)[0]
        A,B=[],[]
        self.data["leaf"]=[]
        for atom in self.leaflets[0]: A.append(self.data["mol"][atom-1])
        for atom in self.leaflets[1]: B.append(self.data["mol"][atom-1])
        for atom in self.data["id"]:
            if self.data["mol"][atom-1] in A:
                self.data["leaf"].append("A")
            elif self.data["mol"][atom-1] in B:
                self.data["leaf"].append("B")
            else:
                exit("Error: Molecule not in A or B")
        return
    def add_data(self, key, value):
        try:
            self.calc_data
        except:
            self.calc_data={}
        self.calc_data[key]=value
        return
    def calc_thick(self, hatoms):
        A,B=[],[]
        for atom in self.data["id"]:
            if self.data["type"][atom-1] in hatoms:
                if self.data["leaf"][atom-1]=="A":
                    A.append(72*self.data["z"][atom-1])
                elif self.data["leaf"][atom-1]=="B":
                    B.append(72*self.data["z"][atom-1])
        AcomZ = np.sum(A)/(len(A)*72)
        BcomZ = np.sum(B)/(len(B)*72)
        dcomZ = AcomZ - BcomZ
        dcomZ = dcomZ - self.L[2]*np.round(dcomZ/self.L[2])
        self.add_data("thickness",np.abs(dcomZ))
        return





def find_leaflet(frame,hatom):
    def PBC(frame, curid, newid):
        # Note the -1's account for difference in counting between lammps + python
        dx=frame.data["x"][curid-1]-frame.data["x"][newid-1]
        dy=frame.data["y"][curid-1]-frame.data["y"][newid-1]
        dz=frame.data["z"][curid-1]-frame.data["z"][newid-1]
        dx = dx - frame.L[0]*np.round(dx/frame.L[0])
        dy = dy - frame.L[1]*np.round(dy/frame.L[1])
        dz = dz - frame.L[2]*np.round(dz/frame.L[2])
        return np.sqrt(dx**2. + dy**2. + dz**2.)
    def calc_neigh(frame,currentid,hmask,rcut):
        partners=[]
        for atom in frame.data["id"][hmask]:
            dr = PBC(frame, currentid, atom)
            if dr < rcut and currentid!=atom: partners.append(atom)
        return partners
    def find_connections(connectivity):
        seen = set()
        for root in connectivity:
            if root not in seen:
                seen.add(root)
                component=[]
                queue = deque([root])
                while queue:
                    node = queue.popleft()
                    component.append(node)
                    for neighbor in connectivity[node]:
                        if neighbor not in seen:
                            seen.add(neighbor)
                            queue.append(neighbor)
                yield component
    connectivity={}
    hmask = frame.data["type"]==hatom
    for atom in frame.data["id"][hmask]:
        connectivity[atom]=calc_neigh(frame,atom,hmask,11)
    leaflets = list(find_connections(connectivity))
    frame.add_leaflets(leaflets)
    return frame,leaflets
    
    


def get_frames(workdir,dumpbase, start, stop, nreps,walkloc):
    wframes, rframes = [],[]
    print("Reading trajectory frames")
    for walker in range(nreps):
        print("Reading walker %d of %d" % (walker, nreps))
        tmp = []
        for fl in range(start,stop):
            with open(workdir+"/"+str(walker)+"/"+dumpbase+"-"+str(fl)+".dat",'r') as f:
                while True:
                    fr=None
                    try:
                        fr = frame("lmps",f)
                    except:
                        break
                    tmp.append(fr)
                tmp.pop() # deletes final frame which is same as first frame of next dump
        wframes.append(tmp)
    print("Read complete")
    print("Beginning dump sorting")
    rat=int(np.shape(walkloc)[0]/np.shape(wframes)[1])
    wframes = np.array(wframes)
    for r in range(nreps):
        tmp=wframes.T[np.where(walkloc[::rat]==r)]
        pickle.dump(tmp,open(dumpbase+"_"+str(r)+".pckl",'wb'))
    return 
            
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-fname', default="dump", type=str, help='dump file name')
    parser.add_argument('-ext', default=".dat", type=str, help='dump file extension')
    args = parser.parse_args()
    fname    = args.fname
    ext      = args.ext
    frames=[]
    with open(fname+ext,'r') as f:
        while True:
            try:
                fr = frame("lmps",f)
            except:
                break
            frames.append(fr)
        print("Read COMPLETE for %s%s, found %d frames" % (fname,ext,len(frames)))
        print("Deleting final frame")
        frames.pop()
        print(len(frames))
    pickle.dump(frames,open(fname+".pckl", 'wb'))
