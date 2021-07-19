#!/usr/bin/env python
from walker import walkers
import numpy as np
import pickle
from collections import deque


class frame:
    """
    This is a pretty large class that contains a lot of info about each frame.

    The nice thing about this is that it allows for iteration on the frames, and lets you work with them
    as objects with class functions if necessary. 

    Class Objects:
    self.T = timestep of frame [int]
    self.N = Number of Atoms [int]
    self._hi, self._lo = {x|y|z} boundary values [float]
    self.L = Dimension of box [array shape (3,) of floats]
    self.fields = Fields in data file [array of strings]
    self.data[field] = array of per-atom values output from data file [float or int]
    self.calc_data[calc_item] = Calculated average of a calculated value [float]
    """
    def __init__(self,dumptype,f):
        # Intializes frame by reading data
        if dumptype=="lmps":
            self.read_lmpsdump(f)
    def read_lmpsdump(self,f):
        # Reads in the frame (assumes traditional lammps custom dump file)
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
        # This next section assigns type based off of common choices
        # Currently no handling of unknown types.
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
        # Sorts by atom ID using the zip function (sorts by id in first column)
        for field in self.fields[1:]:
            _,self.data[field] = zip(*sorted(zip(self.data[self.fields[0]],self.data[field])))
            self.data[field]=np.array(self.data[field])
        self.data[self.fields[0]]=np.array(sorted(self.data[self.fields[0]]))
        return
    def write_lmpsdump(self,f):
        # Writes out the frame in a dump format to file "f"
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
        # Adds leaflet information to the frame based off of calculations
        # There are assumed to be exactly two leaflets.
        self.leaflets=leaflets
        if np.shape(self.leaflets)[0] != 2:
            print("Error: Incorrect number of leaflets")
            print("Leaflet shape: ",np.shape(self.leaflets))
            print(self.leaflets)
            exit()
        A,B=[],[]
        self.data["leaf"]=[]
        # Saves atoms in the two leaflets
        for atom in self.leaflets[0]: A.append(self.data["mol"][atom-1])
        for atom in self.leaflets[1]: B.append(self.data["mol"][atom-1])
        for atom in self.data["id"]:
            if self.data["mol"][atom-1] in A:
                self.data["leaf"].append("A")
            elif self.data["mol"][atom-1] in B:
                self.data["leaf"].append("B")
            else:
                exit("Error: Molecule not in A or B") # should never happen
        return
    def add_data(self, key, value):
        # Adds data of type key to a frame object
        # If this is first time calling object, makes a new dictionary
        try:
            self.calc_data
        except:
            self.calc_data={}
        self.calc_data[key]=value
        return
    def add_mass(self,M):
        #Adds mass to be assumed for all atoms
        self.M = M
        return
    def calc_area(self):
        # Calculates bilayer area
        self.add_data("area",np.abs(self.L[0]*self.L[1]))
        return
    def calc_thick(self, hatoms):
        # Calculates bilayer thickness
        A,B=[],[]
        # sort distances by leaflet A and B
        try:
            self.M
        except:
            self.add_mass(72)
        for atom in self.data["id"]:
            if self.data["type"][atom-1] in hatoms:
                if self.data["leaf"][atom-1]=="A":
                    A.append(72*self.data["z"][atom-1])
                elif self.data["leaf"][atom-1]=="B":
                    B.append(72*self.data["z"][atom-1])
        # Calc z component of COM for both leaflets
        AcomZ = np.sum(A)/(len(A)*self.M)
        BcomZ = np.sum(B)/(len(B)*self.M)
        dcomZ = AcomZ - BcomZ
        # Apply PBC
        dcomZ = dcomZ - self.L[2]*np.round(dcomZ/self.L[2])
        self.add_data("thickness",np.abs(dcomZ))
        return
    def calc_p2(self,hatoms,tatoms):
        try:
            self.M
        except:
            self.add_mass(72)
        def COM(self,atoms):
            comx,comy,comz=0,0,0
            denom = 0
            for atom in atoms:
                comx = comx + self.data["x"][atom]*self.M
                comy = comy + self.data["y"][atom]*self.M
                comz = comz + self.data["z"][atom]*self.M
                denom = denom + self.M
            rcom = np.array([comx/denom, comy/denom, comz/denom])
            return rcom
        def unitvec(r,x,y,z):
            drx = r[0]-x
            dry = r[1]-y
            drz = r[2]-z
            print(drx,dry,drz)
            dr = np.sqrt(drx**2. + dry**2. + drz**2)
            return np.array([drx/dr, dry/dr, drz/dr])
        def c2(vec):
            z = np.zeros(3)
            z[2]=1
            costhta = np.dot(vec,z)
            return 0.5*(3*costhta**2. - 1)

        # This calculates P2=0.5*(3 cos^2 theta -1)
        # where theta = angle between bilayer normal (z axis) and the vector pointing from 
        # second c2 bead to the FIRST headgroup atom
        
        # Calculate atoms per mol and num lipids
        apermol = np.sum((self.data["mol"]==1)*1)
        nlipids = int(self.N/apermol)
        
        header,tail = {},{}
        # Initialize
        for lip in range(nlipids): header[lip], tail[lip]=[],[]
        # Find header atoms, tail atoms
        for atom in range(len(self.data["type"])):
            a_id = atom + 1
            if self.data["type"][atom] in hatoms:
                header[self.data["mol"][atom]-1].append(atom)
            if (a_id-tatoms[0])%apermol==0:
                tail[self.data["mol"][atom]-1].append(atom)
        hcom, tr = [],[]
        for key in header:
            hcom.append(COM(self,header[key]))
        vecs=[]
        # Find all vectors from lipid tail to head group, and note which leaflet is used
        for lip in range(nlipids):
            for atom in tail[lip]:
                vecs.append(unitvec(hcom[lip],self.data["x"][atom],self.data["y"][atom],self.data["z"][atom]))
        c2vals=[]
        # Loop over vectors and calculate P2
        for vec in vecs:
            c2vals.append(c2(vec))
        avc2 = np.average(c2vals)
        self.add_data("c2",avc2)
        return

def find_leaflet(frame,hatom,rcut):
    """
    Probably one of the more complicated functions in this code - this calculates which leaflet a lipid belongs to (A or B)
    """
    def PBC(frame, curid, newid):
        # Note the -1's account for difference in counting between lammps + python
        # Distance calculation
        dx=frame.data["x"][curid-1]-frame.data["x"][newid-1]
        dy=frame.data["y"][curid-1]-frame.data["y"][newid-1]
        dz=frame.data["z"][curid-1]-frame.data["z"][newid-1]
        # Apply PBC correction
        dx = dx - frame.L[0]*np.round(dx/frame.L[0])
        dy = dy - frame.L[1]*np.round(dy/frame.L[1])
        dz = dz - frame.L[2]*np.round(dz/frame.L[2])
        # Return distance
        return np.sqrt(dx**2. + dy**2. + dz**2.)
    def calc_neigh(frame,currentid,hmask,rcut):
        # Find neighbors
        partners=[]
        # Loop over header atoms ids, calculate PBC based off of that
        for atom in frame.data["id"][hmask]:
            dr = PBC(frame, currentid, atom) # gets min image distance
            if dr < rcut and currentid!=atom: partners.append(atom) # if image distance less than cutoff, atom connected
        return partners
    def find_connections(connectivity):
        # This function fins connections using graph theory
        # connectivity is a dictionary of nearest neighbors within cutoff
        seen = set()
        # Loop over dictionary keys (aka root atom ids)
        for root in connectivity:
            # If root has not already been assigned to a leaflet
            if root not in seen:
                #Add root as having been seen
                seen.add(root)
                component=[]
                # Faster type of list
                queue = deque([root])
                # Continues until deque is depleted
                while queue:
                    node = queue.popleft() # removes left element and returns it as node
                    component.append(node) # save node as part of component list
                    # Loop over neighbors of node
                    for neighbor in connectivity[node]:
                        # Checks if neighbor of node have been seen yet
                        if neighbor not in seen:
                            # If not: add neighbor and add to queue
                            seen.add(neighbor)
                            queue.append(neighbor)
                yield component
    connectivity={}
    hmask = frame.data["type"]==hatom
    # For every atom, find connectivity within distance rcut away
    for atom in frame.data["id"][hmask]:
        connectivity[atom]=calc_neigh(frame,atom,hmask,rcut)
    # Sorts into leaflets
    leaflets = list(find_connections(connectivity))
    # Adds to frame
    frame.add_leaflets(leaflets)
    return frame,leaflets
    
    


def get_frames(workdir,dumpbase, start, stop, nreps,walkloc):
    # Function that reads the frames for each walker.
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
    # Double checks that only frames that match dumps in the overall lammps log are included
    rat=int(np.shape(walkloc)[0]/np.shape(wframes)[1])
    # Writes a pickle file that stores the frame object for post-processing.
    wframes = np.array(wframes)
    for r in range(nreps):
        tmp=wframes.T[np.where(walkloc[::rat]==r)]
        pickle.dump(tmp,open(dumpbase+"_"+str(r)+".pckl",'wb'))
    return 
            
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-fname', default="dump", type=str, help='dump file name')
    parser.add_argument('-nbins', default=12, type=int, help='Number of bins')
    parser.add_argument('-opt', default=1, type=int, help='Choose Option')
    args = parser.parse_args()
    fname    = args.fname
    nbins    = args.nbins
    opt      = args.opt
    if opt == 1:
        frames=[]
        with open(fname+".dat",'r') as f:
            while True:
                try:
                    fr = frame("lmps",f)
                except:
                    break
                frames.append(fr)
            print("Read COMPLETE for %s%s, found %d frames" % (fname,".dat",len(frames)))
            print("Deleting final frame")
            frames.pop()
            print(len(frames))
        pickle.dump(frames,open(fname+".pckl", 'wb'))
    else:
        frames=[]
        for i in range(nbins):
            fr = pickle.load(open(fname + "_" +str(i)+".pckl",'rb'))
            frames.append(fr[-1])
            frames[i].T=i
        with open("dump_singlesnapshot.lmpstraj",'w') as f:
            for i in range(nbins):
                frames[i].write_lmpsdump(f)
        print("Dumps have been written")

