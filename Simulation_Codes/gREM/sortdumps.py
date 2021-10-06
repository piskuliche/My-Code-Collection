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
    def clear_data(self):
        for key in self.data:
            self.data[key]=[]
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
    
def unit_vec_with_pbc(r1, r2, L):
    # Simple function that takes two numpy arrays for coordinates and a list of box lengths (assumed to be lx,ly,lz)
    # and then it returns a unit vector pointing from r2 to r1 as a numpy array.
    dr_vec = np.zeros(3)
    for i in range(3): 
        dr_vec[i] = r1[i]-r2[i]
        dr_vec[i] = dr_vec[i] - L[i] * np.round(dr_vec[i]/L[i])
    drsq = np.sqrt(np.sum(np.dot(dr_vec,dr_vec)))
    return dr_vec/drsq

def P2(x):
    return 0.5*(3*np.power(x,2)-1)

def find_vecs(frame,tail_pairs):
        # This function takes the frame class object and a 2d array (shape (X,2) where X == number of pairs you choose)
        # and it pulls coordinates and then calculates unit position vectors.
        atoms_per_lipid_mol = np.sum((frame.data["mol"]==1)*1)
        num_lipids = int(frame.N/atoms_per_lipid_mol)

        vecs = []
        for lipid in range(num_lipids):
            lipid_start_id = lipid * atoms_per_lipid_mol + 1
            for pair in tail_pairs:
                # Note: -1 here accounts for the fact that both count from 1.
                # This allows for non-pythonic input.
                atom1_id = pair[0] - 1 + lipid_start_id
                atom2_id = pair[1] - 1 + lipid_start_id
                atom1_r, atom2_r = [],[]
                cart = ["x","y","z"]
                for i in range(3):
                    atom1_r.append(frame.data[cart[i]][atom1_id-1])
                    atom2_r.append(frame.data[cart[i]][atom2_id-1])
                atom1_r, atom2_r = np.array(atom1_r), np.array(atom2_r)
                vecs.append(unit_vec_with_pbc(atom1_r, atom2_r, frame.L))
        return vecs


def calc_P2(frame, tail_pairs):
    # Fnd unit vectors
    bond_vectors = find_vecs(frame, tail_pairs)
    # Find P2 Values
    z = np.array([0,0,1])
    P2_values = []
    for vec in bond_vectors:
        P2_values.append(P2(np.dot(vec,z)))
    P2_av = np.average(P2_values)
    frame.add_data("P2",P2_av)
    return P2_values



                     





    


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
    def grab_leaf(frame,hatom,rcut):
        connectivity={}
        hmask = frame.data["type"]==hatom
        # For every atom, find connectivity within distance rcut away
        for atom in frame.data["id"][hmask]:
            connectivity[atom]=calc_neigh(frame,atom,hmask,rcut)
        # Sorts into leaflets
        leaflets = list(find_connections(connectivity))
        return leaflets
    leaflets=[]
    c=0
    while len(leaflets) != 2:
        if c < 10:
            leaflets=grab_leaf(frame,hatom,rcut)
            if len(leaflets) != 2:
                rcut = rcut + 0.5
                print("redefining leaflet cutoff to %10.5f" % rcut)
        else:
            exit("Error: too many leaf iterations")



    # Adds to frame
    frame.add_leaflets(leaflets)
    return frame,leaflets
    
    

# Potentially want to make this better so that not everything is stored in memory. 
def get_frames(workdir,dumpbase, start, stop, nreps,walkloc,hatom,rcut):
    # Function that reads the frames for each walker.
    wframes, rframes = [],[]
    print("Reading trajectory frames")
    for walker in range(nreps):
        print("Reading walker %d of %d" % (walker, nreps))
        tmp = []
        for fl in range(start,stop):
            print("fl",fl)
            with open(workdir+"/"+str(walker)+"/"+dumpbase+"-"+str(fl)+".dat",'r') as f:
                while True:
                    fr=None
                    try: # read frame, calculate properties, clear xyz data
                        fr = frame("lmps",f)
                        fr,_ = find_leaflet(fr,hatom,rcut)
                        fr.calc_thick([3,4])
                        fr.calc_area()
                        calc_P2(fr, [[5,6],[9,10]])
                        fr.clear_data()
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
    elif opt == 2:
        fr = pickle.load(open(fname + "_"+str(nbins)+".pckl",'rb'))
        with open("dump_walker_"+str(nbins)+".lmpstraj",'w') as f:
            for frame in fr:
                frame.write_lmpsdump(f)
        print("wrote walker %d" % str(nbins))

    else:
        frames=[]
        for i in range(nbins):
            fr = pickle.load(open(fname + "_" +str(i)+".pckl",'rb'))
            frames.append(fr[-1])
            frames[i].T=i
        with open("dump_singlesnapshot.lmpstraj",'w') as f:
            for i in range(nbins):
                frames[i].write_lmpsdump(f)
                #np.savetxt("temp"+str(i),frames[i].calc_p2([3],[6,10]))
        print("Dumps have been written")

