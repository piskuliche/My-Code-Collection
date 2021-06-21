import numpy as np
import argparse

# Read Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-itp', default="None", type=str, help='Name of the ITP file [default="None"]')
parser.add_argument('-out', default="sys.lt", type=str, help='Name of the output file [default="sys.lt"]')
parser.add_argument('-name', default="TEST", type=str, help='Name of the molecule [default="TEST"]')
parser.add_argument('-coord', default="None", type=str, help='Name of coordinate file [default="None"]')
args= parser.parse_args()

itpfile =   args.itp
outfile =   args.out
name    =   args.name
coordfile = args.coord

if itpfile == "None": exit("Error: Must specify an ITP file for operation")

def read_itp(itpfile):
    flatoms =   "[ atoms ]"
    flbonds =   "[ bonds ]"
    flangls  =  "[ angles ]"
    flagatoms, flagbonds, flagangls = 0, 0, 0
    atoms={"id":[],"type":[],"name":[], "Q":[]}
    bonds={"atom1":[],"atom2":[], "type":[]}
    angls={"atom1":[],"atom2":[],"atom3":[],"type":[]}
    with open(itpfile,'r') as f:
        lines = f.readlines()
        for line in lines:
            if flatoms in line:
                flagatoms = 1 # turn on atoms
                continue
            if flbonds in line:
                flagatoms = 0 # turn of atoms
                flagbonds = 1 # turn on bonds
                continue
            if flangls in line:
                flagatoms, flagbonds = 0, 0 # turn of atoms, bonds
                flagangls = 1 # turn on angles
                continue
            if line == "\n": continue # skips if empty line
            if ";" in line: continue  # skips if commented
            if flagatoms == 1:
                msk=[0,1,4,6]
                tmp = [line.strip().split()[i] for i in msk]
                [x.append(y) for x,y in zip([atoms["id"],atoms["type"],atoms["name"],atoms["Q"]],tmp)]
            if flagbonds == 1:
                msk=[0,1,3]
                tmp = [line.strip().split()[i] for i in msk]
                [x.append(y) for x,y in zip([bonds["atom1"],bonds["atom2"],bonds["type"]],tmp)]
            if flagangls == 1:
                msk=[0,1,2,4]
                tmp = [line.strip().split()[i] for i in msk]
                [x.append(y) for x,y in zip([angls["atom1"],angls["atom2"],angls["atom3"],angls["type"]],tmp)]
    print("There are %d atoms" % len(atoms["id"]))
    print("There are %d bonds" % len(bonds["atom1"]))
    print("There are %d angles" % len(angls["atom1"]))
    print("The ITP file has been read successfully!")
    atoms["id"]=np.array(atoms["id"],dtype=int)
    atoms["type"]=np.array(atoms["type"],dtype=str)
    atoms["name"]=np.array(atoms["name"],dtype=str)
    atoms["Q"]=np.array(atoms["Q"],dtype=float)
    bonds["atom1"]=np.array(bonds["atom1"],dtype=int)
    bonds["atom2"]=np.array(bonds["atom2"],dtype=int)
    bonds["type"]=np.array(bonds["type"],dtype=str)
    angls["atom1"]=np.array(angls["atom1"],dtype=int)
    angls["atom2"]=np.array(angls["atom2"],dtype=int)
    angls["atom3"]=np.array(angls["atom3"],dtype=int)
    angls["type"]=np.array(angls["type"],dtype=str)
    return atoms, bonds, angls

def write_header(outf,name):
    outf.write('import "drymartini.lt"\n')
    outf.write('\n')
    outf.write('%s inherits DRYMARTINI {\n' % name)
    return

def write_atoms(outf,atoms,r):
    if r.size == 0:
        r = np.zeros((3,len(atoms["id"])))
    outf.write("\n")
    outf.write('write("Data Atoms") {\n')
    for atm in range(len(atoms["id"])):
        outf.write('  $atom:%s $mol:. @atom:%s % 6.4f % 10.6f % 10.6f % 10.6f\n' % (atoms["name"][atm], atoms["type"][atm], atoms["Q"][atm], r[0][atm],r[1][atm],r[2][atm]))
    outf.write("}")
    outf.write("\n")
    return

def write_bonds(outf,atoms,bonds):
    outf.write("\n")
    outf.write('write("Data Bond List") {\n')
    for bnd in range(len(bonds["atom1"])):
        outf.write('  $bond:b%d $atom:%s $atom:%s\n' % (bnd,atoms["name"][bonds["atom1"][bnd]-1],atoms["name"][bonds["atom2"][bnd]-1]))
    outf.write("}")
    outf.write("\n")

def write_angls(outf,atoms,angls):
    outf.write("\n")
    outf.write('write("Data Angles By Type") {\n')
    for ang in range(len(angls["atom1"])):
        outf.write('  $angle:%s @atom:%s @atom:%s @atom:%s\n' % (angls["type"][ang], atoms["name"][angls["atom1"][ang]-1],atoms["name"][angls["atom2"][ang]-1],atoms["name"][angls["atom3"][ang]-1]))
    outf.write('}\n')
    outf.write('\n')

def write_footer(outf,name):
    outf.write('\n')
    outf.write('%s.scale(10)\n'%name)
    outf.write('}\n')

def read_coords(coordname):
    r=np.genfromtxt(coordname,usecols=(3,4,5),skip_header=2,skip_footer=1).T
    com=np.zeros(3)
    for i in range(3):
        com[i]=np.sum(r[i]*75)
    com=com/(72*np.shape(r)[1])
    print("Center of mass original", com)
    for i in range(3): r[i]=r[i]-com[i]
    r[2]=r[2]-np.min(r[2])+0.1
    print("Coords are shape ", np.shape(r))
    return r
                
                
if __name__ == "__main__":
    atoms, bonds, angls = read_itp(itpfile)
    r = np.array([])
    with open(outfile,'w') as outf:
        if coordfile != "None": r=read_coords(coordfile)
        write_header(outf,name)
        write_atoms(outf, atoms,r)
        write_bonds(outf,atoms,bonds)
        write_footer(outf,name)
        #write_angls(outf,atoms,angls)
                
                
            
