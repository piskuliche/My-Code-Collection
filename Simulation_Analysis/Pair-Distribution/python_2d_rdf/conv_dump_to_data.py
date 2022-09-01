import numpy as np
def pull_header(lines):
    h_lines = []
    for line in lines:
        if "Masses" in line:
            break
        h_lines.append(line)
    return h_lines

def pull_sections(lines):
    header=pull_header(lines)
    nheader = len(header)
    lines=lines[nheader:]
    flag=0
    sections={}
    sec = None
    for line in lines:
        cleanline = line.strip().split()
        if len(cleanline)==1 and cleanline != "\n":
            sec=cleanline[0]
            sections[sec]=[]
        if sec is not None and len(cleanline) > 1:
                sections[sec].append(line)
    return header,sections

def Read_Data(file):
    """This function takes a data file and saves the important bits, but also clears the non-important bits."""
    with open(file,'r') as f:
        lines=f.readlines()
        header,sections=pull_sections(lines)
        print("Found the following keys:")
        for key in sections:
            print("-"+key+" has %d lines" % len(sections[key]))
        return header, sections


def set_map(sections):
    masslines=sections["Masses"]
    M={}
    for mass in masslines:
        mvals=mass.strip().split()
        M[int(mvals[0])]=float(mvals[1])
    atomlines=sections["Atoms"]
    mid={}
    typ={}
    chg={}
    for atom in atomlines:
        avals=atom.strip().split()
        mid[avals[0]]=int(avals[1])
        typ[avals[0]]=int(avals[2])
        chg[avals[0]]=float(avals[3])
    """
    bondlines=sections["Bonds"]
    btyp=[]
    bat1=[]
    bat2=[]
    for bond in bondlines:
        bvals = bond.strip().split()
        btyp.append(int(bvals[1]))
        bat1.append(int(bvals[2]))
        bat2.append(int(bvals[3]))
    anglines=sections["Angles"]
    atyp=[]
    aat1=[]
    aat2=[]
    aat3=[]
    for ang in anglines:
        angvals=ang.strip().split()
        atyp=angvals[1]
        aat1.append(int(angvals[2]))
        aat2.append(int(angvals[3]))
        aat3.append(int(angvals[4]))
    return M,mid,typ,chg,btyp,bat1,bat2,atyp,aat1,aat2,aat3
    """
    return M,mid

def Read_Config(df,natoms):
    L=[0,0,0]
    r=[]
    molid=[]
    typval=[]
    df.readline()
    ts=int(df.readline().strip())
    df.readline()
    df.readline()
    df.readline()
    line = df.readline().strip().split()
    L[0] = -float(line[0])+float(line[1])
    line = df.readline().strip().split()
    L[1] = -float(line[0])+float(line[1])
    line = df.readline().strip().split()
    L[2] = -float(line[0])+float(line[1])
    df.readline()
    for l in range(natoms):
        line = df.readline().strip().split()
        mol,typ,x,y,z = int(line[1]),int(line[2]),float(line[3]),float(line[4]),float(line[5])
        r.append([x,y,z])
        molid.append(mol)
        typval.append(typ)
    return np.array(L),np.array(r),ts,molid,typval


def Get_Distances(r,L):
    """
    This is a pretty complicated function to do a simple thing.
    This code takes two vectors of size(m,3) and size(n,3)
    Then it uses numpy broadcasting to calculate ALL the pairwise distances and outputs them in a matrix, sized (m,n)
    Then I use it to calculate the distances.
    (described here: https://stackoverflow.com/questions/60039982/numpy-python-vectorize-distance-function-to-calculate-pairwise-distance-of-2-ma/60040269#60040269)
    """
    vecdr = r[:,np.newaxis,:]-r[np.newaxis,:,:]
    vecdr = vecdr - np.multiply(L,np.round(np.divide(vecdr,L)))
    dr = np.linalg.norm(vecdr,axis=-1)
    return dr

def Sort_Distances(dr):
    """
    This function takes an array (dr) and sorts it by the last axis.

    It returns:
    1) the idx map of the sort
    2) the sorted dr array (sdr)
    """
    natoms = np.shape(dr)[0]
    idx = np.argsort(dr,axis=-1,kind='quicksort')
    #print(idx)
    sdr=np.take_along_axis(dr, idx, axis=-1)
    return idx,sdr

def Prep_Dir(mol,kspace=False):
    import os
    if not os.path.exists("./include_dir"):
        os.makedirs("./include_dir")
    if not os.path.exists("./ener_dir"):
        os.makedirs("./ener_dir")
    if not os.path.exists("./calc_dir"):
        os.makedirs("./calc_dir")
    if not os.path.exists("./calc_dir/logs"):
        os.makedirs("./calc_dir/logs")
    nmols=len(np.unique(mol))
    for i in range(nmols):
        fi = open("./include_dir/include.groups-%d"%i,'w')
        fi.write("group solu molecule 0\n")
        fi.write("group close molecule 0\n")
        fi.write("group far molecule 0\n")
        fi.write("compute sske solu ke\n")
        fi.write("compute sspair solu group/group solu pair yes\n")
        if kspace==True: fi.write("compute sskspc solu group/group solu kspace yes\n")
        fi.write("compute ssbnd solu bond\n")
        fi.write("compute ssang solu angle\n")
        #fi.write("compute ssdih solu dihedral\n")
        #fi.write("compute ssimp solu improper\n")
        #fi.write("variable ssintra equal c_ssbnd+c_ssang+c_ssdih+c_ssimp\n")
        fi.write("variable ssintra equal c_ssbnd[1]+c_ssang[1]\n")
        fi.write("compute ccpair close group/group close pair yes\n")
        if kspace==True: fi.write("compute cckspc close group/group close kspace yes\n")
        fi.write("compute ccke close ke\n")
        fi.write("compute ccbnd close bond\n")
        fi.write("compute ccang close angle\n")
        #fi.write("compute ccdih close dihedral\n")
        #fi.write("compute ccimp close improper\n")
        fi.write("variable ccintra equal c_ccbnd[1]+c_ccang[1]\n")
        fi.write("compute ffpair far group/group far pair yes\n")
        if kspace==True: fi.write("compute ffkspc far group/group far kspace yes\n")
        fi.write("compute ffke far ke\n")
        fi.write("compute ffbnd far bond\n")
        fi.write("compute ffang far angle\n")
        #fi.write("compute ffdih far dihedral\n")
        #fi.write("compute ffimp far improper\n")
        #fi.write("variable ffintra equal c_ffbnd+c_ffang+c_ffdih+c_ffimp\n")
        fi.write("variable ffintra equal c_ffbnd[1]+c_ffang[1]\n")
        fi.write("compute sfpair solu group/group far pair yes\n")
        if kspace==True: fi.write("compute sfkspc solu group/group far kspace yes\n")
        fi.write("compute scpair solu group/group close pair yes\n")
        if kspace==True: fi.write("compute sckspc solu group/group close kspace yes\n")
        fi.write("compute cfpair close group/group far pair yes\n")
        if kspace==True: fi.write("compute cfkspc close group/group far kspace yes\n")
        fi.write("fix 1 all ave/time 1000 1 1000 c_sske c_ccke c_ffke file ../ener_dir/ke.compute%d\n"%i)
        fi.write("fix 2 all ave/time 1000 1 1000 c_sspair c_ccpair c_ffpair c_scpair c_sfpair c_cfpair file ../ener_dir/pair.compute%d\n"%i)
        if kspace==True: fi.write("fix 3 all ave/time 1000 1 1000 c_sskspc c_cckspc c_ffkspc c_sckspc c_sfkspc c_cfkspc file ../ener_dir/kspc.compute%d\n"%i)
        fi.write("fix 4 all ave/time 1000 1 1000 v_ssintra v_ccintra v_ffintra file ../ener_dir/intra.compute%d\n"%i)
        fi.close()






def Write_Groups(dr,ts):
    nmols=np.shape(dr)[0]
    for i in range(nmols):
        fi = open("./include_dir/include.groups-%d"%i,'a')
        close=np.where(dr[0]<20)[-1]
        fi.write("group lt molecule ")
        for i in close:
            fi.write("%d "% (i+1))
        fi.write("\n")
        fi.write("group solu molecule %d\n"%i)
        fi.write("group close subtract lt solu\n")
        fi.write("group far subtract all lt\n")
        fi.write('print "beginning rerun %d"\n'%ts)
        fi.write("rerun ../dump.lammpsdump first %d last %d dump x y z\n" % (ts,ts))
        fi.write('print "finishing rerun"\n')
        fi.close()

def Calc_Com(r,mol,typ,M):
    natoms=len(mol)
    nmols=len(np.unique(mol))
    comr=np.zeros((nmols,3))
    Mtot=np.zeros(nmols)
    for i in range(natoms):
        mid=mol[i]-1
        atype=typ[i]
        for j in range(3):
            comr[mid][j] += r[i][j]*M[atype]
        Mtot[mid]+=M[atype]
    for i in range(nmols):
        comr[i]=comr[i]/Mtot[i]
    return comr 

def Gen_Configs(datafile,nconfigs,sections,sysfile):
    print("Beginning Mapping Process")
    M,mid = set_map(sections)
    natoms=len(mid)
    nmols=None
    with open(datafile, 'r') as df:
        for i in range(nconfigs):
            if i%100 == 0:
                print("Reached config %d of %d" % (i,nconfigs))
            L,r,ts,mol,typid=Read_Config(df,natoms)
            if i == 0:
                Prep_Dir(mol)
            comr=Calc_Com(r,mol,typid,M)
            nmols=len(np.unique(mol))
            dr = Get_Distances(comr,L)
            #idx,sdr=Sort_Distances(dr)
            Write_Groups(dr,ts)
    Write_System(sysfile,nmols)
    return

def Write_System(sysfile,nmols):
    print("Writing System Files")
    lines=None
    with open(sysfile,'r') as f:
        lines = f.readlines()
    for i in range(nmols):
        with open("calc_dir/"+sysfile+"-%d"%i,'w') as g:
            for line in lines:
                g.write(line)
            g.write("log logs/log.%d.out\n"%i)
            g.write("include ../include_dir/include.groups-%d\n"%i)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='This code creates a series of files for lammps simulations that quickly rerun and re-calculate energies')
    parser.add_argument('-data',default="equil.data",type=str,help='Name of data file [default equil.data]')
    parser.add_arugment('-dump',default="dump.lammpsdump",type=str,help='Name of dump file [default dump.lammpsdump]')
    parser.add_argument('-nc',default=10000, type=int, help='Number of configurations [default 10000]')
    parser.add_arguemtn('-sf', default="system.in", type=str, help='System input file [default system.in]')
    args = parser.parse_args()

    datafile = args.data
    dumpfile = args.dump
    nconfigs = args.nc
    sysfile  = args.sf

    header,sections=Read_Data(datafile)
    Gen_Configs(dumpfile,nc,sections,sysfile)
