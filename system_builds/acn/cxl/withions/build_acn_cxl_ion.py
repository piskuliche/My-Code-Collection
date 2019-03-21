#!/usr/bin/env python

import numpy as np
import sys, subprocess
import os
"""
This is a python program which creates a box of acn for use in lammps simulations.
"""

if len(sys.argv) != 6:
    print("Usage python build_acn.py numacns numco2s numions type boxlength(F or #)")
    exit()

numacns=int(sys.argv[1])
numco2s=int(sys.argv[2])
numions=int(sys.argv[3])
type=int(sys.argv[4])
boxparam=0.0
if sys.argv[5] != "F":
    boxparam = float(sys.argv[5])


def gen_packmol(nacn, nco2, nions):
    n=nacn+nco2
    if boxparam == 0.0:
        L=(n/0.034)**(1/3.)-2.
    else:
        L=boxparam-2.
    pmolfile = "acncx.pmol"
    pm = open(pmolfile,'w')
    pm.write("tolerance 2.0\n")
    pm.write("filetype xyz\n")
    pm.write("output system.xyz\n")
    pm.write("\n")
    pm.write("structure acn.xyz\n")
    pm.write("  number %s\n" % nacn)
    pm.write("  inside box 2. 2. 2. %s %s %s\n" % (L, L, L))
    pm.write("end structure\n")
    pm.write("\n")
    pm.write("structure co2.xyz\n")
    pm.write("  number %s\n" % nco2)
    pm.write("  inside box 2. 2. 2. %s %s %s\n" % (L, L, L))
    pm.write("end structure\n")
    pm.write("\n")
    pm.write("structure clo4.xyz\n")
    pm.write("  number %s\n" % nions)
    pm.write("  inside box 2. 2. 2. %s %s %s\n" % (L, L, L))
    pm.write("end structure\n")
    pm.write("\n")
    pm.write("structure li.xyz\n")
    pm.write("  number %s\n" % nions)
    pm.write("  inside box 2. 2. 2. %s %s %s\n" % (L, L, L))
    pm.write("end structure\n") 
    pm.close()

    acnfile = "acn.xyz"
    acn = open(acnfile, 'w')
    acn.write("3\n")
    acn.write("acn\n")
    acn.write("ch3          0.000   0.000   0.000\n")
    acn.write("c            1.540   0.000   0.000\n")
    acn.write("n            2.697   0.000   0.000\n")
    acn.close()

    co2file = "co2.xyz"
    co2 = open(co2file, 'w')
    co2.write("3\n")
    co2.write("co2\n")
    co2.write("C             0.000   0.000   0.000\n")
    co2.write("O             1.160   0.000   0.000\n")
    co2.write("O            -1.160   0.000   0.000\n")
    co2.close()

    clo4file = "clo4.xyz"
    per = open(clo4file, 'w')
    per.write("5\n")
    per.write("clo4\n")
    per.write("Cl            0.245   0.253   0.064\n")
    per.write("O             0.889   1.166  -0.802\n")
    per.write("O            -0.996   0.811   0.437\n")
    per.write("O             0.970   0.081   1.265\n")
    per.write("O            -0.009  -0.987  -0.567\n")
    per.close()

    lifile = "li.xyz"
    li = open(lifile, 'w')
    li.write("1\n")
    li.write("li\n")
    li.write("Li             0.000   0.000   0.000\n")
    li.close()

    print("Please run the following command: packmol < acncx.pmol")


def tolammps(nacn,nco2,nions):
    n=nacn+nco2+nions
    datafile = "data.acncx"
    dat = open(datafile, 'w')
    x,y,z = np.genfromtxt('system.xyz', usecols=(1,2,3), unpack=True, skip_header=2)
    q=[0.269, 0.129, -0.398,0.700,-0.350,-0.350,0.974754, -0.4936885, -0.4936885, -0.4936885, -0.4936885,1.00]
    type=[1,2,3,4,5,5,6,7,7,7,7,8]
    if boxparam == 0.0:
        L = (n/0.034)**(1/3.) 
    else:
        L = boxparam
    dat.write("LAMMPS Description\n")
    dat.write("\n")
    dat.write("%s atoms\n" % (nacn*3+nco2*3+nions*6))
    dat.write("%s bonds\n" % (nacn*2+nco2*2+nions*4))
    dat.write("%s angles\n" % (nacn+nco2+nions*4))
    dat.write("0 dihedrals\n")
    dat.write("0 impropers\n")
    dat.write("\n")
    dat.write("8 atom types\n")
    dat.write("4 bond types\n")
    dat.write("3 angle types\n")
    dat.write("\n")
    dat.write("0.0 %s xlo xhi\n" % L)
    dat.write("0.0 %s  ylo yhi\n" % L)
    dat.write("0.0 %s  zlo zhi\n" % L)
    dat.write("\n")
    dat.write("Masses\n")
    dat.write("\n")
    # Acetonitrile
    dat.write("1 15.035\n")
    dat.write("2 12.011\n")
    dat.write("3 14.007\n")
    # Carbon Dioxide
    dat.write("4 12.011\n")
    dat.write("5 15.999\n")
    # Lithium Perchlorate
    dat.write("6 35.45\n")
    dat.write("7 15.999\n")
    dat.write("8 6.94\n")
    dat.write("\n")
    dat.write("Atoms\n")
    dat.write("\n")
    for i in range(nacn):
        for j in range(0,3):
            aindex=i*3+1+j
            mindex=i+1
            dat.write("%s %s %s %s %s %s %s\n" % (aindex, mindex, type[j], q[j], x[aindex-1], y[aindex-1], z[aindex-1]))
    for i in range(nco2):
        for k in range(3,6):
            j=k-3
            aindex=nacn*3+i*3+1+j
            mindex=nacn+i+1
            dat.write("%s %s %s %s %s %s %s\n" % (aindex, mindex, type[k], q[k], x[aindex-1], y[aindex-1], z[aindex-1]))
    for i in range(nions):
        for k in range(6,11):
            j=k-5
            aindex=nacn*3+nco2*3+j+i*5
            mindex=nacn+nco2+i+1
            dat.write("%s %s %s %s %s %s %s\n" % (aindex,mindex,type[k],q[k],x[aindex-1],y[aindex-1],z[aindex-1]))
    for i in range(nions):
        k = 11
        aindex = nacn*3 + nco2*3 + nions*5 + i+1
        mindex = nacn + nco2 + nions + i + 1
        dat.write("%s %s %s %s %s %s %s\n" % (aindex,mindex,type[k],q[k],x[aindex-1],y[aindex-1],z[aindex-1]))
    dat.write("\n")
    dat.write("Bonds\n")
    dat.write("\n")
    for i in range(nacn):
        b1index=i*2+1
        b2index=i*2+2
        btype = [1,2]
        ch3index =i*3+1
        cindex=i*3+2
        nindex=i*3+3
        dat.write("%s %s %s %s\n" % (b1index,btype[0],cindex,ch3index))
        dat.write("%s %s %s %s\n" % (b2index,btype[1],cindex,nindex))
    for i in range(nco2):
        b1index=nacn*2+i*2+1
        b2index=nacn*2+i*2+2
        btype=3
        cindex=nacn*3+i*3+1
        o1index=nacn*3+i*3+2
        o2index=nacn*3+i*3+3
        dat.write("%s %s %s %s\n" % (b1index,btype,cindex,o1index))
        dat.write("%s %s %s %s\n" % (b2index,btype,cindex,o2index))
    for i in range(nions):
        b1index=nacn*2+nco2*2+1+i*4
        b2index=b1index+1
        b3index=b2index+1
        b4index=b3index+1
        clindex=nacn*3+nco2*3+1+i*5
        o1index=nacn*3+nco2*3+2+i*5
        o2index=nacn*3+nco2*3+3+i*5
        o3index=nacn*3+nco2*3+4+i*5
        o4index=nacn*3+nco2*3+5+i*5
        btype = 4
        dat.write("%s %s %s %s\n" % (b1index,btype,clindex,o1index))
        dat.write("%s %s %s %s\n" % (b2index,btype,clindex,o2index))
        dat.write("%s %s %s %s\n" % (b3index,btype,clindex,o3index))
        dat.write("%s %s %s %s\n" % (b4index,btype,clindex,o4index))

    dat.write("\n")
    dat.write("Angles\n")
    dat.write("\n")
    for i in range(nacn):
        aindex=i+1
        atype=1
        CH3index =i*3+1
        Cindex=i*3+2
        Nindex=i*3+3
        dat.write("%s %s %s %s %s\n" % (aindex, atype, CH3index, Cindex, Nindex))
    for i in range(nco2):
        aindex=nacn+i+1
        atype=2
        cindex=nacn*3+i*3+1
        o1index=nacn*3+i*3+2
        o2index=nacn*3+i*3+3
        dat.write("%s %s %s %s %s\n" % (aindex, atype, o1index, cindex, o2index))
    for i in range(nions):
        a1index=nacn+nco2+1+i*4
        a2index=a1index+1
        a3index=a2index+1
        a4index=a3index+1
        atype=3
        clindex=nacn*3+nco2*3+i*5+1
        oindex=nacn*3+nco2*3+i*5+2
        dat.write("%s %s %s %s %s\n" %(a1index,atype,oindex,clindex,oindex+1))
        dat.write("%s %s %s %s %s\n" %(a2index,atype,oindex+1,clindex,oindex+2))
        dat.write("%s %s %s %s %s\n" %(a3index,atype,oindex+2,clindex,oindex+3))
        dat.write("%s %s %s %s %s\n" %(a4index,atype,oindex+3,clindex,oindex))
    dat.close()

if type == 0:
    gen_packmol(numacns,numco2s,numions)
elif type == 1:
    tolammps(numacns,numco2s,numions)
elif type == 3:
    gen_packmol(numacns,numco2s,numions)
    subprocess.call(["/panfs/pfs.local/work/laird/e924p726/thompsonwork/Programs/Executables/packmol < acncx.pmol"],shell=True)
    tolammps(numacns,numco2s,numions)
else:
    print("Error: Incorrect chosen type")
