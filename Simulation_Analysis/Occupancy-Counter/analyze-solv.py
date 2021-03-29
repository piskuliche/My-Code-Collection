#!/usr/bin/env python
import numpy as np
from scipy import stats

nblocks=5
t_val=stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)

names = np.genfromtxt("inp_vals", usecols=(0),dtype=str, skip_header=1)
data = np.genfromtxt("solv_data.dat",dtype=int, unpack=True)
data2 = np.genfromtxt("hbonds.dat", dtype=int, unpack=True)
nwater=0.0
nacryl=0.0
# read input file
with open("solv.in",'r') as f:
    lines=f.readlines()
    nwater,nacryl=lines[-1].strip().split()
    nacryl=float(nacryl)
    nwater=float(nwater)


av={}
err={}
name={}
avhbnd={}
erhbnd={}

dv = []
dv2 = []
for i in range(data.shape[0]):
    name[names[i]]=names[i]
    av[names[i]]=np.average(data[i])
    avhbnd[names[i]]=np.average(data2[i])
    dv.append(av[names[i]])
    dv2.append(avhbnd[names[i]])
    bl_data = np.average(np.split(data[i],nblocks),axis=0)
    bl_data2 = np.average(np.split(data2[i],nblocks),axis=0)
    err[names[i]]=np.std(bl_data)*t_val
    erhbnd[names[i]]=np.std(bl_data2)*t_val

f=open("acryl_solvation_values.dat",'w')
g=open("acryl_solv_scale.dat",'w')
h=open("acryl_hbnd_values.dat",'w')
h2=open("acryl_hbnd_scale.dat",'w')
Nsum=0.0
Nerr=0.0
Osum=0.0
Oerr=0.0
Ncount=0
Ocount=0
for key in av:
    f.write("%s %s %s\n" % (name[key],av[key]/nacryl,err[key]/nacryl))
    g.write("%s %s\n" % (name[key],av[key]/nacryl/(np.max(dv)/nacryl)))
    h.write("%s %s %s\n" % (name[key],avhbnd[key]/nacryl,erhbnd[key]/nacryl))
    h2.write("%s %s\n" % (name[key],avhbnd[key]/nacryl/(np.max(dv2)/nacryl)))
    if "N" in key:
        Nsum += avhbnd[key]/nacryl
        Nerr += erhbnd[key]/nacryl
        Ncount+=1
    if "O" in key:
        Osum += avhbnd[key]/nacryl
        Oerr += erhbnd[key]/nacryl
        Ocount+=1
chk=open("N-hbondvals.dat",'w')
chk.write("%s %s %s\n" % (Ncount, Nsum/Ncount, Nerr/Ncount))
chk.close()
chkO=open("O-hbondvals.dat",'w')
chkO.write("%s %s %s\n" % (Ocount, Osum/Ocount, Oerr/Ocount))
chkO.close()

f.close()
g.close()
h.close()
h2.close()
