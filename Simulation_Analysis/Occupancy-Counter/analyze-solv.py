import numpy as np
from scipy import stats

nblocks=5
t_val=stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)

names = np.genfromtxt("inp_vals", usecols=(0),dtype=str, skip_header=1)
data = np.genfromtxt("solv_data.dat",dtype=int, unpack=True)
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

dv = []
for i in range(data.shape[0]):
    name[names[i]]=names[i]
    av[names[i]]=np.average(data[i])
    dv.append(av[names[i]])
    bl_data = np.average(np.split(data[i],nblocks),axis=0)
    err[names[i]]=np.std(bl_data)*t_val

f=open("acryl_solvation_values.dat",'w')
g=open("acryl_solv_scale.dat",'w')
for key in av:
    f.write("%s %s %s\n" % (name[key],av[key]/nacryl,err[key]/nacryl))
    g.write("%s %s\n" % (name[key],av[key]/nacryl/(np.max(dv)/nacryl)))
f.close()
g.close()

