import numpy as np

ev=2
vv=11

data={}
with open("log.lammps") as f:
    lines=f.readlines()
    flag=0
    keys=[]
    for line in lines:
        if "Loop" in line:
            flag=0
            print("stop")
        if flag == 1:
            for key in keys:
                data[key].append(float(line.strip().split()[loc[key]]))
        if "Step" in line:
            flag=1
            data={}
            loc={}
            keys=line.strip().split()
            count = 0
            for key in keys:
                data[key]=[]
                loc[key]=count
                count+=1
            print("start")

for key in data:
    data[key].pop()
    np.savetxt("%s_init.out"%key,np.c_[data[key]])
    if key == "Volume":
        L = np.array(data[key])**(1./3.)
        np.savetxt("L.dat",np.c_[L])

