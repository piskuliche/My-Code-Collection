#!/usr/bin/env python


import numpy as np
import argparse
from scipy import stats

parser = argparse.ArgumentParser()
parser.add_argument('-f', default="log.lammps", help='Log file to parse')
parser.add_argument('-nb', default=5, help='Number of blocks')
parser.add_argument('-info', default="value", help='Information to append to output files')
parser.add_argument('-ew', default=0, help='Information to append to output files')


args = parser.parse_args()
fname = str(args.f)
nblocks=int(args.nb)
info=str(args.info)
ew=int(args.ew)


convp = 1.439E-5

readflag=0
runs = []
runstep = -1
items = []
data={}
with open(fname, 'r') as f:
    lines=f.readlines()
    for line in lines:
        if "Step" in line:
            readflag = 1
            runstep += 1
            data = {}
            if runstep == 0:
                for item in line.split():
                    items.append(item)
            for item in items:
                data[item] = []
        if "Loop time" in line:
            readflag = 0
            runs.append(data)
        if readflag == 1:
            if line.split()[0] != items[0]:
                for i in range(len(items)):
                    if "WARNING" not in line:
                        data[items[i]].append(float(line.split()[i]))

print("There are %s run commands" % len(runs))
for r in range(len(runs)):
    print("*****************************************************************")
    print("Run %d has %d data points" % (r, len(runs[r][items[0]])))
    if runs[r]["Volume"][0]==runs[r]["Volume"][1]:
        print("Run Type is NVT")
    else:
        print("Run Type is NPT")

    nitems = len(runs[r][items[0]])-1
    nperb = int(nitems/nblocks)
    print(nperb)
    print(nitems)
    t_val=stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)
    dE=[]
    d2E=[]
    bl_dE=[]
    bl_d2E=[]
    dV=[]
    d2V=[]
    bl_dV=[]
    bl_d2V=[]
    for item in items:
        print(item)
        np.savetxt(str(item)+'_init.out', np.c_[runs[r][item]])
        bl = []
        for b in range(nblocks):
            bstart = b*nperb
            bend = bstart + nperb
            bl.append(np.average(runs[r][item][bstart:bend]))
        std = np.std(bl)*t_val

        print("The average of item %s is %.5f +/- %.5f" %(item,np.average(runs[r][item]),std))
        f = open("run"+str(r+1)+"_"+str(item)+'.avg','w')
        f.write('%s %.5f %.5f\n' % (info,np.average(runs[r][item]),std))
        f.close()
        if item=="Volume":
            L=np.asarray(runs[r][item])**(1/3.)
            np.savetxt('L.dat', np.c_[L])
        if (r == len(runs)-1 and not any(x in item for x in ["Px","Py","Pz","Time","Step","Press","E_"]) and ew==1):
            if item == "TotEng":
                av_en = np.average(runs[r][item])
                for i in range(len(runs[r][item])):
                    dE.append(runs[r][item][i]-av_en)
                    d2E.append(dE[i]**2.)
                for b in range(nblocks):
                    tmp1 = []
                    tmp2 = []
                    bstart=b*nperb
                    bend=(b+1)*nperb
                    av_en = np.average(runs[r][item][bstart:bend])
                    for i in range(bstart,bend):
                        tmp1.append(runs[r][item][i]-av_en)
                        tmp2.append((runs[r][item][i]-av_en)**2.)
                    bl_dE.append(tmp1)
                    bl_d2E.append(tmp2)
            elif item == "Volume":
                av_V = np.average(runs[r][item])
                for i in range(len(runs[r][item])):
                    dV.append(runs[r][item][i]-av_V)
                    d2V.append(dV[i]**2.)
                for b in range(nblocks):
                    tmp1 = []
                    tmp2 = []
                    bstart=b*nperb
                    bend=(b+1)*nperb
                    av_V = np.average(runs[r][item][bstart:bend])
                    for i in range(bstart,bend):
                        tmp1.append(runs[r][item][i]-av_V)
                        tmp2.append((runs[r][item][i]-av_V)**2.)
                    bl_dV.append(tmp1)
                    bl_d2V.append(tmp2)
    if (ew==1): 
        avd2E=np.average(d2E)
        avd2V=np.average(d2V)
        bl_avd2E=[]
        bl_avd2V=[]
        for b in range(nblocks):
            bl_avd2E.append(np.average(bl_d2E[b]))
            bl_avd2V.append(np.average(bl_d2V[b]))
        for item in items:
            ew=[]
            e2w=[]
            vw=[]
            v2w=[]
            if (r == len(runs)-1 and not any(x in item for x in ["Px","Py","Pz","Time","Step","Press","E_"]) and ew==1):
                print("****STEP: %s****" % item)
                for i in range(len(dE)):
                    ew.append(-runs[r][item][i]*dE[i])
                    e2w.append((d2E[i]-avd2E)*runs[r][item][i])
                    vw.append(-runs[r][item][i]*dV[i])
                    v2w.append((d2V[i]-avd2V)*runs[r][item][i])
                ed= np.average(ew)
                e2d=np.average(e2w)
                vd= np.average(vw)
                v2d=np.average(v2w)
                tderiv=ed + convp*vd
                t2deriv=e2d + convp**2.*v2d
                bl_ed = []
                bl_e2d = []
                bl_vd = []
                bl_v2d = []
                bl_tderiv=[]
                bl_t2deriv=[]
                for b in range(nblocks):
                    bstart = b*nperb
                    bend = (b+1)*nperb
                    tmp1=[]
                    tmp2=[]
                    tmp3=[]
                    tmp4=[]
                    for i in range(bstart,bend):
                        tmp1.append(-runs[r][item][i]*bl_dE[b][i-bstart])
                        tmp2.append((bl_d2E[b][i-bstart]-bl_avd2E[b])*runs[r][item][i])
                        tmp3.append(-runs[r][item][i]*bl_dV[b][i-bstart])
                        tmp4.append((bl_d2V[b][i-bstart]-bl_avd2V[b])*runs[r][item][i])
                    bl_ed.append(np.average(tmp1))
                    bl_e2d.append(np.average(tmp2))
                    bl_vd.append(np.average(tmp3))
                    bl_v2d.append(np.average(tmp4))
                    bl_tderiv.append(bl_ed[b]+convp*bl_vd[b])
                    bl_t2deriv.append(bl_e2d[b]+convp**2.*bl_v2d[b])
                ed_err = np.std(bl_ed)*t_val
                e2d_err = np.std(bl_e2d)*t_val
                vd_err = np.std(bl_vd)*t_val
                v2d_err = np.std(bl_v2d)*t_val
                tder_err = np.std(bl_tderiv)*t_val
                t2der_err = np.std(bl_t2deriv)*t_val
                print("e 1st derivative of %s is %s +/- %s" % (item,ed,ed_err))
                print("e 2nd derivative of %s is %s +/- %s" % (item,e2d,e2d_err))
                print("v 1st derivative of %s is %s +/- %s" % (item,vd,vd_err))
                print("v 2nd derivative of %s is %s +/- %s" % (item,v2d,v2d_err))
                print("Total 1st beta deriv of %s is %s +/- %s" %(item, tderiv,tder_err))
                print("Total 2nd beta deriv of %s is %s +/- %s" % (item, t2deriv,t2der_err))

