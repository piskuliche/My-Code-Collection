#!/usr/bin/env python

import numpy as np
from scipy import interpolate,stats
import argparse, sys
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('-fname', default="TandS_STWHAM.out", type=str, help='dump file name')
parser.add_argument('-Tmin', default=250, type=float, help='Minimum temperature')
parser.add_argument('-Tmax', default=350, type=float, help='Maximum temperature')
parser.add_argument('-tol', default=1e-5, type=float, help='Tolerance for search')
parser.add_argument('-binstep', default=0.000001, type=float, help='Step for temperature bins')
parser.add_argument('-surface_area', default=900, type=float, help='Surface area')
parser.add_argument('-nlipids', default=32, type=int, help='Number of lipids')
parser.add_argument('-nblocks', default=5, type=int, help='Number of blocks')
args = parser.parse_args()


fname   = args.fname
Tmin    = args.Tmin
Tmax    = args.Tmax
binstep = args.binstep
tol     = args.tol
surface_area = args.surface_area
nlipids = args.nlipids
nblocks = args.nblocks

def read_data(filname,hcol,tcol,scol):
    H, T, S = np.genfromtxt(filname, usecols=(hcol,tcol,scol),unpack=True)
    msk = np.where(S!=0)
    H,T,S=H[msk],T[msk],S[msk]
    return H,T,S

def hull_lin(x1,x2,y1,y2,x):
    a = (y2-y1)/(x2-x1)
    b = y2-x2*a
    return a*x+b

def calc_gibbs_hull(roots,H,T,S,brange,verbose):
    fitS = interpolate.splrep(H,S)
    print("Entropy values at roots", interpolate.splev(roots,fitS))
    #define hull
    hx = [roots[0],roots[-1]]
    hy = [interpolate.splev(roots[0],fitS),interpolate.splev(roots[-1],fitS)]
    # hull coeffs
    a=(hy[1]-hy[0])/(hx[1]-hx[0])
    b=hy[1]-hx[1]*(hy[1]-hy[0])/(hx[1]-hx[0])
    msk1 = H>roots[0]
    msk2 = H<roots[-1]
    msk  = msk1 * msk2
    trunc_T = T[msk]
    trunc_H = H[msk]
    trunc_S = S[msk]
    hull    = [hull_lin(hx[0],hx[1],hy[0],hy[1],x) for x in trunc_H]
    ent     = [interpolate.splev(x,fitS) for x in trunc_H]
    diff    = [(hull[i]-ent[i]) for i in range(len(hull))]
    S_Surf  = np.max(diff)
    latent_heat = roots[-1]-roots[0]
    entropy_fusion = interpolate.splev(roots,fitS)[-1]-interpolate.splev(roots,fitS)[0]
    surface_tension = max(diff)/brange/(2*surface_area)
    fen = trunc_H - 1/brange*trunc_S
    data={}
    data["S_Surf"]=S_Surf
    data["S_Fus"]=entropy_fusion
    data["latent_heat"]=latent_heat
    data["surface_tension"]=surface_tension
    data["Gbar"]=max(fen)-min(fen)
    if verbose == 1:
        np.savetxt("thermo_data", np.c_[trunc_T, trunc_H, trunc_S,fen,hull])
        print("Surface entropy is %10.5f kcal/mol/K" % S_Surf)
        print("in cal/mol/K = %10.5f" % (S_Surf*1000))
        print("Located at enthalpy %10.5f kcal/mol" % trunc_H[diff.index(max(diff))])
        print("Entropy at this point %10.5f kcal/mol/K" % trunc_S[diff.index(max(diff))])
        print("Gibbs hull at this point %10.5f kcal/mol" % hull[diff.index(max(diff))])
        print("Surface Area chosen as %10.5f A^2" % surface_area)
        print("Entropy of fusion %10.5f kcal/mol/K" % entropy_fusion)
        print("...in cal/mol/K = %10.5f" % (entropy_fusion*1000))
        print("...per lipid %10.5f" % (entropy_fusion*1000/nlipids))
        print("Latent Heat %10.5f kcal/mol" % (latent_heat))
        print("...per lipid %10.5f" % (latent_heat/nlipids))
        print("Surface tension in kcal/mol A^2 %10.5f" % surface_tension)
        print("...in cal/m^2 = %10.5f" % (surface_tension*nlipids*1000/6.02214129))
        print("Free Energy Barrier %10.5f kcal/mol" % (max(fen)-min(fen)))
        print("Limit of stability is %10.5f K" % max(trunc_T))

    return data

def construct_maxwell(H,T,S,binstep,verbose):
    import matplotlib.pyplot as plt
    beta = np.divide(1,T)
    beta_range = np.arange(1/Tmax,1/Tmin, binstep)
    count, found = 0, 0
    data={}
    while found == 0:
        count += 1
        if count>100: sys.exit("Error: count too high")
        for t in range(len(beta_range)):
            dT = beta - beta_range[t]
            # Represent dT as a spline interpolation
            dTfit = interpolate.splrep(H,dT)
            """ 
            hmin, hmax = H.min(), H.max()
            hh = np.linspace(hmin,hmax,200)
            spline = interpolate.BSpline(dTfit[0],dTfit[1],dTfit[2],extrapolate=False)
            plt.plot(H,dT,'bo')
            plt.plot(hh, spline(hh),color="red")
            plt.show()
            """
            
            # Find the roots of dT
            roots = interpolate.sproot(dTfit, mest=200)
            if len(roots)<3:
                print(roots,np.shape(dTfit))
                sys.exit("Error: Not enough crossing in interpolation")
            area = []
            # Calculate area of roots by integrating between them
            for r in range(len(roots)-1):
                area.append(interpolate.splint(roots[r],roots[r+1],dTfit))
            area1 = 0
            area2 = 0
            # Sum areas
            for a in area:
                if a>0: area1 += a
                else: area2 += a
            area_sum = area1 + area2

            # Check whether total area is zero
            if np.abs(area_sum)<tol:
                found = 1
                Tm = 1/beta_range[t]
                print("Transition temperature is %10.5f K" % (Tm))
                data=calc_gibbs_hull(roots,H,T,S,beta_range[t],verbose)
                data["Tm"] = Tm
            # If area is less than zero, already passed Tm
            # Reduce the window, and try again
            if area_sum<0:
                bmin = beta_range[t-1]
                bmax = beta_range[t+1]
                binstep /= 10
                if binstep<1e-12:
                    sys.exit("Error: bin too small")
                beta_range = np.arange(bmin,bmax,binstep)
                if len(beta_range) == 0:
                    sys.exit("Error: ran out of temperatures")
                break
    return data

def Error(total, block):
    print("*****Final Data*****")
    t_val = stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)
    units={"Tm":"K","latent_heat":"kcal/mol", "Gbar":"kcal/mol","S_Surf":"kcal/mol/K", "S_Fus":"kcal/mol/K", "surface_tension":"kcal/mol/A^2"}
    for key in total_data:
        blockvals = []
        for b in range(nblocks): blockvals.append(block[b][key])
        pickle.dump(blockvals,open(key + "_blockvals.pckl",'wb'),protocol=pickle.HIGHEST_PROTOCOL)
        err = np.std(blockvals)*t_val
        print("%.20s: %10.5f +/- %10.5f %s" % (key,total[key],err,units[key]))
        if key not in ["Tm","surface_tension","Gbar"]:
            print("...%.20s %10.5f +/- %10.5f %s" % ("per lipid",total[key]/nlipids, err/nlipids, units[key]))
    return
if __name__ == "__main__":
    H,T,S = read_data(fname,0,1,3)
    bl_data=[]
    for b in range(nblocks):
        bname = "TandS_bl_"+str(b)+".out"
        h, t, s = read_data(bname,0,1,2)
        bl_data.append(construct_maxwell(h,t,s,binstep,0))
    total_data = construct_maxwell(H,T,S,binstep,1)
    print("test")
    Error(total_data,bl_data)
