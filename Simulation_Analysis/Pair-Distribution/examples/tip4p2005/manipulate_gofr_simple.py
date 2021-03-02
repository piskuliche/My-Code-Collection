import numpy as np
from scipy import stats,integrate

# This is a more streamlined code to calculate the parameters for g(r)

kb=0.0019872041
T=298.15
nblocks=5
t_val=stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)
cut=125


def calc_PMF(gofr): 
    return -kb*T*np.log(gofr,where=gofr!=0,out=np.zeros_like(gofr))

def calc_A(x,pmf):
    return pmf-2*kb*T*np.log(x,where=x!=0,out=np.zeros_like(x))

def calc_U(gofr,egofr):
    return np.divide(egofr,gofr,where=gofr!=0,out=np.zeros_like(gofr))

def calc_S(U,A):
    return (U-A)/T

def Error(arr):
    return np.std(arr,axis=0)*t_val

x,data = np.genfromtxt("pairdist_1_1.dat",usecols=(0,1),unpack=True,skip_header=cut)
edata = np.genfromtxt("epairdist_1_1.dat",usecols=(1),unpack=True,skip_header=cut)
PMF = calc_PMF(data)
A   = calc_A(x,PMF)
U   = calc_U(data,edata)
S   = calc_S(U,A)
mTdS = S*-T

bl_data,bl_edata=[],[]
bl_PMF, bl_A, bl_U,bl_S = [], [], [],[]
for block in range(nblocks):
    bl_data.append(np.genfromtxt("bl_%d_pairdist_1_1.dat"%(block+1),usecols=1,unpack=True,skip_header=cut))
    bl_edata.append(np.genfromtxt("bl_%d_epairdist_1_1.dat"%(block+1), usecols=1,unpack=True,skip_header=cut))
    bl_PMF.append(calc_PMF(bl_data[block]))
    bl_A.append(calc_A(x,bl_PMF[block]))
    bl_U.append(calc_U(bl_data[block],bl_edata[block]))
    bl_S.append(calc_S(bl_U[block],bl_A[block]))

A_err = Error(bl_A)
PMF_err = Error(bl_PMF)
RDF_err = Error(bl_data)
eRDF_err = Error(bl_edata)
U_err = Error(bl_U)
S_err = Error(bl_S)
mTdS_err = S_err*T



np.savetxt("A.dat",np.c_[x,A,A_err])
np.savetxt("PMF.dat",np.c_[x,PMF,PMF_err])
np.savetxt("RDF.dat",np.c_[x,data,RDF_err])
np.savetxt("eRDF.dat",np.c_[x,edata,eRDF_err])
np.savetxt("U.dat", np.c_[x,U,U_err])
np.savetxt("mTdS.dat",np.c_[x,mTdS,mTdS_err])
np.savetxt("S.dat",np.c_[x,S,S_err])

            
