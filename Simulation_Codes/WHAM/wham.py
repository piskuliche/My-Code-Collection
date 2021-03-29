import numpy as np
from numpy import inf,nan
import argparse
import matplotlib.pyplot as plt

kb=0.0019872041
T=298.15
kbT=kb*T

np.set_printoptions(threshold=np.inf)

def calc_bias(x,k):
    return 0.5*k*x**2.

def calc_p(F):
    numerator=np.sum(ni,axis=0)
    denominator=np.sum(np.multiply(cnt,np.exp(np.divide(-np.subtract(U.T,F-np.min(F)),kbT))),axis=1)
    #numerator=np.zeros(nbins)
    #denominator=np.zeros(nbins)
    #for i in range(Nwindows):
    #    for j in range(nbins):
    #        numerator[j]+=ni[i][j]
    #        denominator[j]+=cnt[i]*np.exp(F[i]/kbT)*np.exp(-U[i][j]/kbT)
    return numerator/denominator

def calc_f(P):
    F=-kbT*np.log(np.sum(np.multiply(P,np.exp(-U/kbT)),axis=1))
    #F=F-F[0]
    F[F==inf] = 0.0
    return F




def wham_iteration(F):
    P=calc_p(F)
    F=calc_f(P)
    if np.isnan(np.array(P).any()):
        exit("Nan")
    dF = np.sum(np.abs(np.subtract(F,F_old)))
    print(dF,np.shape(P))
    isconverged=False
    if (dF < tolerance):
        isconverged=True
        np.savetxt("wham_pmf.out", np.c_[xvals,-kbT*np.log(P)-np.min(-kbT*np.log(P)),P])
        np.savetxt("wham_F.out", np.c_[F-np.min(F)])
    return F, isconverged


parser=argparse.ArgumentParser()
parser.add_argument('-Nw', default=30,type=int, help='Number of Windows')
parser.add_argument('-rlow', default=1.5,type=float, help='Low cutoff')
parser.add_argument('-rhigh', default=8,type=float, help='High cutoff')
parser.add_argument('-nbin', default=100, type=int, help='Histogram bins')
parser.add_argument('-k', default=11.0, type=float, help='Force constant')
parser.add_argument('-plot', default=0, type=int, help='If 1, plots histograms in matplotlib')
parser.add_argument('-subfile', default="lif.distance", type=str, help="File name for colvar in subdirectory")
parser.add_argument('-unit', default=1, type=int, help="[0] Angstrom [1] bohr (output always angstrom, kcal/mol)")
args = parser.parse_args()

Nwindows=args.Nw
rlow=args.rlow
rhi=args.rhigh
nbins=args.nbin
k=args.k
unit=args.unit
subfile=args.subfile
shouldplot=args.plot

# Convert Units
BohrToAng=0.529177
convdist=1.0
convk=1.0
if unit==1:
    print("Converting r,k to units of angstroms, kcal/mol/ang^2 from bohr, hartree/bohr^2")
    convdist=BohrToAng
    convk=627.509/(BohrToAng**2.)

k=k*convk
rlow=rlow*convdist
rhi=rhi*convdist
print("New low: %s New hi: %s" % (rlow, rhi))
print("new k: %s" % k)

ni=[]
cnt=[]
center=[]
F_old=np.zeros(Nwindows)
xvals=np.linspace(rlow,rhi,num=nbins)
xF=np.linspace(rlow,rhi,num=Nwindows)
U=[]
xc=np.genfromtxt("wham_metadata.info",usecols=1,unpack=True)
for window in range(Nwindows):
    data=np.genfromtxt(str(window)+"/"+subfile, usecols=1,unpack=True)
    data=data*convdist
    hist,bins=np.histogram(data,bins=nbins,range=(rlow,rhi),density=False)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
    ni.append(hist)
    cnt.append(np.sum(hist))
    dx=np.subtract(xc[window],xvals)
    U.append(calc_bias(dx,k))

U=np.array(U)

if shouldplot==1: plt.show()
isconverged = False
tolerance = 1e-8
maxiter=10000
iteration=0
#does the wham iterations
while ( isconverged == False):
    F_old,isconverged=wham_iteration(F_old)
    iteration+=1
    if iteration > maxiter:
        exit("Too many iterations")
    
    

     



