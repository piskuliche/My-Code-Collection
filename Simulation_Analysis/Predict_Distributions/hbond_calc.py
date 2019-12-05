import numpy as np
import pickle
import os,time, argparse
from scipy import stats
from scipy.spatial import ConvexHull, Voronoi
from scipy.spatial.distance import pdist, squareform

def find_closest(dr_arr, drop):
    """
    This finds the minimum values
    """
    m,n=dr_arr.shape
    mask=np.ones((m,n),dtype=bool)
    mask[range(m),drop]=False
    finarray=dr_arr[mask].reshape(m,n-1)
    minval=finarray.min(axis=1)
    minh,mino=np.where(dr_arr==finarray.min(axis=1,keepdims=True))
    return minval, mino



def distance_wrap(r):
    """
    This calculates the distances with periodic boundaries
    """
    dr=np.sqrt(squareform(pdist(r[:,0])-np.round(pdist(r[:,0])))**2. + squareform(pdist(r[:,1])-np.round(pdist(r[:,1])))**2. + squareform(pdist(r[:,2])-np.round(pdist(r[:,2])))**2.)
    return dr

def calc_hbonds(frame):
    """
    This calculates the oo, oh, and angle distributions for hydrogen bond exchanges.
    """
    roo,roh,cosang=[],[],[]
    natoms = len(frame["type"])

    # Masks dr to be based on hatoms and all others
    hdr=np.array(frame["dr"][frame["h"]])
    hodr=np.array(hdr[:,frame["o"]])
    # Masks dr ot be based on oatoms 
    odr=np.array(frame["dr"][frame["o"]])

    # This is an array that is 686 long and tells what  each hatom is closest to
    oatms=(np.array(frame["co"])/3)[frame["h"]].astype(int)

    # Calculates distances and angle.
    sides_roh,closest=find_closest(hodr,oatms)
    sides_rho=hodr.min(axis=1)
    oarr=odr[:,np.array(frame["o"])]
    sides_roo=oarr[oatms,closest]
    roh = np.multiply(sides_roh,frame["L"])
    rho = np.multiply(sides_rho,frame["L"])
    roo = np.multiply(sides_roo,frame["L"])
    #Law of cosines to get the hbond angle
    #cos_theta=(side_rho**2.+side_roo**2.-side_roh**2.)/(2*side_rho*side_roo)
    cosang = np.divide(np.power(rho,2)+np.power(roo,2)-np.power(roh,2),np.multiply(2,np.multiply(rho,roo)))
    # Does histogramming of distances and angles
    histoh,bins=np.histogram(roh,bins=200,range=(1.0,4.0),density=False)
    histoh = histoh/len(roh)
    histoo,bins=np.histogram(roo,bins=200,range=(1.0,4.0),density=False)
    histoo = histoo/len(roo)
    histang, bins=np.histogram(np.arccos(cosang), bins=200, range=(0.0,1.0), density=False)
    histang = histang/len(cosang)
    return np.asarray(histoo),np.asarray(histoh),np.asarray(histang)
