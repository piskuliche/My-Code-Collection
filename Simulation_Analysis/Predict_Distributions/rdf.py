import numpy as np

def init_rdf(params,n):
    const=4.0*np.pi*n/3.0
    h_id=[]
    r = 0.0
    while r < params["rmax"]-params["dr"]:
        h_id.append(const*(r+params["dr"])**3.-const*r**3.)
        r += params["dr"]
    return h_id

def remove_likes(dr_arr, drop):
    m,n=dr_arr.shape
    mask=np.ones((m,n),dtype=bool)
    mask[range(m),drop]=False
    finarray=dr_arr[mask].reshape(m,n-1)
    return finarray



def radial_distribution(params,frame):
    volume=frame["L"]**3.
    oodr=frame["dr"][frame["o"]][:,frame["o"]]
    ohdr=frame["dr"][frame["o"]][:,frame["h"]]
    hhdr=frame["dr"][frame["h"]][:,frame["h"]]
    goo,goh,ghh=0,0,0
    goo=np.histogram(oodr,bins=params["rdfbins"],range=(0.0,params["rmax"]/frame["L"]),density=False)[0]
    goh=np.histogram(ohdr,bins=params["rdfbins"],range=(0.0,params["rmax"]/frame["L"]),density=False)[0]
    ghh=np.histogram(hhdr,bins=params["rdfbins"],range=(0.0,params["rmax"]/frame["L"]),density=False)[0]
    goo = np.divide(goo,len(oodr))
    goo = np.divide(goo,params["h_o"]/volume)
    goh = np.divide(goh,len(ohdr))
    goh = np.divide(goh,params["h_h"]/volume)
    ghh = np.divide(ghh,len(hhdr))
    ghh = np.divide(ghh,params["h_h"]/volume)
    goo[0], goh[0], ghh[0]=0.0,0.0,0.0
    return np.asarray(goo), np.asarray(goh), np.asarray(ghh)





