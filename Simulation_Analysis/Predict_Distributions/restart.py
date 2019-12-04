import numpy as np
import pickle


def read_pkl(params,resno):
    """
    This reads the pickled binary restart file from a subdirectory
    """
    with open(str(resno)+"/restart_distribution.pkl","rb") as g:
        Res_Data=pickle.load(g)
    return Res_Data

def read_restart(params):
    """
    This loops over subdirectories, and then for every key, adds that to a data directory.
    """
    data = read_pkl(params,0)
    for r in range(1,params["numR"]):
        tmp = read_pkl(params,r)
        for key in data:
            for item in tmp[key]:
                data[key].append(item)
    return data
        
