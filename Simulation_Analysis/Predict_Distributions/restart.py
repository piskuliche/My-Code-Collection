import numpy as np
import pickle
from post_calculation import manipulate_data


def read_pkl(params,resno,ftype):
    """
    This reads the pickled binary restart file from a subdirectory
    """
    with open(str(resno)+"/restart_distribution."+str(ftype)+".pkl","rb") as g:
        Res_Data=pickle.load(g)
    return Res_Data

def read_restart(params,ftype):
    """
    This loops over subdirectories, and then for every key, adds that to a data directory.
    """
    # Reads in the original data
    data = read_pkl(params,0,ftype)
    # Loops over the number of different directories the calculation has been subdivided over.
    for r in range(1,params["numR"]):
        tmp = read_pkl(params,r,ftype)
        for key in data:
            tmp[key]=np.array(tmp[key])
            # loops over items within each key, with options depending on what the type of thing is being
            # read back in.
            for item in tmp[key]:
                if (ftype=="data"): data[key].append(item)
                else: data[key]=np.append(data[key],item)
    print("After restart: there are %d total %s points" % (int(len(data[key])),ftype))
    return data

def call_calculation(params):
    """
    This calls the file read of the pkl files, and calculates output properties.
    """
    compiled_data=read_restart(params,"data")
    compiled_ener=read_restart(params,"energy")
    manipulate_data(params,compiled_data,compiled_ener)
    return 


        
