import numpy as np


def lj(eps,sig,r):
    rfact=(sig/r)**6
    U = 4*eps*(rfact**2 - rfact)
    return U




