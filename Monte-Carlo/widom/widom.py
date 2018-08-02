import numpy as np

class system:
    def __init__(self,nmols, rc):
        self.nmols = nmols
        self.x = []
        self.y = []
        self.z = []
        self.q = []
        self.type = []
        self.id = []
        self.molid = []
        self.eps = np.zeros((nmols, nmols))
        self.sig = np.zeros((nmols, nmols))
        self.L = L
        self.rc = rc

    def init_frame(self, file):
        self.id, self.molid, self.type, self.q, self.x, self.y, self.z = np.genfromtxt(file, usecols=(0,1,2,3,4,5,6), unpack=True)

    # Mixing Rules
    def geo_mix(a,b):
        out = sqrt(a*b)
        return out

    def ari_mix(a,b):
        out = (a+b)/2.0
        return out

    def pair_mix(self, eps, sig):
        for i in range(len(eps)):
            for j in range(len(eps)):
                self.eps[i][j] = self.geo_mix(eps[i],eps[j])
                self.sig[i][j] = self.ari_mix(sig[i],sig[j])

    # Read in parameters from file [type] [eps] [sig]
    def read_params(self, paramfile):
        eps = []
        sig = []
        eps, sig = np.genfromtxt(paramfile, usecols = (1,2), unpack=True)
        self.pair_mix(eps,sig)
    # Function for periodic boundary conditions
    def pbc(self,dxyz):
        dxyz = dxyz - self.L*int(round(dxyz/L))

    def calc_dist(self, atomi, atomj):
        dx = x[atomi]-x[atomj]
        dy = y[atomi]-y[atomj]
        dz = z[atomi]-z[atomj]
        self.pbc(dx)
        self.pbc(dy)
        self.pbc(dz)
        r = sqrt(dx**2 + dy**2 + dz**2)
        return r

    def calc_lj(self,dr,atomi,atomj):
        epsilon = self.eps[self.type[atomi]][self.type[atomj]]
        sigma = self.sig[self.type[atomi]][self.type[atomj]]
        rfact = (sigma/dr)**6
        U = 4*epsilon*(rfact**2 - rfact)
        return U
    # Ewald Sum Section
    def calc_elec():
        """

        """
        Vr = self.elec_realsum()
        Vf = self.elec_recipsum()
        Vs = self.selfsum()

    def elec_realsum():
        Vr = 0.0
    def elec_recipsum():
        Vf = 0.0
    def elec_selfsum():
        Vs = 0.0
    def calc_energy(self):
        U = 0.0
        for atomi in range(nmols):
            atomj = 0
            while atomj < atomi:
                rij= self.calc_dist(atomi, atomj)
                if rij > self.rc:
                    U += 0.0
                else:
                    U += calc_lj(rij, atomi, atomj)

                atomj += 1

def
