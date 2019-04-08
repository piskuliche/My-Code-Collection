import numpy as np
import itertools
import matplotlib.pyplot as plt

class distances:
    def __init__(self,nframes):
        self.r = np.zeros((3,2,nframes))
        self.dr = np.zeros((3,nframes))
        self.drmag = np.zeros(nframes)
    def read_frames(self, nframes, infile):
        rx, ry, rz = np.genfromtxt(infile, usecols = (1,2,3), unpack=True,max_rows=nframes*2)
        self.r[0,0] = rx[::2]
        self.r[1,0] = ry[::2]
        self.r[2,0] = rz[::2]
        self.r[0,1] = rx[1::2]
        self.r[1,1] = ry[1::2]
        self.r[2,1] = rz[1::2]
    def calc_dist(self,nframes, L):
        for frame in range(nframes):
            self.dr[0,frame]=self.r[0,0,frame]-self.r[0,1,frame]
            self.dr[1,frame]=self.r[1,0,frame]-self.r[1,1,frame]
            self.dr[2,frame]=self.r[2,0,frame]-self.r[2,1,frame]
            self.dr[0,frame]=self.dr[0,frame]-L*int(round(self.dr[0,frame]/L))
            self.dr[1,frame]=self.dr[1,frame]-L*int(round(self.dr[1,frame]/L))
            self.dr[2,frame]=self.dr[2,frame]-L*int(round(self.dr[2,frame]/L))
            self.drmag[frame] = np.sqrt(self.dr[0,frame]**2 + self.dr[1,frame]**2 + self.dr[2,frame]**2)
            print frame, self.drmag[frame], self.dr[0,frame], self.dr[1,frame], self.dr[2,frame]

    def bin_dist(self,nframes, L):
        dr = 0.1
        rcut = L/2.0
        self.hist = np.histogram(self.drmag,range=(0.0,rcut), bins = int(rcut)*100,density=True)
        np.savetxt('out.hist', np.c_[self.hist[1][1:], self.hist[0]*(L**3/(4*np.pi*self.hist[1][1:]**2))])
        

nframes = 100001
infile = 'lif.xyz'
v = distances(nframes)
L = 14.219066

v.read_frames(nframes,infile)
v.calc_dist(nframes,L)
v.bin_dist(nframes, L)
