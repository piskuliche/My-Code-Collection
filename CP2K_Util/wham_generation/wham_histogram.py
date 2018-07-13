import numpy as np
import matplotlib.pyplot as plt
import sys

basename=str(sys.argv[1])
nbins = int(sys.argv[2])

print("For %s bins, histogramming %s" % (nbins, basename))

x = []

for bin in range(nbins):
    print("Bin number %s completed" % bin)
    x.append([])
    x[bin]=np.genfromtxt(str(bin)+'/lif.distance',unpack=True, usecols=(1),comments='#')
    plt.hist(x[bin], 20, normed=1,histtype='step')

plt.xlabel('distance (bohr)')
plt.ylabel('Probability')
plt.show()
