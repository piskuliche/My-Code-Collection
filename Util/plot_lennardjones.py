def lennard_jones(r, epsilon, sigma):
    rat = pow(sigma/r,6)
    lj = 4 * epsilon * (rat**2. - rat)
    return lj

if __name__ == "__main__":
    import sys
    import numpy as np
    import matplotlib.pyplot as plt
    if "-h" in sys.argv or len(sys.argv) == 1:
        print("Usage python plot_lj.py -n numplots -eps epsval1 .. -sig sigval1 .. -rmin minr")
        sys.exit()
    if "-eps" not in sys.argv or "-sig" not in sys.argv:
        print("Usage directions with -h")
        sys.exit()
    if "-n" in sys.argv:
        index = int(sys.argv.index("-n")+1)
        n = int(sys.argv[index])
    else:
        n=1
    eps, sig = [], []
    if "-eps" in sys.argv:
        for i in range(n):
            index = int(sys.argv.index("-eps")+i+1)
            eps.append(float(sys.argv[index]))
    if "-sig" in sys.argv:
        for i in range(n):
            index = int(sys.argv.index("-sig")+i+1)
            sig.append(float(sys.argv[index]))
    if "-rmin" in sys.argv:
        index = int(sys.argv.index("-rmin")+1)
        rmin = float(sys.argv[index])
    else:
        rmin = sig[0]*0.95
    assert len(eps) == len(sig), "Must have same number of parameters"

    r = np.arange(rmin,2.0*max(sig),0.1)
    plt.ylim((-1.5*max(eps),0.2))
    for i in range(len(eps)):
        plt.plot(r,lennard_jones(r,eps[i],sig[i]))
    plt.show()
    

