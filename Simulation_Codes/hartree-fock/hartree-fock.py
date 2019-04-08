import numpy as np
"""
Calculates all the integrals for the scf calculation
"""
def integrl(iop,n,r,zeta1,zeta2,za,zb):
    coef  = np.zeros(3,3)
    expon = np.zeros(3,3)
    d1    = np.zeros(3)
    a1    = np.zeros(3)
    d2    = np.zeros(3)
    a2    = np.zeros(3)
    pi    = math.pi
    """
    These are the contraction coefficients and exponents for a normalized
    slater orbital with exponent 1.0 in terms of normalized 1s primative gaussians
    """
    coef =  {1:     {1: 1.0,        2: 0.678914,    3: 0.444635},
             2:     {1: None,       2: 0.430129,    3: 0.535328},
             3:     {1: None,       2: None,        3: 0.154329}}
    expon = {1:     {1: 0.270950,   2: 0.151623,    3: 0.109818},
             2:     {1: None,       2: 0.851819,    3: 0.405771},
             3:     {1: None,       2: None,        3: 2.22766}} 
    for i in range(0,n):
        a1[i] = expon[i][n] * (zeta1 **2)
        d1[i] = coef[i][n] * ((2.0 * a1[i] / pi) ** 0.75)
        a2[i] = expon[i][n] * (zeta2 **2)
        d2[i] = coef[i][n] * ((2.0 * a2[i] / pi) ** 0.75)
    r2 = r*r
    # Scale and normalize
    s12 = 0.0
    t11 = 0.0
    t12 = 0.0
    t22 = 0.0
    v11A = 0.0
    v12A = 0.0
    v22A = 0.0
    v11B = 0.0
    v12B = 0.0
    v22B = 0.0
    v1111 = 0.0
    v2111 = 0.0
    v2121 = 0.0
    v2211 = 0.0
    v2221 = 0.0
    v2222 = 0.0
    # Calcualte one electron integrals
    # Center A is atom1, B atom2
    # Origin on A
    # V12A off-diagonal nuclear attraction to A
    for i in range(0,n):
        for j in range(0,n):
            #rap2 sqared distance between a and p
            rap = a2[j] *r/(a1[i]+a2[j])
            rap2 = rap**2
            rbp22 = (r - rap)**2
            s12 = s12+s_func(a1[i],a2[j],r2)*d1[i]*d1[j]
            t11 = t11+t_func(a1[i],a1[j],0.0)*d1[i]*d2[j]
            t12 = t12+t_func(a1[i],a2[j],r2)*d1[i]*d2[j]
            t22 = t22+t_func(a2[i],a2[j],0.0)*d2[i]*d2[j]
            v11A = v11A+v_func(a1[i],a1[j],0.0,0.0,za)*d1[i]*d1[j]
            v12A = v12A+v_func(a1[i],a2[j],r2,rap2,za)*d1[i]*d2[j]
            v22A = v22A+v_func(a2[i],a2[j],0.0,0.0,za)*d2[i]*d2[j]
            v11B = v11B+v_func(a1[i],a1[j],0.0,0.0,zb)*d1[i]*d1[j]
            v12B = v12B+v_func(a1[i],a2[j],r2,rap2,zb)*d1[i]*d2[j]
            v22B = v22B+v_func(a2[i],a2[j],0.0,0.0,zb)*d2[i]*d2[j]
    # Calculate two-electron integrals
    for i in range(0,n):
        for j in range(0,n):
            for k in range(0,n):
                for l in range(0,n):



"""
Does a Hartree-Fock Calculationf or a two-electron diatomic using the 1s minimal
STO-NG Basis Set
Minimal basis set has basis functions 1 and 2 on nuclei a and b

iop = 0 no printing
iop = 1 print only converged results
iop = 2 print every iteration
n       STO-NG Caclualtion (1,2 or 3)
r       Bond Length (au)
zeta1   Slater orbital exponent (function1)
zeta2   Slater orbital exponent (function2)
za      Atomic number atomA
zb      Atomic Number atomB
"""
def HFCALC(iop,n,r,zeta1,zeta2,za,zb):
    if iop != 0:
        print("STO-%sG FOR ATOMIC NUMBERS %s AND %s" %(n,za,zb))
    # Calculate all the one and two-electron integrals
    integrl(iop, n, r, zeta1, zeta2, za, zb)
    # Be inefficient and put all the integrals in pretty arrays
    colect(iop,n,r,zeta1,zeta2,za,zb)
    # Perform the SCF calculation
    scf(iop,n,r,zeta1,zeta2,za,zb)
    


"""
Minimal BASIS STO-3G CALCULATION ON HEH+
adapted from Szabo and Ostlund
"""
iop     = 2
n       = 3 
r       = 1.4632
zeta1   = 2.0925
zeta2   = 1.24
za      = 2.0
zb      = 1.0
HFCALC(iop,n,r,zeta1,zeta2,za,zb)



