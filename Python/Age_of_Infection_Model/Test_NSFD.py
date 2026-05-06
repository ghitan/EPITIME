import numpy as np
import matplotlib.pyplot as plt
from NSFD_AoI import NSFD_AoI 

# Test_NSFD - Prototype script to run age-of-infection model esperiments with 
#             the NSFD scheme (part of EPITIME software)
#
# ==========================================================================
# Version:    1.0
# Date:       April, 2026
# Authors:    M. Pezzella, G. Zanghirati
# See:        https://arxiv.org/abs/2605.00067v1
# License:    GNU GPL v.3
# Repository: https://github.com/ghitan/EPITIME


def run_test():
    problem = {
        'S0': 99/100*1000,
        'N': 1000,
        'T': 15.0,
        'h': 1e-2,
        'beta': 1e-3,
        'A': lambda t: 2/(1+1*t**2),
        'phi0': lambda x: (N - S0)*A(x),
        'verbosity': 1
    }
    
    T = 15.0
    h = 1e-2
    A = lambda t: 2/(1+1*t**2)
    N = 1e3
    beta = 1e-3
    S0 = 99/100*N
    phi0 = lambda x: (N - S0)*A(x)
    verbosity = 1

    t, Y, P = NSFD_AoI(S0=S0, T=T, h=h, A=A, phi0=phi0, beta=beta, verbosity=verbosity)
    
    plt.figure()
    plt.subplot(1, 2, 1)
    plt.plot(t, Y[0, :], linewidth = '2')
    plt.title("Susceptibles")
    plt.xlabel("t")
    plt.ylabel("S(t)")
    plt.subplot(1, 2, 2)
    plt.plot(t, Y[1, :], linewidth = '2')
    plt.title("Total Infectivity")
    plt.xlabel("t")
    plt.ylabel("φ(t)")
    plt.show()
    
    return t, Y, P

if __name__ == "__main__":
    t, Y, P = run_test()