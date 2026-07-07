import numpy as np
import matplotlib.pyplot as plt
from NSFD_AoI import NSFD_AoI

# TEST 1 of the Paper: EXPERIMENTAL ORDER OF CONVERGENCE (part of EPITIME software)
#
# E. Messina, M. Pezzella, A. Vecchio.
# "A non-standard numerical scheme for an age-of-infection epidemic model"
# Journal of Computational Dynamics, 2022, 9(2): 239–252.
# DOI: 10.3934/jcd.2021029
#
# ==========================================================================
# Version:    1.0
# Date:       April, 2026
# Authors:    M. Pezzella, G. Zanghirati
# See:        https://arxiv.org/abs/2605.00067v1
# License:    GNU GPL v.3
# Repository: https://github.com/ghitan/EPITIME

print("\nEPITIME AoI: test 1: experimental order of convergence\n")

#====================|
# Problem definition |
#====================|
AoI_problem = {
    'T': 1,
    'h': 2**-16,
    'A': lambda t: (1 + t) ** -2,
    'N': 10.0,
    'beta': 1e-3,
    'S0': 0.9 * 10.0,
    'verbosity': 1
}
AoI_problem['phi0'] = lambda t: (AoI_problem['N'] - AoI_problem['S0']) * AoI_problem['A'](t)

h_ref = AoI_problem['h']

#====================|
# Reference Solution |
#====================|
print("Computing the reference solution...", end = "")
_, Y_ref, P_ref = NSFD_AoI(problem=AoI_problem)
print(f"done.\nThe reference solution was computed in {P_ref['elapsed_time']:.3f} seconds with a stepsize h = {h_ref:.2e}.")

#================================|
# Experimental Convergence Study |
#================================|
start_expn = 6
TestN = int(-np.log2(h_ref) - 2 - start_expn)
H = []
S_err = []
Phi_err = []

# Prepare the running echo to console
Ndisp = max(1., TestN / 10.) * np.arange( min(TestN, 10) + 1. )
k = 1
print("\nConvergence study started...\n0%..", end = "")

for j in range(1, TestN + 1):
    h_test = 2**(- j - start_expn)
    H.append(h_test)
    AoI_test = AoI_problem.copy()
    AoI_test['h'] = h_test
    _, Y_test, _ = NSFD_AoI(problem=AoI_test)
    zooming = np.arange(0, np.size(Y_ref,1), int(h_test / h_ref))
    S_err.append(np.linalg.norm(Y_ref[0, zooming] - Y_test[0, :]) / np.linalg.norm(Y_test[0, :]))
    Phi_err.append(np.linalg.norm(Y_ref[1, zooming] - Y_test[1, :]) / np.linalg.norm(Y_test[1, :]))
    del Y_test
    # Display running iterations percentage
    if ( j >= Ndisp[k] ):
        print(f"{np.int32(100. * Ndisp[k] / TestN):d}%..", end = "")
        k += 1

print("done.\nConvergence study finished.")

H = np.array(H)
S_err = np.array(S_err)
Phi_err = np.array(Phi_err)

OrderS = np.zeros(TestN)
OrderPhi = np.zeros(TestN)
OrderS[1:] = np.log2(S_err[:-1] / S_err[1:])
OrderPhi[1:] = np.log2(Phi_err[:-1] / Phi_err[1:])

#==============================================|
# Relative Errors and Experimental Order Table |
#==============================================|
print(f"\n{'h':<10} {'RelError_S':<15} {'Order_S':<10} {'RelError_Phi':<15} {'Order_Phi':<10}")
for j in range(TestN):
    print(f"2^-{j+1:<7} {S_err[j]:<15.3e} {OrderS[j]:<10.2f} {Phi_err[j]:<15.3e} {OrderPhi[j]:<10.2f}")
    
#==============================|
# Plots of the Relative Errors |
#==============================|
shift = 1 / (H[-1] * 2)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
ax1.loglog(H, S_err, 'o-b', linewidth=2)
ax1.loglog(H, H * S_err[-1] * shift, ':b', linewidth=1.5)
ax1.set_xlabel('h')
ax1.set_ylabel('E(h)')
ax1.set_title('Relative Error on S(t)')
ax1.set_aspect('equal', adjustable='box')
ax2.loglog(H, Phi_err, 's-c', linewidth=2)
ax2.loglog(H, H * Phi_err[-1] * shift, ':c', linewidth=1.5)
ax2.set_xlabel('h')
ax2.set_ylabel('E(h)')
ax2.set_title('Relative Error on φ(t)')
ax2.set_aspect('equal', adjustable='box')
plt.show()

print("\nDone.\n\n")

# ==============================================================================
# End of Test1_Paper_JCD.py
# ==============================================================================
