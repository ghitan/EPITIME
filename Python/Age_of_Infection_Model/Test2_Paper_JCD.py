from NSFD_AoI import NSFD_AoI
import numpy as np
from scipy.stats import norm
from scipy.integrate import quad
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# TEST 2 of the Paper: LONG TIME BEHAVIOUR INVESTIGATION (part of EPITIME software)
#
# E. Messina, M. Pezzella, A. Vecchio.
# "A non-standard numerical scheme for an age-of-infection epidemic model"
# Journal of Computational Dynamics, 2022, 9(2): 239–252.
# DOI: 10.3934/jcd.2021029
#
# ==============================================================================
# Version:    1.0
# Date:       April, 2026
# Authors:    M. Pezzella, G. Zanghirati
# See:        https://arxiv.org/abs/2605.00067v1
# License:    GNU GPL v.3
# Repository: https://github.com/ghitan/EPITIME

#====================|
# Problem definition |
#====================|
mu = 0.2
sigma = 2 * mu
A_func = lambda t: norm.pdf(t, mu, sigma)

AoI_problem = {
    'T': 30,
    'h': 1e-1,
    'A': A_func,
    'N': 1e5,
    'beta': 3e-5,
    'S0': 99950,
    'phi0': lambda t: (1e5 - 99950) * A_func(t),
    'verbosity': 0
}

TestN = 3
S_final = np.zeros(TestN)
Phi_final = np.zeros(TestN)
H = np.zeros(TestN)
R0_num = np.zeros(TestN)

betaN = AoI_problem['beta'] * AoI_problem['N']
for j in range(1, TestN + 1):
    h_test = 10**-j
    H[j - 1] = h_test
    AoI_test = AoI_problem.copy()
    AoI_test['h'] = h_test
    _, Y_test, _ = NSFD_AoI(problem=AoI_test)
    S_final[j - 1] = Y_test[0, -1]
    Phi_final[j - 1] = Y_test[1, -1]
    A_vals = AoI_problem['A'](np.arange(h_test, 1e4 + h_test, h_test))
    R0_num[j - 1] = betaN * h_test * np.sum(A_vals)
    del Y_test

#==============================|
# Continuous R₀ and Final Size |
#==============================|
R0_cont, _ = quad(AoI_problem['A'], 0, np.inf)
R0_cont *= betaN

FS_function = lambda x: np.log(AoI_problem['S0'] / x) - R0_cont * (1 - x / AoI_problem['N'])
S_infty = fsolve(FS_function, 0.1, xtol=1e-14)[0]
err_infty = np.abs(S_final - S_infty) / S_infty

# Print results
print(f"\nThe Continuous Basic Reproduction Number is R0 = {R0_cont:.2f}")
print(f"The Continuous Final Size is S(∞) = {S_infty:.5e}\n")

print(f"{'h':<10} {'R₀(h)':<15} {'S_∞(h)':<15} {'Rel. S_∞ Err.':<20} {'φ_∞(h)':<15}")
for j in range(TestN):
    print(f"10^-{j+1:<6} {R0_num[j]:<15.5e} {S_final[j]:<15.5e} {err_infty[j]:<20.5e} {Phi_final[j]:<15.5e}")

#============================================|
# Plot of the numerical solution  for h=1e-3 |
#============================================|
AoI_plot = AoI_problem.copy()
AoI_plot['h'] = 1e-3
t_num, Sol_num, _ = NSFD_AoI(problem=AoI_plot)

plt.figure()
plt.plot(t_num, Sol_num[0], linewidth=2, color='b', label='S(t)')
plt.plot(t_num, Sol_num[1], linewidth=2, color='c', label='φ(t)')
plt.scatter(t_num[-1], S_infty, linewidth=2, edgecolor='b', facecolor='b', label='S(∞)')
plt.xlabel('t')
plt.legend()
plt.title(f'NSFD Numerical Solution with h = {AoI_plot["h"]:.2e}')
plt.show()