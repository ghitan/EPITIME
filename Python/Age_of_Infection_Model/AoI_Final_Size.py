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

from NSFD_AoI import NSFD_AoI
import numpy as np
from scipy.stats import norm
from scipy.integrate import quad
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

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
plt.hlines(y=S_infty,xmin=t_num[0], xmax=t_num[-1], colors='darkcyan', linestyles='--', linewidth=2, label='S(∞)')
plt.xlabel('t')
plt.legend()
plt.title(f'NSFD Numerical Solution with h = {AoI_plot["h"]:.2e}')
plt.show()


#========================================================================|
# Comparison with Trapezoidal Direct Quadrature (Structure-Preservation) |
#========================================================================|
from scipy.optimize import fsolve

def Trap_AoI(S0, T, h, A, Phi0, beta):

    N = int(np.ceil(T / h))
    t = np.linspace(0.0, T, N + 1)

    Y = np.zeros((2, N + 1))

    phi0 = Phi0(0.0)
    sp0 = S0 * phi0
    A0 = A(0.0)

    Y[:, 0] = [S0, phi0]

    for n in range(N):
        tmp1 = 0.0
        tmp2 = 0.0
        if n >= 2:
            for j in range(1, n):

                Sj = Y[0, j]
                Phij = Y[1, j]
                prod = Sj * Phij

                tmp1 -= 2.0 * prod
                tmp2 += 2.0 * A(t[n-j-1]) * prod
        tmp = np.array([
            S0,
            Phi0(t[n])
        ]) + (beta * h / 2.0) * np.array([
            -sp0 + tmp1,
             A(t[n]) * sp0 + tmp2
        ])
        def FUN(X):
            return np.array([
                X[0] + (beta * h / 2.0) * X[0] * X[1] - tmp[0],
                X[1] - (beta * h / 2.0) * A0 * X[0] * X[1] - tmp[1]
            ])
        sol = fsolve(
            FUN,
            Y[:, n],
            xtol=1e-12,
            maxfev=2000
        )
        Y[:, n + 1] = sol
    return t, Y

AoI_problem['T'] = 8
AoI_problem['h'] = 2.5e-1
AoI_problem['beta'] = 6e-5

#====================================|
# NSFD and Trapezoidal DQ solutions  |
#====================================|
t, Y_nsfd, _ = NSFD_AoI(problem=AoI_problem)

t_trap, Y_trap = Trap_AoI(
    AoI_problem['S0'],
    AoI_problem['T'],
    AoI_problem['h'],
    AoI_problem['A'],
    AoI_problem['phi0'],
    AoI_problem['beta']
)

#=============================|
# Comparison plots            |
#=============================|
fig, ax = plt.subplots(1, 2, figsize=(10, 4))

# S(t)
ax[0].plot(
    t,
    Y_nsfd[0],
    linewidth=2,
    color=(0.96, 0.47, 0.16),
    label='NSFD'
)

ax[0].plot(
    t_trap,
    Y_trap[0],
    linewidth=2,
    color=(0.51, 0.25, 0.09),
    label='TrapDQ'
)

ax[0].axhline(0.0, linestyle=':')
ax[0].set_xlabel('t')
ax[0].set_ylabel('S(t)')
ax[0].legend(loc='best')
ax[0].set_box_aspect(1)

# φ(t)
ax[1].plot(
    t,
    Y_nsfd[1],
    linewidth=2,
    color='g',
    label='NSFD'
)

ax[1].plot(
    t_trap,
    Y_trap[1],
    linewidth=2,
    color='lightseagreen',
    label='TrapDQ'
)

ax[1].axhline(0.0, linestyle=':')
ax[1].set_xlabel('t')
ax[1].set_ylabel('φ(t)')
ax[1].legend(loc='best')
ax[1].set_box_aspect(1)

fig.suptitle('Numerical Solutions for h = 0.25')

plt.tight_layout()
plt.show()

# ==============================================================================
# End of AoI_Final_Size.py
# ==============================================================================

