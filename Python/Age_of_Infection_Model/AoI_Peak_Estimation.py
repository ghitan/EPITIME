# EPIDEMIC PEAK INVESTIGATION (part of EPITIME software)
#
# ==========================================================================
# Version:    1.0
# Date:       June, 2026
# Authors:    M. Pezzella, G. Zanghirati
# See:        https://arxiv.org/abs/2605.00067v1
# License:    GNU GPL v.3
# Repository: https://github.com/ghitan/EPITIME

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from NSFD_AoI import NSFD_AoI

# ================================================================
# Parameters
# ================================================================
N = 100

i0 = np.array([0.0268, 1.27e-6, 0.05, 1.33e-4, 5e-6, 0.001])
s0 = 1 - i0

betaN = np.array([4.6291, 0.5, 10, 0.13905, 0.3333, 0.5])
gamma = np.array([2.82, 0.3, 1, 0.018379, 0.1111, 0.3])

N_cases = len(i0)

peak_Continuous = np.zeros(N_cases)
t_peak = np.zeros(N_cases)
E_peak = np.zeros(N_cases)

# Prepare the running echo to console
Ndisp = max(1., N_cases / 10.) * np.arange( min(N_cases, 10) + 1. )
k = 1
print( "\nAoI peak estimation started...\n0%..", end = '' )


# ================================================================
# Main loop (peak computation)
# ================================================================
for j in range(N_cases):

    def A(t):
        return np.exp(-gamma[j] * t)

    AoI_problem = {
        "S0": s0[j] * N,
        "N": N,
        "T": 100,
        "h": 1e-3,
        "beta": betaN[j] / N,
        "A": A,
        "phi0": lambda t: N * i0[j] * A(t),
        "verbosity": 0
    }

    # R0
    R0 = AoI_problem["beta"] * N * quad(A, 0, np.inf)[0]

    s_tmp = s0[j]
    i_tmp = i0[j]
    sm = 1 / R0

    def peak_integr(s):
        return (1 / betaN[j]) * 1.0 / (
            s * (s - sm * np.log(s / s_tmp) - i_tmp - s_tmp)
        )

    peak_Continuous[j] = quad(
        peak_integr, s_tmp, sm,
        epsabs=1e-12, epsrel=1e-12
    )[0]

    # NSFD solve
    t, Y, P = NSFD_AoI( problem = AoI_problem )

    incidence = AoI_problem["beta"] * Y[0, :] * Y[1, :]

    idx = np.argmin(np.abs(Y[0, :] / N - 1 / R0))
    t_peak[j] = t[idx]
    E_peak[j] = AoI_problem["beta"] * Y[1, idx]

    # Display running iterations percentage
    if ( j >= Ndisp[k] ):
        print( f"{ np.int32(100. * Ndisp[k] / N_cases) :d}%..", end = '' ); k += 1

print( "100%\nPeak estimation finished.\n" )

# ================================================================
# Convergence study
# ================================================================
h_vec = np.array([
    1e-1, 5e-2, 2.5e-2, 1.25e-2,
    6.25e-3, 3.125e-3, 1.5625e-3
])

Nh = len(h_vec)

t_peak_h = np.zeros((Nh, N_cases))
err_tp = np.zeros((Nh, N_cases))
eoc_tp = np.full((Nh, N_cases), np.nan)

# Prepare the running echo to console
# Ndisp = max(1., N_cases / 10.) * np.arange( min(N_cases, 10) + 1. )
kk = 1
print( "\nAoI convergence study started...\n0%..", end = '' )

for j in range(N_cases):

    def A(t):
        return np.exp(-gamma[j] * t)

    R0 = betaN[j] * quad(A, 0, np.inf)[0]

    tp_exact = peak_Continuous[j]

    for k in range(Nh):

        AoI_problem = {
            "S0": s0[j] * N,
            "N": N,
            "T": 100,
            "h": h_vec[k],
            "beta": betaN[j] / N,
            "A": A,
            "phi0": lambda t: N * i0[j] * A(t),
            "verbosity": 0
        }

        t, Y, P = NSFD_AoI(problem=AoI_problem)

        idx = np.argmin(np.abs(Y[0, :] / N - 1 / R0))

        t_peak_h[k, j] = t[idx]
        err_tp[k, j] = abs(t_peak_h[k, j] - tp_exact) / tp_exact

        if k > 0:
            eoc_tp[k, j] = np.log(err_tp[k-1, j] / err_tp[k, j]) / np.log(2)

    # Display running iterations percentage
    if ( j >= Ndisp[kk] ):
        print( f"{ np.int32(100. * Ndisp[kk] / N_cases) :d}%..", end = '' ); kk += 1

print( "100%\nConvergence study finished.\n" )

# ================================================================
# Plot (log-log convergence)
# ================================================================
plt.figure()

colors = plt.cm.jet(np.linspace(0, 1, N_cases))

for j in range(N_cases):
    plt.loglog(
        h_vec,
        err_tp[:, j],
        '-o',
        linewidth=1.5,
        markersize=5,
        color=colors[j]
    )

# reference slope O(h)
C = 0.5 * err_tp[-1, 0] / h_vec[-1]
plt.loglog(h_vec, C * h_vec, 'k--', linewidth=1.5)

plt.xlabel('h')
plt.ylabel('Relative error of $t_p$')
plt.title('Convergence of epidemic peak time')

plt.grid(True, which='both')
plt.legend([f'Case {j+1}' for j in range(N_cases)] + [r'$\mathcal{O}(h)$'])

plt.tight_layout()
plt.show()

# ==============================================================================
# End of AoI_Peak_Estimation.py
# ==============================================================================
