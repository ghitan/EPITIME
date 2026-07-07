import numpy as np
import matplotlib.pyplot as plt
from NSFD_behavioural import NSFD_behavioural

# ==== Behavioural epidemic peaks investigation (part of EPITIME software) ====
#
# B. Buonomo, E. Messina, C. Panico, A. Vecchio.
# "An integral renewal equation approach to behavioural epidemic models with 
#  information index"
# Journal of Mathematical Biology, 2025, 90(8)
# DOI : 10.1007/s00285-024-02172-y
#
# ==============================================================================
# Version:    1.0
# Date:       April, 2026
# Authors:    C. Panico, G. Zanghirati
# See:        https://arxiv.org/abs/2605.00067v1
# License:    GNU GPL v.3
# Repository: https://github.com/ghitan/EPITIME

print("\nEPITIME behavioural: epidemic peaks investigation\n")

# =========================================================
# Input Parameters 
# =========================================================
T  = 1000
h  = 0.1
N  = 5 * (10**7)
mu = 1.0 / (75 * 365)

# =========================================================
# Functions
# =========================================================
ni = 1 / 7
R0 = 20
bb = (R0 / N) * (mu + ni)**2

A = lambda x: bb * x * np.exp(-(mu + ni) * x)
a = 1 / 30
K = lambda x: a * np.exp(-a * x)
g = lambda x: x
alpha = 8 * (10**3)
betaM = lambda x: 1 / (1 + alpha * x)

# =========================================================
# Equilibria 
# =========================================================
Feq = (mu / N) * (R0 - 1) / (1 + alpha * mu)
lam = mu
Seq = lam * (1 / N + alpha * Feq) / (mu * (R0 / N))

# =========================================================
# Initial data
# =========================================================
S0  = 0.99 * Seq * N

# =========================================================   
# Solver execution
# =========================================================
t, Y, P = NSFD_behavioural(
    T, h, N, mu, S0, A, K,
    betaM, g, True, 1
    )

# =========================================================
# Variable Extraction 
# =========================================================
S_N = Y[0, :] * N 
F_N = Y[1, :] * N 
Nt = len(t)
incidence_NS = np.zeros(Nt - 1)

for n in range(Nt - 1): 
    D = h
    if (T - np.floor(T / h) * h != 0) and (n == Nt - 2):
        D = T - np.floor(T / h) * h
    
    incidence_NS[n] = N * ((Y[0, n] - Y[0, n + 1]) / h + mu * (1 - Y[0, n + 1]))

incidence_NS = np.append(incidence_NS, incidence_NS[-1])

# =========================================================
# NUMERICAL PEAK DETECTION 
# =========================================================
# 1. Peaks of the Driving Factor F(t)
f_diff1 = F_N[1:-1] > F_N[:-2]
f_diff2 = F_N[1:-1] > F_N[2:]
idx_f = np.where(f_diff1 & f_diff2)[0] + 1
pks_f = F_N[idx_f]
locs_f = t[idx_f]

valid_f = pks_f > (np.max(F_N) * 0.01)
pks_f = pks_f[valid_f]
locs_f = locs_f[valid_f]

# 2. Peaks of the Incidence
inc_diff1 = incidence_NS[1:-1] > incidence_NS[:-2]
inc_diff2 = incidence_NS[1:-1] > incidence_NS[2:]
idx_inc = np.where(inc_diff1 & inc_diff2)[0] + 1
pks_inc = incidence_NS[idx_inc]
locs_inc = t[idx_inc]

valid_inc = pks_inc > (np.max(incidence_NS) * 0.01)
pks_inc = pks_inc[valid_inc]
locs_inc = locs_inc[valid_inc]

# =========================================================
# Output
# =========================================================
print('\n=========================================================')
print('  RESULTS ')
print('=========================================================')

print('\n--- PEAKS OF DRIVING FACTOR F(t) ---')
if len(locs_f) == 0:
    print('No peaks detected.')
else:
    for i in range(len(locs_f)):
        print(f'Peak {i+1}: Time = {locs_f[i]:6.2f}, F(t) = {pks_f[i]:12.4e}')

print('\n--- PEAKS OF INCIDENCE ---')
if len(locs_inc) == 0:
    print('No peaks detected.')
else:
    for i in range(len(locs_inc)):
        print(f'Peak {i+1}: Time = {locs_inc[i]:6.2f}, Incidence(t) = {pks_inc[i]:12.2f}')
print('=========================================================')

# =========================================================
# Data Visualization
# =========================================================
fig = plt.figure('Epidemic Peaks')

# Subplot 1: Driving Factor F(t)
plt.subplot(2, 1, 1)
plt.plot(t, F_N, 'b-', linewidth=1.5, label='F(t)')
plt.plot(locs_f, pks_f, 'ro', markerfacecolor='r', markersize=6, label='Detected Peaks')
plt.grid(True)
plt.xlabel('t')
plt.ylabel('F(t)')
plt.legend(loc='best')

# Subplot 2: Incidence
plt.subplot(2, 1, 2)
plt.plot(t, incidence_NS, 'g-', linewidth=1.5, label='Inc(t)')
plt.plot(locs_inc, pks_inc, 'ko', markerfacecolor='k', markersize=6, label='Detected Peaks')
plt.grid(True)
plt.xlabel('t')
plt.ylabel('Inc(t)')
plt.legend(loc='best')

plt.tight_layout() 

# -----------------------------
# Vector graphics export
# -----------------------------
plt.savefig('peaks_beh.pdf', format='pdf', transparent=True, dpi=1200)

plt.show()

print("\nDone.\n\n")

# ==============================================================================
# End of Behavioural_Trapezoidal_IF.py
# ==============================================================================


