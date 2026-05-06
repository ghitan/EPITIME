import numpy as np
from NSFD_behavioural import NSFD_behavioural
from plotData import plotData

# =============================== Input Parameters ===============================
T = 3000
h = 0.1
N = 5 * 10**7
mu = 1 / (75 * 365)

# =============================== Infectivity function (trapezoidal SARS) ===============================

def kernel_sars(tau):
    ta, tb, tc, td = 4, 7, 11, 14

    tau = np.array(tau)
    A = np.zeros_like(tau, dtype=float)

    mask1 = (tau > ta) & (tau < tb)
    mask2 = (tau >= tb) & (tau <= tc)
    mask3 = (tau > tc) & (tau < td)

    A[mask1] = np.exp(-mu * tau[mask1]) * (tau[mask1] - ta) / (tb - ta)
    A[mask2] = np.exp(-mu * tau[mask2])
    A[mask3] = np.exp(-mu * tau[mask3]) * (td - tau[mask3]) / (td - tc)

    return A

# =============================== Infectivity function ===============================
R0 = 3.3
R0_scaled = R0 / N
tp = np.linspace(4, 14, 100)
A_vals = kernel_sars(tp)
Ahat0 = np.trapezoid(A_vals, tp)
p0 = R0_scaled / Ahat0
A = lambda x: p0 * kernel_sars(x)

# =============================== Memory kernel ===============================
a = 1 / 30
K = lambda x: (a**2) * x * np.exp(-a * x)

# =============================== Message function ===============================
g = lambda x: x

# =============================== Inhibition function ===============================
alpha = 5 * 10**4
betaM = lambda x: 1 / (1 + alpha * x)

# =============================== Equilibria ===============================

Feq = mu / N * (R0 - 1) / (1 + alpha * mu)
lam = mu
Seq = lam * (1 / N + alpha * Feq) / (mu * (R0 / N))

# =============================== Initial data ===============================
S0 = 0.7 * N

# =============================== Run NSFD_behavioural ===============================
t, Y, P = NSFD_behavioural(
    T, h, N, mu, S0, A, K,
    betaM, g, 1
)

# =============================== Extract Variables ===============================
S = Y[0, :]
F = Y[1, :]

# =============================== Plot ===============================
plotData(t, F, Feq)
