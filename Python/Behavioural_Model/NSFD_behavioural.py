import numpy as np
import time
import warnings
def NSFD_behavioural(*args, **kwargs):
    """
NSFD_behavioural
____________________________________________________________________________________________
Weighted Non Standard Finite Difference (NSFD) scheme for the solution of the behavioural
integro-differential epidemic model with demographic turnover and information index

Version:    1.0  (part of EPITIME software)
Date:       April, 2026
Authors:    B. Buonomo,  E. Messina, C. Panico, M. Pezzella, G. Zanghirati
See:        https://arxiv.org/abs/2605.00067v1
License:    GNU GPL v.3
Repository: https://github.com/ghitan/EPITIME

********************************************************************************************
SYNOPSIS:
t, Y, P = NSFD_behavioural(*args, **kwargs)

********************************************************************************************
INPUT:
The function supports two alternative input modes:

1) The function has a structured input BEH_prob with the following optional fields:

   T             : Final time of the integration                    [Default: 1000]
   h             : Time step used for discretization                [Default: 1]
   mu            : Natural death rate                               [Default: 1/(75*365)]
   N             : Total population size                            [Default: 5e7]
   S0            : Initial susceptible population                   [Default: 0.2*6e7]
   A             : Infectivity function                             [Default: bb* t* exp(-(mu +ni).*t)]
   betaM         : Inhibition function                              [Default: 1/((1 / N) + alpha* x)]
   g             : Message function                                 [Default: x)]
   K             : Memory kernel                                    [Default: a * exp(-a * x)]
   verbosity     : 0 to suppress warnings, 1 to enable              [Default: 1]

   EXAMPLE: NSFD_behavioural(BEH_prob)

 2) The function has positional inputs:

   EXAMPLE: NSFD_behavioural(T, h, N, mu, S0, A, K, betaM, g, Verbosity)

OUTPUT:
t               : vector of time
S =Y[0,:]       : numerical approximation of the solution S(t)
F =Y[1,:]       : numerical approximation of the solution F(t)
P               : (optional) dict containing performance metrics and resource usage:
                  • elapsed_time   - total computation time [seconds]
                  • steps_number   - number of time steps or iterations (Nt)
                  • flops_number   - estimated number of floating-point operations
                  • memory_bytes   - total memory used by key variables [bytes]

_________________________________________________________________________________________
The Python implementation is organized into logical Code Blocks,
labelled PCB:A to PCB:F, each corresponding to a specific stage of the
solver. These labels are used consistently in the accompanying paper
to reference individual parts of the implementation.

PCB:A : Default parameters and input processing
PCB:B : Time grid construction and memory allocation
PCB:C : Input consistency and admissibility checks
PCB:D : NSFD time-stepping scheme
PCB:E : Performance metrics and diagnostics
PCB:F : Auxiliary routine 1
PCB:G : Auxiliary routine 2
_________________________________________________________________________________________

For further details we refer to the following papers
[1]  B. Buonomo, E. Messina, C. Panico, A. Vecchio.
      A stable numerical method for integral epidemic models with behavioral changes in
      contact patterns.
      Electronic Transactions on Numerical Analysis, 2024, 61, 137–156.
      DOI : 10.1553/etna_vol61s137

[2]  B. Buonomo, E. Messina, C. Panico, A. Vecchio.
      An integral renewal equation approach to behavioural epidemic models with
      information index
      Journal of Mathematical Biology, 2025, 90(8).
      DOI : 10.1007/s00285-024-02172-y
_________________________________________________________________________________________

DISCLAIMER

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
_________________________________________________________________________________________

"""

    # ------------------- PCB:A | Default parameters and input processing -------------------
    defaults = {
        'T': 1000,
        'h': 1,
        'N': 5e7,
        'mu': 1/(75*365),
        'S0': 0.2*5e7,
        'A': None,
        'K': None,
        'betaM': None,
        'g': None,
        'verbosity': 1
    }
    # ----- Input parsing new -----
    BEH_prob = parse_input(defaults, *args, **kwargs)

    # ------------------- PCB:B | Time grid construction and memory allocation -------------------
    start_time = time.time()
    T = BEH_prob['T']
    h = BEH_prob['h']
    N = BEH_prob['N']
    mu = BEH_prob['mu']
    S0 = BEH_prob['S0'] / N

    ni_t = int(np.floor(T / h))
    D = h
    D_tilde = 0
    t = np.arange(0, T+h, h)
    if t[-1] != T:
        D_tilde = T - ni_t * h
        t = np.append(t, T)
    Nt = len(t)

    # Memory allocation
    Y = np.zeros((2, Nt))
    M = np.zeros(Nt)
    betaM_vec = np.zeros(Nt)
    g_vec = np.zeros(Nt)

    # Initial conditions
    Y[0,0] = S0
    A_vec = BEH_prob['A'](t)
    K_vec = BEH_prob['K'](t)
    f_0_cost = (1 - S0/N) * mu
    F0_vec = f_0_cost * A_vec
    Y[1,0] = F0_vec[0]
    M0_vec = (1/N)*BEH_prob['g'](N*Y[1, 0]) * K_vec
    M[0] = M0_vec[0]

    allocation_time = time.time() - start_time

    # ------------------- PCB:C | Input consistency check -------------------
    _check_input(BEH_prob, defaults, A_vec)

    # ------------------- PCB:D | NSFD weight definition and time-stepping ------------
    start_weight = time.time()
    Sum1 = np.sum(K_vec)
    T_incr = T + h
    incr = h * BEH_prob['K'](T_incr)
    it_max = 0
    Nt_2 = Nt**2
    K_eval = Nt
    while incr > 1e-11 and it_max < Nt_2:
        Sum1 += incr
        T_incr += h
        incr = h * BEH_prob['K'](T_incr)
        it_max += 1
        K_eval += 1
    gammah = 1 / (h * Sum1)
    weight_definition_time = time.time() - start_weight

# =============================== Time iteration ===============================
    start_step = time.time()

    for n in range(Nt-1):
        if D_tilde != 0 and n == Nt-2:
            D = D_tilde

        betaM_val = BEH_prob['betaM'](N * M[n])
        g_val = BEH_prob['g'](N * Y[1, n])

        Y[0, n+1] = (Y[0, n] + mu*D) / (1 + D * (mu + N * betaM_val * Y[1, n]))

        betaM_vec[n] = N * betaM_val
        g_vec[n] = g_val

        if n >= 1:
          var1 = Y[0,1:n+1] * betaM_vec[0:n] * Y[1,0:n]
          int1 = np.dot(A_vec[n:0:-1], var1)

          int2 = np.dot(K_vec[n-1::-1], g_vec[0:n])
        else:
          int1 = 0.0
          int2 = 0.0

        Y[1, n+1] = F0_vec[n+1] + D * int1

        M[n+1] = M0_vec[n+1] + gammah * (1/N) * D * int2

    stepping_time = time.time() - start_step
    elapsed_time = allocation_time + weight_definition_time + stepping_time

    # ------------------- PCB:E | Performance metrics and diagnostics -------------------
    flops_counter = np.sum(9 + 6*np.arange(1,Nt)) + 18*D
    memory_bytes = Y.nbytes + t.nbytes + A_vec.nbytes + F0_vec.nbytes
    P = {}
    P['elapsed_time'] = elapsed_time
    P['steps_number'] = Nt
    P['flop_number'] = flops_counter
    P['Kern_number'] = Nt + K_eval
    P['memory_bytes'] = memory_bytes

    return t, Y, P

# =============================== AUXILIARY MODULES ===========================
# ------------------- PCB:F | Auxiliary routine 1------------------
def parse_input(defaults, *args, **kwargs):
    """Parse and complete the BEH problem input."""
    BEH_prob = defaults.copy()

    # Update from dict or positional arguments
    if len(args) == 1 and isinstance(args[0], dict):
        BEH_prob.update(args[0])
    elif len(args) > 0:
        keys = list(defaults.keys())
        for k,arg in enumerate(args):
            if k < len(keys) and arg is not None:
                BEH_prob[keys[k]] = arg
    BEH_prob.update(kwargs)

    verbosity = BEH_prob['verbosity']

    # Default function definitions
    ni = 1/7
    R0 = 20
    bb = R0 * (BEH_prob['mu'] + ni)**2
    # ----- Fill in missing fields with defaults -----
    if BEH_prob['A'] is None:
        BEH_prob['A'] = lambda x: bb * x * np.exp(-(BEH_prob['mu'] + ni)*x)
    if BEH_prob['K'] is None:
        a = 1/30
        BEH_prob['K'] = lambda x: a * np.exp(-a*x)
    if BEH_prob['betaM'] is None:
        alpha = 8e3
        BEH_prob['betaM'] = lambda x: 1 / ((1 / BEH_prob['N']) + alpha*x)
    if BEH_prob['g'] is None:
        BEH_prob['g'] = lambda x: x
    return BEH_prob

# ------------------- PCB:G | Auxiliary routine 2------------------
""" It performs validation and consistency checks on scalar, vector and
 matrix parameters of the BEH problem definition
 INPUT:
   => BEH_prob          : struct containing problem definition
   => defaults          : struct with default values for each field
   => A_vec             : vector of A(t_n) values 

 CHECKS PERFORMED:
   (1) Positivity and finiteness of scalar fields
   (2) Non-negativity of A(t)
   (3) Logical relation: N > S0

 NOTES:
   Errors stop execution.
"""
def _check_input(BEH_prob, defaults, A_vec):
    """Validate BEH_prob structure fields and consistency."""
    scalar_fields = ['S0','N','T','h','mu']
    for f in scalar_fields:
        value = BEH_prob[f]
        if not (np.isscalar(value) and value > 0 and np.isfinite(value)):
            raise ValueError(f"Invalid scalar field: {f}")
    if np.any(A_vec < 0):
        raise ValueError("A_vec must be non-negative")
    if BEH_prob['N'] <= BEH_prob['S0']:
        raise ValueError("Invalid values: N must be greater than S0")