import numpy as np
import time

def NSFD_AoI(*, problem=None, S0=None, N=None, T=None, h=None, beta=None, A=None, phi0=None, verbosity=None):
    """
    NSFD Method
    Unconditionally Positive, Monotonicity-preserving Non Standard Finite
    Difference (NSFD) scheme for the solution of the Age of Infection (AoI)
    epidemic model:

      S'(t) = -β S(t) φ(t)
      φ(t)  = φ₀(t) + β ∫₀ᵗ A(t-τ) S(τ)φ(τ) dτ
    
    where S(t) is the number of susceptibles and φ(t) is the total infectivity.

    Version:    1.0 (part of EPITIME software)
    Date:       April, 2026
    Authors:    M. Pezzella, G. Zanghirati
    See:        https://arxiv.org/abs/2605.00067v1
    License:    GNU GPL v.3
    Repository: https://github.com/ghitan/EPITIME

    ============================================================================
    INPUT USAGE:
    The function supports two alternative input modes:

    (1) Structured input:
        => NSFD_AoI(problem=dict(...))
        where problem is a dictionary with optional keys:
            • S0       : Initial number of susceptibles S(0)     [Default: 990]
            • N        : Population Size                         [Default: 1000]
            • T        : Final time T (integration over [0,T])   [Default: 100]
            • h        : Time step size                          [Default: 0.1]
            • beta     : Effective contact rate                  [Default: 0.001]
            • A        : Infectivity kernel A(t)                 [Default: exp(-t)]
            • phi0     : Initial infectivity φ₀(t)               [Default: (N-S0)·A(t)]
            • verbosity: 0 = silent, 1 = warnings                [Default: 1]

    (2) Positional input (via keyword arguments):
        => NSFD_AoI(S0=..., N=..., T=..., h=..., beta=..., A=..., phi0=..., verbosity=...)

    OUTPUTS:
    <= t : Time vector (uniform grid from 0 to T with step h)
    <= Y : Solution array. Each row is [Sₙ; φₙ] ≈ [S(tₙ); φ(tₙ)]
    <= P : (Optional) Structure with performance metrics:
            • elapsed_time: Total execution time (seconds) excluding checks
            • steps_number: Total number of time steps performed
            • flop_number : Estimated number of floating-point operations
            • Kern_number : Number of kernel function evaluations
            • memory_bytes: Estimated memory used by core variables (bytes)
    ==========================================================================
    
    For further details see:
    E. Messina, M. Pezzella, A. Vecchio.
    "A non-standard numerical scheme for an age-of-infection epidemic model"
    Journal of Computational Dynamics, 2022, 9(2): 239-252.
    DOI: 10.3934/jcd.2021029
    
    ==========================================================================
    
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
    
    ==========================================================================
    The Python implementation is organized into logical Python  Code Blocks,
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
    ==========================================================================
    """
    #==========================================================================
    # PCB:A | Default parameters and input processing
    #--------------------------------------------------------------------------
    # Parsing of structured or keyword-based input, completion with default
    # values, and verbosity management.
    #==========================================================================
    
    #======================================================|
    # Parse input: either a struct or positional arguments |
    #       Verbosity and Warnings definitions             |
    #======================================================|
    params = _parse_input(problem, S0, N, T, h, beta, A, phi0, verbosity)
    S0, N, T, h, beta, A, phi0, verbosity = [params[k] for k in ['S0','N','T','h','beta','A','phi0','verbosity']]

    #==========================================================================
    # PCB:B | Time grid construction and memory allocation
    #--------------------------------------------------------------------------
    # Construction of the temporal grid, handling the case T/h not integer,
    # and allocation of all core arrays.
    #==========================================================================

    #=====================|
    # Time discretization |
    #=====================|
    tic_alloc = time.time()
    t = np.arange(0, T + h, h)
    flag_small = False
    if t[-1] > T:
        t = t[:-1]
    if t[-1] != T:
        h_small = T - t[-1]
        t = np.append(t, T)
        flag_small = True
    Nt = len(t)
    H = beta * h
    A_val = A(t.reshape(-1, 1))
    P0 = phi0(t)
    Y = np.zeros((2, Nt))
    Y[:, 0] = [S0, P0[0]]
    alloc_time = time.time() - tic_alloc
    
    #==========================================================================
    # PCB:C | Input consistency and admissibility checks
    #--------------------------------------------------------------------------
    # Verification of positivity, logical consistency and admissibility of
    # model parameters and functions.
    #==========================================================================

    #=========================|
    # Input consistency check |
    #=========================|
    _check_input(params, A_val, P0, t)
    
    #==========================================================================
    # PCB:D | NSFD time-stepping scheme
    #--------------------------------------------------------------------------
    # Core numerical loop implementing the NSFD discretization of the AoI
    # model, ensuring unconditional positivity and monotonicity preservation.
    #==========================================================================

    #=====================|
    # Time-stepping loop  |
    #=====================|
    tic_step = time.time()
    H = beta * h
    for n in range(Nt - 1):
        if flag_small and n == Nt - 2: H = beta * h_small
        Y[0, n + 1] = Y[0, n] / (1 + H * Y[1, n])
        conv = np.dot(Y[1, :n + 1] * Y[0, 1:n + 2], A_val[n + 1:0:-1, 0])
        Y[1, n + 1] = P0[n + 1] + H * conv
    step_time = time.time() - tic_step
    elapsed_time = alloc_time + step_time
    
    #==========================================================================
    # PCB:E | Performance metrics and diagnostics
    #--------------------------------------------------------------------------
    # Computation of timing, FLOP count estimate, kernel evaluations and
    # memory usage.
    #==========================================================================

    #=======================================|
    # Performance metrics (optional output) |
    #=======================================|
    flop_count = 1 + sum([5 + 3 * k for k in range(1, Nt)]) + 2*int(flag_small)
    memory_bytes = Y.nbytes + t.nbytes + A_val.nbytes + P0.nbytes
    P = {
        'elapsed_time': elapsed_time,
        'steps_number': Nt - 1,
        'flop_number' : flop_count,
        'Kern_number' : Nt,
        'memory_bytes': memory_bytes
    }

    return t, Y, P

#==========================================================================
# PCB:F | Auxiliary routine 1
#--------------------------------------------------------------------------
# Input parsing, default completion and verbosity handling.
#==========================================================================

# ✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦|
# AUXILIARY MODULE 1                         |
# Parse and complete the AoI problem input   |
# ✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦|
def _parse_input(problem, S0, N, T, h, beta, A, phi0, verbosity):
    defaults = {
        'S0': 990, 'N': 1000, 'T': 10, 'h': 0.1, 'beta': 2.5e-3,
        'A': lambda t: np.exp(-3.0*t), 'phi0': None, 'verbosity': 1
    }
    # Collect input values into a dictionary
    user_input = {
        'S0': S0, 'N': N, 'T': T, 'h': h, 'beta': beta,
        'A': A, 'phi0': phi0, 'verbosity': verbosity
    }
    # If structured problem dict is provided, update user_input
    if problem is not None:
        user_input.update(problem)
    # Fill missing keys with defaults
    params = {k: user_input[k] if user_input[k] is not None else defaults[k] for k in defaults}
    # Default phi0 if not provided
    if params['phi0'] is None:
        params['phi0'] = lambda t: (params['N'] - params['S0']) * params['A'](t)
    # Verbosity: show warnings about missing inputs
    if params['verbosity'] > 0:
        missing = [k for k in defaults if user_input[k] is None]
        if missing:
            print(f"Warning: Missing inputs. Default values used for: {', '.join(missing)}")
    return params


#==========================================================================
# PCB:G | Auxiliary routine 2
#--------------------------------------------------------------------------
# Validation of AoI problem structure, admissibility and consistency checks.
#==========================================================================

# ✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦|
# AUXILIARY MODULE 2                                     |
# Validate AoI_problem structure fields and consistency  |
# ✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦✦|
def _check_input(params, A_val, P0, t):
    S0 = params['S0']
    N = params['N']
    T = params['T']
    h = params['h']
    beta = params['beta']
    verbosity = params['verbosity']
    # Check that numerical inputs (if user-defined) are positive and consistent
    # Scalars positive
    if any(p <= 0 for p in [S0, N, T, h, beta]):
        raise ValueError("Invalid input: all parameters must be positive.")
    # N > S0
    if N <= S0:
        raise ValueError("Invalid input: N must be greater than S0.")
    # phi0 consistency
    upper_bound = (N - S0) * A_val.flatten()                             
    bound_sensitivity = 2e-14
    if np.any(P0.flatten() > upper_bound + bound_sensitivity) and verbosity > 0:
        print("Warning: Some phi0 values exceed theoretical upper bound.")
    # Dimensions
    if A_val.shape[0] != len(t):
        raise ValueError("A(t) must return array with same length as t.")
    if P0.shape[0] != len(t):
        raise ValueError("phi0(t) must return array with same length as t.")