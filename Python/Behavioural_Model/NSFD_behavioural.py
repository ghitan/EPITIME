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

   T          : Final time of the integration                    [Default: 1000]
   h          : Time step used for discretization                [Default: 1]
   mu         : Natural death rate                               [Default: 1/(75*365)]
   N          : Total population size                            [Default: 5e7]
   S0         : Initial susceptible population                   [Default: 0.2*6e7]
   A          : Infectivity function                             [Default: bb* t* exp(-(mu +ni).*t)]
   betaM      : Inhibition function                              [Default: 1/((1 / N) + alpha* x)]
   g          : Message function                                 [Default: x)]
   K          : Memory kernel                                    [Default: a * exp(-a * x)]
   compGamh   : Flag for automatic computation of cutoff index   [Default: 1]
   verbosity  : 0 to suppress warnings, 1 to enable              [Default: 1]

   EXAMPLE: NSFD_behavioural(BEH_prob)

 2) The function has positional inputs:

   EXAMPLE: NSFD_behavioural(T, h, N, mu, S0, A, K, betaM, g, Verbosity)

OUTPUT:
t             : vector of time
S = Y[0,:]    : numerical approximation of the solution S(t)
F = Y[1,:]    : numerical approximation of the solution F(t)
P             : (optional) dict containing performance metrics and resource usage:
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
        'T':         1000,
        'h':         1,
        'N':         5e7,
        'mu':        1/(75*365),
        'S0':        0.2*5e7,
        'A':         None,
        'K':         None,
        'betaM':     None,
        'g':         None,
        'compGamh':  True,
        'verbosity': 1
    }
    # ----- Input parsing new -----
    BEH_prob = parse_input(defaults, *args, **kwargs)

    # ------------------- PCB:B | Time grid construction and memory allocation -------------------
    allocation_time = -time.time()
    T  = BEH_prob['T']
    h  = float(BEH_prob['h'])
    N  = BEH_prob['N']
    mu = BEH_prob['mu']
    S0 = BEH_prob['S0'] / N

    ni_t = int(np.floor(T / h))
    t = np.array([i * h for i in range(ni_t + 1)])
    
    D = h
    D_tilde = 0.0
    
    if abs(t[-1] - T) > 1e-12: 
        D_tilde = T - t[-1]
        t = np.append(t, T)
    else:
        t[-1] = T  
        
    Nt = len(t)

    # Memory allocation
    Y = np.zeros((2, Nt))
    M = np.zeros(Nt)
    betaM_vec = np.zeros(Nt)
    g_vec     = np.zeros(Nt)

    # Initial conditions
    Y[0,0]   = S0
    A_vec    = BEH_prob['A'](t)
    K_vec    = BEH_prob['K'](t)
    f_0_cost = (1 - S0/N) * mu
    F0_vec   = f_0_cost * A_vec
    Y[1,0]   = F0_vec[0]
    M0_vec   = (1/N)*BEH_prob['g'](N*Y[1, 0]) * K_vec
    M[0]     = M0_vec[0]

    allocation_time += time.time()

    # ------------------- PCB:C | Input consistency check -------------------
    _check_input(BEH_prob, defaults, A_vec)

    # ------------------- PCB:D | NSFD weight definition and time-stepping ------------
    weight_definition_time = -time.time()

    cutoffIndexComput_time = 0.0
    K_eval = Nt

    # ***** Added after referee report, round 1 ************************************
    if BEH_prob['compGamh']:
        cutoffIndexComput_time = -time.time()
        epsAbs     = 1.0E-10; # absolute threshold for tail tolerance
        Ctail      = 1.0E-04; # proportionality constant for discrete tail negligibility
        epsTailh   = min(epsAbs, Ctail * h) # discrete tail negligibility threshold
        blkLen     = int(max(50, np.ceil(1.0 / h))) # block length for discrete weights-sequence splitting
        blkLenOrig = blkLen
        blkSeqLen  = 3 # number of consecutive tail blocks satisfying tail negl. cond.

        # Start by taking the last blkSeqLen blocks available in K_vec.
        # If the prescribed block length is too large compared with the current
        # history length Nt, a temporary reduced block size is used only for the
        # initial tail-negligibility assessment.
    
        # Handle short histories: ensure that at least blkSeqLen blocks can be
        # extracted from the available kernel samples.
        if (blkLen * blkSeqLen > Nt):
            blkLen = int(np.floor(Nt / blkSeqLen))
            if BEH_prob['verbosity']:
                warnings.warn(f"Insufficient history length Nt={Nt}. "
                              f"Using temporary blkLen={blkLen} instead of {blkLenOrig} "
                              "for the initial tail-negligibility check.")
            reducedBlkLen = True
        else:
            reducedBlkLen = False
    
        extraTimeSteps = h * np.arange(1, blkLen + 1)
        # compute these first blkSeqLen partial sums, because the values are 
        # already available in K_vec
        blkSum = K_vec[-blkLen * blkSeqLen:].reshape((blkLen, blkSeqLen), order='F').sum(axis=0)
        # block sequence acceptance criteria:
        # 1) blkSeqLen consecutive tail blocks satisfying blkSum(j) <= epsTailh
        # 2) their sums are non-increasing, namely blkSum(j) <= blkSum(j-1)
        found = (np.all(blkSum <= epsTailh) and np.all(np.diff(blkSum) <= 0.0))
    
        # Restore the prescribed block length for all subsequent tail extensions.
        # The temporary reduction is used only to initialize the block sequence. 
        if reducedBlkLen:
            blkLen = blkLenOrig
            extraTimeSteps = h * np.arange(1, blkLen + 1)
            
        Tstart = T
        # The first blkSeqLen block sums have already been initialized from K_vec.
        nBlks  = blkSeqLen
        
        while not found:
            # Tail negligibility has not yet been detected.
            # Extend the kernel support by one additional block.
            K_vec_extra = BEH_prob['K'](Tstart + extraTimeSteps)
            nBlks += 1
            blkSum = np.append(blkSum, np.sum(K_vec_extra))
            Tstart = Tstart + extraTimeSteps[-1]
            # check both the negligibility and the non-increasing conditions of the
            # last blkSeqLen blocks
            found = np.all(blkSum[-blkSeqLen:] <= epsTailh) and \
                    np.all(np.diff(blkSum[-blkSeqLen:]) <= 0.0)
                    
        if not found: 
            warnings.warn("Unable to satisfy tail negligibility conditions")
            
        K_eval = K_eval + (nBlks - blkSeqLen) * blkLen
        Sum1 = np.sum(blkSum) + np.sum(K_vec[:-blkSeqLen])
        cutoffIndexComput_time += time.time()
        gammah = 1.0 / (h * Sum1)
    else:
        gammah = 1.0
    # ******************************************************************************

    weight_definition_time += time.time()

# =============================== Time iteration ===============================
    stepping_time = -time.time()

    for n in range(Nt - 1):
        if D_tilde != 0 and n == (Nt - 2):
            D = D_tilde

        betaM_val = BEH_prob['betaM'](N * M[n])
        g_val = BEH_prob['g'](N * Y[1, n])

        Y[0, n+1] = (Y[0, n] + mu*D) / (1 + D * (mu + N * betaM_val * Y[1, n]))

        betaM_vec[n] = N * betaM_val
        g_vec[n] = g_val

        var1 = Y[0, 1:n+2] * betaM_vec[0:n+1] * Y[1, 0:n+1]
        
        int1 = np.dot(A_vec[n+1:0:-1], var1)
        int2 = np.dot(K_vec[n::-1], g_vec[0:n+1])

        Y[1, n+1] = F0_vec[n+1] + D * int1
        M[n+1] = M0_vec[n+1] + gammah * (1/N) * D * int2

    stepping_time += time.time()
    elapsed_time = allocation_time + weight_definition_time + stepping_time

    # ------------------- PCB:E | Performance metrics and diagnostics -------------------
    Nt_D = int(Nt - D) 
    flops_counter = np.sum(9 + 6 * np.arange(1, Nt_D + 1)) + 18 * D
    memory_bytes = Y.nbytes + t.nbytes + A_vec.nbytes + F0_vec.nbytes
    
    P = {}
    P['elapsed_time'] = elapsed_time
    P['cutoff_time']  = cutoffIndexComput_time
    P['steps_number'] = Nt
    P['flops_number'] = flops_counter
    P['Kern_number']  = K_eval
    P['memory_bytes'] = memory_bytes

    return t, Y, P

# =============================== AUXILIARY MODULES ===========================
# ------------------- PCB:F | Auxiliary routine 1------------------
def parse_input(defaults, *args, **kwargs):
    BEH_prob = defaults.copy()
    
    param_order = ['T', 'h', 'N', 'mu', 'S0', 'A', 'K', 'betaM', 'g', 'compGamh', 'verbosity']

    if len(args) > 0:
        if isinstance(args[0], dict):
            user_dict = args[0]
            for key, val in user_dict.items():
                if val is not None:
                    BEH_prob[key] = val
        else:
            for i, val in enumerate(args):
                if i < len(param_order) and val is not None:
                    BEH_prob[param_order[i]] = val

    for key, val in kwargs.items():
        if val is not None:
            BEH_prob[key] = val

    for key, default_val in defaults.items():
        val = BEH_prob.get(key)
        if val is None or (isinstance(val, (list, np.ndarray)) and len(val) == 0):
            BEH_prob[key] = default_val

    ni = 1/7
    bb = (20 / 5e7) * (BEH_prob['mu'] + ni)**2 
    
    if not callable(BEH_prob['A']):
        BEH_prob['A'] = lambda t: bb * t * np.exp(-(BEH_prob['mu'] + ni) * t)
    if not callable(BEH_prob['K']):
        a = 1/30
        BEH_prob['K'] = lambda t: a * np.exp(-a * t)
    if not callable(BEH_prob['betaM']):
        alpha = 8e3
        BEH_prob['betaM'] = lambda x: 1.0 / (1.0 + alpha * x)
    if not callable(BEH_prob['g']):
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

# ==============================================================================
# End of NSFD_behavioural.py
# ==============================================================================
