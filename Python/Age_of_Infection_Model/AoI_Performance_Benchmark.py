# ==========================================================================
#  COMPUTATIONAL COMPLEXITY ANALYSIS OF NSFD_AoI  (part of EPITIME software)
#
#  This script evaluates the computational cost of the NSFD_AoI solver by
#  measuring:
#
#   1) CPU execution time using Python's timeit module.
#   2) Memory consumption of the main variables returned by the solver.
#
#  The analysis is repeated for increasing numbers of time steps Nt.
#
#  Output:
#   - Semilog plots of memory and CPU time versus Nt.
#   - Summary table printed for Nt multiples of 100.
#
# ==========================================================================

import numpy as np
import matplotlib.pyplot as plt
import timeit
from scipy.io import savemat

from NSFD_AoI import NSFD_AoI

print("\nEPITIME AoI: computational complexity analysis\n")

# ==========================================================================
# TEST CONFIGURATION
# ==========================================================================
N_tried = np.arange(1000, 10001, 5, dtype=int)
AoI_Problem = {
    "T": 100,
    "verbosity": 0
}

# ==========================================================================
# PREALLOCATION
# ==========================================================================
execution_time   = np.zeros(len(N_tried))
memory_allocated = np.zeros(len(N_tried))

# ==========================================================================
# WARM-UP EXECUTION
#
# The first execution may include one-time overheads associated with
# imports, cache initialization and memory allocation.
# ==========================================================================
AoI_Problem["h"] = AoI_Problem["T"] / 5
_, _, _ = NSFD_AoI(problem=AoI_Problem)

# ==========================================================================
# MAIN BENCHMARK LOOP
# ==========================================================================

# Prepare the running echo to console
Nsize = len(N_tried)
Ndisp = max(1., Nsize / 10.) * np.arange( min(Nsize, 10) + 1. )
k = 1
print("\nPerformance benchmark started...\n0%..", end = "")

for i, nt in enumerate(N_tried):
    # ----------------------------------------------------------------------
    # Current number of time steps
    # ----------------------------------------------------------------------
    AoI_Problem["h"] = AoI_Problem["T"] / nt

    # ----------------------------------------------------------------------
    # CPU TIME MEASUREMENT
    #
    # Repeat the timing experiment several times and retain the minimum
    # average runtime, similarly to MATLAB's timeit philosophy.
    # ----------------------------------------------------------------------
    stmt = lambda: NSFD_AoI( problem = AoI_Problem )

    execution_time[i] = min(
        timeit.repeat(
            stmt,
            repeat=5,
            number=10
        )
    ) / 10.0

    # ----------------------------------------------------------------------
    # MEMORY MEASUREMENT
    # ----------------------------------------------------------------------
    _, _, P = NSFD_AoI(problem=AoI_Problem)
    memory_allocated[i] = P["memory_bytes"]

    # Display running iterations percentage
    if ( i >= Ndisp[k] ):
        print(f"{np.int32(100. * Ndisp[k] / Nsize):d}%..", end = "")
        k += 1

print("100%\nBenchmark finished.")

# ==========================================================================
# GRAPHICS SETTINGS
# ==========================================================================
from platform import system
fontname = None
match system()[0]:
    case 'W': fontname = "Palatino Linotype"  # Windows
    case 'D': fontname = "Palatino"           # Mac OSX
    case 'L': fontname = "URW Palladio L"     # Linux
    case _:   print(f"WARNING: unknown underlying operating system {system()}")
if ( fontname != None ): plt.rcParams["font.family"] = fontname

# ==========================================================================
# COMPLEXITY PLOTS
# ==========================================================================
fig, ax = plt.subplots(1, 2, figsize=(10, 5))

# ==========================================================================
# MEMORY COMPLEXITY
# ==========================================================================
ax[0].semilogy(
    N_tried,
    memory_allocated/(1024**2),
    linewidth=2,
    color='b',
    label='Core Variables Memory'
)

ax[0].semilogy(
    N_tried,
    0.8 * N_tried *
    ((memory_allocated[-1]/(1024**2)) / N_tried[-1]),
    '--k',
    label=r'$O(N_t)$'
)

ax[0].set_xlabel(r'$N_t$')
ax[0].set_ylabel('Memory (MB)')
ax[0].set_title('Memory complexity analysis')
ax[0].legend(loc='lower right')
ax[0].set_box_aspect(1)

# ==========================================================================
# TIME COMPLEXITY
# ==========================================================================

ax[1].semilogy(
    N_tried,
    execution_time,
    linewidth=2,
    color='g',
    label='Measured CPU time'
)

ax[1].semilogy(
    N_tried,
    0.5 * (N_tried**2) *
    (execution_time[-1] / (N_tried[-1]**2)),
    '-.k',
    label=r'$O(N_t^2)$'
)

ax[1].set_xlabel(r'$N_t$')
ax[1].set_ylabel('CPU Time (s)')
ax[1].set_title('Time complexity analysis')
ax[1].legend(loc='lower right')
ax[1].set_box_aspect(1)

plt.tight_layout()
plt.show()

# ==========================================================================
# SUMMARY TABLE
# ==========================================================================

idx = (N_tried % 100) == 0

print()
print("--------------------------------------------------------")
print("   Nt         CPU Time (s)         Memory (MB)")
print("--------------------------------------------------------")

for nt, cpu, mem in zip(
        N_tried[idx],
        execution_time[idx],
        memory_allocated[idx]):

    print(
        f"{int(nt):6d}    "
        f"{cpu:14.4e}    "
        f"{mem/(1024**2):14.4f}"
    )

print("--------------------------------------------------------")

# ==========================================================================
# SAVE RESULTS
# ==========================================================================

savemat(
    "AOI_Python_Benchmark_Results.mat",
    {
        "N_tried": N_tried,
        "execution_time": execution_time,
        "memory_allocated": memory_allocated
    }
)

print("\nDone.\n\n")

# ==============================================================================
# End of AoI_Performance_Benchmark.py
# ==============================================================================
