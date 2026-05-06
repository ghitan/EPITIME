import matplotlib.pyplot as plt

def plotData(t, F, Feq):
    """
    Plot of Test1_JMB for the paper experiments (part of EPITIME software)
    
    B. Buonomo, E. Messina, C. Panico, A. Vecchio.
    "An integral renewal equation approach to behavioural epidemic models with 
     information index"
    Journal of Mathematical Biology, 2025, 90(8)
    DOI : 10.1007/s00285-024-02172-y
    
    ============================================================================
    Version:    1.0
    Date:       April, 2026
    Authors:    C. Panico, G. Zanghirati
    See:        https://arxiv.org/abs/2605.00067v1
    License:    GNU GPL v.3
    Repository: https://github.com/ghitan/EPITIME
    """

    # -----------------------------
    # Default
    # -----------------------------
    axisfontsize   = 11
    labelfontsize  = 11
    titlefontsize  = 13
    legendfontsize = labelfontsize

    # -----------------------------
    # Figure setup
    # -----------------------------
    fig, ax = plt.subplots(figsize=(8,4))

    # Plot F/Feq
    ax.plot(t, F / Feq, 'g', linewidth=2, label=r'$F/F_e$')

    ax.axhline(1, linestyle='--', color=[0.4, 0.8, 0.4], linewidth=1.5, label=r'$F/F_e=1$')

    # Labels and title
    ax.set_xlabel('t', fontsize=labelfontsize)
    ax.set_ylabel(r'$F(t)\,/\,F_e$', fontsize=labelfontsize)
    ax.set_title('NSFD Numerical Solution', fontsize=titlefontsize)

    # Legend
    ax.legend(fontsize=legendfontsize, loc='upper right')

    # Axes style
    ax.tick_params(axis='both', labelsize=axisfontsize)

    width_cm = 16.51
    width_in = width_cm / 2.54
    fig_width = width_in
    bbox = fig.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    aspect = bbox.height / bbox.width
    fig_height = fig_width * aspect
    fig.set_size_inches(fig_width, fig_height)

    plt.tight_layout()

    plt.show()