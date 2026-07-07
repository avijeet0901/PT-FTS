import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# ── CONFIGURATION ─────────────────────────────────────────────────────────────
ND   = 2    # Number of variables (generalize this to any value)
num  = 50   # Number of iterations
rows = 27   # Number of points on the string

# Variable names and plot styles -- add sequentially as given in the phi-psi-ite column information
var_names  = ['phi', 'psi', 'Chi1', 'theta', 'zeta', 'alpha', 'beta', 'gamma', 'delta']
var_labels = [r'$\phi$', r'$\psi$', r'$\theta$', r'$\zeta$',r'\chi1', 
              r'$\alpha$', r'$\beta$', r'$\gamma$', r'$\delta$']
var_styles = ['ko-', 'ro-', 'bs-', 'gD-', 'm^-', 'cP-', 'yh-', 'kX-']

# Trim to ND
var_names  = var_names[:ND]
var_labels = var_labels[:ND]
var_styles = var_styles[:ND]

# ── READ DATA ─────────────────────────────────────────────────────────────────
# data[i] holds all values for variable i across all iterations
data = [[] for _ in range(ND)]

for count in range(1, num + 1):
    with open("phi-psi-ite-" + str(count), 'r') as f:
        for x in f:
            if not x.strip():
                continue
            cols = x.split()
            for i in range(ND):
                data[i].append(cols[i])

data = [np.array(d, dtype=float) for d in data]

# ── COMPUTE CONVERGENCE ───────────────────────────────────────────────────────
sigma = [[] for _ in range(ND)]

for j in range(rows, num * rows, rows):
    s = [0.0] * ND
    for k in range(j, j + rows):
        for i in range(ND):
            s[i] += np.abs(data[i][k] - data[i][k - rows])
    for i in range(ND):
        sigma[i].append(s[i] / rows)

# ── WRITE OUTPUT ──────────────────────────────────────────────────────────────
fmt = "\t".join(["%.11f"] * ND) + "\n"
with open("convergence.dat", "a") as myf:
    for j in range(len(sigma[0])):
        myf.write(fmt % tuple(sigma[i][j] for i in range(ND)))

# ── PLOT ──────────────────────────────────────────────────────────────────────
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

fig, axes = plt.subplots(ND, 1, figsize=(6,6), sharex=True)

# Make axes iterable if only one subplot
if ND == 1:
    axes = [axes]

for i, ax in enumerate(axes):

    # Plot
    ax.plot(sigma[i], var_styles[i], label=var_labels[i])

    # Border styling
    ax.patch.set_edgecolor('black')
    ax.patch.set_linewidth(3.0)

    # X label ONLY for bottom subplot
    if i == ND - 1:
        ax.set_xlabel(r'Iteration Number',
                      fontsize=30,
                      fontname='Arial')

        ax.tick_params(axis='x', labelsize=30)

    else:
        # Remove x ticks and labels from upper plots
        ax.tick_params(axis='x',
                       which='both',
                       bottom=False,
                       top=False,
                       labelbottom=False)

    # Y ticks
    ax.tick_params(axis='y', labelsize=25)

    # Tick locators
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))

    ax.tick_params(axis="x",
                   direction="out",
                   length=10,
                   width=2,
                   color="k")

    ax.tick_params(axis="x",
                   which='minor',
                   length=6,
                   width=2,
                   color='k')

    ax.tick_params(axis="y",
                   direction="out",
                   length=10,
                   width=2,
                   color="k")

    ax.tick_params(axis="y",
                   which='minor',
                   length=6,
                   width=2,
                   color='k')

    # Number of ticks
    ax.locator_params(axis='x', nbins=6)
    ax.locator_params(axis='y', nbins=6)

    # Legend
    ax.legend(frameon=True,
              prop={'family': 'Arial', 'size': 20},
              loc=0)

# Shared Y label for all plots
fig.supylabel(r'Convergence',
              fontsize=30,
              fontname='Arial',
              x=0.02)

plt.tight_layout()

plt.savefig('convergence.png',
            bbox_inches='tight',
            pad_inches=0.1,
            dpi=300)

plt.show()
