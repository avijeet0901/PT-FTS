import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import AutoMinorLocator

# ==================================================
# PARAMETERS
# ==================================================

base_dir = "./output"

ntraj = 100

# fitting region (ps)
xmin = 0.0
xmax = 1.2

# ==================================================
# READ AVERAGE MSD
# ==================================================

avg_data = np.loadtxt("avg_msd_initial.dat")

time_avg = avg_data[:, 0]
msd_avg  = avg_data[:, 1]

# ==================================================
# PLOT ALL MSD CURVES
# ==================================================

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
ax.patch.set_edgecolor('black')  ##Border color
ax.patch.set_linewidth(3.0)  #Border linewidt
ax.set_xlabel(r'Time (ps)', fontsize=30, fontname='Arial')
ax.set_ylabel(r'MSD', fontsize=30, fontname='Arial')
plt.xticks(fontsize=30,fontname='Arial')
plt.yticks(fontsize=30, fontname='Arial')

plt.locator_params(axis='x', nbins=4)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.tick_params(axis="x", direction="out", length=10, width=2, color="k")
ax.tick_params(axis="x", which='minor', length=6, width=2, color='k')
plt.locator_params(axis='y', nbins=6)
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.tick_params(axis="y", direction="out", length=10, width=2, color="k")
ax.tick_params(axis="y", which='minor', length=6, width=2, color='k')

for traj_id in range(1, ntraj + 1):

    msd_file = f"{base_dir}/{traj_id}/msd_initial.dat"

    if not os.path.exists(msd_file):
        print("NOT FOUND!")
        continue

    data = np.loadtxt(msd_file)

    time = data[:, 0]
    msd  = data[:, 1]

    #plt.plot(time,msd,color='k',alpha=0.5)

# ==================================================
# PLOT AVERAGE MSD
# ==================================================

plt.plot(
    time_avg,
    msd_avg,
    'k-',
    linewidth=3,
    label='Average MSD'
)

# ==================================================
# FIT LINEAR REGION
# ==================================================

mask = (time_avg >= xmin) & (time_avg <= xmax)

fit_x = time_avg[mask]
fit_y = msd_avg[mask]

# linear fit WITH covariance matrix
coeff, cov = np.polyfit(
    fit_x,
    fit_y,
    1,
    cov=True
)

print("Ceofficient: ", coeff)
slope = coeff[0]
intercept = coeff[1]

# uncertainty in slope
slope_err = np.sqrt(cov[0,0])

# fitted line
fit_line = slope * fit_x + intercept

plt.plot(
    fit_x,
    fit_line,
    'r--',
    linewidth=3,
    label=f'Linear Fit'  #Linear fit ({xmin}-{xmax} ps)
)

# ==================================================
# DIFFUSIVITY
# ==================================================

# MSD = 2Dt
D = slope / 2.0

# propagated uncertainty
D_err = slope_err / 2.0

print("\n====================================")
print("DIFFUSIVITY ANALYSIS")
print("====================================")

print(f"Fit region : {xmin:.3f} to {xmax:.3f} ps")

print(
    f"Slope      : "
    f"{slope:.6e} ± {slope_err:.6e}"
)

print("\nMSD = 2Dt")

print(
    f"D           : "
    f"{D:.6e} ± {D_err:.6e}"
)

print("\n------------------------------------")
print("COPY-PASTE")
print("------------------------------------")

print(f"D     = {D:.6e}")
print(f"D_err = {D_err:.6e}")

# ==================================================
# FINAL PLOT
# ==================================================

plt.legend(frameon=True, prop={'family': 'Arial', 'size': 20}, loc=0)

plt.tight_layout()

plt.xlim(0, 1.2)
plt.ylim(0, 0.3)

plt.savefig(
    "msd_fit.png",
    dpi=300
)

plt.show()
