
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
ax.patch.set_edgecolor('black')  ##Border color
ax.patch.set_linewidth(3.0)  #Border linewidt
ax.set_xlabel(r'Arc length', fontsize=30, fontname='Arial')
ax.set_ylabel(r'Free energy (kJ/mol)', fontsize=30, fontname='Arial')
plt.xticks(fontsize=30,fontname='Arial')
plt.yticks(fontsize=30, fontname='Arial')

plt.locator_params(axis='x', nbins=6)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.tick_params(axis="x", direction="out", length=10, width=2, color="k")
ax.tick_params(axis="x", which='minor', length=6, width=2, color='k')
plt.locator_params(axis='y', nbins=6)
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.tick_params(axis="y", direction="out", length=10, width=2, color="k")
ax.tick_params(axis="y", which='minor', length=6, width=2, color='k')

data = np.loadtxt("2D_string.dat", comments="#")
x   = data[:, 0]
y   = data[:, 1]
err = data[:, 2]/2
plt.errorbar(x, y, yerr=err, fmt="k-", alpha=0.8, label="2D")  

data = np.loadtxt("3D_string.dat", comments="#")
x   = data[:, 0]
y   = data[:, 1]
err = data[:, 2]/2
plt.errorbar(x, y, yerr=err, fmt="r-", alpha=0.8, label="3D")  

data = np.loadtxt("1d_fes_with_error.dat", comments="#")
x   = data[:, 0]
y   = data[:, 1]
err = data[:, 2]/2
#plt.plot(x, y, "r-", alpha=0.8, label="MBAR")

#err = data[:, 3]

#x = (x - x.min()) / (x.max() - x.min())
#y = y - np.min(y)

plt.errorbar(x, y, yerr=err, fmt="g-", alpha=0.8, label="2D-TASM")  # plt.error → plt.errorbar, "r-o" → fmt=

plt.legend(frameon=True, prop={'family': 'Arial', 'size': 20}, loc=0)
plt.xlim(0,1); plt.ylim(0,40)
plt.savefig("compare.png",bbox_inches = 'tight', pad_inches = 0.1, dpi=300)
plt.show()



