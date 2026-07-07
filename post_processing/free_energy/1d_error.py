import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

# ---------- SETTINGS ----------
n_path_points = 50 # You can interpolate more point on the path

# ---------- STEP 1: LOAD PATH ----------
f2 = open("../../phi-psi-ite-50", "r") # This is MFEP
lines2 = f2.readlines()
f2.close()
phi_c_raw = np.array([float(x.split('\t')[0]) for x in lines2])
psi_c_raw = np.array([float(x.split('\t')[1]) for x in lines2])

# ---------- STEP 2: ARC-LENGTH PARAMETRIZE PATH ----------
dx = np.diff(phi_c_raw);  dy = np.diff(psi_c_raw)
lxy = np.concatenate([[0], np.cumsum(np.sqrt(dx**2 + dy**2))])
lxy /= lxy[-1]
inter = np.linspace(0, 1, n_path_points)
phi_path = np.interp(inter, lxy, phi_c_raw)
psi_path = np.interp(inter, lxy, psi_c_raw)

# ---------- STEP 3: LOAD fes_2D_error.dat ----------
data = np.loadtxt("fes_2D_error.dat", comments="#") 
phi_2d = data[:, 0]
psi_2d = data[:, 1]
fes_2d = data[:, 2]
err_2d = data[:, 3]

points = np.column_stack([phi_2d, psi_2d])

# ---------- STEP 4: INTERPOLATE FES AND ERROR ONTO PATH ----------
fes_path = griddata(points, fes_2d, (phi_path, psi_path), method='linear')
err_path = griddata(points, err_2d, (phi_path, psi_path), method='linear')

fes_path -= np.nanmin(fes_path)   # shift minimum to zero

# ---------- STEP 5: SAVE ----------
with open("1d_fes_with_error.dat", "w") as f:
    f.write("# arc_length  FreeEnergy(kT)  Error_2sigma(kT)\n")
    for i in range(n_path_points):
        f.write(f"{inter[i]:.11f}\t{fes_path[i]:.11f}\t{err_path[i]:.11f}\n")
print("Saved 1d_fes_with_error.dat")

# ---------- STEP 6: PLOT ----------
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(inter, fes_path, 'o--', color='royalblue', lw=1.5, ms=4, label='FES')
ax.fill_between(inter,
                fes_path - err_path,
                fes_path + err_path,
                alpha=0.3, color='royalblue', label=r'$\pm 2\sigma$')
ax.set_xlabel('Arc length (normalized)', fontsize=16)
ax.set_ylabel('Free Energy (kT)', fontsize=16)
ax.tick_params(labelsize=13)
ax.legend(fontsize=13)
fig.tight_layout()
fig.savefig('1d_fes_with_error.png', dpi=300, bbox_inches='tight')
plt.show()
