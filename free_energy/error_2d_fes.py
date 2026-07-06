import numpy as np
import matplotlib.pyplot as plt
import os

# ---------- SETTINGS ----------
base_path        = "../output/"
weights_filename = "weights.dat"

colvar_filename  = "weights.dat"

stride_A         = 1
stride_W         = 1
kB               = 0.0083144621       # set to 0.0083144621 for kJ/mol
T                = 272.0       # set to 300.0 (or your T) for kJ/mol
N_iter           = 3
nbins            = 100

# ---------- READ PATH ----------
f2 = open("../../phi-psi-ite-50", "r")
lines2 = f2.readlines()
f2.close()
phi_e = [float(x.split('\t')[0]) for x in lines2]
psi_e = [float(x.split('\t')[1]) for x in lines2]
phie  = np.array(phi_e)
psie  = np.array(psi_e)

# ---------- STEP 1: READ WINDOWS ----------
phi_list = []
psi_list = []
w_list   = []

n_windows = 0
for i in range(1, 251):
    if os.path.isdir(os.path.join(base_path, str(i))):
        n_windows = i
    else:
        break
print(f"Found {n_windows} windows")

for i in range(1, n_windows + 1):
    colvar_path  = os.path.join(base_path, str(i), colvar_filename)
    weights_path = os.path.join(base_path, str(i), weights_filename)

    if not os.path.exists(colvar_path):
        print(f"Warning: Missing {colvar_path}")
        continue
    if not os.path.exists(weights_path):
        print(f"Warning: Missing {weights_path}")
        continue

    colvar_data = np.loadtxt(colvar_path, comments="#")
    if colvar_data.ndim == 1:
        colvar_data = colvar_data.reshape(1, -1)

    phi_vals = colvar_data[::stride_A, 1]
    psi_vals = colvar_data[::stride_A, 2]

    #w_vals = np.loadtxt(weights_path, comments="#")
    w_vals = colvar_data[::stride_W, 3]

    if len(phi_vals) != len(w_vals):
        print(f"Warning: Window {i} length mismatch — skipping")
        continue

    phi_list.append(phi_vals)
    psi_list.append(psi_vals)
    w_list.append(w_vals)

n_windows = len(phi_list)

# ---------- STEP 2: FULL-DATA 2D FES ----------
phi_all = np.concatenate(phi_list)
psi_all = np.concatenate(psi_list)
W_all   = np.concatenate(w_list)
print(f"Total samples : {len(phi_all)}")
print(f"Sum of weights: {W_all.sum():.6f}")

# Fixed bin edges from full data — reused in bootstrap for consistency
phi_edges = np.linspace(phi_all.min(), phi_all.max(), nbins)
psi_edges = np.linspace(psi_all.min(), psi_all.max(), nbins)

H, xedges, yedges = np.histogram2d(
    phi_all, psi_all,
    bins=[phi_edges, psi_edges],
    weights=W_all,
    density=False
)

P = H / np.sum(H)
with np.errstate(divide='ignore', invalid='ignore'):
    F = np.where(P > 0, -np.log(P), np.nan)
F -= np.nanmin(F)

phi_centers = 0.5 * (xedges[:-1] + xedges[1:])
psi_centers = 0.5 * (yedges[:-1] + yedges[1:])

kBT = kB * T
Z   = F * kBT

# Save mean 2D FES
with open("fes_2D.dat", "w") as f:
    for i in range(len(phi_centers)):
        for j in range(len(psi_centers)):
            f.write(f"{phi_centers[i]:.6f} {psi_centers[j]:.6f} {Z[i, j]:.6f}\n")
print("Saved fes_2D.dat")

# ---------- STEP 3: BOOTSTRAP 2D FES ----------
F_boot_2d = np.full((N_iter, nbins - 1, nbins - 1), np.nan)
print(f"Running {N_iter} iterations...")

for b in range(N_iter):
    phi_b_list = []
    psi_b_list = []
    w_b_list   = []

    for j in range(n_windows):
        N_j    = len(phi_list[j])
        n_draw = max(1, int(N_j / N_iter))      # frames to draw from this window
        idx    = np.random.choice(N_j, size=n_draw, replace=True)
        phi_b_list.append(phi_list[j][idx])
        psi_b_list.append(psi_list[j][idx])
        w_b_list.append(w_list[j][idx])

    phi_b = np.concatenate(phi_b_list)
    psi_b = np.concatenate(psi_b_list)
    w_b   = np.concatenate(w_b_list)

    H_b, _, _ = np.histogram2d(
        phi_b, psi_b,
        bins=[phi_edges, psi_edges],
        weights=w_b,
        density=False
    )

    P_b = H_b / np.sum(H_b)
    with np.errstate(divide='ignore', invalid='ignore'):
        F_b = np.where(P_b > 0, -np.log(P_b), np.nan)
    F_b -= np.nanmin(F_b)

    F_boot_2d[b] = F_b * kBT

    if (b + 1) % 50 == 0:
        print(f"  Iteration {b+1}/{N_iter}")

F_err_2d = 2 * np.nanstd(F_boot_2d, axis=0)

# Save 2D error map
with open("fes_2D_error.dat", "w") as f:
    f.write("# phi psi FreeEnergy(kT) Error_2sigma(kT)\n")
    for i in range(len(phi_centers)):
        for j in range(len(psi_centers)):
            f.write(f"{phi_centers[i]:.6f} {psi_centers[j]:.6f} "
                    f"{Z[i,j]:.6f} {F_err_2d[i,j]:.6f}\n")
print("Saved fes_2D_error.dat")

# Save bootstrap surfaces for use in Script 2
np.save("F_boot_2d.npy",  F_boot_2d)
np.save("phi_centers.npy", phi_centers)
np.save("psi_centers.npy", psi_centers)
np.save("Z_mean_2d.npy",   Z)
print("Saved bootstrap arrays (.npy) for 1D path analysis")

# ---------- STEP 4: PLOT ----------
X, Y = np.meshgrid(phi_centers, psi_centers)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

cont = axes[0].contourf(X, Y, Z.T, levels=30, cmap='gnuplot')
cbar = fig.colorbar(cont, ax=axes[0])
cbar.ax.set_ylabel('Free Energy (kT)', fontsize=14)
axes[0].plot(phie, psie, 'ko-', lw=0.5, markersize=1.0)
axes[0].set_xlabel(r'$\Phi$', fontsize=16)
axes[0].set_ylabel(r'$\Psi$', fontsize=16)
axes[0].set_title('Mean 2D FES', fontsize=14)

cont2 = axes[1].contourf(X, Y, F_err_2d.T, levels=30, cmap='Reds')
cbar2 = fig.colorbar(cont2, ax=axes[1])
cbar2.ax.set_ylabel('2σ Error (kT)', fontsize=14)
axes[1].plot(phie, psie, 'ko-', lw=0.5, markersize=1.0)
axes[1].set_xlabel(r'$\Phi$', fontsize=16)
axes[1].set_ylabel(r'$\Psi$', fontsize=16)
axes[1].set_title('Bootstrap 2σ Error', fontsize=14)

fig.tight_layout()
fig.savefig('fes2D_with_error.png', dpi=300, bbox_inches='tight')
plt.show()
print("Done.")
