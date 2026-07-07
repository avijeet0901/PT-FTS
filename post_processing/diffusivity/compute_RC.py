import numpy as np
import os
from scipy.interpolate import CubicSpline

# =========================================================
# SETTINGS
# =========================================================
stride = 1
input_list_file = "sample_data.dat"
output_base = "./output"
chunk_size = 500

# =========================================================
# PERIODIC DIFFERENCE
# =========================================================
def periodic_diff(a, b):
    d = a - b
    return d - 2*np.pi * np.rint(d / (2*np.pi))

# =========================================================
# 1. LOAD REFERENCE PATH -- THIS IS MFEP
# =========================================================
path_data = np.loadtxt("phi-psi-ite-50")
phi_path = path_data[:, 0]
psi_path = path_data[:, 1]

# =========================================================
# 2. ARC-LENGTH PARAMETERIZATION
# =========================================================
ds = np.sqrt(np.diff(phi_path)**2 + np.diff(psi_path)**2)
s = np.concatenate(([0], np.cumsum(ds)))
t = s / s[-1]

# =========================================================
# 3. SPLINE FIT
# =========================================================
cs_phi = CubicSpline(t, phi_path)
cs_psi = CubicSpline(t, psi_path)
cs_s   = CubicSpline(t, s)

# =========================================================
# 4. DENSE PATH
# =========================================================
t_dense = np.linspace(0, 1, 5000)
phi_dense = cs_phi(t_dense)
psi_dense = cs_psi(t_dense)
s_dense   = cs_s(t_dense)

# =========================================================
# PERIODIC PROJECTION (BATCHED)
# =========================================================
def project_to_curve_batch(phi_q, psi_q):
    n = len(phi_q)
    t_out   = np.zeros(n)
    phi_out = np.zeros(n)
    psi_out = np.zeros(n)
    s_out   = np.zeros(n)

    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        phi_chunk = phi_q[start:end][:, None]
        psi_chunk = psi_q[start:end][:, None]
        dphi = periodic_diff(phi_chunk, phi_dense[None, :])
        dpsi = periodic_diff(psi_chunk, psi_dense[None, :])
        dist2 = dphi**2 + dpsi**2
        idx = np.argmin(dist2, axis=1)
        t_out[start:end]   = t_dense[idx]
        phi_out[start:end] = phi_dense[idx]
        psi_out[start:end] = psi_dense[idx]
        s_out[start:end]   = s_dense[idx]

    return t_out, phi_out, psi_out, s_out

# =========================================================
# 5. READ INPUT LIST
# =========================================================


with open(input_list_file, "r") as f:
    window_idx = 0                              # only increments on valid lines
    for line in f:
        if not line.strip() or line.strip().startswith("#"):  # skip blank lines and comments
            continue

        window_idx += 1                         # increment only for real data lines
        cols = line.split()
        colvar_path = cols[0]
        print(f"\nProcessing window {window_idx}: {colvar_path}")

        # OUTPUT FOLDER: ./output/1/, ./output/2/, ...
        out_dir = os.path.join(output_base, str(window_idx))
        os.makedirs(out_dir, exist_ok=True)
        save_path = os.path.join(out_dir, "COLVAR_updated")

        # =========================================================
        # READ COLVAR
        # =========================================================
        times = []
        phi_data = []
        psi_data = []
        counter = 0

        with open(colvar_path, "r") as cf:
            for line2 in cf:
                if line2.startswith("#") or not line2.strip():
                    continue
                if counter % stride != 0:
                    counter += 1
                    continue
                vals = line2.split()
                try:
                    times.append(float(vals[0]))
                    phi_data.append(float(vals[1]))
                    psi_data.append(float(vals[2]))
                except:
                    pass
                counter += 1

        if len(phi_data) == 0:
            print("No valid data, skipping.")
            continue

        times    = np.array(times)
        phi_data = np.array(phi_data)
        psi_data = np.array(psi_data)

        # =========================================================
        # PROJECT
        # =========================================================
        t_vals, phi_proj, psi_proj, s_vals = project_to_curve_batch(
            phi_data, psi_data
        )

        # =========================================================
        # SAVE OUTPUT
        # =========================================================
        updated_data = np.column_stack((
            times,
            phi_data,
            psi_data,
            phi_proj,
            psi_proj,
            t_vals,
            s_vals
        ))

        np.savetxt(
            save_path,
            updated_data,
            header="time phi psi phi_proj psi_proj t s",
            fmt="%.6f"
        )

        print(f"Saved: {save_path}")
