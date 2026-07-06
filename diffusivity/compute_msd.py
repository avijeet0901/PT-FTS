import numpy as np
import os
import matplotlib.pyplot as plt

# ==================================================
# PARAMETERS
# ==================================================

base_dir = "./output"
file_name = "COLVAR_updated"

ntraj = 100

# ==================================================
# STORAGE
# ==================================================

all_msds = []
common_lag = None

# ==================================================
# LOOP OVER TRAJECTORIES
# ==================================================

for traj_id in range(1, ntraj + 1):

    colvar_path = f"{base_dir}/{traj_id}/{file_name}"

    if not os.path.exists(colvar_path):
        print(f"Missing: {colvar_path}")
        continue

    print(f"Reading {colvar_path}")

    # --------------------------------------------------
    # LOAD DATA
    # --------------------------------------------------

    data = np.loadtxt(colvar_path, comments="#")

    time = data[:, 0]
    s    = data[:, 5]

    nframes = len(s)

    # --------------------------------------------------
    # COMPUTE MSD
    # --------------------------------------------------

    msd = np.zeros(nframes)

    for lag in range(nframes):

        displacements = s[lag:] - s[:nframes-lag]

        msd[lag] = np.mean(displacements**2)

    # lag times
    lag_time = time - time[0]

    # --------------------------------------------------
    # SAVE MSD INTO EACH FOLDER
    # --------------------------------------------------

    out_data = np.column_stack((lag_time, msd))

    out_file = f"{base_dir}/{traj_id}/msd_lag.dat"

    np.savetxt(
        out_file,
        out_data,
        fmt="%.8f",
        header="lag_time(ps)   MSD"
    )

    print(f"Saved: {out_file}")

    # --------------------------------------------------
    # STORE FOR AVERAGING
    # --------------------------------------------------

    all_msds.append(msd)

    if common_lag is None:
        common_lag = lag_time

# ==================================================
# ENSEMBLE AVERAGE MSD
# ==================================================

all_msds = np.array(all_msds)

avg_msd = np.mean(all_msds, axis=0)

std_msd = np.std(all_msds, axis=0)

# ==================================================
# SAVE AVERAGE MSD
# ==================================================

avg_data = np.column_stack((common_lag, avg_msd, std_msd))

np.savetxt(
    "avg_msd_lag.dat",
    avg_data,
    fmt="%.8f",
    header="lag_time(ps)   avg_MSD   std_MSD"
)

print("\nSaved average MSD: avg_msd_lag.dat")

# ==================================================
# PLOT
# ==================================================

plt.figure(figsize=(7,5))

plt.plot(common_lag, avg_msd)

plt.xlabel("Lag time (ps)")
plt.ylabel("MSD")

plt.title("Average MSD at Barrier")

plt.tight_layout()

plt.savefig("avg_msd_lag.png", dpi=300)

plt.show()
