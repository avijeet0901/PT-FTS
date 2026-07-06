import numpy as np
import re
from collections import defaultdict

# =========================================================
# USER INPUT
# =========================================================

check_file = "phe-data.dat"

ncv = int(input("Enter number of CVs: "))

switch_freq = int(input("Enter exchange attempt frequency: "))

kB = 0.0083144621  # kj/mol/K

# =========================================================
# PARSE check.dat
# =========================================================

# store all independent runs for each umbrella window
window_runs = defaultdict(list)

with open(check_file, 'r') as f:

    for line in f:

        if not line.strip():
            continue

        parts = line.split()

        colvar_file = parts[0]

        values = list(map(float, parts[1:]))

        centers = np.array(values[:ncv])

        kappa = np.array(values[ncv:2*ncv])

        temperature = values[-1]

        beta = 1.0 / (kB * temperature)

        # extract md number from path
        match = re.search(r'/md(\d+)/', colvar_file)

        if match is None:
            raise ValueError(f"Cannot determine window from path:\n{colvar_file}")

        window_id = int(match.group(1))

        window_runs[window_id].append({
            "file": colvar_file,
            "center": centers,
            "kappa": kappa,
            "temp": temperature,
            "beta": beta
        })

# =========================================================
# SORT WINDOWS
# =========================================================

window_ids = sorted(window_runs.keys())

windows = []

# =========================================================
# READ + ACCUMULATE TRAJECTORIES
# =========================================================

for wid in window_ids:

    runs = window_runs[wid]

    all_traj = []

    for run in runs:

        data = np.loadtxt(run["file"])

        # assuming:
        # col0 = time
        # next ncv columns = CVs

        cvs = data[:, 1:1+ncv]

        # subsample according to exchange frequency
        cvs = cvs[::switch_freq]

        all_traj.append(cvs)

    all_traj = np.vstack(all_traj)

    windows.append({
        "id": wid,
        "center": runs[0]["center"],
        "kappa": runs[0]["kappa"],
        "beta": runs[0]["beta"],
        "traj": all_traj
    })

print(f"\nTotal umbrella windows: {len(windows)}")

# =========================================================
# HARMONIC BIAS ENERGY
# =========================================================

def bias_energy(x, center, kappa):

    dx = x - center
    dx = dx - np.pi * np.round(dx / np.pi)

    return 0.5 * np.sum(kappa * dx * dx)

# =========================================================
# EXCHANGE PROBABILITY CALCULATION
# =========================================================

exchange_probs = {}

for i in range(len(windows)-1):

    wi = windows[i]
    wj = windows[i+1]

    traj_i = wi["traj"]
    traj_j = wj["traj"]

    ni = len(traj_i)
    nj = len(traj_j)

    nswap = min(ni, nj)

    probs = []

    # odd-even style pairing is not relevant anymore
    # because we are estimating ensemble probabilities
    # from accumulated independent simulations

    for step in range(nswap):

        xi = traj_i[step]
        xj = traj_j[step]

        # energies before exchange
        Ui_xi = bias_energy(xi, wi["center"], wi["kappa"])
        Uj_xj = bias_energy(xj, wj["center"], wj["kappa"])

        # energies after exchange
        Ui_xj = bias_energy(xj, wi["center"], wi["kappa"])
        Uj_xi = bias_energy(xi, wj["center"], wj["kappa"])

        delta = (Ui_xj + Uj_xi) - (Ui_xi + Uj_xj)

        beta_avg = 0.5 * (wi["beta"] + wj["beta"])

        prob = min(1.0, np.exp(-beta_avg * delta))

        probs.append(prob)

    exchange_probs[(wi["id"], wj["id"])] = np.array(probs)

# =========================================================
# PRINT RESULTS
# =========================================================

print("\n========================================")
print("Replica Exchange Probabilities")
print("========================================\n")

for pair in exchange_probs:

    probs = exchange_probs[pair]

    avg_prob = np.mean(probs)

    print(f"md{pair[0]} <-> md{pair[1]} : {avg_prob:.4f}")
