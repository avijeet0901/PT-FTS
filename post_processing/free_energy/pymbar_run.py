import os
os.environ["JAX_PLATFORM_NAME"] = "cpu"
os.environ["JAX_ENABLE_X64"] = "False"
import numpy as np
from pymbar import MBAR
kB = 0.0083144621  # kJ/mol/K
target_temp = 272.0 
# ---------- STEP 0: USER INPUT ----------
D = int(input("Enter number of biased CVs: "))
stride = 10
# ---------- STEP 1: READ INPUT FILE ----------
input_file = "sample_data.dat"
paths = []
centers = []
ksprings = []
temps = []
with open(input_file) as f:
    for line in f:
        if not line.strip() or line.strip().startswith("#"):
            continue
        tokens = line.split()
        path = tokens[0]
        c = list(map(float, tokens[1:1+D]))
        k = list(map(float, tokens[1+D:1+2*D]))
        T = float(tokens[-1])
        paths.append(path)
        centers.append(c)
        ksprings.append(k)
        temps.append(T)
paths = np.array(paths)
centers = np.array(centers)
ksprings = np.array(ksprings)
temps = np.array(temps)
K = len(paths)
print(f"Detected {K} umbrella windows with {D} CVs")
# ---------- STEP 2: LOAD DATA ----------
cv_all = []
U_all = []
N_k = []
for k in range(K):
    data = np.loadtxt(paths[k])
    CV = data[::stride, 1:1+D]
    U  = data[::stride, 1+D]
    cv_all.append(CV)
    U_all.append(U)
    N_k.append(len(U))
cv_all = np.concatenate(cv_all)
U_all = np.concatenate(U_all)
N_k = np.array(N_k)
N = len(U_all)
print(f"Total samples: {N}")
# ---------- STEP 3: BUILD u_kn ----------
u_kn = np.zeros((K, N), dtype=np.float32)
L = 2 * np.pi
periodic = np.array([True] * D)  # EDIT if needed
for k in range(K):
    beta_k = 1.0 / (kB * temps[k])
    diff = cv_all - centers[k]
    for d in range(D):
        if periodic[d]:
            diff[:, d] -= L * np.round(diff[:, d] / L)
    bias = 0.5 * np.sum(ksprings[k] * diff**2, axis=1)
    u_kn[k, :] = beta_k * (U_all + bias)
    
# ---------- STEP 4: RUN PYMBAR ----------
print("Running PyMBAR...")
mbar = MBAR(u_kn, N_k, verbose=True)
f_k = mbar.f_k
# ---------- STEP 5: COMPUTE UNBIASED WEIGHTS ----------
beta0 = 1.0 / (kB * target_temp)
u0_n = beta0 * U_all   # shape (N,)
# denominator: sum_k N_k * exp(f_k - u_kn)
log_terms = f_k[:, None] - u_kn + np.log(N_k[:, None])
# log-sum-exp (stable)
max_log = np.max(log_terms, axis=0)
log_denom = max_log + np.log(np.sum(np.exp(log_terms - max_log), axis=0))
# weights
log_w_n = -u0_n - log_denom
w_n = np.exp(log_w_n)
# normalize
w_n /= np.sum(w_n)
print(f"Check normalization of weights: {np.sum(w_n)}")
# ---------- STEP 6: SPLIT WEIGHTS AND SAVE PER-WINDOW ----------
# Output format per window: time | CV1 | CV2 | ... | CVD | weight
output_base = "./output"
starts = np.concatenate([[0], np.cumsum(N_k[:-1])])  # start index for each window

# Build header: "time  CV1  CV2  ...  CVD  weight"
cv_header = "  ".join([f"CV{d+1}" for d in range(D)])
col_header = f"time  {cv_header}  weight"

for k in range(K):
    start = starts[k]
    end   = start + N_k[k]
    w_k   = w_n[start:end]   # weights for window k
    cv_k  = cv_all[start:end] # CVs for window k (already strided)

    # Reload original data to recover the time column (apply same stride)
    raw_data = np.loadtxt(paths[k])
    time_k   = raw_data[::stride, 0]

    # Stack columns: time | CV1 ... CVD | weight
    out_array = np.column_stack([time_k, cv_k, w_k])

    folder = os.path.join(output_base, str(k + 1))   # ./output/1, ./output/2, ...
    os.makedirs(folder, exist_ok=True)
    out_path = os.path.join(folder, "weights.dat")
    np.savetxt(
        out_path,
        out_array,
        header=col_header,
        fmt=["%.6f"] + ["%.6f"] * D + ["%.10e"],
    )
    print(f"Window {k+1}: {N_k[k]} samples → {out_path}  [time + {D} CVs + weight]")

# Save global free energies
np.savetxt("fk.dat", f_k, header="dimensionless free energies")
print("Done. Weights saved to per-window folders.")
