import numpy as np

# ---------- SETTINGS ----------
ND       = 2          # number of dimensions (change to any N)
filename = "phi-psi-ite-2"
outfile  = "phi-psi-ite-2"

# ---------- READ DATA ----------
f = open(filename, "r")
lines = f.readlines()
f.close()

# Store all dimensions as rows in a 2D array: shape (ND, count)
coords = []
for x in lines:
    cols = x.strip().split('\t')
    coords.append([float(cols[d]) for d in range(ND)])

coords = np.array(coords).T      # shape (ND, count)
count  = coords.shape[1]

# ---------- SMOOTHING (endpoint-preserving) ----------
coords_n = np.zeros_like(coords)

# Fix endpoints for all dimensions
coords_n[:, 0]  = coords[:, 0]
coords_n[:, -1] = coords[:, -1]

# 3-point moving average for interior points
for j in range(1, count - 1):
    coords_n[:, j] = (coords[:, j-1] + coords[:, j] + coords[:, j+1]) / 3.0

coords = coords_n

# ---------- NORMALIZE EACH DIMENSION TO [0,1] FOR ARC-LENGTH ----------
# Prevents dimensions with large ranges from dominating the arc-length
ranges = coords.max(axis=1) - coords.min(axis=1)
ranges[ranges == 0] = 1.0          # avoid division by zero for flat dimensions
coords_norm = (coords - coords.min(axis=1, keepdims=True)) / ranges[:, None]

# ---------- ARC-LENGTH (using normalized coords) ----------
diffs = np.diff(coords_norm, axis=1)           # shape (ND, count-1)
dl    = np.sqrt(np.sum(diffs**2, axis=0))      # segment lengths, shape (count-1,)

lxy = np.concatenate([[0], np.cumsum(dl)])
lxy /= lxy[-1]                                 # normalize to [0, 1]

# ---------- EQUIDISTANT RESAMPLING (fix endpoints explicitly) ----------
inter = np.linspace(0, 1, count)

resampled = np.zeros_like(coords)
for d in range(ND):
    resampled[d] = np.interp(inter, lxy, coords[d])

# Fix endpoints exactly (interp can have tiny floating point drift)
resampled[:, 0]  = coords_n[:, 0]
resampled[:, -1] = coords_n[:, -1]

# ---------- WRITE OUTPUT ----------
fmt = '\t'.join(['%.8f'] * ND) + '\n'
with open(outfile, "w") as myfile:
    for i in range(count):
        myfile.write(fmt % tuple(resampled[d, i] for d in range(ND)))

#print(f"Done. Written {count} points ({ND}D) to {outfile}")
