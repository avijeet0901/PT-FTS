from __future__ import division
import numpy as np

# Number of collective variables (dimensions)
ND = 2

# Read COLVAR file
data = np.loadtxt("COLVAR")

# Extract the ND columns starting from column index 2
cv_data = data[:, 1:1+ND]

# Compute mean of each column
cv_means = np.mean(cv_data, axis=0)

# Write the means
with open(f"../../phi-psi-ite-2", "a") as myfile: 
    for value in cv_means[:-1]:
        myfile.write(f"{value:.6f}\t")
    myfile.write(f"{cv_means[-1]:.6f}\n")
