import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

# ==================================================
# INPUT
# ==================================================

fes_file = "2D_FE_TASM.dat"

# temperature (K)
T = 272.0

# Boltzmann constant in kJ/mol/K
kB = 0.0083144621

# beta
beta = 1.0 / (kB * T)

# ==================================================
# DIFFUSIVITY - From MSD fit
# ==================================================
D = 9.219830e-02 #8.387353e-02
D_err = 2.052373e-03 #2.272732e-03

# ==================================================
# approximate locations of minima and barrier
# ==================================================

s0_guess = 0.32653061224
sb_guess = 0.81632653061
s1_guess = 1.0 #0.98

# ==================================================
# LOAD FES
# ==================================================
data = np.loadtxt(fes_file)

s    = data[:, 0]
F    = data[:, 1]
Ferr = np.zeros_like(s) #You can also read the error from the file

# ==================================================
# SPLINE FIT
# ==================================================

spline = UnivariateSpline(s, F, s=0)

s_dense = np.linspace(s.min(), s.max(), 5000)

F_dense = spline(s_dense)

d2F = spline.derivative(n=2)(s_dense)

# ==================================================
# FIND CLOSEST POINTS
# ==================================================

idx0 = np.argmin(np.abs(s_dense - s0_guess))
idxb = np.argmin(np.abs(s_dense - sb_guess))
idx1 = np.argmin(np.abs(s_dense - s1_guess))

s0 = s_dense[idx0]
sb = s_dense[idxb]
s1 = s_dense[idx1]

F0 = F_dense[idx0]
Fb = F_dense[idxb]
F1 = F_dense[idx1]

# ==================================================
# GET FREE ENERGY ERRORS
# ==================================================

# nearest points in ORIGINAL data

idx0_raw = np.argmin(np.abs(s - s0))
idxb_raw = np.argmin(np.abs(s - sb))
idx1_raw = np.argmin(np.abs(s - s1))

F0_err = Ferr[idx0_raw]
Fb_err = Ferr[idxb_raw]
F1_err = Ferr[idx1_raw]

# ==================================================
# BARRIER HEIGHTS
# ==================================================

dF_0_to_1 = Fb - F0
dF_1_to_0 = Fb - F1

# propagated barrier uncertainty --> I have zero
dF_0_to_1_err = np.sqrt(Fb_err**2 + F0_err**2)
dF_1_to_0_err = np.sqrt(Fb_err**2 + F1_err**2)

# ==================================================
# CURVATURES
# ==================================================

Fpp_0 = d2F[idx0]
Fpp_b = abs(d2F[idxb])
Fpp_1 = d2F[idx1]

# ==================================================
# KRAMERS PREFACTORS
# ==================================================

prefactor_0_to_1 = (
    np.sqrt(beta * Fpp_0 * beta * Fpp_b)
    / (2.0 * np.pi)
)

prefactor_1_to_0 = (
    np.sqrt(beta * Fpp_1 * beta * Fpp_b)
    / (2.0 * np.pi)
)

# ==================================================
# KRAMERS RATES
# ==================================================

K01 = (
    D
    * prefactor_0_to_1
    * np.exp(-beta * dF_0_to_1)
)

K10 = (
    D
    * prefactor_1_to_0
    * np.exp(-beta * dF_1_to_0)
)

# ==================================================
# ERROR PROPAGATION -- ignoring curvature uncertainty
# ==================================================

rel_err_01 = np.sqrt(
    (D_err / D)**2 +
    (beta * dF_0_to_1_err)**2
)

rel_err_10 = np.sqrt(
    (D_err / D)**2 +
    (beta * dF_1_to_0_err)**2
)

K01_err = K01 * rel_err_01
K10_err = K10 * rel_err_10

# ==================================================
# PRINT RESULTS
# ==================================================

print("\n====================================")
print("FREE ENERGY ANALYSIS")
print("====================================")

print(f"\nTemperature          : {T:.2f} K")
print(f"beta                 : {beta:.6f} mol/kJ")

print("\n------------------------------------")
print("STATE LOCATIONS")
print("------------------------------------")

print(f"s0                    : {s0:.6f}")
print(f"sb                    : {sb:.6f}")
print(f"s1                    : {s1:.6f}")

print("\n------------------------------------")
print("FREE ENERGIES")
print("------------------------------------")

print(f"F(s0)                 : {F0:.6f} ± {F0_err:.6f}")
print(f"F(sb)                 : {Fb:.6f} ± {Fb_err:.6f}")
print(f"F(s1)                 : {F1:.6f} ± {F1_err:.6f}")

print("\n------------------------------------")
print("BARRIER HEIGHTS")
print("------------------------------------")

print(
    f"ΔF(0→1)               : "
    f"{dF_0_to_1:.6f} ± "
    f"{dF_0_to_1_err:.6f} kJ/mol"
)

print(
    f"ΔF(1→0)               : "
    f"{dF_1_to_0:.6f} ± "
    f"{dF_1_to_0_err:.6f} kJ/mol"
)

print("\n------------------------------------")
print("CURVATURES")
print("------------------------------------")

print(f"F''(s0)               : {Fpp_0:.6f}")
print(f"|F''(sb)|             : {Fpp_b:.6f}")
print(f"F''(s1)               : {Fpp_1:.6f}")

print("\n------------------------------------")
print("PREFACTORS")
print("------------------------------------")

print(f"P01                   : {prefactor_0_to_1:.6e}")
print(f"P10                   : {prefactor_1_to_0:.6e}")

print("\n------------------------------------")
print("DIFFUSIVITY")
print("------------------------------------")

print(
    f"D                     : "
    f"{D:.6e} ± {D_err:.6e}"
)

print("\n------------------------------------")
print("KRAMERS RATES")
print("------------------------------------")

print(
    f"K01                   : "
    f"{K01:.6e} ± {K01_err:.6e} ps^-1"
)

print(
    f"K10                   : "
    f"{K10:.6e} ± {K10_err:.6e} ps^-1"
)

# ==================================================
# PLOT
# ==================================================

plt.figure(figsize=(7,5))

plt.plot(s, F, 'o', label='data')
plt.plot(s_dense, F_dense, '-', label='spline')

plt.axvline(s0, ls='--', label='state 0')
plt.axvline(sb, color='r', ls='--', label='barrier')
plt.axvline(s1, ls='--', label='state 1')

plt.xlabel("Arc-length")
plt.ylabel("Free energy (kJ/mol)")

plt.legend()

plt.savefig(
    "locations.png",
    bbox_inches='tight',
    pad_inches=0.1,
    dpi=300
)

plt.show()
