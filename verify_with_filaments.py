#!/usr/bin/env python3
"""
Verification: Scaling law WITH filament data (as in Ver 9.2 Figure 3)
"""
import numpy as np
import pandas as pd
from scipy.stats import linregress

df = pd.read_csv('4b_QIC_S_Result_N170.csv')

# Filament data from Tudorache et al. 2025 (as cited in Ver 9.2 Table 4)
filament_R = np.array([50, 1700, 15000])      # kpc
filament_v = np.array([110, 110, 110])         # km/s
filament_Deff = filament_R * filament_v        # D_eff = R * v

# Galaxy data
gal_R = df['R'].values
gal_D = df['D_eff'].values

# Combined (170 galaxies + 3 filaments = 173 data points)
all_R = np.concatenate([gal_R, filament_R])
all_D = np.concatenate([gal_D, filament_Deff])

log_R = np.log10(all_R)
log_D = np.log10(all_D)

slope, intercept, r_value, p_value, std_err = linregress(log_R, log_D)
R2 = r_value ** 2

print("=" * 70)
print("SCALING LAW WITH FILAMENT DATA (173 points = 170 galaxies + 3 filaments)")
print("=" * 70)
print(f"  Slope (alpha):  {slope:.4f} ± {std_err:.4f}")
print(f"  R²:             {R2:.4f}")
print(f"  Predicted:      alpha=1.38, R²=0.920")

v_alpha = "PASS" if abs(slope - 1.38) < 0.02 else "FAIL"
v_R2 = "PASS" if abs(R2 - 0.920) < 0.005 else "FAIL"
print(f"  Alpha match:    [{v_alpha}]  (diff = {abs(slope-1.38):.4f})")
print(f"  R² match:       [{v_R2}]  (diff = {abs(R2-0.920):.4f})")

# Bootstrap with filaments
print(f"\n  Running Bootstrap (N=10,000) with filament data...")
N_BS = 10000
np.random.seed(42)
n = len(log_R)
bs_slopes = np.zeros(N_BS)
for i in range(N_BS):
    idx = np.random.choice(n, size=n, replace=True)
    s, _, r, _, _ = linregress(log_R[idx], log_D[idx])
    bs_slopes[i] = s

ci_lo = np.percentile(bs_slopes, 2.5)
ci_hi = np.percentile(bs_slopes, 97.5)
bs_mean = np.mean(bs_slopes)
bs_std = np.std(bs_slopes)

print(f"  Bootstrap mean:  {bs_mean:.3f} ± {bs_std:.3f}")
print(f"  95% CI:          [{ci_lo:.3f}, {ci_hi:.3f}]")
print(f"  Predicted CI:    [1.24, 1.59]")
print(f"  Excludes α=1.0:  {ci_lo > 1.0}")

v_ci = "PASS" if (abs(ci_lo - 1.24) < 0.05 and abs(ci_hi - 1.59) < 0.10) else "FAIL"
print(f"  CI match:        [{v_ci}]")

# Also report galaxy-only for comparison
print(f"\n  --- For reference: Galaxy-only (N=170) ---")
log_Rg = np.log10(gal_R)
log_Dg = np.log10(gal_D)
sg, ig, rg, pg, seg = linregress(log_Rg, log_Dg)
print(f"  Galaxy-only alpha: {sg:.4f}, R²: {rg**2:.4f}")
print(f"  This is DIFFERENT from the published 1.38 because")
print(f"  Ver 9.2 includes filament data in the regression.")

print("\n" + "=" * 70)
