#!/usr/bin/env python3
"""
QIC-S Theory: Complete Numerical Verification (Phase 1)
========================================================
Mission: Independently verify ALL statistical claims in Chapter 5
against the raw data (4b_QIC_S_Result_N170.csv).

Author: Claude (Verification role)
Date: 2026-04-13
"""

import numpy as np
import pandas as pd
from scipy.stats import linregress

# ============================================================
# Load Data
# ============================================================
CSV_PATH = '4b_QIC_S_Result_N170.csv'
df = pd.read_csv(CSV_PATH)

print("=" * 70)
print("QIC-S NUMERICAL VERIFICATION: PHASE 1")
print("=" * 70)
print(f"\nData loaded: {len(df)} galaxies")
print(f"Columns: {list(df.columns)}")
print(f"Sample:\n{df.head(3)}\n")

# ============================================================
# VERIFICATION 1: Phase Classification (Order vs Chaos)
# ============================================================
print("=" * 70)
print("VERIFICATION 1: Phase Classification (Order vs Chaos)")
print("=" * 70)

M_threshold = 0.5
order_mask = df['M'] < M_threshold
chaos_mask = df['M'] >= M_threshold

n_total = len(df)
n_order = order_mask.sum()
n_chaos = chaos_mask.sum()
pct_order = n_order / n_total * 100
pct_chaos = n_chaos / n_total * 100

# Statistical properties of M
M_mean = df['M'].mean()
M_median = df['M'].median()
M_min = df['M'].min()
M_max = df['M'].max()
M_min_galaxy = df.loc[df['M'].idxmin(), 'Galaxy']
M_max_galaxy = df.loc[df['M'].idxmax(), 'Galaxy']

print(f"\n  Total galaxies:    {n_total}")
print(f"  Order (M < 0.5):   {n_order}  ({pct_order:.1f}%)")
print(f"  Chaos (M >= 0.5):  {n_chaos}  ({pct_chaos:.1f}%)")
print(f"\n  M statistics:")
print(f"    Mean:    {M_mean:.3f}")
print(f"    Median:  {M_median:.3f}")
print(f"    Min:     {M_min:.3f}  ({M_min_galaxy})")
print(f"    Max:     {M_max:.3f}  ({M_max_galaxy})")

# Compare with Ver 9.2 predictions
pred_order_pct = 78.2
pred_chaos_pct = 21.8
pred_M_mean = 0.330
pred_M_median = 0.178

v1_phase = "PASS" if (abs(pct_order - pred_order_pct) < 0.1 and 
                       abs(pct_chaos - pred_chaos_pct) < 0.1) else "FAIL"
v1_mean = "PASS" if abs(M_mean - pred_M_mean) < 0.005 else "FAIL"
v1_median = "PASS" if abs(M_median - pred_M_median) < 0.005 else "FAIL"

print(f"\n  --- JUDGMENT ---")
print(f"  Phase distribution: Computed {pct_order:.1f}%/{pct_chaos:.1f}%  "
      f"vs Predicted 78.2%/21.8%  -> [{v1_phase}]")
print(f"  M mean:   Computed {M_mean:.3f} vs Predicted {pred_M_mean:.3f}  -> [{v1_mean}]")
print(f"  M median: Computed {M_median:.3f} vs Predicted {pred_M_median:.3f}  -> [{v1_median}]")

# ============================================================
# VERIFICATION 2: Universal Scaling Law (alpha from linear regression)
# ============================================================
print("\n" + "=" * 70)
print("VERIFICATION 2: Universal Scaling Law (D_eff ∝ R^alpha)")
print("=" * 70)

# Use only positive values for log
valid = (df['R'] > 0) & (df['D_eff'] > 0)
log_R = np.log10(df.loc[valid, 'R'].values)
log_D = np.log10(df.loc[valid, 'D_eff'].values)

slope, intercept, r_value, p_value, std_err = linregress(log_R, log_D)
R_squared = r_value ** 2

print(f"\n  Valid data points: {valid.sum()} / {n_total}")
print(f"\n  Linear regression on log10(D_eff) vs log10(R):")
print(f"    Slope (alpha):       {slope:.4f} ± {std_err:.4f}")
print(f"    Intercept:           {intercept:.4f}")
print(f"    R²:                  {R_squared:.4f}")
print(f"    p-value:             {p_value:.2e}")

pred_alpha = 1.38
pred_R2 = 0.920

v2_alpha = "PASS" if abs(slope - pred_alpha) < 0.02 else "FAIL"
v2_R2 = "PASS" if abs(R_squared - pred_R2) < 0.005 else "FAIL"

print(f"\n  --- JUDGMENT ---")
print(f"  Alpha: Computed {slope:.3f} vs Predicted {pred_alpha}  -> [{v2_alpha}]")
print(f"  R²:    Computed {R_squared:.3f} vs Predicted {pred_R2}  -> [{v2_R2}]")

# ============================================================
# VERIFICATION 3: Bootstrap Analysis (95% CI for alpha)
# ============================================================
print("\n" + "=" * 70)
print("VERIFICATION 3: Bootstrap Analysis (N=10,000)")
print("=" * 70)

N_BOOTSTRAP = 10000
np.random.seed(42)  # Reproducibility

n_data = len(log_R)
bootstrap_slopes = np.zeros(N_BOOTSTRAP)
bootstrap_R2 = np.zeros(N_BOOTSTRAP)

for i in range(N_BOOTSTRAP):
    idx = np.random.choice(n_data, size=n_data, replace=True)
    s, ic, r, p, se = linregress(log_R[idx], log_D[idx])
    bootstrap_slopes[i] = s
    bootstrap_R2[i] = r ** 2

# 95% CI
ci_lower = np.percentile(bootstrap_slopes, 2.5)
ci_upper = np.percentile(bootstrap_slopes, 97.5)
bs_mean = np.mean(bootstrap_slopes)
bs_std = np.std(bootstrap_slopes)

R2_ci_lower = np.percentile(bootstrap_R2, 2.5)
R2_ci_upper = np.percentile(bootstrap_R2, 97.5)
R2_mean = np.mean(bootstrap_R2)
R2_std = np.std(bootstrap_R2)

# Does the 95% CI exclude alpha = 1.0?
excludes_1 = ci_lower > 1.0

print(f"\n  Bootstrap samples: {N_BOOTSTRAP}")
print(f"\n  Scaling Exponent (alpha):")
print(f"    Original estimate:   {slope:.3f}")
print(f"    Bootstrap mean:      {bs_mean:.3f}")
print(f"    Standard error:      {bs_std:.3f}")
print(f"    95% CI:              [{ci_lower:.3f}, {ci_upper:.3f}]")
print(f"    Bias:                {bs_mean - slope:.4f}")
print(f"\n  Coefficient of Determination (R²):")
print(f"    Original estimate:   {R_squared:.3f}")
print(f"    Bootstrap mean:      {R2_mean:.3f}")
print(f"    Standard error:      {R2_std:.3f}")
print(f"    95% CI:              [{R2_ci_lower:.3f}, {R2_ci_upper:.3f}]")
print(f"\n  Excludes alpha=1.0?:   {excludes_1}")

# Compare with Ver 9.2 predictions
pred_bs_alpha = 1.40
pred_bs_std = 0.10
pred_ci_lower = 1.24
pred_ci_upper = 1.59

v3_ci = "PASS" if (abs(ci_lower - pred_ci_lower) < 0.05 and 
                    abs(ci_upper - pred_ci_upper) < 0.05) else "FAIL"
v3_exclude = "PASS" if excludes_1 else "FAIL"
v3_mean = "PASS" if abs(bs_mean - pred_bs_alpha) < 0.05 else "FAIL"

print(f"\n  --- JUDGMENT ---")
print(f"  Bootstrap mean alpha: Computed {bs_mean:.2f} vs Predicted {pred_bs_alpha}  -> [{v3_mean}]")
print(f"  95% CI: Computed [{ci_lower:.2f}, {ci_upper:.2f}] "
      f"vs Predicted [{pred_ci_lower}, {pred_ci_upper}]  -> [{v3_ci}]")
print(f"  Excludes α=1.0: {excludes_1}  -> [{v3_exclude}]")

# ============================================================
# VERIFICATION 4: Representative galaxies spot check
# ============================================================
print("\n" + "=" * 70)
print("VERIFICATION 4: Representative Galaxy Spot Checks")
print("=" * 70)

spot_checks = {
    'NGC0100': (0.16, 'Order'),
    'UGC00128': (0.25, 'Order'),
    'NGC2403': (0.40, 'Order'),
    'NGC6503': (0.57, 'Chaos'),
}

for name, (pred_M, pred_phase) in spot_checks.items():
    row = df[df['Galaxy'].str.contains(name, case=False)]
    if len(row) > 0:
        actual_M = row['M'].values[0]
        actual_phase = 'Order' if actual_M < 0.5 else 'Chaos'
        match_M = abs(actual_M - pred_M) < 0.02
        match_phase = (actual_phase == pred_phase)
        status = "PASS" if (match_M and match_phase) else "FAIL"
        print(f"  {name:12s}: M={actual_M:.3f} (pred {pred_M:.2f}), "
              f"Phase={actual_phase} (pred {pred_phase})  -> [{status}]")
    else:
        print(f"  {name:12s}: NOT FOUND in dataset  -> [SKIP]")

# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("FINAL VERIFICATION SUMMARY")
print("=" * 70)

all_results = {
    'Phase distribution (78.2%/21.8%)': v1_phase,
    'M mean (0.330)': v1_mean,
    'M median (0.178)': v1_median,
    'Scaling exponent alpha (1.38)': v2_alpha,
    'R² (0.920)': v2_R2,
    'Bootstrap mean alpha (~1.40)': v3_mean,
    'Bootstrap 95% CI ([1.24, 1.59])': v3_ci,
    'CI excludes alpha=1.0': v3_exclude,
}

n_pass = sum(1 for v in all_results.values() if v == 'PASS')
n_fail = sum(1 for v in all_results.values() if v == 'FAIL')

for item, result in all_results.items():
    print(f"  [{result}] {item}")

print(f"\n  Total: {n_pass} PASS / {n_fail} FAIL / {len(all_results)} TOTAL")

if n_fail == 0:
    print("\n  ★★★ ALL VERIFICATIONS PASSED ★★★")
    print("  The numerical claims in Chapter 5 are CONFIRMED by independent computation.")
else:
    print(f"\n  ⚠ {n_fail} verification(s) did NOT match predictions.")
    print("  Review required before Chapter 5 can be finalized.")

print("\n" + "=" * 70)
