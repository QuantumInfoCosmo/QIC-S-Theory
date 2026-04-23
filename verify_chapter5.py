#!/usr/bin/env python3
"""
QIC-S Theory: Complete Numerical Verification (Phase 1) - Ver 1.2 Updated
========================================================
Mission: Independently verify ALL statistical claims in Chapter 5 (Ver 1.2)
against the raw data (4b_QIC_S_Result_N170.csv).

Updates: Scaling exponent and R2 values synchronized with Sasada_QICS_Theory_v1_2.pdf.
"""

import numpy as np
import pandas as pd
from scipy.stats import linregress

# ============================================================
# 1. Load Data
# ============================================================
CSV_PATH = '4b_QIC_S_Result_N170.csv'
try:
    df = pd.read_csv(CSV_PATH)
except FileNotFoundError:
    print(f"Error: {CSV_PATH} not found.")
    exit()

print("=" * 70)
print("QIC-S NUMERICAL VERIFICATION: PHASE 1 (Ver 1.2)")
print("=" * 70)
print(f"\nData loaded: {len(df)} galaxies")

# ============================================================
# VERIFICATION 1: Phase Classification (Order vs Chaos)
# ============================================================
print("\n" + "=" * 70)
print("VERIFICATION 1: Phase Classification (Order vs Chaos)")
print("=" * 70)

M_threshold = 0.5
n_total = len(df)
n_order = (df['M'] < M_threshold).sum()
n_chaos = (df['M'] >= M_threshold).sum()
pct_order = n_order / n_total * 100
pct_chaos = n_chaos / n_total * 100

M_mean = df['M'].mean()
M_median = df['M'].median()

# Ver 1.2 Expectations
pred_pct_order = 78.2
pred_pct_chaos = 21.8
pred_M_mean = 0.330
pred_M_median = 0.178

v1_phase = "PASS" if abs(pct_order - pred_pct_order) < 0.1 else "FAIL"
v1_mean = "PASS" if abs(M_mean - pred_M_mean) < 0.01 else "FAIL"
v1_median = "PASS" if abs(M_median - pred_M_median) < 0.01 else "FAIL"

print(f"  Order Phase: {n_order} ({pct_order:.1f}%) [Pred: {pred_pct_order}%] -> {v1_phase}")
print(f"  Chaos Phase: {n_chaos} ({pct_chaos:.1f}%) [Pred: {pred_pct_chaos}%] -> {v1_phase}")
print(f"  M Mean:      {M_mean:.3f} [Pred: {pred_M_mean:.3f}] -> {v1_mean}")
print(f"  M Median:    {M_median:.3f} [Pred: {pred_M_median:.3f}] -> {v1_median}")

# ============================================================
# VERIFICATION 2: Universal Scaling Law (N=170)
# ============================================================
print("\n" + "=" * 70)
print("VERIFICATION 2: Universal Scaling Law (N=170)")
print("=" * 70)

log_R = np.log10(df['R'])
log_D = np.log10(df['D_eff'])
slope, intercept, r_value, p_value, std_err = linregress(log_R, log_D)
R_squared = r_value**2

# Ver 1.2 Expectations
pred_alpha = 1.573
pred_R2 = 0.945

v2_alpha = "PASS" if abs(slope - pred_alpha) < 0.005 else "FAIL"
v2_R2 = "PASS" if abs(R_squared - pred_R2) < 0.005 else "FAIL"

print(f"  Scaling Exponent (alpha): {slope:.3f} [Pred: {pred_alpha}] -> {v2_alpha}")
print(f"  Coefficient (R2):         {R_squared:.3f} [Pred: {pred_R2}] -> {v2_R2}")

# ============================================================
# VERIFICATION 3: Bootstrap Analysis (N=10,000)
# ============================================================
print("\n" + "=" * 70)
print("VERIFICATION 3: Bootstrap Analysis (N=10,000)")
print("=" * 70)

np.random.seed(42)
N_BS = 10000
bs_slopes = np.zeros(N_BS)
n = len(log_R)

for i in range(N_BS):
    idx = np.random.choice(n, size=n, replace=True)
    s, _, _, _, _ = linregress(log_R.iloc[idx], log_D.iloc[idx])
    bs_slopes[i] = s

ci_lower = np.percentile(bs_slopes, 2.5)
ci_upper = np.percentile(bs_slopes, 97.5)
bs_mean = np.mean(bs_slopes)

# Ver 1.2 Expectations (95% CI: [1.518, 1.627])
pred_ci_lower = 1.518
pred_ci_upper = 1.627

v3_mean = "PASS" if abs(bs_mean - pred_alpha) < 0.01 else "FAIL"
v3_ci = "PASS" if (abs(ci_lower - pred_ci_lower) < 0.01 and abs(ci_upper - pred_ci_upper) < 0.01) else "FAIL"
v3_exclude = "PASS" if ci_lower > 1.0 else "FAIL"

print(f"  Bootstrap Mean alpha: {bs_mean:.3f}")
print(f"  95% CI: [{ci_lower:.3f}, {ci_upper:.3f}] [Pred: [{pred_ci_lower}, {pred_ci_upper}]] -> {v3_ci}")
print(f"  Exclude alpha=1.0:    {'Yes' if v3_exclude == 'PASS' else 'No'} -> {v3_exclude}")

# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("FINAL VERIFICATION SUMMARY (Ver 1.2 Compliance)")
print("=" * 70)

all_results = {
    'Phase distribution (78.2%/21.8%)': v1_phase,
    'M mean (0.330)': v1_mean,
    'M median (0.178)': v1_median,
    'Scaling exponent alpha (1.573)': v2_alpha,
    'R2 (0.945)': v2_R2,
    'Bootstrap 95% CI ([1.518, 1.627])': v3_ci,
    'CI strictly excludes alpha=1.0': v3_exclude,
}

for item, result in all_results.items():
    print(f"  {item:40s}: [{result}]")

n_pass = sum(1 for v in all_results.values() if v == 'PASS')
print(f"\nScore: {n_pass}/{len(all_results)} passed.")

if n_pass == len(all_results):
    print(">>> VERDICT: SUCCESS. All numerical claims are consistent with Ver 1.2 Theory.")
else:
    print(">>> VERDICT: INCONSISTENCY DETECTED. Please check the data/theory alignment.")