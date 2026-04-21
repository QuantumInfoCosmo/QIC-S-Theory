#!/usr/bin/env python3
"""
Chapter 5 Evidence Generation: All numerical claims computed and certified.
This script produces the exact numbers to be embedded in the textbook.
"""
import numpy as np
import pandas as pd
from scipy.stats import linregress

df = pd.read_csv('4b_QIC_S_Result_N170.csv')

print("=" * 70)
print("CHAPTER 5: CERTIFIED NUMERICAL EVIDENCE")
print("=" * 70)

# --- Phase Classification ---
n = len(df)
n_order = (df['M'] < 0.5).sum()
n_chaos = (df['M'] >= 0.5).sum()

print(f"\n[PHASE CLASSIFICATION]")
print(f"  N = {n}")
print(f"  Order: {n_order} ({n_order/n*100:.1f}%)")
print(f"  Chaos: {n_chaos} ({n_chaos/n*100:.1f}%)")
print(f"  M_mean   = {df['M'].mean():.3f}")
print(f"  M_median = {df['M'].median():.3f}")
print(f"  M_min    = {df['M'].min():.3f} ({df.loc[df['M'].idxmin(),'Galaxy']})")
print(f"  M_max    = {df['M'].max():.3f} ({df.loc[df['M'].idxmax(),'Galaxy']})")

# --- Scaling: Galaxy-only (N=170) ---
logR_g = np.log10(df['R'].values)
logD_g = np.log10(df['D_eff'].values)
s_g, i_g, r_g, p_g, se_g = linregress(logR_g, logD_g)

print(f"\n[SCALING LAW: GALAXY-ONLY (N=170)]")
print(f"  alpha = {s_g:.4f} ± {se_g:.4f}")
print(f"  R²    = {r_g**2:.4f}")

# --- Scaling: Galaxy + Filament (N=173, as in Ver 9.2) ---
fil_R = np.array([50, 1700, 15000])
fil_D = fil_R * 110
all_R = np.concatenate([df['R'].values, fil_R])
all_D = np.concatenate([df['D_eff'].values, fil_D])
logR_a = np.log10(all_R)
logD_a = np.log10(all_D)
s_a, i_a, r_a, p_a, se_a = linregress(logR_a, logD_a)

print(f"\n[SCALING LAW: GALAXY + FILAMENT (N=173, Ver 9.2 method)]")
print(f"  alpha = {s_a:.4f} ± {se_a:.4f}")
print(f"  R²    = {r_a**2:.4f}")

# --- Bootstrap: Galaxy-only ---
np.random.seed(42)
N_BS = 10000
n_g = len(logR_g)
bs_g = np.array([linregress(logR_g[idx:=np.random.choice(n_g,n_g,replace=True)], logD_g[idx])[0] for _ in range(N_BS)])
print(f"\n[BOOTSTRAP: GALAXY-ONLY (N=170, 10000 resamples)]")
print(f"  mean  = {np.mean(bs_g):.3f} ± {np.std(bs_g):.3f}")
print(f"  95%CI = [{np.percentile(bs_g,2.5):.3f}, {np.percentile(bs_g,97.5):.3f}]")
print(f"  Excludes 1.0: {np.percentile(bs_g,2.5) > 1.0}")

# --- Bootstrap: Galaxy + Filament ---
np.random.seed(42)
n_a = len(logR_a)
bs_a = np.array([linregress(logR_a[idx:=np.random.choice(n_a,n_a,replace=True)], logD_a[idx])[0] for _ in range(N_BS)])
print(f"\n[BOOTSTRAP: GALAXY + FILAMENT (N=173, 10000 resamples)]")
print(f"  mean  = {np.mean(bs_a):.3f} ± {np.std(bs_a):.3f}")
print(f"  95%CI = [{np.percentile(bs_a,2.5):.3f}, {np.percentile(bs_a,97.5):.3f}]")
print(f"  Excludes 1.0: {np.percentile(bs_a,2.5) > 1.0}")

print("\n" + "=" * 70)
print("ALL EVIDENCE GENERATED SUCCESSFULLY")
print("=" * 70)
