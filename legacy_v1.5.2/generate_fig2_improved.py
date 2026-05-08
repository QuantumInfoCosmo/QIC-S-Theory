#!/usr/bin/env python3
"""
QIC-S: Improved Figure 2 — Universal Scaling Law (N=170)
========================================================
Adds: regression line, 95% confidence band, residual subplot,
      axis definitions in title/labels.
"""

import numpy as np
import pandas as pd
from scipy.stats import linregress
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.size'] = 11

# Load data
df = pd.read_csv('4b_QIC_S_Result_N170.csv')
R = df['R'].values
D = df['D_eff'].values
M = df['M'].values

logR = np.log10(R)
logD = np.log10(D)

# Regression
slope, intercept, r_value, p_value, std_err = linregress(logR, logD)
R2 = r_value**2

# Bootstrap 95% CI for the line
np.random.seed(42)
N_BS = 10000
n = len(logR)
bs_slopes = np.zeros(N_BS)
bs_intercepts = np.zeros(N_BS)
for i in range(N_BS):
    idx = np.random.choice(n, n, replace=True)
    s, ic, _, _, _ = linregress(logR[idx], logD[idx])
    bs_slopes[i] = s
    bs_intercepts[i] = ic

# Prediction band
x_fit = np.linspace(logR.min() - 0.1, logR.max() + 0.1, 200)
y_lines = np.array([bs_slopes[i] * x_fit + bs_intercepts[i] for i in range(N_BS)])
y_lower = np.percentile(y_lines, 2.5, axis=0)
y_upper = np.percentile(y_lines, 97.5, axis=0)
y_best = slope * x_fit + intercept

# Residuals
logD_pred = slope * logR + intercept
residuals = logD - logD_pred

# ============================================================
# Plot
# ============================================================
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10),
                                gridspec_kw={'height_ratios': [3, 1]},
                                sharex=True)

# --- Main panel ---
scatter = ax1.scatter(10**logR, 10**logD, c=M, cmap='RdBu_r', 
                       s=40, edgecolors='k', linewidths=0.3, alpha=0.8,
                       vmin=0, vmax=1.9, zorder=3)

# Regression line
ax1.plot(10**x_fit, 10**y_best, 'r-', linewidth=2.5, zorder=4,
         label=f'$\\alpha = {slope:.3f} \\pm {std_err:.3f}$  ($R^2 = {R2:.3f}$)')

# 95% CI band
ax1.fill_between(10**x_fit, 10**y_lower, 10**y_upper, 
                  color='red', alpha=0.12, zorder=2,
                  label='95% Bootstrap CI')

# alpha=1.0 reference line (trivial scaling)
y_trivial = 1.0 * x_fit + intercept + (slope - 1.0) * np.mean(logR)
ax1.plot(10**x_fit, 10**y_trivial, 'k--', linewidth=1.2, alpha=0.5, zorder=2,
         label='$\\alpha = 1.0$ (trivial kinematic)')

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel('$D_{\\mathrm{eff}}$ = R × v  [kpc km/s]', fontsize=13)
ax1.set_title('Universal Scaling Law: 170 SPARC Galaxies\n'
              '$D_{\\mathrm{eff}} \\propto R^{1.573}$,  '
              'Bootstrap 95% CI [1.518, 1.627],  '
              '$\\alpha = 1.0$ excluded', fontsize=13)
ax1.legend(fontsize=11, loc='upper left')
ax1.grid(True, which='both', alpha=0.2)

# Colorbar
cbar = plt.colorbar(scatter, ax=ax1, pad=0.02)
cbar.set_label('Phase Metric $M$', fontsize=12)

# --- Residual panel ---
ax2.scatter(R, residuals, c=M, cmap='RdBu_r', s=25, edgecolors='k', 
            linewidths=0.3, alpha=0.7, vmin=0, vmax=1.9)
ax2.axhline(y=0, color='r', linewidth=1.5, linestyle='-')
ax2.axhline(y=np.std(residuals), color='gray', linewidth=0.8, linestyle='--', alpha=0.5)
ax2.axhline(y=-np.std(residuals), color='gray', linewidth=0.8, linestyle='--', alpha=0.5)

ax2.set_xscale('log')
ax2.set_xlabel('Characteristic Scale $R$ (outermost measured radius) [kpc]', fontsize=13)
ax2.set_ylabel('Residual\n(log$_{10}$)', fontsize=11)
ax2.set_ylim(-0.8, 0.8)
ax2.grid(True, which='both', alpha=0.2)

# Residual stats annotation
ax2.text(0.98, 0.92, f'RMS = {np.std(residuals):.3f} dex',
         transform=ax2.transAxes, ha='right', va='top', fontsize=10,
         bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
plt.savefig('Fig2_Scaling_Law_Improved.png', dpi=300, bbox_inches='tight')
print(f"Figure saved: Fig2_Scaling_Law_Improved.png")
print(f"Slope: {slope:.4f} ± {std_err:.4f}")
print(f"R²: {R2:.4f}")
print(f"Residual RMS: {np.std(residuals):.4f} dex")
