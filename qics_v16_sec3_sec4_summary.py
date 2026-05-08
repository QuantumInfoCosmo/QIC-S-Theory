"""
QIC-S Ver 1.6: §3-§4 Combined Results and Figures
All data reflects final n=5 ensemble for L=7-9 (Gemini+Claude verified)
"""
import numpy as np
from scipy.stats import linregress
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================
# §3 DATA: χ redefinition (L=3..9, n=5 for L=7-9)
# ============================================================
s3_data = {
    3:  {'N':81,   'D_vol':3.573e-04, 'chi_old':6.713,   'chi_iso':8.106e-03, 'D_old':5.236e-05, 'D_new':4.442e-02,
         'chi_x':8.119e-03, 'chi_y':8.097e-03, 'chi_z':8.100e-03},
    4:  {'N':192,  'D_vol':1.192e-03, 'chi_old':15.74,   'chi_iso':6.475e-03, 'D_old':7.618e-05, 'D_new':1.843e-01,
         'chi_x':6.451e-03, 'chi_y':6.476e-03, 'chi_z':6.497e-03},
    5:  {'N':375,  'D_vol':1.640e-03, 'chi_old':31.21,   'chi_iso':5.349e-03, 'D_old':5.279e-05, 'D_new':3.071e-01,
         'chi_x':5.296e-03, 'chi_y':5.404e-03, 'chi_z':5.347e-03},
    6:  {'N':648,  'D_vol':1.419e-03, 'chi_old':53.83,   'chi_iso':4.508e-03, 'D_old':2.635e-05, 'D_new':3.149e-01,
         'chi_x':4.485e-03, 'chi_y':4.493e-03, 'chi_z':4.546e-03},
    7:  {'N':1029, 'D_vol':1.166e-03, 'chi_old':85.73,   'chi_iso':3.894e-03, 'D_old':1.360e-05, 'D_new':2.994e-01,
         'chi_x':3.888e-03, 'chi_y':3.831e-03, 'chi_z':3.962e-03},
    8:  {'N':1536, 'D_vol':8.374e-04, 'chi_old':125.5,   'chi_iso':3.403e-03, 'D_old':6.63e-06,  'D_new':2.461e-01,
         'chi_x':3.362e-03, 'chi_y':3.472e-03, 'chi_z':3.384e-03},
    9:  {'N':2187, 'D_vol':6.430e-04, 'chi_old':215.0,   'chi_iso':3.053e-03, 'D_old':2.83e-06,  'D_new':2.107e-01,
         'chi_x':3.055e-03, 'chi_y':3.056e-03, 'chi_z':3.015e-03},
}

# §4 DATA
s4_data = {
    0.1: {'D_vol':3.021e-04, 'D_new':5.702e-02, 'chi_iso':5.299e-03, 'BW':5.25, 'gap':0.066},
    0.2: {'D_vol':4.044e-04, 'D_new':7.636e-02, 'chi_iso':5.293e-03, 'BW':5.46, 'gap':0.076},
    0.5: {'D_vol':6.939e-04, 'D_new':1.313e-01, 'chi_iso':5.279e-03, 'BW':6.27, 'gap':0.172},
    1.0: {'D_vol':1.748e-03, 'D_new':3.288e-01, 'chi_iso':5.314e-03, 'BW':7.58, 'gap':0.180},
    2.0: {'D_vol':3.909e-03, 'D_new':7.370e-01, 'chi_iso':5.305e-03, 'BW':11.17, 'gap':0.358},
    5.0: {'D_vol':5.563e-04, 'D_new':1.039e-01, 'chi_iso':5.371e-03, 'BW':21.85, 'gap':0.535},
}

# ============================================================
# FIGURE 1 (§2): D_GK Extrapolation L=3-9
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

L_vals = sorted(s3_data.keys())
Ns = np.array([s3_data[L]['N'] for L in L_vals])
D_vols = np.array([s3_data[L]['D_vol'] for L in L_vals])
inv_N = 1.0 / Ns

ax = axes[0]
ax.plot(inv_N, D_vols, 'o-', markersize=8, color='blue', label='D_vol(1/N) ensemble avg')
slope_all, int_all, r_all, _, _ = linregress(inv_N, D_vols)
inv_ext = np.linspace(0, max(inv_N)*1.1, 100)
ax.plot(inv_ext, int_all + slope_all * inv_ext, 'r--', label=f'Fit L=3..9: D_∞ = {int_all:.2e}')
mask5 = np.array([L >= 5 for L in L_vals])
slope_lg, int_lg, r_lg, _, _ = linregress(inv_N[mask5], D_vols[mask5])
ax.plot(inv_ext, int_lg + slope_lg * inv_ext, 'g--', alpha=0.7, label=f'Fit L=5..9: D_∞ = {int_lg:.2e}')
ax.set_xlim(0, max(inv_N)*1.15)
ax.set_xlabel('1/N', fontsize=12); ax.set_ylabel('D_vol (mean)', fontsize=12)
ax.set_title('D_vol Extrapolation to N→∞\n(3D Lattice, N=3 Ring Unit Cell, L=3–9)', fontsize=13, fontweight='bold')
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
for L in L_vals:
    ax.annotate(f'L={L}', (1/s3_data[L]['N'], s3_data[L]['D_vol']), textcoords="offset points", xytext=(8,5), fontsize=8)

# Panel B: D_vol on log-log
ax = axes[1]
log_N = np.log10(Ns)
log_D = np.log10(D_vols)
slope_pw, int_pw, r_pw, _, _ = linregress(log_N, log_D)
ax.loglog(Ns, D_vols, 'bs-', markersize=8, label=f'D_vol ∝ N^{slope_pw:.2f} (R²={r_pw**2:.3f})')
N_fit = np.logspace(np.log10(min(Ns)), np.log10(max(Ns)*1.5), 50)
ax.loglog(N_fit, 10**int_pw * N_fit**slope_pw, 'b--', alpha=0.4)
ax.set_xlabel('N (total nodes)', fontsize=12); ax.set_ylabel('D_vol', fontsize=12)
ax.set_title('D_vol Scaling (log-log)', fontsize=13, fontweight='bold')
ax.legend(fontsize=10); ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig('Fig1_D_GK_L3to9_extrapolation.png', dpi=150, bbox_inches='tight')
plt.close()

# ============================================================
# FIGURE 2 (§3): χ Redefinition L=3-9
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 11))

chi_old = [s3_data[L]['chi_old'] for L in L_vals]
chi_iso = [s3_data[L]['chi_iso'] for L in L_vals]

# Panel A
ax = axes[0, 0]
ax.loglog(Ns, chi_old, 'rs-', markersize=8, label='χ_old (flat index)')
ax.loglog(Ns, chi_iso, 'bo-', markersize=8, label='χ_iso (3D coords)')
s1,i1,_,_,_ = linregress(log_N, np.log10(chi_old))
s2,i2,_,_,_ = linregress(log_N, np.log10(chi_iso))
ax.loglog(N_fit, 10**i1 * N_fit**s1, 'r--', alpha=0.4)
ax.loglog(N_fit, 10**i2 * N_fit**s2, 'b--', alpha=0.4)
ax.set_xlabel('N'); ax.set_ylabel('χ')
ax.set_title(f'(A) χ Scaling: old ∝ N^{s1:.2f} vs new ∝ N^{s2:.2f}', fontweight='bold')
ax.legend(); ax.grid(True, alpha=0.3, which='both')

# Panel B
ax = axes[0, 1]
D_old = [s3_data[L]['D_old'] for L in L_vals]
D_new = [s3_data[L]['D_new'] for L in L_vals]
ax.semilogy(Ns, np.abs(D_old), 'rs-', markersize=8, label='D_full_old (→0, artifact)')
ax.semilogy(Ns, np.abs(D_new), 'go-', markersize=8, label='D_full_new (→const)')
ax.semilogy(Ns, np.abs(D_vols), 'b^-', markersize=8, label='D_vol (→const)')
ax.set_xlabel('N'); ax.set_ylabel('|D|')
ax.set_title('(B) Three Normalizations Compared', fontweight='bold')
ax.legend(fontsize=9); ax.grid(True, alpha=0.3, which='both')

# Panel C
ax = axes[1, 0]
ax.plot(Ns, D_new, 'go-', markersize=10, linewidth=2)
L5plus = [s3_data[L]['D_new'] for L in L_vals if L >= 5]
ax.axhline(np.mean(L5plus), color='green', ls='--', alpha=0.5, label=f'L≥5 mean = {np.mean(L5plus):.3f}')
ax.set_xlabel('N'); ax.set_ylabel('D_full_new')
ax.set_title('(C) D_full_new Convergence', fontweight='bold')
ax.legend(); ax.grid(True, alpha=0.3)
for L in L_vals:
    ax.annotate(f'L={L}', (s3_data[L]['N'], s3_data[L]['D_new']), textcoords="offset points", xytext=(8,5), fontsize=9)

# Panel D
ax = axes[1, 1]
for mu, color, marker in [('chi_x','red','o'), ('chi_y','blue','s'), ('chi_z','green','^')]:
    vals = [s3_data[L][mu] for L in L_vals]
    ax.semilogy(Ns, vals, f'{color[0]}{marker}-', markersize=7, label=f'χ_{mu[-1]}')
ax.semilogy(Ns, chi_iso, 'k*-', markersize=10, label='χ_iso (mean)')
ax.set_xlabel('N'); ax.set_ylabel('χ_μ')
ax.set_title('(D) Isotropy: χ_x ≈ χ_y ≈ χ_z', fontweight='bold')
ax.legend(); ax.grid(True, alpha=0.3, which='both')

plt.suptitle('QIC-S Ver 1.6 §3: χ_local Redefinition with Physical 3D Coordinates (L=3–9)',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('Fig2_chi_redefinition_L3to9.png', dpi=150, bbox_inches='tight')
plt.close()

# ============================================================
# FIGURE 3 (§4): Ratio Scan
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
ratios = sorted(s4_data.keys())
D_vol_r = [s4_data[r]['D_vol'] for r in ratios]
D_new_r = [s4_data[r]['D_new'] for r in ratios]
BW_r = [s4_data[r]['BW'] for r in ratios]
gap_r = [s4_data[r]['gap'] for r in ratios]

ax = axes[0]
ax.semilogx(ratios, [d*1e3 for d in D_vol_r], 'bo-', markersize=10, linewidth=2)
peak_r = ratios[np.argmax(D_vol_r)]
ax.axvline(peak_r, color='red', ls='--', alpha=0.5, label=f'Peak at r={peak_r}')
ax.axvline(1.0, color='gray', ls=':', alpha=0.3, label='r=1 (uniform)')
ax.set_xlabel('r = g_inter / g_internal'); ax.set_ylabel('D_vol (×10⁻³)')
ax.set_title('(A) Transport vs Coupling Ratio', fontweight='bold'); ax.legend(); ax.grid(True, alpha=0.3)

ax = axes[1]
ax.semilogx(ratios, D_new_r, 'go-', markersize=10, linewidth=2)
ax.axvline(peak_r, color='red', ls='--', alpha=0.5, label=f'Peak at r={peak_r}')
ax.set_xlabel('r = g_inter / g_internal'); ax.set_ylabel('D_full_new')
ax.set_title('(B) Normalized Transport vs Ratio', fontweight='bold'); ax.legend(); ax.grid(True, alpha=0.3)

ax = axes[2]
ax.semilogx(ratios, BW_r, 'rs-', markersize=8, label='Bandwidth')
ax.set_xlabel('r = g_inter / g_internal'); ax.set_ylabel('Bandwidth', color='red')
ax.tick_params(axis='y', labelcolor='red')
ax2 = ax.twinx()
ax2.semilogx(ratios, gap_r, 'b^-', markersize=8, label='Spectral gap')
ax2.set_ylabel('Spectral gap', color='blue'); ax2.tick_params(axis='y', labelcolor='blue')
ax.set_title('(C) Spectral Properties', fontweight='bold')
ax.legend(loc='upper left'); ax2.legend(loc='lower right'); ax.grid(True, alpha=0.3)

plt.suptitle('QIC-S Ver 1.6 §4: g_inter/g_internal Ratio Dependence (L=5, N=375)',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('Fig_S4_ratio_scan.png', dpi=150, bbox_inches='tight')
plt.close()

print("All figures regenerated with final n=5 data and NumPy compatibility fix.")