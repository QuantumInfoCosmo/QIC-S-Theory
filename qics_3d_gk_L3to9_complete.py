"""
===============================================================================
QIC-S 3D Lattice D_GK Scaling: L=3〜9 完全再現パッケージ
===============================================================================
目的: N=3 Ring 単位セルの3次元立方格子における Green-Kubo 輸送係数の
      N→∞ 外挿精度の検証。
===============================================================================
実行環境
===============================================================================
Python: 3.12.3
numpy:  2.4.4
scipy:  1.17.1
CPU: シングルスレッド実行
メモリ: L=9 (N=2187) で行列サイズ 2187×2187 complex128 → ~73MB per matrix
===============================================================================
パラメータ
===============================================================================
J_mean  = 1.0      # 結合定数の平均値
J_std   = 0.3      # 結合定数の乱れ（標準偏差）
beta    = 1.0      # 逆温度 (kT = 1)
t_max   = 150.0    # GK時間積分の上限
t_points= 1000     # 時間積分の離散点数
threshold= 1e-12   # rho_J2 の寄与カットオフ（微小項の除外）
アンサンブル数:
  L=3 (N=81):   20 samples
  L=4 (N=192):  20 samples
  L=5 (N=375):  10 samples
  L=6 (N=648):   5 samples
  L=7 (N=1029):  5 samples
  L=8 (N=1536):  5 samples
  L=9 (N=2187):  5 samples
乱数シード: seed = i * 1000 + L  (i = サンプル番号 0, 1, ...)
===============================================================================
"""

import numpy as np
from scipy.linalg import eigh
from scipy.stats import linregress
import time
import sys

def calc_3d_disorder_GK_fast(L, seed, beta=1.0, J_mean=1.0, J_std=0.3, t_max=150.0):
    N_cells = L**3
    N_total = 3 * N_cells
    def idx(x, y, z, alpha):
        return 3 * ((x % L) * (L**2) + (y % L) * L + (z % L)) + alpha
    np.random.seed(seed)
    g_internal = np.clip(np.random.normal(J_mean, J_std, N_total), 0.1, None)
    g_inter = np.clip(np.random.normal(J_mean, J_std, N_total), 0.1, None)
    H = np.zeros((N_total, N_total), dtype=complex)
    J_op = np.zeros((N_total, N_total), dtype=complex)
    for x in range(L):
        for y in range(L):
            for z in range(L):
                for alpha in range(3):
                    i = idx(x, y, z, alpha)
                    j = idx(x, y, z, (alpha + 1) % 3)
                    H[i, j] = -g_internal[i]; H[j, i] = -g_internal[i]
                    J_op[i, j] = 1j * g_internal[i]; J_op[j, i] = -1j * g_internal[i]
                i0, j0 = idx(x, y, z, 0), idx(x+1, y, z, 0)
                H[i0, j0] = -g_inter[i0]; H[j0, i0] = -g_inter[i0]
                J_op[i0, j0] = 1j * g_inter[i0]; J_op[j0, i0] = -1j * g_inter[i0]
                i1, j1 = idx(x, y, z, 1), idx(x, y+1, z, 1)
                H[i1, j1] = -g_inter[i1]; H[j1, i1] = -g_inter[i1]
                J_op[i1, j1] = 1j * g_inter[i1]; J_op[j1, i1] = -1j * g_inter[i1]
                i2, j2 = idx(x, y, z, 2), idx(x, y, z+1, 2)
                H[i2, j2] = -g_inter[i2]; H[j2, i2] = -g_inter[i2]
                J_op[i2, j2] = 1j * g_inter[i2]; J_op[j2, i2] = -1j * g_inter[i2]
    evals, evecs = eigh(H)
    rho_eig = np.exp(-beta * evals); rho_eig /= np.sum(rho_eig)
    J_en = evecs.conj().T @ J_op @ evecs
    rho_J2 = (rho_eig[:, None] * np.abs(J_en)**2).flatten()
    E_diff = (evals[:, None] - evals[None, :]).flatten()
    mask = rho_J2 > 1e-12
    rho_J2_f = rho_J2[mask]; E_diff_f = E_diff[mask]
    t_eval = np.linspace(0, t_max, 1000)
    corr_vals = np.array([np.sum(rho_J2_f * np.cos(E_diff_f * t)) for t in t_eval])
    try:
        from numpy import trapezoid
        integrals = [trapezoid(corr_vals[:i+1], t_eval[:i+1]) for i in range(1, len(t_eval))]
    except ImportError:
        integrals = [np.trapz(corr_vals[:i+1], t_eval[:i+1]) for i in range(1, len(t_eval))]
    D_raw = np.mean(integrals[-200:])
    D_vol = D_raw / N_total
    rho_eq_site = evecs @ np.diag(rho_eig) @ evecs.conj().T
    n_op = np.diag(np.arange(N_total, dtype=float))
    n2_avg = np.real(np.trace(rho_eq_site @ n_op @ n_op))
    n_avg = np.real(np.trace(rho_eq_site @ n_op))
    chi_local = beta * (n2_avg - n_avg**2) / N_total
    D_full = D_vol / chi_local if chi_local > 1e-15 else 0
    return {'D_raw': D_raw, 'D_vol': D_vol, 'D_full': D_full, 'chi_local': chi_local}

if __name__ == "__main__":
    print("=" * 70)
    print("QIC-S 3D LATTICE: D_GK Scaling (L=3-9)")
    print("=" * 70)
    configs = [
        (3, 20), (4, 20), (5, 10), (6, 5), (7, 5), (8, 5), (9, 5),
    ]
    all_results = {}
    for L, n_samples in configs:
        N_total = 3 * L**3
        print(f"\nL={L} (N={N_total}), ensemble={n_samples}")
        sys.stdout.flush()
        raw_list, vol_list, full_list, chi_list = [], [], [], []
        t0 = time.time()
        for i in range(n_samples):
            seed = i * 1000 + L
            r = calc_3d_disorder_GK_fast(L, seed)
            raw_list.append(r['D_raw']); vol_list.append(r['D_vol'])
            full_list.append(r['D_full']); chi_list.append(r['chi_local'])
        all_results[L] = {
            'N': N_total, 'n': n_samples,
            'D_vol_mean': np.mean(vol_list),
            'D_vol_std': np.std(vol_list, ddof=1) if len(vol_list)>1 else 0,
        }
        print(f"  D_vol = {all_results[L]['D_vol_mean']:.4e} ± {all_results[L]['D_vol_std']:.4e}  ({time.time()-t0:.1f}s)")
