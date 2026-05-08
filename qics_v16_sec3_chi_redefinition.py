"""
QIC-S Ver 1.6 §3: χ_local の3D格子に適した再定義
=================================================

問題: 従来の χ_local は n_op = diag(0, 1, ..., N-1) を使っていたが、
      3D格子ではノード番号と物理座標に対応がなく、χ_local ∝ N となる。

解決: 物理座標 (x, y, z) に基づく方向別感受率を定義する。

  X̂_μ = Σ_i  μ_i |i><i|   (μ = x, y, z)

  χ_μ = (β/N) [ <X̂_μ²> - <X̂_μ>² ]

これなら χ_μ は intensive (N 非依存) になるはず。

検証:
  1. χ_old (従来) と χ_new (座標ベース) の N スケーリングを比較
  2. D_full_new = D_vol / χ_new が N 非依存の定数に収束するか確認
  3. 等方性チェック: χ_x ≈ χ_y ≈ χ_z か
"""

import numpy as np
from scipy.linalg import eigh
from scipy.stats import linregress
import time
import sys


def calc_3d_with_proper_chi(L, seed, beta=1.0, J_mean=1.0, J_std=0.3, t_max=150.0):
    """
    3D格子GK計算 + 座標ベースの χ_μ を同時計算。
    """
    N_cells = L**3
    N_total = 3 * N_cells
    
    def idx(x, y, z, alpha):
        return 3 * ((x % L) * (L**2) + (y % L) * L + (z % L)) + alpha
    
    # 物理座標テーブルの構築
    coords = np.zeros((N_total, 3))  # (x, y, z) for each node
    for x in range(L):
        for y in range(L):
            for z in range(L):
                for alpha in range(3):
                    i = idx(x, y, z, alpha)
                    coords[i, 0] = x  # x座標
                    coords[i, 1] = y  # y座標
                    coords[i, 2] = z  # z座標
    
    # Hamiltonian & Current 構築
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
    
    # 対角化
    evals, evecs = eigh(H)
    rho_eig = np.exp(-beta * evals)
    rho_eig /= np.sum(rho_eig)
    
    # GK積分
    J_en = evecs.conj().T @ J_op @ evecs
    J_en_abs2 = np.abs(J_en)**2
    rho_J2 = (rho_eig[:, None] * J_en_abs2).flatten()
    E_diff = (evals[:, None] - evals[None, :]).flatten()
    mask = rho_J2 > 1e-12
    rho_J2_f = rho_J2[mask]
    E_diff_f = E_diff[mask]
    
    t_eval = np.linspace(0, t_max, 1000)
    corr_vals = np.array([np.sum(rho_J2_f * np.cos(E_diff_f * t)) for t in t_eval])
    
    # 台形積分
    try:
        from numpy import trapezoid
        integrals = [trapezoid(corr_vals[:i+1], t_eval[:i+1]) for i in range(1, len(t_eval))]
    except ImportError:
        integrals = [np.trapz(corr_vals[:i+1], t_eval[:i+1]) for i in range(1, len(t_eval))]
        
    D_raw = np.mean(integrals[-200:])
    D_vol = D_raw / N_total
    
    # 密度行列（サイト基底）
    rho_eq = evecs @ np.diag(rho_eig) @ evecs.conj().T
    
    # ── 従来の χ_old (flat index) ──
    n_op_old = np.diag(np.arange(N_total, dtype=float))
    n2_old = np.real(np.trace(rho_eq @ n_op_old @ n_op_old))
    n1_old = np.real(np.trace(rho_eq @ n_op_old))
    chi_old = beta * (n2_old - n1_old**2) / N_total
    
    # ── 新しい χ_μ (物理座標ベース) ──
    chi_new = {}
    for mu, label in enumerate(['x', 'y', 'z']):
        X_mu = np.diag(coords[:, mu])
        X2_avg = np.real(np.trace(rho_eq @ X_mu @ X_mu))
        X1_avg = np.real(np.trace(rho_eq @ X_mu))
        chi_mu = beta * (X2_avg - X1_avg**2) / N_total
        chi_new[label] = chi_mu
    
    # 等方平均
    chi_iso = np.mean([chi_new['x'], chi_new['y'], chi_new['z']])
    
    # D_full の各バージョン
    D_full_old = D_vol / chi_old if chi_old > 1e-15 else 0
    D_full_new = D_vol / chi_iso if chi_iso > 1e-15 else 0
    
    return {
        'D_raw': D_raw,
        'D_vol': D_vol,
        'D_full_old': D_full_old,
        'D_full_new': D_full_new,
        'chi_old': chi_old,
        'chi_x': chi_new['x'],
        'chi_y': chi_new['y'],
        'chi_z': chi_new['z'],
        'chi_iso': chi_iso,
    }


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("=" * 75)
    print("QIC-S Ver 1.6 §3: χ_local REDEFINITION")
    print("  OLD: n_op = diag(0,1,...,N-1)  → χ_old ∝ N")
    print("  NEW: X̂_μ = diag(x_i) etc.    → χ_μ should be O(1)")
    print("=" * 75)
    
    # 論文 Ver 1.6 Final Draft と整合させた計算設定
    configs = [
        (3, 20),
        (4, 20),
        (5, 10),
        (6, 5),
        (7, 5),
        (8, 5),
        (9, 5),
    ]
    
    all_results = {}
    
    for L, n_samples in configs:
        N = 3 * L**3
        print(f"\n{'─'*75}")
        print(f"L={L} (N={N}), ensemble={n_samples}")
        print(f"{'─'*75}")
        sys.stdout.flush()
        
        keys = ['D_raw', 'D_vol', 'D_full_old', 'D_full_new',
                'chi_old', 'chi_x', 'chi_y', 'chi_z', 'chi_iso']
        accum = {k: [] for k in keys}
        t0 = time.time()
        
        for i in range(n_samples):
            seed = i * 1000 + L
            ts = time.time()
            r = calc_3d_with_proper_chi(L, seed)
            for k in keys:
                accum[k].append(r[k])
            dt = time.time() - ts
            print(f"  [{i+1}/{n_samples}] D_vol={r['D_vol']:.3e}  "
                  f"χ_old={r['chi_old']:.3e}  χ_iso={r['chi_iso']:.3e}  "
                  f"D_old={r['D_full_old']:.3e}  D_new={r['D_full_new']:.3e}  "
                  f"({dt:.1f}s)")
            sys.stdout.flush()
        
        means = {k: np.mean(accum[k]) for k in keys}
        stds = {k: (np.std(accum[k], ddof=1) if len(accum[k]) > 1 else 0) for k in keys}
        means['N'] = N
        means['L'] = L
        means['n'] = n_samples
        all_results[L] = means
        
        print(f"  → D_vol      = {means['D_vol']:.4e}")
        print(f"  → χ_old      = {means['chi_old']:.4e}")
        print(f"  → χ_x,y,z    = {means['chi_x']:.4e}, {means['chi_y']:.4e}, {means['chi_z']:.4e}")
        print(f"  → χ_iso      = {means['chi_iso']:.4e}")
        print(f"  → D_full_old = {means['D_full_old']:.4e}")
        print(f"  → D_full_new = {means['D_full_new']:.4e}")
        print(f"  → time: {time.time()-t0:.1f}s")
    
    # スケーリング解析
    print(f"\n{'='*75}")
    print("SCALING ANALYSIS: χ_old vs χ_new")
    print(f"{'='*75}")
    
    L_vals = sorted(all_results.keys())
    Ns = np.array([all_results[L]['N'] for L in L_vals])
    log_N = np.log10(Ns)
    
    quantities = {
        'χ_old':      [all_results[L]['chi_old'] for L in L_vals],
        'χ_iso(new)': [all_results[L]['chi_iso'] for L in L_vals],
        'D_vol':      [all_results[L]['D_vol'] for L in L_vals],
        'D_full_old': [all_results[L]['D_full_old'] for L in L_vals],
        'D_full_new': [all_results[L]['D_full_new'] for L in L_vals],
    }
    
    print(f"\n{'Quantity':>15} {'∝ N^β':>10} {'β':>8} {'R²':>8} {'結論':>20}")
    print(f"{'─'*15} {'─'*10} {'─'*8} {'─'*8} {'─'*20}")
    
    for name, vals in quantities.items():
        arr = np.array(vals)
        if np.all(arr > 0):
            slope, _, r, _, stderr = linregress(log_N, np.log10(arr))
            if abs(slope) < 0.1:
                conclusion = "≈ 定数 (intensive)"
            elif abs(slope - 1) < 0.15:
                conclusion = "∝ N¹ (extensive/N)"
            elif abs(slope + 1) < 0.15:
                conclusion = "∝ N⁻¹ (artifact)"
            else:
                conclusion = f"∝ N^{slope:.2f}"
            print(f"{name:>15} {'N^'+f'{slope:.3f}':>10} {slope:8.3f} {r**2:8.4f} {conclusion:>20}")
    
    # 等方性チェック
    print(f"\n{'='*75}")
    print("ISOTROPY CHECK: χ_x ≈ χ_y ≈ χ_z ?")
    print(f"{'='*75}")
    print(f"{'L':>3} {'N':>6} {'χ_x':>12} {'χ_y':>12} {'χ_z':>12} {'max/min':>10}")
    for L in L_vals:
        r = all_results[L]
        vals = [r['chi_x'], r['chi_y'], r['chi_z']]
        ratio = max(vals) / min(vals) if min(vals) > 0 else float('inf')
        print(f"{L:3d} {r['N']:6d} {r['chi_x']:12.4e} {r['chi_y']:12.4e} {r['chi_z']:12.4e} {ratio:10.4f}")
    
    # 結果テーブル
    print(f"\n{'='*75}")
    print("FULL RESULTS TABLE")
    print(f"{'='*75}")
    print(f"{'L':>3} {'N':>6} {'D_vol':>11} {'χ_old':>11} {'χ_iso':>11} {'D_old':>11} {'D_new':>11}")
    print(f"{'─'*3} {'─'*6} {'─'*11} {'─'*11} {'─'*11} {'─'*11} {'─'*11}")
    for L in L_vals:
        r = all_results[L]
        print(f"{L:3d} {r['N']:6d} {r['D_vol']:11.4e} {r['chi_old']:11.4e} "
              f"{r['chi_iso']:11.4e} {r['D_full_old']:11.4e} {r['D_full_new']:11.4e}")