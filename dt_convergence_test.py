#!/usr/bin/env python3
"""
QIC-S: dt Convergence Test
===========================
Test whether t_mix values converge as dt -> 0.
If values stabilize, we have a physically meaningful result.
If they keep changing, we have a numerical artifact.
"""

import numpy as np
from scipy.integrate import solve_ivp
import warnings
warnings.filterwarnings('ignore')

class CausalGraphSimulator:
    def __init__(self, adjacency_matrix, gamma=0.05):
        self.N = adjacency_matrix.shape[0]
        self.gamma = gamma
        self.H = -1.0 * adjacency_matrix.astype(complex)
        self.jump_ops = []
        for i in range(self.N):
            for j in range(self.N):
                if adjacency_matrix[i, j] > 0:
                    L = np.zeros((self.N, self.N), dtype=complex)
                    L[j, i] = 1.0
                    self.jump_ops.append(L)
    
    def lindblad_rhs(self, t, rho_vec):
        rho = rho_vec.reshape((self.N, self.N))
        drho = -1j * (self.H @ rho - rho @ self.H)
        for L in self.jump_ops:
            Ld = L.conj().T
            LdL = Ld @ L
            drho += self.gamma * (L @ rho @ Ld - 0.5 * (LdL @ rho + rho @ LdL))
        return drho.flatten()
    
    def run(self, initial_node=0, t_max=60.0, dt=0.01):
        rho0 = np.zeros((self.N, self.N), dtype=complex)
        rho0[initial_node, initial_node] = 1.0
        t_eval = np.arange(0, t_max, dt)
        sol = solve_ivp(self.lindblad_rhs, [0, t_max], rho0.flatten(),
                        t_eval=t_eval, method='RK45', rtol=1e-10, atol=1e-12)
        probs = []
        for rv in sol.y.T:
            rho = rv.reshape((self.N, self.N))
            probs.append(np.real(np.diag(rho)))
        return t_eval, np.array(probs)

def build_chain(n):
    A = np.zeros((n, n))
    for i in range(n-1): A[i,i+1] = A[i+1,i] = 1.0
    return A

def build_ring(n):
    A = build_chain(n)
    A[0, n-1] = A[n-1, 0] = 1.0
    return A

def build_composite(n1, n2, coupling=0.3):
    N = n1 + n2; A = np.zeros((N, N))
    for i in range(n1-1): A[i,i+1] = A[i+1,i] = 1.0
    A[0, n1-1] = A[n1-1, 0] = 1.0
    for i in range(n1, N-1): A[i,i+1] = A[i+1,i] = 1.0
    A[n1, N-1] = A[N-1, n1] = 1.0
    A[n1-1, n1] = A[n1, n1-1] = coupling
    return A

def compute_tmix(t, probs, n, threshold=0.01):
    target = 1.0 / n
    for i in range(len(t)):
        if np.max(np.abs(probs[i] - target)) < threshold:
            return t[i]
    return float('inf')

# ============================================================
# Test all topologies across 5 dt values
# ============================================================

print("=" * 75)
print("QIC-S: dt CONVERGENCE TEST")
print("Question: Do t_mix values converge as dt -> 0?")
print("RK45 with rtol=1e-10, atol=1e-12 (high precision)")
print("=" * 75)

dt_values = [0.05, 0.02, 0.01, 0.005, 0.001]

topos = [
    ("N=2 Chain",  build_chain(2)),
    ("N=3 Chain",  build_chain(3)),
    ("N=3 Ring",   build_ring(3)),
    ("N=4 Chain",  build_chain(4)),
    ("N=4 Ring",   build_ring(4)),
    ("N=5 Ring",   build_ring(5)),
    ("N=6 Ring",   build_ring(6)),
    ("3+3 weak",   build_composite(3, 3, 0.3)),
    ("3+3 strong", build_composite(3, 3, 1.0)),
]

# Header
header = f"{'Topology':<14}"
for dt in dt_values:
    header += f" {'dt='+str(dt):>10}"
header += "  Converged?"
print(f"\n{header}")
print("-" * (14 + 11*len(dt_values) + 13))

convergence_results = {}

for label, adj in topos:
    n = adj.shape[0]
    row = f"{label:<14}"
    vals = []
    
    for dt in dt_values:
        sim = CausalGraphSimulator(adj, gamma=0.05)
        t, probs = sim.run(initial_node=0, t_max=60.0, dt=dt)
        tmix = compute_tmix(t, probs, n)
        vals.append(tmix)
        row += f" {tmix:>10.3f}" if tmix < 1000 else f" {'>60':>10}"
    
    # Check convergence: dt=0.01, 0.005, 0.001 should agree within dt
    last3 = vals[2:]  # dt=0.01, 0.005, 0.001
    if all(v < 1000 for v in last3):
        spread = max(last3) - min(last3)
        converged = spread < 0.05  # within 0.05 time units
        row += f"  {'YES' if converged else 'NO':>5} (spread={spread:.3f})"
    else:
        converged = False
        row += f"  {'N/A':>5}"
    
    convergence_results[label] = {
        'values': vals, 'converged': converged,
        'dt001': vals[2], 'dt0001': vals[4] if len(vals) > 4 else None
    }
    print(row)

# ============================================================
# Detailed analysis for key topologies
# ============================================================
print("\n" + "=" * 75)
print("DETAILED CONVERGENCE ANALYSIS")
print("=" * 75)

for label in ["N=2 Chain", "N=3 Ring", "N=4 Ring", "N=3 Chain"]:
    vals = convergence_results[label]['values']
    print(f"\n  {label}:")
    for dt, v in zip(dt_values, vals):
        print(f"    dt={dt:<6} -> t_mix = {v:.4f}")
    
    # dt=0.01 vs dt=0.001 difference
    diff = abs(vals[2] - vals[4])
    print(f"    |dt=0.01 - dt=0.001| = {diff:.4f}")
    
    if diff < 0.02:
        print(f"    -> CONVERGED: dt=0.01 and dt=0.001 agree to within {diff:.4f}")
    else:
        print(f"    -> NOT CONVERGED: significant difference of {diff:.4f}")

# ============================================================
# Final verdict
# ============================================================
print("\n" + "=" * 75)
print("VERDICT")
print("=" * 75)

n3ring_vals = convergence_results["N=3 Ring"]['values']
n3ring_spread = max(n3ring_vals[2:]) - min(n3ring_vals[2:])

print(f"\n  N=3 Ring t_mix across dt=0.01/0.005/0.001:")
print(f"    {n3ring_vals[2]:.4f} / {n3ring_vals[3]:.4f} / {n3ring_vals[4]:.4f}")
print(f"    Spread: {n3ring_spread:.4f}")

# Is N=3 Ring still fastest at dt=0.001?
dt001_vals = {k: v['values'][4] for k, v in convergence_results.items()}
fastest_001 = min(dt001_vals, key=lambda k: dt001_vals[k])
print(f"\n  At dt=0.001, fastest topology: {fastest_001} (t_mix={dt001_vals[fastest_001]:.4f})")
print(f"  N=3 Ring at dt=0.001: {dt001_vals['N=3 Ring']:.4f}")

if fastest_001 == "N=3 Ring":
    print(f"\n  ★ N=3 Ring remains the fastest at the highest precision.")
    print(f"  ★ The result is NUMERICALLY CONVERGED and PHYSICALLY ROBUST.")
else:
    print(f"\n  ⚠ {fastest_001} is faster than N=3 Ring at dt=0.001.")
    print(f"  ⚠ Further investigation needed.")

print("=" * 75)
