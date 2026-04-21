#!/usr/bin/env python3
"""
QIC-S Chapter 6: Numerical Verification v3 (DEFINITIVE)
========================================================
dt = 0.01 (refined), t_first (first crossing) definition.

This is the definitive version for GitHub/OSF publication.
All values are computed with dt=0.01 for consistency.
The mixing time t_mix is defined as the FIRST time the
total variation distance drops below threshold epsilon=0.01.

Author: Claude (Verification role)
Date: 2026-04-19
"""

import numpy as np
from scipy.integrate import solve_ivp
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Edge-based Lindblad Simulator
# ============================================================

class CausalGraphSimulator:
    """Edge-based Lindblad dissipation on causal graph."""
    def __init__(self, adjacency_matrix, gamma=0.05):
        self.N = adjacency_matrix.shape[0]
        self.gamma = gamma
        self.H = -1.0 * adjacency_matrix.astype(complex)
        self.jump_ops = []
        for i in range(self.N):
            for j in range(self.N):
                if adjacency_matrix[i, j] > 0:
                    L = np.zeros((self.N, self.N), dtype=complex)
                    L[j, i] = 1.0  # |j><i|
                    self.jump_ops.append(L)
    
    def lindblad_rhs(self, t, rho_vec):
        rho = rho_vec.reshape((self.N, self.N))
        drho = -1j * (self.H @ rho - rho @ self.H)
        for L in self.jump_ops:
            Ld = L.conj().T
            LdL = Ld @ L
            drho += self.gamma * (L @ rho @ Ld - 0.5 * (LdL @ rho + rho @ LdL))
        return drho.flatten()
    
    def run(self, initial_node=0, t_max=120.0, dt=0.01):
        rho0 = np.zeros((self.N, self.N), dtype=complex)
        rho0[initial_node, initial_node] = 1.0
        t_eval = np.arange(0, t_max, dt)
        sol = solve_ivp(self.lindblad_rhs, [0, t_max], rho0.flatten(),
                        t_eval=t_eval, method='RK45', rtol=1e-8, atol=1e-10)
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

def compute_mixing_time(t, probs, n, threshold=0.01):
    """First time max|P_i - 1/N| < threshold."""
    target = 1.0 / n
    for i in range(len(t)):
        if np.max(np.abs(probs[i] - target)) < threshold:
            return t[i]
    return float('inf')


# ============================================================
# Run all topologies with dt=0.01
# ============================================================

print("=" * 70)
print("QIC-S CHAPTER 6: NUMERICAL VERIFICATION v3 (DEFINITIVE)")
print("dt = 0.01, t_first definition, gamma = 0.05, epsilon = 0.01")
print("=" * 70)

gamma = 0.05
DT = 0.01

configs = [
    ("N=2 Chain",  build_chain(2)),
    ("N=3 Chain",  build_chain(3)),
    ("N=3 Ring",   build_ring(3)),
    ("N=4 Chain",  build_chain(4)),
    ("N=4 Ring",   build_ring(4)),
    ("N=5 Ring",   build_ring(5)),
    ("N=6 Ring",   build_ring(6)),
]

results = {}
print(f"\n{'Topology':<15} {'N':>3} {'Edges':>5} {'t_mix':>8} {'t/edge':>8}")
print("-" * 45)

for label, adj in configs:
    n = adj.shape[0]
    n_edges = int(np.sum(adj) / 2)
    sim = CausalGraphSimulator(adj, gamma=gamma)
    t, probs = sim.run(initial_node=0, t_max=120.0, dt=DT)
    t_mix = compute_mixing_time(t, probs, n, threshold=0.01)
    t_per_edge = t_mix / n_edges if n_edges > 0 and t_mix < 1000 else float('inf')
    results[label] = {'t_mix': t_mix, 'n': n, 'edges': n_edges, 't_per_edge': t_per_edge}
    
    t_str = f"{t_mix:.2f}" if t_mix < 1000 else ">120"
    e_str = f"{t_per_edge:.2f}" if t_per_edge < 1000 else "N/A"
    print(f"{label:<15} {n:>3} {n_edges:>5} {t_str:>8} {e_str:>8}")

print("\n--- Composite Topologies ---")
print(f"{'Topology':<20} {'N':>3} {'t_mix':>8}")
print("-" * 35)

composites = [
    ("3+3 weak (c=0.3)",  build_composite(3, 3, 0.3)),
    ("3+3 strong (c=1.0)", build_composite(3, 3, 1.0)),
    ("3+4 weak (c=0.3)",  build_composite(3, 4, 0.3)),
]

for label, adj in composites:
    n = adj.shape[0]
    sim = CausalGraphSimulator(adj, gamma=gamma)
    t, probs = sim.run(initial_node=0, t_max=120.0, dt=DT)
    t_mix = compute_mixing_time(t, probs, n, threshold=0.01)
    results[label] = {'t_mix': t_mix, 'n': n}
    t_str = f"{t_mix:.2f}" if t_mix < 1000 else ">120"
    print(f"{label:<20} {n:>3} {t_str:>8}")

# ============================================================
# Qualitative Checks
# ============================================================
print("\n" + "=" * 70)
print("QUALITATIVE CHECKS")
print("=" * 70)

prime = {k: v for k, v in results.items() if 'weak' not in k and 'strong' not in k}
fastest = min(prime, key=lambda k: prime[k]['t_mix'])
q1 = fastest == "N=3 Ring"
print(f"  [{'PASS' if q1 else 'FAIL'}] N=3 Ring is fastest: "
      f"{fastest} ({prime[fastest]['t_mix']:.2f})")

q2 = results['N=3 Ring']['t_mix'] < results['N=3 Chain']['t_mix']
ratio = results['N=3 Chain']['t_mix'] / results['N=3 Ring']['t_mix']
print(f"  [{'PASS' if q2 else 'FAIL'}] N=3 Ring << N=3 Chain: "
      f"{results['N=3 Ring']['t_mix']:.2f} vs {results['N=3 Chain']['t_mix']:.2f} "
      f"(ratio={ratio:.0f}x)")

comp_slower = all(results[k]['t_mix'] > results['N=3 Ring']['t_mix']
                   for k in results if 'weak' in k or 'strong' in k)
q3 = comp_slower
print(f"  [{'PASS' if q3 else 'FAIL'}] All composites slower than N=3 Ring")

prime_eff = {k: v for k, v in prime.items() if 'edges' in v and v['t_per_edge'] < 1000}
most_eff = min(prime_eff, key=lambda k: prime_eff[k]['t_per_edge'])
q4 = most_eff == "N=3 Ring"
print(f"  [{'PASS' if q4 else 'FAIL'}] N=3 Ring best efficiency: "
      f"{most_eff} ({prime_eff[most_eff]['t_per_edge']:.2f})")

print("\n" + "=" * 70)
if q1 and q2 and q3 and q4:
    print("★★★ ALL CHECKS PASSED ★★★")
else:
    n_pass = sum([q1,q2,q3,q4])
    print(f"⚠ {n_pass}/4 checks passed")
print("=" * 70)
