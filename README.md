# QIC-S Theory: Universal Scaling Law of 170 Galaxies and the Minimum Stable Unit of Causal Graphs (N=3)

**Version 1.5.2 (Revised May 2026)** | **DOI:** [10.17605/OSF.IO/UJBPW](https://doi.org/10.17605/OSF.IO/UJBPW)

---

## ⚠️ Update Notice (May 2, 2026)

**Paper updated from Ver. 1.2 to Ver. 1.5.2.** This is a major revision addressing structural deficiencies identified through multi-AI peer review simulations (Claude/Anthropic, Gemini/Google, Le Chat/Mistral). All three systems independently flagged the same four issues; all have been resolved. See [Revision History](#revision-history) for details.

The previous `verify_chapter5.py` fix (April 24, 2026) remains in effect. No changes to data files or computational scripts.

---

## Overview

This repository contains the paper, data, verification codes, and figures for the Quantum Information Cosmology – Sasada (QIC-S) theory.

QIC-S proposes that the universe is a collection of discrete causal graphs, where each galaxy possesses a unique effective Hamiltonian (Multi-Hamiltonian hypothesis). The theory explains galaxy rotation curves without invoking particle dark matter, instead attributing the observed "missing mass" to the gradient energy of a scalar information field $D_{\text{eff}}$.

### Key Results

| Result | Value | Method |
|:---|:---|:---|
| Parameter-free agreement | 99.46% ± 2.53% | 7 SPARC galaxies, zero free parameters |
| SLACS lens verification | Δ < 3% (all 3 systems) | κ(θ) vs D_eff(r) cross-validation (Table 1) |
| Universal scaling law | $D_{\text{eff}} \propto R^{1.573}$ ($R^2 = 0.945$) | 170 SPARC galaxies, Tier 1 only |
| Trivial scaling excluded | α = 1.0 outside 95% CI [1.518, 1.627] | Bootstrap 10,000 resamples |
| Phase classification | 78.2% Order / 21.8% Chaos | Empirical threshold M = 0.5 |
| Minimum stable unit | N=3 Ring ($t_{\text{mix}}$ = 0.684) | Edge-based Lindblad, dt → 0 |
| Post-mixing stability | $t_{\text{stable}}$ = 27.73 (N=3 Ring only) | Among tested topologies |

---

## Repository Structure

```
.
├── Sasada_QICS_Theory_v1_5_2.pdf      # Main paper (English, 7 pages)
├── README.md                           # This file
├── 4b_QIC_S_Result_N170.csv           # 170-galaxy dataset (Galaxy, M, R, D_eff)
├── Comparison_Histogram_N170.pdf      # Figure 1: Mass agreement distribution
├── Fig2_Scaling_Law_Improved.png      # Figure 2: Universal scaling law
├── verify_chapter5.py                 # Phase classification & scaling verification
├── verify_with_filaments.py           # Galaxy + filament combined regression
├── chapter5_evidence.py               # All Chapter 5 statistics (one-shot)
├── verify_chapter6_v3_dt001.py        # Edge-based Lindblad mixing time (dt=0.01)
├── dt_convergence_test.py             # dt convergence test (5 stages)
└── generate_fig2_improved.py          # Figure 2 generation script
```

---

## Reproducing the Results

### Requirements
- Python >= 3.9
- numpy, pandas, scipy, matplotlib

### Section 2: Phase Classification & Scaling Law
```bash
python verify_chapter5.py          # 7/7 PASS expected
python verify_with_filaments.py    # Galaxy+filament combined regression
python chapter5_evidence.py        # Full numerical evidence output
```

### Section 3: Edge-based Lindblad Mixing Times
```bash
python verify_chapter6_v3_dt001.py # All topologies (dt=0.01)
python dt_convergence_test.py      # 5-stage convergence (dt=0.05→0.001)
```

### Figure 2
```bash
python generate_fig2_improved.py   # Requires 4b_QIC_S_Result_N170.csv
```

---

## Data Description

### 4b_QIC_S_Result_N170.csv

| Column | Type | Description |
|:---|:---|:---|
| Galaxy | string | Galaxy identifier (SPARC naming convention) |
| M | float | Phase Metric $M = \text{Var}(\log(\|\nabla H\| + \varepsilon))$ |
| R | float | Characteristic scale [kpc] (outermost measured radius) |
| D_eff | float | Effective transport coefficient [kpc km/s] ($= R \times v$) |

**Source:** Derived from SPARC database (Lelli, McGaugh, & Schombert, 2016, AJ, 152, 157).

**Physical basis of $D_{\text{eff}} = R \times v$:** Grounded in Green-Kubo formalism. The correlation time $\tau \sim R/v$ yields $D_{\text{eff}} \sim \langle v \cdot v \rangle \tau = v^2 \cdot (R/v) = R \times v$. See paper §3.3.

---

## Scope and Limitations

- All analyses are restricted to the **galactic scale (Tier 1)**.
- Application to galaxy cluster scales (Tier 2) is **suspended** to avoid circular reasoning from ΛCDM-biased mass estimates (PSZ2 R₅₀₀/M₅₀₀).
- The N=3 Ring is the most stable topology **among those tested**; absolute universality across all possible topologies remains unproven.
- The micro-macro connection (N=3 → hydrodynamic limit) is **numerically confirmed** ($R^2 = 0.9843$ at N=32) but **not rigorously proven** (requires $N \to \infty$).
- The M = 0.5 phase threshold is an **empirical criterion**; derivation from bifurcation theory is an open problem.

---

## AI Disclosure

All theoretical ideas, physical interpretations, and research decisions are the sole intellectual responsibility of the author. AI systems were employed as auxiliary tools in accordance with their functional roles:

| AI System | Role |
|:---|:---|
| **Claude (Anthropic)** | Theoretical synthesis, independent numerical verification, document preparation |
| **Gemini (Google)** | Code development, simulation execution, LaTeX typesetting, cross-validation |
| **ChatGPT (OpenAI)** | Initial theoretical construction and drafting |
| **Le Chat (Mistral)** | Peer review simulation (v1.3 onwards) |

Multi-AI peer review simulations were conducted at each revision stage. All three independent reviewers (Claude, Gemini, Le Chat) identified the same four structural deficiencies in v1.2, which have been resolved in v1.5.2.

---

## Revision History

| Version | Date | Score | Key Changes |
|:---|:---|:---|:---|
| Ver. 1.2 | Apr 2026 | 61/100 | Initial integrated paper |
| Ver. 1.3 | May 2026 | 78/100 | M=0.5 flagged as empirical; t_stable defined; SLACS method described; §3.3 (D_eff justification) added |
| Ver. 1.4.1 | May 2026 | 86/100 | Figures embedded; Table 3 three-category notation; "among tested topologies" qualification |
| Ver. 1.5 | May 2026 | 89/100 | Variance/log-scale rationale restored; D_eff derivation expanded; Figure 2 cross-ref; data sources clarified |
| Ver. 1.5.1 | May 2026 | 90/100 | Figure 2 reference added to §2.3 |
| **Ver. 1.5.2** | **May 2026** | **91/100** | **§3.7 ADM correspondence expanded; GitHub URL added; OSF DOI added to all self-references** |

Scores are from Claude (Anthropic) strict peer review mode. Gemini and Le Chat independently scored v1.3 revisions at 85/100.

### Prior Research Versions

| Version | Date | Key Contribution |
|:---|:---|:---|
| Ver. 4.4 | Sep 2025 | Emergence and fragility of laws (2–8 node simulations) |
| Ver. 6.2.1 | Oct 2025 | Green-Kubo normalization (volume-normalized susceptibility) |
| Ver. 2.1 | Nov 2025 | 6 galaxies + SLACS gravitational lens verification |
| Ver. 3.9.11 | Dec 2025 | Analytic derivation of BTFR and $a_0$ via entropic matching |
| Ver. 5.1 | Jan 2026 | 7-galaxy parameter-free direct inversion |
| Ver. 7.0 | Jan 2026 | Conformal interface mathematical foundation |
| Ver. 8.1 | Jan 2026 | Two-Tier steady-state cosmology and Phase Metric (M) |
| Ver. 9.2 | Feb 2026 | 170 galaxies + universal scaling law ($D_{\text{eff}} \propto R^{1.38}$, combined Tier 1+2) |

All prior versions are available within this OSF project: [10.17605/OSF.IO/UJBPW](https://doi.org/10.17605/OSF.IO/UJBPW)

---

## Citation

```
Sasada, Y. (2026). QIC-S Theory: Universal Scaling Law of 170 Galaxies
and the Minimum Stable Unit of Causal Graphs (N=3). Version 1.5.2.
DOI: 10.17605/OSF.IO/UJBPW
```

---

## License

Paper and figures: CC BY 4.0 | Code: MIT
