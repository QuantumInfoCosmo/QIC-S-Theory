# QIC-S Theory: Universal Scaling Law of 170 Galaxies and the Minimum Stable Unit of Causal Graphs (N=3)

**Version 1.2 (Revised April 2026)** **DOI:** [10.17605/OSF.IO/UJBPW](https://doi.org/10.17605/OSF.IO/UJBPW)


---

## Overview

This repository contains the paper, data, verification codes, and figures for the Quantum Information Cosmology – Sasada (QIC-S) theory, Version 1.2.

QIC-S proposes that the universe is a collection of discrete causal graphs, where each galaxy possesses a unique effective Hamiltonian (Multi-Hamiltonian hypothesis). The theory explains galaxy rotation curves without invoking particle dark matter, instead attributing the observed "missing mass" to the gradient energy of a scalar information field $D_{\text{eff}}$.

### Key Results

- **Parameter-free verification:** 99.46% ± 2.53% agreement with 7 SPARC galaxies (zero free parameters).
- **Independent validation:** Consistency with SLACS strong gravitational lens convergence profiles.
- **Universal scaling law:** $D_{\text{eff}} \propto R^{1.573}$ ($R^2 = 0.945$), with trivial scaling $\alpha = 1.0$ strictly excluded (Bootstrap 95% CI [1.518, 1.627]).
- **Phase classification:** 78.2% Order Phase / 21.8% Chaos Phase across 170 SPARC galaxies.
- **Minimum stable causal unit:** N=3 Ring ($t_{\text{mix}} = 0.684$, converged at dt→0), proven via Edge-based Lindblad equation.

---

## Repository Structure

```
.
├── Sasada_QICS_Theory_v1_2.pdf        # Main paper (English, 6 pages)
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

### Chapter 5: Phase Classification & Scaling Law
```bash
# Verify all 170-galaxy statistics
python verify_chapter5.py

# Verify galaxy+filament combined regression (Ver. 9.2 comparison)
python verify_with_filaments.py

# Generate all certified numerical evidence
python chapter5_evidence.py
```
**Expected output:** 13/13 verification items PASS.

### Chapter 6: Edge-based Lindblad Mixing Times
```bash
# Run mixing time analysis for all topologies (dt=0.01)
python verify_chapter6_v3_dt001.py

# Run dt convergence test (dt = 0.05, 0.02, 0.01, 0.005, 0.001)
python dt_convergence_test.py
```
**Expected output:** N=3 Ring $t_{\text{mix}} = 0.69$ (dt=0.01), converging to 0.684 (dt=0.001). All qualitative checks PASS.

### Figure 2: Scaling Law Plot
```bash
# Generate improved Figure 2 (requires 4b_QIC_S_Result_N170.csv)
python generate_fig2_improved.py
```

---

## Data Description

### 4b_QIC_S_Result_N170.csv
| Column | Type | Description |
| :--- | :--- | :--- |
| Galaxy | string | Galaxy identifier (SPARC naming convention) |
| M | float | Phase Metric (log-variance of Hamiltonian gradient) |
| R | float | Characteristic scale [kpc] (outermost measured radius) |
| D_eff | float | Effective transport coefficient [kpc km/s] ($= R \times v$) |

**Source:** Derived from SPARC database (Lelli, McGaugh, & Schombert, 2016, AJ, 152, 157).

---

## Scope and Limitations
- All analyses are restricted to the **galactic scale (Tier 1)**.
- Application to galaxy cluster scales (Tier 2) is currently suspended to avoid circular reasoning from ΛCDM-biased mass estimates.
- The 2D-to-4D extension of the conformal interface mechanism remains an **unproven working hypothesis**.
- The mathematical bridge from discrete N=3 causal graphs to the macroscopic hydrodynamic limit has been **numerically confirmed** ($R^2 = 0.9843$ at N=32) but **not rigorously proven**.

---

## AI Disclosure
All theoretical ideas, physical interpretations, and theoretical frameworks are the sole responsibility of the author. AI systems were employed as auxiliary tools:
- **Claude (Anthropic):** Theoretical synthesis, independent numerical verification, and document preparation.
- **Gemini (Google):** Code development, simulation execution, and document preparation (LaTeX typesetting).
- **ChatGPT (OpenAI):** Initial theoretical construction and drafting.

---

### Prior Versions

| Version | Date | Key Contribution |
| :--- | :--- | :--- |
| Ver. 4.4 | Sep 2025 | Emergence and fragility of laws (2–8 node simulations) |
| Ver. 6.2.1 | Oct 2025 | Green-Kubo normalization correction |
| Ver. 2.1 | Nov 2025 | 6 galaxies + gravitational lens verification |
| Ver. 3.9.11 | Dec 2025 | Analytic derivation of BTFR and $a_0$ via entropic matching |
| Ver. 5.1 | Jan 2026 | 7-galaxy parameter-free derivation |
| Ver. 7.0 | Jan 2026 | Conformal interface mathematical foundation |
| Ver. 8.1 | Jan 2026 | Two-Tier steady-state cosmology |
| Ver. 9.2 | Feb 2026 | 170 galaxies + universal scaling law |
| **Ver. 1.2** | **Apr 2026** | **Integrated paper (this work)** |

## Prior Versions

| Version | Date | Key Contribution |
|---------|------|------------------|
| Ver. 4.4 | Sep 2025 | Emergence and fragility of laws (2–8 node simulations) |
| Ver. 2.1 | Nov 2025 | 6 galaxies + gravitational lens verification |
| Ver. 6.2.1 | Oct 2025 | Green-Kubo normalization correction |
| Ver. 5.1 | Jan 2026 | 7-galaxy parameter-free derivation |
| Ver. 7.0 | Jan 2026 | Conformal interface mathematical foundation |
| Ver. 8.1 | Jan 2026 | Two-Tier steady-state cosmology |
| Ver. 9.2 | Feb 2026 | 170 galaxies + universal scaling law |
| **Ver. 1.2** | **Apr 2026** | **Integrated paper (this work)** |

---

## Citation
```
Sasada, Y. (2026). QIC-S Theory: Universal Scaling Law of 170 Galaxies
and the Minimum Stable Unit of Causal Graphs (N=3). Version 1.2.
```

---

## License
Paper and figures: CC BY 4.0 / Code: MIT
