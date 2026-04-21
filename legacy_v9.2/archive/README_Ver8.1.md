# QIC-S Theory: Quantum Information Cosmology

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/Version-8.1-blue.svg)](https://github.com/QuantumInfoCosmo/QuantumInfoCosmo_NGC2403)
[![DOI](https://img.shields.io/badge/DOI-10.17605%2FOSF.IO%2FMBJN9-blue)](https://doi.org/10.17605/OSF.IO/MBJN9)

## Overview

**QIC-S (Quantum Information Cosmology - Sasada)** is a theoretical framework that explains galactic rotation curves without invoking Dark Matter particles. The theory treats spacetime as an information processing medium and proposes a **Two-Tier System**:

- **Tier 1 (Local)**: Galaxies as information processing units with local entropy production
- **Tier 2 (Global)**: Universal steady-state maintained by Interface Energy supply

### Key Features (Ver. 8.1 Scientific Update)

- **Zero-Parameter Prediction**: Uses only the theoretical acceleration scale a₀ = cH₀/2π ≈ 1.23×10⁻¹⁰ m/s²
- **Rigorous Data-Driven Analysis**: Removed randomized visualization. Now uses strictly physical data.
- **Phase Metric**: Quantitatively distinguishes mature (Order) vs. germinating (Chaos) galaxies using Log-Variance Metric.

---

## 📂 Repository Structure

The repository is organized as follows. Please ensure you have the `data/` directory populated with SPARC `.dat` files.

```
.
├── qics_analyzer.py          # [MAIN SCRIPT] The Python analysis tool (Ver. 2.1)
├── data/                     # [DATA FOLDER] Place SPARC galaxy files here
│   ├── NGC2403_rotmod.dat
│   ├── NGC6503_rotmod.dat
│   ├── UGC00128_rotmod.dat
│   ├── ID830_rotmod.dat
│   └── ...
├── QICS_Universality_NGC6503_UGC00128.png   # Figure 1: Order vs Order
├── QICS_Scientific_Result.png               # Figure 2: Order vs Chaos
├── QIC-S_Ver8_1_Sasada_2026.pdf             # Full Scientific Paper
└── README.md                                # This document
```

---

## Latest Results (Ver. 8.1)

### Phase Metric Analysis

| Galaxy | Type | Phase | Metric (M) | Interpretation |
|--------|------|-------|------------|----------------|
| NGC 6503 | Standard Spiral | 5 | 0.17 | Stable Order |
| UGC 128 | LSB Galaxy | 5 | 0.26 | Stable Order |
| NGC 2403 | SABcd Spiral | 5 | 0.30 | Stable Order |
| ID 830 | Transitional | 4 | 1.91 | Chaotic (Entropic Release) |

**Threshold**: M = 0.5 (Order/Chaos boundary)

### Key Figures

#### Figure 1: Universality of QIC-S across galaxy morphologies
![Universality](QICS_Universality_NGC6503_UGC00128.png)

*NGC 6503 (Standard Spiral, Metric: 0.17) vs UGC 128 (LSB Galaxy, Metric: 0.26). Both exhibit stable Hamiltonian landscapes despite different surface brightness.*

#### Figure 2: Phase transition from Order to Chaos
![Phase Transition](QICS_Scientific_Result.png)

*NGC 2403 (Phase 5, Metric: 0.30) vs ID 830 (Phase 4, Metric: 1.91). The chaotic landscape in ID 830 visualizes Entropic Release during the germinating phase.*

---

## 🚀 Quick Start

### 1. Requirements

Install the necessary Python libraries:

```bash
pip install numpy matplotlib
```

### 2. Prepare Data

Ensure your galaxy rotation curve data files (SPARC format `.dat`) are placed in the `data/` folder.

You can download SPARC data from: [http://astroweb.cwru.edu/SPARC/](http://astroweb.cwru.edu/SPARC/)

### 3. Run Analysis

Use the **`qics_analyzer.py`** script to perform the analysis.

**Example 1: Single Galaxy Analysis**

```bash
python qics_analyzer.py --file1 data/NGC2403_rotmod.dat
```

**Example 2: Comparison Analysis (Order vs Order)**

```bash
python qics_analyzer.py --file1 data/NGC6503_rotmod.dat --file2 data/UGC00128_rotmod.dat
```

**Example 3: Phase Transition Analysis (Order vs Chaos)**

```bash
python qics_analyzer.py --file1 data/NGC2403_rotmod.dat --file2 data/ID830_rotmod.dat
```

### Output

The script generates:

1. **Rotation curve plot**: Observed vs. QIC-S prediction (Zero-Parameter)
2. **Hamiltonian Landscape**: Polar visualization of H_eff gradient (No random noise)
3. **Phase Metric**: Quantitative classification displayed on the plot

---

## Statistical Validation (Ver. 5.1)

QIC-S achieves **99.46% ± 2.53% agreement** with observed rotation curves across 7 galaxies:

| Galaxy | Type | Agreement |
|--------|------|-----------|
| DDO 154 | Dwarf | 100.21% |
| NGC 2403 | SABcd | 98.46% |
| NGC 3198 | SBc | 100.49% |
| UGC 128 | LSB | 104.28% |
| NGC 5055 | SAbc | 98.25% |
| NGC 6503 | SAcd | 98.11% |
| IC 2574 | Irr | 96.43% |

---

## Theory Summary

### Fundamental Equation

The characteristic acceleration scale is derived from cosmological parameters:

```
a₀ = cH₀ / 2π ≈ 1.23 × 10⁻¹⁰ m/s²
```

### Phase Metric Definition

```python
M = Var(log(|∇H_eff| + ε))
```

- **M < 0.5**: Phase 5 (Mature, Ordered)
- **M > 0.5**: Phase 4 (Germinating, Chaotic)

---

## Version History

| Version | Date | Description |
|---------|------|-------------|
| 8.1 | 2026-01-26 | Scientific update: removed randomized visualization |
| 8.0.5 | 2026-01-25 | Two-Tier System formulation |
| 7.0 | 2026-01 | Conformal Interface theory integration |
| 5.1 | 2026-01-16 | 7-galaxy validation (99.46% agreement) |
| 3.9.11 | 2025-12-30 | CSH derivation of a₀ |

---

## Citation & Resources

### OSF Preprint

The full manuscript and supplementary materials are available on OSF:

**DOI:** [10.17605/OSF.IO/MBJN9](https://doi.org/10.17605/OSF.IO/MBJN9)

### BibTeX

If you use QIC-S in your research, please cite:

```bibtex
@article{Sasada2026QICS,
  author  = {Sasada, Yoshiaki},
  title   = {Quantum Information Cosmology (QIC-S): A Two-Tier Steady-State Universe},
  year    = {2026},
  version = {8.1},
  url     = {https://github.com/QuantumInfoCosmo/QuantumInfoCosmo_NGC2403},
  doi     = {10.17605/OSF.IO/MBJN9}
}
```

---

## Data Sources

- **SPARC Database**: [Lelli, McGaugh & Schombert (2016)](http://astroweb.cwru.edu/SPARC/)
- **ID 830**: Obuchi et al. (2026), ApJ, 997, 156

---

## References

1. Kokubo, M. & Harikane, Y. (2025). ApJ, 995, 24.
2. Obuchi, S. et al. (2026). ApJ, 997, 156.
3. Tudorache, M. N. et al. (2025). MNRAS, 544, 4306–4319.
4. Komatsu, S. et al. (2025). arXiv:2512.11045.
5. Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). SPARC Database. AJ, 152, 157.

---

## Acknowledgments

This research was assisted by AI systems (Claude for theoretical articulation and Gemini for numerical analysis). All physical interpretations and theoretical frameworks are the sole responsibility of the author.

---

## License

MIT License - See [LICENSE](LICENSE) for details.

---

## Contact

- **Author**: Yoshiaki Sasada
- **GitHub**: [@QuantumInfoCosmo](https://github.com/QuantumInfoCosmo)
