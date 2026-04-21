# Quantum Information Cosmology with Scalar Field (QIC-S)
**Parameter-Free Derivation of Galaxy Rotation Curves**

[![GitHub](https://img.shields.io/badge/GitHub-Repo-black.svg?logo=github)](https://github.com/QuantumInfoCosmo/QuantumInfoCosmo_NGC2403)
[![OSF](https://img.shields.io/badge/OSF-10.17605%2FOSF.IO%2FYADKR-blue.svg)](https://osf.io/yadkr/)
![Version](https://img.shields.io/badge/Version-5.1-green.svg)

* **Version:** 5.1
* **Author:** Yoshiaki Sasada (Independent Researcher)
* **Status:** Complete / Open for Review

---

## 🌌 Overview

**Does the universe really need Dark Matter?**

This project presents **Quantum Information Cosmology with Scalar Field (QIC-S)**, a novel theoretical framework that resolves the galactic "missing mass" problem without postulating non-baryonic dark matter.

By treating gravity as an emergent phenomenon driven by the conservation of **Information Field Energy**, QIC-S naturally reproduces the flat rotation curves of galaxies.

### 🚀 Key Achievements

1.  **Parameter-Free**: Unlike $\Lambda$CDM (NFW profiles) or MOND ($a_0$), QIC-S requires **no tuning parameters**. The halo mass is mathematically derived solely from the visible baryonic distribution.

2.  **Universal Validity**: Tested against 7 diverse galaxies from the **SPARC database**, covering 2 orders of magnitude in mass ($10^{9.4} - 10^{11.3} M_\odot$).

3.  **High Precision**: Mean statistical agreement of **99.46%** between theoretical predictions and dynamically required masses.

---

## 📊 Results Summary

| Galaxy | Type | M_Newton ($M_\odot$) | M_QIC ($M_\odot$) | Agreement |
| :--- | :--- | :--- | :--- | :--- |
| **DDO 154** | Dwarf | $2.46 \times 10^9$ | $2.46 \times 10^9$ | **100.21%** |
| **NGC 2403** | Spiral | $6.92 \times 10^{10}$ | $6.81 \times 10^{10}$ | **98.46%** |
| **NGC 3198** | Spiral | $1.65 \times 10^{11}$ | $1.66 \times 10^{11}$ | **100.49%** |
| **UGC 128** | LSB | $1.67 \times 10^{11}$ | $1.75 \times 10^{11}$ | **104.28%** |
| **NGC 5055** | Spiral | $2.02 \times 10^{11}$ | $1.98 \times 10^{11}$ | **98.25%** |
| **NGC 6503** | Spiral | $5.57 \times 10^{10}$ | $5.47 \times 10^{10}$ | **98.11%** |
| **IC 2574** | Irregular | $7.72 \times 10^9$ | $7.44 \times 10^9$ | **96.43%** |

> **Global Mean Agreement: 99.46% ± 2.53%**

---

## 🧠 Theoretical Framework

QIC-S posits that the "missing mass" is a manifestation of field energy associated with a scalar information field, $D(x,t)$, which emerges from the information entropy of spacetime.

### Core Equations

* **Field-Density Correspondence:**
    $\rho_{\rm eff} = (\nabla D)^2$

* **Energy Conservation:**
    $M_{\rm eff}(r) = \int 4\pi r^2 (\nabla D)^2 dr$

* **Local Response:**
    The field $D$ is induced locally by the baryonic matter distribution.

---

## 📂 Repository Contents

### Paper
* **`Sasada_QIC-S_Paper_Ver5.1.pdf`** : The full scientific paper (Main Text).

### Analysis Code
* **`QICS_ZeroParam_Analysis.py`** : True zero-parameter analysis (Definitive Edition, Ver 4.0).
* **`QIC_S_MultiGalaxy_Analysis.py`** : Multi-galaxy analysis code (8 galaxies).
* **`QIC_S_NGC2403_Landscape.py`** : Hamiltonian landscape visualization for NGC 2403.

### Data
* **`data/`** : Rotation curve data from SPARC database.
* **`figures/`** : Analysis plots and figures.

---

## 🔗 Citation

If you use this theory or code in your research, please cite as:

```bibtex
@misc{Sasada2026QICS,
  title   = {Quantum Information Cosmology with Scalar Field (QIC-S): 
             Parameter-Free Derivation of Galaxy Rotation Curves 
             from Information-Theoretic Field Energy},
  author  = {Sasada, Yoshiaki},
  year    = {2026},
  publisher = {OSF},
  doi     = {10.17605/OSF.IO/YADKR},
  url     = {https://doi.org/10.17605/OSF.IO/YADKR}
}
```
