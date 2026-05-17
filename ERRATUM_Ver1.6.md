# ERRATUM: QIC-S Theory Ver 1.6 (May 2026)

**Issue date:** May 16, 2026
**Author:** Yoshiaki Sasada
**Affected document:** *QIC-S Theory Ver 1.6 — Hydrodynamic Limit of Causal Networks: Convergence of Transport Coefficients and the Origin of Scaling Laws in Three-Dimensional Lattices*
**OSF DOI:** 10.17605/OSF.IO/KB75P (Ver 1.6)
**Status:** Section 2, parts of Sections 3 and 7 are **retracted** pending Ver 1.7. The Ver 1.6 document remains available as a frozen scientific record.

---

## Summary of the issue

After publication of Ver 1.6, an extension of the principal numerical experiment to L = 10 (N = 3000) was carried out using the methodology and code released alongside Ver 1.6. The result is incompatible with the central claim of §2 and §7 of Ver 1.6.

**Observed at L = 10:** D_vol = 4.92 × 10⁻⁴ (n = 5 ensemble).

**Combined with the L = 3–9 data of Ver 1.6:**

- Power-law fit (L = 5–10) gives D_vol ∝ L⁻¹·⁷⁷⁶, R² = 0.957 (compared with L⁻¹·⁶⁰⁹ for L = 5–9 alone). The decay exponent **worsens** with the addition of L = 10 rather than flattening.
- The L = 10 measurement is within 9% of the power-law-decay prediction (~4.5 × 10⁻⁴) and ~33% below the 1/N convergence prediction (~7.4 × 10⁻⁴).
- The intercept of the 1/N linear extrapolation has a t-statistic of 0.01, meaning **D_∞ = 0 cannot be statistically distinguished from D_∞ > 0** within the available data.

**Conclusion:** Within the specific model studied in Ver 1.6 — a three-dimensional cubic lattice with the N = 3 Ring as the unit cell, uniform-degree connectivity, periodic boundary conditions, disorder J_std = 0.3, and **pure unitary** time evolution — D_vol decays monotonically with system size and the existence of a strictly positive hydrodynamic limit D_∞ > 0 is **not supported** by the present data.

## What is retracted

The following statements in Ver 1.6 are retracted:

1. **Abstract:** "the volume-normalized transport coefficient D_vol = D_raw/N converges to a constant independent of system size N, establishing D_∞ > 0."
2. **§2.5 conclusion** insofar as it asserts that D_∞ > 0 is supported beyond the L = 3–9 fit range.
3. **§7 Conclusion (1):** "The hydrodynamic limit D_∞ > 0 is established."

The phrasing of §5 ("α = 1.573 follows naturally from the BTFR and the galaxy mass–size relation") is also weakened, because the bridge formula D_eff = D_GK × τ_R/τ_c requires D_GK to be finite in the limit; this requirement is not satisfied by the current model.

## What is NOT affected

The following results of Ver 1.5.2 and prior versions are **independent of the §2 calculation** and remain valid:

- The 99.46% ± 2.53% agreement with **7 SPARC galaxy** rotation curves via the Direct Inversion Method (Ver 5.1).
- The derivation of the Baryonic Tully–Fisher Relation from Cauchy Slice Holography and CNMG (Ver 3.9.11).
- The theoretical derivation of a_0 = cH_0/2π.
- The N = 3 Ring as the minimal stable causal unit via Edge-based Lindblad dynamics (t_mix = 0.684).

The following results introduced in Ver 1.6 also remain valid:

- §3: The redefinition of the local susceptibility χ_iso based on physical coordinates (the diagnosis that χ_old ∝ N is a 1D-residual artifact stands).
- §4: The peak in transport efficiency near r = g_inter / g_internal ≈ 2 is a real numerical observation in the model, although n = 5 ensembles warrant a higher-statistics rerun.
- §5.5: The absence of correlation between scaling residuals and the Phase Metric M (p = 0.66).
- §6: The period matrix analysis on N = 3 Ring.

## Why this matters and how Ver 1.7 will address it

The L = 10 result indicates that **the specific microscopic model of Ver 1.6 — uniform 3D lattice + pure unitary dynamics + N = 3 Ring unit cell — exhibits localization-like behavior even in three dimensions** under modest disorder (J_std = 0.3). This is reminiscent of, but not formally identified with, Anderson localization; further diagnostics are required to distinguish disorder-induced from topology-induced mechanisms.

Ver 1.7 will:

1. Replace the §2 D_∞ > 0 claim with the observed power-law decay and explicitly characterize it as a model limitation, not a theory failure.
2. Add diagnostic experiments separating (i) disorder-induced localization (testable by J_std = 0), (ii) Ring-topology-induced localization (testable with alternative unit cells), and (iii) PBC-induced artifacts (testable with open boundaries).
3. Re-examine the microscopic model requirements for a viable bridge to the macroscopic D_eff = R × v formula. The candidates under consideration are (a) introducing local Lindblad dissipation (already used elsewhere in QIC-S, e.g. Edge-based Lindblad for N = 3 Ring), and (b) modifying the network topology (e.g. adding long-range couplings).
4. Retain §3, §4, §5.5, and §6 with appropriate updates.

## On the integrity of this record

The L = 10 calculation was carried out using exactly the same code, parameters, and random-seed convention released with Ver 1.6, without modification. The discrepancy with the Ver 1.6 conclusion was therefore reproducible and not the result of a coding error or a parameter change. This erratum is issued as soon as the result was confirmed, in line with standard practice for promptly correcting the scientific record when an author can no longer endorse a published claim.

---

**Yoshiaki Sasada**
QIC-S Project, Machida, Tokyo
May 16, 2026
