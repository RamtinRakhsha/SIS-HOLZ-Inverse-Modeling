# SIS-HOLZ Inverse Modeling

MATLAB code, experimental phase data, and supplementary robustness analyses supporting the manuscript:

**Generalized Inverse Modeling of SIS-HOLZ Strain Profiles**

**Authors:** Mana Norouzpour and Ramtin Rakhsha

## Overview

This repository contains the MATLAB implementation used for the generalized forward and inverse analysis of self-interference of split higher order Laue zone (SIS-HOLZ) data.

The framework represents the through-thickness projected displacement using surface-relaxed sine basis functions rather than prescribing a single displacement profile in advance.

The repository includes:

- MATLAB scripts used to generate the main manuscript figures;
- the experimental phase dataset used for the reanalysis;
- supplementary robustness analyses of the representative experimental \(N=5\) inversion.

## Repository structure

```text
SIS-HOLZ-Inverse-Modeling/
│
├── README.md
│
├── Main_Figures/
│   ├── Figure1_Forward_Simulations.m
│   ├── Figure3_Synthetic_Validation.m
│   ├── Figure4a_Noise_Basis_Recoverability.m
│   ├── Figure4b_Thickness_Sensitivity.m
│   ├── Figure4c_Instrumental_Phase_Sensitivity.m
│   └── Figure5_Experimental_Reanalysis_N5.m
│
├── Data/
│   └── ExPhase.txt
│
└── Supplementary/
    ├── Supplementary_False_Double_Hump_Test_N5.m
    └── Experimental_N5_Thickness_Sensitivity.m
```

## Main figure scripts

### `Main_Figures/Figure1_Forward_Simulations.m`

Generates the generalized SIS-HOLZ forward simulations for representative symmetric and asymmetric through-thickness displacement profiles.

The figure includes:

- prescribed projected displacement profiles;
- simulated normalized SIS-HOLZ interferograms;
- corresponding relative phase traces.

### `Main_Figures/Figure3_Synthetic_Validation.m`

Performs the controlled synthetic inverse validation using three surface-relaxed sine basis functions.

The prescribed coefficients are:

\[
a_1=32~\mathrm{pm}, \qquad
a_2=-4~\mathrm{pm}, \qquad
a_3=2~\mathrm{pm}.
\]

The script evaluates recovery of the known displacement profile after addition of controlled complex-wave noise.

### `Main_Figures/Figure4a_Noise_Basis_Recoverability.m`

Evaluates the recoverability of individual basis functions as a function of basis-function order and complex-wave noise level.

### `Main_Figures/Figure4b_Thickness_Sensitivity.m`

Evaluates the sensitivity of the recovered displacement to uncertainty in specimen thickness in the synthetic robustness analysis.

### `Main_Figures/Figure4c_Instrumental_Phase_Sensitivity.m`

Evaluates the sensitivity of the recovered displacement to an incorrect instrumental phase-curvature parameter in the synthetic robustness analysis.

### `Main_Figures/Figure5_Experimental_Reanalysis_N5.m`

Performs the generalized phase-based inverse reanalysis of the previously reported SIS-HOLZ experimental dataset using five surface-relaxed sine basis functions and the published specimen thickness of 166 nm.

The fitted parameter set includes:

- five displacement basis coefficients;
- a constant phase offset;
- a residual linear phase ramp.

## Experimental data

### `Data/ExPhase.txt`

Contains the reconstructed experimental phase profile used for the generalized experimental reanalysis.

### Data-path note

The current experimental MATLAB scripts expect `ExPhase.txt` to be accessible from the script location.

If the data file is retained only in the `Data/` folder, either:

1. copy `ExPhase.txt` temporarily into the same folder as the script before running it; or
2. update the data-file path in the script to point to:

```matlab
../Data/ExPhase.txt
```

## Supplementary analyses

### `Supplementary/Supplementary_False_Double_Hump_Test_N5.m`

Performs an additional synthetic robustness test of the experimental \(N=5\) inversion.

The test asks whether the same \(N=5\) inversion procedure can artificially generate a double-maximum displacement profile when the true displacement is the previously reported single-maximum 32.5 pm surface-relaxed bell-shaped profile.

The reported calculation uses:

- 200 independent Gaussian phase perturbations;
- perturbation standard deviation \(\sigma_\phi=0.015\) rad;
- specimen thickness \(t=166\) nm;
- five displacement basis coefficients;
- a constant phase offset;
- a residual linear phase ramp;
- three optimization starts per realization.

The value \(\sigma_\phi=0.015\) rad is used as a perturbation scale comparable to the RMSE of the representative experimental \(N=5\) fit. It is not interpreted as an independent estimate of experimental measurement noise.

For the fixed random-number seed supplied with the script, the reported result was:

- **false double-maximum recoveries:** 0/200;
- **median recovered peak displacement:** 33.00 pm;
- **16th–84th percentile interval:** 31.56–34.50 pm;
- **median phase RMSE:** 0.0141 rad.

This test does not establish uniqueness or prove that the double-maximum structure obtained from the experimental phase is uniquely physical. It specifically tests whether random phase perturbations of this magnitude cause the selected \(N=5\) inversion itself to generate a spurious double maximum from a known single-maximum truth.

### `Supplementary/Experimental_N5_Thickness_Sensitivity.m`

Tests the sensitivity of the representative experimental \(N=5\) inverse solution to the assumed specimen thickness.

The full experimental inversion is repeated for:

\[
t/t_0 = 0.90,\;0.95,\;1.00,\;1.05,\;1.10,
\]

where:

\[
t_0 = 166~\mathrm{nm}.
\]

These values correspond to assumed specimen thicknesses of approximately:

\[
149.4,\;157.7,\;166.0,\;174.3,\;182.6~\mathrm{nm}.
\]

At each thickness, the script refits:

- \(a_1,\ldots,a_5\);
- the constant phase offset \(\phi_0\);
- the residual linear phase ramp \(c_1\).

The two-maximum displacement structure is retained across the tested thickness range. The recovered peak amplitude and depth of the central minimum vary with assumed thickness, while the phase residual increases substantially for assumed thicknesses above the published 166 nm value.

This test indicates that the representative double-maximum structure is not generated solely by fixing the specimen thickness at 166 nm. It does not establish uniqueness of the recovered displacement field.

## Requirements

The scripts were developed in MATLAB.

**Optimization Toolbox** is required for scripts that use `lsqnonlin`.

Some scripts use percentile calculations such as `prctile`; depending on the MATLAB release, this may require **Statistics and Machine Learning Toolbox**.

The supplementary false-double-maximum test runs serially and does not require Parallel Computing Toolbox.

## Running the scripts

1. Clone or download the repository.
2. Open MATLAB.
3. Navigate to the script of interest.
4. For scripts that use the experimental dataset, make sure `ExPhase.txt` is accessible as described above.
5. Run the selected M-file.

Generated figures, CSV files, and MAT files are saved according to the output names specified within each script.

## Supplementary reproducibility

The supplementary robustness analyses associated with the manuscript are reproduced by:

```text
Supplementary/Supplementary_False_Double_Hump_Test_N5.m
Supplementary/Experimental_N5_Thickness_Sensitivity.m
```

The first script generates the false-double-maximum robustness result and recovered-profile envelope.

The second script generates the experimental thickness-sensitivity profiles and the corresponding phase-RMSE dependence on assumed specimen thickness.

## Code availability

The repository is available at:

https://github.com/RamtinRakhsha/SIS-HOLZ-Inverse-Modeling

The repository provides the MATLAB scripts and experimental phase data associated with the analyses reported in the manuscript and its supplementary material.

## Citation

If you use this repository, please cite the associated manuscript once the final journal citation becomes available:

> Mana Norouzpour and Ramtin Rakhsha, *Generalized Inverse Modeling of SIS-HOLZ Strain Profiles*.

The final journal citation and DOI can be added here after publication.

## Interpretation note

A single HOLZ reflection constrains only the displacement component projected along its reciprocal-lattice vector. The recovered quantity in these scripts is therefore a **through-thickness projected displacement profile**, not a complete three-dimensional displacement vector.

The representative experimental \(N=5\) profile should not be interpreted as a uniquely determined true displacement field. The manuscript and supplementary analyses explicitly examine sensitivity to noise or perturbation level, specimen thickness, instrumental phase, and basis representation.
