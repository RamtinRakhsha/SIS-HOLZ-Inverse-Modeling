# SIS-HOLZ Inverse Modeling

MATLAB code, experimental phase data, and supplementary analysis supporting the manuscript:

**Generalized Inverse Modeling of SIS-HOLZ Strain Profiles**

**Authors:** Mana Norouzpour and Ramtin Rakhsha

## Overview

This repository contains the MATLAB implementation used for the generalized forward and inverse analysis of self-interference of split higher order Laue zone (SIS-HOLZ) data.

The framework represents the through-thickness projected displacement using surface-relaxed sine basis functions rather than prescribing a single displacement profile in advance.

The repository includes:

- MATLAB scripts used to generate the main manuscript figures;
- the experimental phase dataset used for the reanalysis;
- a supplementary synthetic robustness test for the experimental \(N=5\) inversion.

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
    └── Supplementary_False_Double_Hump_Test_N5.m
```

## Main figure scripts

### `Main_Figures/Figure1_Forward_Simulations.m`

Generates the generalized SIS-HOLZ forward simulations for representative symmetric and asymmetric through-thickness displacement profiles.

The figure includes prescribed projected displacement profiles, simulated normalized SIS-HOLZ interferograms, and corresponding relative phase traces.

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

Evaluates the sensitivity of the recovered displacement to uncertainty in the specimen thickness used in the inverse model.

### `Main_Figures/Figure4c_Instrumental_Phase_Sensitivity.m`

Evaluates the sensitivity of the recovered displacement to an incorrect instrumental phase-curvature parameter.

### `Main_Figures/Figure5_Experimental_Reanalysis_N5.m`

Performs the generalized phase-based inverse reanalysis of the previously reported SIS-HOLZ experimental dataset using five surface-relaxed sine basis functions and the published specimen thickness of 166 nm.

The fitted parameter set includes five displacement basis coefficients, a constant phase offset, and a residual linear phase ramp.

## Experimental data

### `Data/ExPhase.txt`

Contains the experimental phase profile used for the generalized experimental reanalysis.

### Data-path note

The current experimental MATLAB script expects `ExPhase.txt` to be accessible from the script location.

If the data file is retained only in the `Data/` folder, either copy `ExPhase.txt` temporarily into the same folder as the experimental script before running it, or update the data-file path in the script to point to:

```matlab
../Data/ExPhase.txt
```

## Supplementary analysis

### `Supplementary/Supplementary_False_Double_Hump_Test_N5.m`

Performs an additional synthetic robustness test of the \(N=5\) experimental inversion.

The test asks whether the same \(N=5\) inversion procedure can artificially generate a double-maximum displacement profile when the true displacement is the previously reported single-maximum 32.5 pm bell-shaped profile.

The default test uses 50 independent phase-noise realizations, a phase-noise standard deviation of \(\sigma_\phi=0.015\) rad, specimen thickness \(t=166\) nm, five displacement basis coefficients, a constant phase offset, a residual linear phase ramp, and three optimization starts per realization.

For the fixed random-number seed supplied with the script, the reported result was:

- **False double-maximum recoveries:** 0/50
- **Median recovered peak displacement:** 33.07 pm
- **16th–84th percentile interval:** 31.63–34.48 pm
- **Median phase RMSE:** 0.0139 rad

This supplementary test does not establish uniqueness of the experimental inverse solution. It specifically tests whether the selected \(N=5\) inversion procedure itself shows a tendency to generate a double-maximum profile from a known single-maximum truth under a phase-noise level comparable to the experimental residual.

## Requirements

The scripts were developed in MATLAB.

**Optimization Toolbox** is required for scripts that use `lsqnonlin`.

Some scripts use percentile calculations such as `prctile`; depending on the MATLAB release, this may require **Statistics and Machine Learning Toolbox**.

The supplementary false-double-maximum test runs serially and does not require Parallel Computing Toolbox.

## Running the scripts

1. Clone or download the repository.
2. Open MATLAB.
3. Navigate to the script of interest.
4. For the experimental reanalysis or supplementary test, make sure `ExPhase.txt` is accessible as described above.
5. Run the selected M-file.

Generated figures and numerical outputs are saved according to the output names specified within each script.

## Code availability

The repository is available at:

[https://github.com/RamtinRakhsha/SIS-HOLZ-Inverse-Modeling](https://github.com/RamtinRakhsha/SIS-HOLZ-Inverse-Modeling)

The repository provides the MATLAB scripts and experimental phase data associated with the analyses reported in the manuscript and its supplementary validation.

## Citation

If you use this repository, please cite the associated manuscript once the final journal citation becomes available:

> Mana Norouzpour and Ramtin Rakhsha, *Generalized Inverse Modeling of SIS-HOLZ Strain Profiles*.

The final journal citation and DOI can be added here after publication.

## Interpretation note

A single HOLZ reflection constrains only the displacement component projected along its reciprocal-lattice vector. The recovered quantity in these scripts is therefore a **through-thickness projected displacement profile**, not a complete three-dimensional displacement vector.
