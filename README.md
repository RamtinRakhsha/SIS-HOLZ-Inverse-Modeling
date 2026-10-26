# SIS-HOLZ Inverse Modeling

MATLAB code, experimental phase data, and supplementary analyses supporting the manuscript:

**Generalized Inverse Modeling of SIS-HOLZ Strain Profiles**

**Authors:** Mana Norouzpour and Ramtin Rakhsha

## Overview

This repository contains the MATLAB implementation of a generalized forward and inverse framework for self-interference of split higher order Laue zone (SIS-HOLZ) analysis.

The main objective is to recover the through-thickness projected displacement profile without prescribing a single displacement shape in advance. The unknown projected displacement is represented using surface-relaxed sine basis functions and estimated by nonlinear least-squares inversion.

The repository includes:

- generalized SIS-HOLZ forward simulations;
- synthetic inverse validation;
- robustness tests for noise, specimen thickness, and instrumental phase;
- experimental reanalysis of previously reported SIS-HOLZ phase data;
- experimental model-order comparison for \(N=1\) through \(8\);
- a supplementary synthetic test evaluating whether an \(N=5\) inversion can generate a false double-maximum profile from a known single-maximum bell-shaped truth.

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
├── Model_Order/
│   └── Experimental_Model_Order_Sweep_N1_N8.m
│
├── Supplementary/
│   └── Supplementary_False_Double_Hump_Test_N5.m
│
├── Data/
│   └── ExPhase.txt
│
└── Supplementary_Material/
    └── SIS_HOLZ_Supplementary_Material.pdf
```

## MATLAB files

### Main figures

**`Main_Figures/Figure1_Forward_Simulations.m`**  
Generates the generalized forward-model comparison for representative symmetric and asymmetric displacement profiles, including projected displacement, simulated SIS-HOLZ interferograms, and relative phase traces.

**`Main_Figures/Figure3_Synthetic_Validation.m`**  
Performs the controlled synthetic inverse validation using three surface-relaxed sine basis functions and 1% complex-wave RMS noise. The prescribed coefficients are \(32\), \(-4\), and \(2\) pm.

**`Main_Figures/Figure4a_Noise_Basis_Recoverability.m`**  
Tests recovery of individual basis functions as a function of basis-function order and complex-wave noise level.

**`Main_Figures/Figure4b_Thickness_Sensitivity.m`**  
Quantifies the effect of assumed specimen-thickness error on the recovered displacement.

**`Main_Figures/Figure4c_Instrumental_Phase_Sensitivity.m`**  
Quantifies sensitivity of the recovered displacement to an incorrect instrumental phase-curvature parameter.

**`Main_Figures/Figure5_Experimental_Reanalysis_N5.m`**  
Performs the generalized \(N=5\) phase-based inverse reanalysis of the original experimental SIS-HOLZ dataset using the published specimen thickness of 166 nm.

### Model-order analysis

**`Model_Order/Experimental_Model_Order_Sweep_N1_N8.m`**  
Repeats the experimental phase inversion for \(N=1,\ldots,8\) surface-relaxed sine basis functions and reports phase RMSE, RSS, AIC, BIC, recovered peak displacement, number of meaningful local maxima, and fitted coefficients.

### Supplementary validation

**`Supplementary/Supplementary_False_Double_Hump_Test_N5.m`**  
Tests whether the experimental \(N=5\) inversion tends to create a false double-maximum displacement profile when the true profile is the previously reported 32.5 pm single-maximum bell.

The test uses:

- 50 independent phase-noise realizations;
- \(\sigma_\phi=0.015\) rad;
- the experimental diffraction-coordinate grid;
- specimen thickness \(t=166\) nm;
- five displacement basis coefficients;
- a constant phase offset and residual linear phase ramp;
- three optimization starts per realization.

For the fixed random-number seed supplied with the script, the test produced:

- **false double-maximum recoveries:** 0/50;
- **median recovered peak displacement:** 33.07 pm;
- **16th–84th percentile interval:** 31.63–34.48 pm;
- **median phase RMSE:** 0.0139 rad.

This test does not establish uniqueness of the experimental solution. Its purpose is specifically to determine whether the selected \(N=5\) inversion procedure itself shows a tendency to generate a double-maximum profile from a known single-maximum truth under a noise level comparable to the experimental residual.

## Experimental data

`Data/ExPhase.txt` contains the original experimental phase profile used for the experimental reanalysis and model-order comparison.

### Important data-path note

The current experimental MATLAB scripts were written to find `ExPhase.txt` in the **same directory as the script**. If the repository keeps the data only in the separate `Data/` folder shown above, either:

1. copy `Data/ExPhase.txt` into the folder containing the script before running it; or
2. update the `dataFile` / `load` line in the relevant script to point to `../Data/ExPhase.txt`.

This applies to:

- `Figure5_Experimental_Reanalysis_N5.m`
- `Experimental_Model_Order_Sweep_N1_N8.m`
- `Supplementary_False_Double_Hump_Test_N5.m`

Keeping the original scripts unchanged preserves the exact calculations used for the manuscript.

## Requirements

The scripts were developed in MATLAB.

**Optimization Toolbox** is required for scripts that use `lsqnonlin`, including the inverse-fitting and sensitivity analyses.

Some scripts use percentile calculations such as `prctile`; depending on the MATLAB release, this may require **Statistics and Machine Learning Toolbox**.

The supplementary false-double-maximum test runs serially and does **not** require Parallel Computing Toolbox.

## Reproducing the figures and analyses

1. Clone or download this repository.
2. Open MATLAB and navigate to the script of interest.
3. For experimental-data scripts, ensure `ExPhase.txt` is accessible as described in the data-path note above.
4. Run the desired M-file.
5. Generated figures, CSV files, and MAT files are written by the individual scripts according to the output names documented in each file.

### Publication settings for the robustness scripts

The Figure 4 robustness scripts include a `fastMode` control.

- `fastMode = true` uses a reduced grid and/or fewer realizations for rapid testing.
- `fastMode = false` uses the higher-resolution settings intended for the final publication calculation.

For reproduction of the final robustness calculations, use:

```matlab
fastMode = false;
```

in the applicable Figure 4 scripts.

## Supplementary material

The supplementary PDF accompanying the manuscript is provided under:

`Supplementary_Material/SIS_HOLZ_Supplementary_Material.pdf`

It describes the additional \(N=5\) false-double-maximum test and refers to the corresponding MATLAB code in the `Supplementary/` folder.

## Code availability

The repository is available at:

https://github.com/RamtinRakhsha/SIS-HOLZ-Inverse-Modeling

The repository is intended to provide the MATLAB implementation and supporting data needed to reproduce the numerical analyses associated with the manuscript.

## Citation

If you use this repository, please cite the associated manuscript once its final bibliographic information is available:

> Mana Norouzpour and Ramtin Rakhsha, *Generalized Inverse Modeling of SIS-HOLZ Strain Profiles*.

The final journal citation and DOI can be added here after publication.

## Notes on interpretation

A single HOLZ reflection constrains the displacement component projected along its reciprocal-lattice vector. The recovered quantity in these scripts is therefore a **through-thickness projected displacement profile**, not a complete three-dimensional displacement vector.

Higher-order basis functions provide greater flexibility but can also introduce increasingly fine or oscillatory structure. The model-order and robustness analyses should therefore be considered when interpreting detailed features of an inverse solution.
