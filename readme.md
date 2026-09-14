# ESMF — Bulk RNA-seq Deconvolution

**ESMF** is a deconvolution method for estimating cell-type proportions from bulk RNA-seq data using a single-cell-derived reference profile matrix. It formulates deconvolution as a non-negative matrix factorization (NMF) problem that jointly infers cell-type fractions and a refined reference signature. The method is benchmarked against a panel of bulk and single-cell deconvolution approaches.

This repository organizes the ESMF implementation together with its benchmarking, ablation, robustness, parameter-tuning, and visualization scripts into functional modules for reproducibility.

---

## Directory Structure

```
esmf_codes_organized/
├── 01_algorithms_and_utilities/        # Deconvolution engine & shared utilities (used by all other scripts)
├── 02_benchmark_scenarios/             # Benchmark drivers (3 scenario classes × multiple tissues)
├── 03_ablation_study/                  # Ablation study (latest versions only)
├── 04_expression_shift_e4/             # Expression-shift robustness experiment
├── 05_parameter_sweep/                 # Parameter sweeps / tuning / sensitivity
├── 06_parameter_sensitivity/           # Parameter-sensitivity benchmark experiments
├── 07_plotting/                        # Publication-grade plotting scripts
└── 08_environment_misc/                # Environment & miscellaneous
```

---

## Module Descriptions

### `01_algorithms_and_utilities/` — Core algorithms & utilities
The deconvolution engine and shared helper functions, called by every other script.
- `esmf_deconvolution.R` — Core ESMF algorithm (NMF-based cell-proportion deconvolution)
- `bulk_deconvolution.R` — Orchestrator (`deconv_lsc`) that dispatches the project's deconvolution methods (see source for the full supported list)
- `esmf_helper_functions.R` — ESMF preprocessing helpers (normalization, etc.)
- `basic_functions.R` — Basic utilities (sparse-matrix operations, etc.)
- `normalize_col_in_matrix.R` — Column-normalization utility for matrices
- `cibersort.R` — Third-party CIBERSORT method (v1.03)

### `02_benchmark_scenarios/` — Benchmark scenarios
Drivers for the paper's **3 scenario classes × multiple tissues** benchmark. Each script embeds the full pipeline: data loading → QC → marker selection → pseudo-bulk → multi-method evaluation → result saving.
- **Methods compared (13 in total)**: the proposed **ESMF** plus 10 bulk methods (`nnls`, `FARDEEP`, `RLR`, `DCQ`, `elastic_net`, `lasso`, `ridge`, `OLS`, `EPIC`, `DSA`) and 2 scRNA-seq methods (`MuSiC`, `CIBERSORTx`). The bulk set is defined by `methods_to_use` in each script; the two scRNA-seq methods are evaluated separately.
- **Cross-platform / database consistency**: `crossPlatformConsistency_*`
- **Generalization from single reference**: `GeneralizationFromSingleReference_*`
- **Multi-source reference integration**: `MultiSourceReferenceIntegration_*`
- Extended-cell-type variants: `*_Mm_Hypothalamus_6CT.R`, `*_Mm_Hypothalamus_7CT.R`
- Clinical validation: `clinical_validation_with_flow_cytometry.R` (real whole-blood RNA-seq + flow-cytometry ground truth)

### `03_ablation_study/` — Ablation study
Module ablations (regularization, prior, literature, etc.) and parameter-sensitivity analyses for each benchmark scenario.

### `04_expression_shift_e4/` — Expression-shift experiment
Evaluates ESMF robustness under systematic gene-expression shifts.
- `E4_*_ExpressionShift.R` — Expression-shift analysis (3 scenarios)
- `E4_*_CIBERSORTx_backfill.R` — CIBERSORTx result backfill (3 scenarios)
- `Export_CIBERSORTx_Input_E4_*.R` — Export CIBERSORTx inputs (Adipose / Brain / Pancreas)

### `05_parameter_sweep/` — Parameter sweeps / tuning / sensitivity
Sweeps and tuning of key ESMF hyperparameters (manuscript "Parameter Tuning Guidelines").
- **Tuning engine**: `tune_esmf.R` — 5-fold cross-validation, 128 parameter sets, selection by highest average Pearson correlation, no ground-truth leakage.
- **Sweeps**: `sweep_esmf_hm_adipose.R`, `sweep_esmf_hm_pancreas.R`, `sweep_esmf_mm_brain.R`

### `06_parameter_sensitivity/` — Parameter-sensitivity benchmark experiments
Full benchmark-pipeline scripts with an appended parameter-sensitivity section, evaluated against ground truth (the "real" sensitivity analyses; distinct from the older `Parameter_Sensitivity_Hm_*` copies).
- `multi_source_reference_integration_hm_adipose_parameter_sensitivity.R`
- `generalization_from_single_reference_hm_pancreas_parameter_sensitivity.R`
- `cross_platform_consistency_mm_brain_parameter_sensitivity.R`

### `07_plotting/` — Publication-grade plotting
Heatmaps, violin plots, ablation/tuning/sweep plots, cell-type-extension comparisons, E4-related figures, etc.
- Notable: `esmf-plots.R`, `plot_ablation.R`, `plot_ablation_contribution.R`, `plot_tune_esmf_heatmap.R`, `Plot_PerSample_Violin*.R`, `plot_cell_type_extension_comparison.R`, `plot_cell_type_expansion.R`, `e4_merge_cibersortx_and_plot.R`, `plot_fig4_fig5_no_stars.R`, `fix_figure3_violin_lasso_ridge.R`, `plot_esmf_sweep_2col_optima.R`

### `08_environment_misc/` — Environment & miscellaneous
- `list_package_versions.R` — Records the versions of dependency packages in the runtime environment

---

## Usage

1. Before running any analysis script, `source()` the engine and utilities in `01_algorithms_and_utilities/`.
2. Benchmark reproduction: run the scenario scripts in `02_benchmark_scenarios/`.
3. Ablation / robustness / tuning / sensitivity: see `03`, `04`, `05`, `06` respectively.
4. Result plotting: see `07_plotting/`.

## Dependencies

- R (≥ 4.0)
- Key packages: NMF, nnls, Matrix, data.table, ggplot2, etc. (exact versions in `08_environment_misc/list_package_versions.R`)
