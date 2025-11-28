# Tumor Mechanics ROM — Optimization & Tradeoff Analysis
This repository provides performance profiling, subsampling strategies, and tradeoff analysis between full‑order models (FOM), subsampled FOM, and ROM.  
It is an updated and streamlined version of the original *PODforTumors* repository.

## Repository Structure
- **`scripts/optimization_tradeoff.m`** — main analysis script (runs FOM, subsampled FOM, ROM, generates plots).
- **`src/`** — helper functions:
  - `get_damper*.m`, `getMechanicsMaps_2D.m`, `mech_matrix_build_2D.m`, `grad_matrix.m`
  - `buildBoundaries_2D.m`, `get_damper_reduced.m`, `augmentCellMaps_2D.m`
  - `getProjectionMatrix.m`, `getDisplacementProjection_2D.m`, `getStrainProjection_2D.m`
  - `getStressProjection_2D.m`, `buildStrainMat.m`, `run_tests.m`
  - `getMechanicsMaps_2D_LUonce.m`, `log_debug.m`
- **`data/Ex5_patient.mat`** — example dataset
- **`tests/`** — unit and smoke tests for reproducibility.
- **`profiling/`** — profiler outputs (PDFs, traces).
- **`results/`** — generated outputs:
  - `full_FOM_results.mat` — von Mises maps for FOM and subsampled runs
  - `optimization_log.txt` — run log with timings
  - `per_step_error.png`, `tradeoff_plot.png` — figures
  - `crash_dumps/` — saved crash dumps on exceptions
  - `run_env.mat` — environment snapshot (MATLAB version, OS, CPU cores)

## Inputs
- **`data/Ex5_patient.mat`** — example dataset containing:
  - `image_data`:  
    - `NTC1`, `NTC2`, `NTC3` — cell maps  
    - `Tissues` — tissue property map  
    - `BreastMask` — binary mask for boundary conditions
  - `schedule_info`:  
    - `times` — imaging time points  
    - `imagedims` — grid dimensions (used to derive spacing)
- **Derived parameters**:  
  - `h` — grid spacing (from `schedule_info.imagedims`)  
  - `bcs` — boundary conditions (from `BreastMask`)

## Outputs
- **`results/full_FOM_results.mat`** — von Mises maps for FOM and subsampled runs.
- **`results/optimization_log.txt`** — run log with timings.
- **`results/per_step_error.png`** — per‑step relative error plot.
- **`results/tradeoff_plot.png`** — tradeoff plot (average von Mises error vs. speedup).
- **`results/run_env.mat`** — environment snapshot (MATLAB version, OS, CPU cores).
- **`results/crash_dumps/`** — saved crash dumps on exceptions.
- **`profiling/`** — profiler outputs (PDFs, traces).

## How to Run
1. Add helper functions to MATLAB path. From repo root, run: `addpath(genpath(pwd))`.
2. Execute the main analysis script: `scripts/optimization_tradeoff.m`.
3. Inspect outputs in the results/ folder:
   - `optimization_log.txt` — run log with timings
   - `full_FOM_results.mat` — saved von Mises maps
   - Figures (per_step_error.png, tradeoff_plot.png)
4. Optional checks
   - Run `run_tests()` for sanity/unit tests and quick validation of core functionality (`smoke_test()`)

## Contact
Code contributors: Chase Christenson, Graham Pash, Casey Stowers, Tarini Thiagarajan
