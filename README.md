# Tumor Mechanics ROM — Optimization & Tradeoff Analysis
MATLAB Reduced Order Model (ROM) for tumor mechanics solver with performance profiling and optimization. (Updated version of the PODforTumors repo with improvements.)

## What this repo contains
- `optimization_tradeoff.m` — main analysis script (FOM, subsampled FOM, ROM, plots).
- `get_damper*.m`, `getMechanicsMaps_2D.m`, `mech_matrix_build_2D.m`, ... — helper functions.
- `data/Ex5_patient.mat` — example dataset (not included if large).
- `tests/` — unit and smoke tests.
- `profiling/results/` — profiler pdf.
- `crash_dumps/` — saved crash dumps on failure.

## Inputs
- `Ex5_patient.mat` — fields: `image_data` (NTC1, NTC2, NTC3, Tissues, BreastMask), `schedule_info` (times, imagedims).
- `h` (grid spacing) and `bcs` are derived from `schedule_info` and `BreastMask`.

## Outputs
- `full_FOM_results.mat` — saved VM maps for FOM and subsampled runs.
- `optimization_log.txt` — run log with timings.
- `profiling/results/` — profiler pdf.

## How to run
1. Add helper functions to MATLAB path: `addpath(genpath(pwd))`.
2. Run `optimization_tradeoff.m`.
3. Check `optimization_log.txt`

## Reproducibility
- Use `smoke_test()` to validate core functionality quickly.

## Tests & verification
See `tests/` for unit tests and verification cases.

## Contact
Code authors: Chase Christenson, Graham Pash, Casey Stowers, Tarini Thiagarajan
