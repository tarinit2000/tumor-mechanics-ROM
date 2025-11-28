# tumor-mechanics-ROM
MATLAB Reduced Order Model (ROM) for tumor mechanics solver with performance profiling and optimization. (Updated version of the PODforTumors repo with improvements.)

## Optimization + Tradeoff Analysis

How to run:
1. Add helper functions to MATLAB path.
2. Run `optimization_analysis.m` (script).
3. Outputs:
   - `optimization_log.txt` (run log)
   - `profiler_report/` (HTML profiler)
   - `crash_dumps/` (if any crashes)
   - `full_FOM_results.mat` (results)

Key numbers to report:
- Full solve time: <t_full_mech>
- Subsample (s=4) time: <t_sub_mech>
- ROM time: <t_ROM_mech>
- Final save I/O time: <t_io_save>
- Memory snapshot: <total_mem_MB>

Notes:
- The script runs a smoke test at startup.
- Profiler HTML saved to `profiler_report/` for hotspot analysis.
