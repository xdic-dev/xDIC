
Parallel Setting
================

Why Use Parallel Computing?
---------------------------
multiDIC can be computationally intensive. Parallel computing speeds up processes like image correlation and 3D reconstruction.

Setting Up Parallel Processing
------------------------------
1. Ensure MATLAB’s Parallel Computing Toolbox is installed.
2. Use `parpool` to start a parallel pool:
```matlab
parpool('local', 4);
```
3. Enable parallel execution by modifying `config.m`:
```matlab
use_parallel = true;
num_workers = 4;
```
