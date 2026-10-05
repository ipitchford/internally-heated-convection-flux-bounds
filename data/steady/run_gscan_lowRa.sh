#!/bin/bash
# Continuations in the period at fixed Ra (Figure 4c,d), as run; paths relative to this file.
cd "$(dirname "$0")/../../code"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 gamma_scan.py --mode rb --Ra 1e5 --start ../data/steady/rb_ns_X_1.000e+05.npy --G0 0.7315 --Gmin 0.55 --Gmax 1.8 --Nz 48 --modes_per_unit 100 --fac 1.04 --out ../data/steady/gscan_rb_1e5.json > ../data/steady/gscan_rb_1e5.log 2>&1
python3 gamma_scan.py --mode rb --Ra 177827.94100389228 --start ../data/steady/rb_ns_X_1.778e+05.npy --G0 0.5964 --Gmin 0.45 --Gmax 1.8 --Nz 48 --modes_per_unit 100 --fac 1.04 --out ../data/steady/gscan_rb_1.78e5.json > ../data/steady/gscan_rb_1.78e5.log 2>&1
python3 gamma_scan.py --mode rb --Ra 316227.7660168379 --start ../data/steady/rb_ns_X_3.162e+05.npy --G0 0.5050 --Gmin 0.40 --Gmax 1.8 --Nz 56 --modes_per_unit 110 --fac 1.04 --out ../data/steady/gscan_rb_3.16e5.json > ../data/steady/gscan_rb_3.16e5.log 2>&1
python3 gamma_scan.py --mode rb --Ra 562341.3251903491 --start ../data/steady/rb_ns_X_5.623e+05.npy --G0 0.4332 --Gmin 0.35 --Gmax 1.8 --Nz 64 --modes_per_unit 130 --fac 1.04 --out ../data/steady/gscan_rb_5.62e5.json > ../data/steady/gscan_rb_5.62e5.log 2>&1
