#!/bin/bash
# The simulations of Section 5, in the order they were run. Each command writes data/dns/<name>.json (the summary that
# run_dns.py also prints), _hist.npz, _prof.npz, _state.npz and _snap.npz.
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1
# 1. Spin-up runs from small random perturbations of the conductive state (version 1 of the paper).
#    Their stored summaries (data/dns/ih_*_256x64.json) were written by an earlier version of IPC2D.run that took
#    equal-weight samples every ten steps; check_dns.py re-averages their windows in physical time from the stored
#    histories. Re-running these commands now uses the corrected averaging and overwrites those summaries.
# python3 run_dns.py --R 1e5 --bc ns --tend 3.0 --tavg 1.0
# python3 run_dns.py --R 1e5 --bc sf --tend 3.0 --tavg 1.0
# python3 run_dns.py --R 1e6 --bc ns --tend 2.0 --tavg 0.5
# python3 run_dns.py --R 1e6 --bc sf --tend 2.0 --tavg 0.5
# 2. Continuation runs used in the paper: restarted from the final states of the spin-up runs, averaged in physical
#    time over every step of the window [t_avg, t_end], with eight equal sub-blocks. Wall-clock times on one core:
#    3, 7, 10 and 13 minutes.
python3 run_dns.py --R 1e5 --bc ns --tend 7.0 --tavg 3.0 --restart ../data/dns/ih_ns_R1e+05_G4_256x64_state.npz --tag _cont
python3 run_dns.py --R 1e5 --bc sf --tend 7.0 --tavg 3.0 --restart ../data/dns/ih_sf_R1e+05_G4_256x64_state.npz --tag _cont
python3 run_dns.py --R 1e6 --bc ns --tend 5.0 --tavg 2.0 --restart ../data/dns/ih_ns_R1e+06_G4_256x64_state.npz --tag _cont
python3 run_dns.py --R 1e6 --bc sf --tend 5.0 --tavg 2.0 --restart ../data/dns/ih_sf_R1e+06_G4_256x64_state.npz --tag _cont
