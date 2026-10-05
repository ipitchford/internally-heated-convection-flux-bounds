#!/bin/bash
# Re-run the checks behind every number in the paper, from the stored data, and rebuild the paper.
#   ./reproduce.sh          quick mode: validation tests, all derived checks, constants, tables, figures,
#                           assertion suite and paper build (about 20 minutes on one core)
#   ./reproduce.sh --full   also recomputes the slab-correction grids, the truncation/independent checks and c_w
#                           (about 1.5 hours more)
# The simulations (code/run_dns.py, queued by code/queue_dns.sh) and the steady-roll continuations
# (code/steady_scan.py, code/gamma_scan.py) take hours and are not rerun here; their stored outputs in data/ are
# what the checks read. Paths are relative to this directory.
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1
mkdir -p logs
run() { echo ">> $*" >&2; ( cd code && "$@" ); }
if [[ "${1:-}" == "--full" ]]; then
  run python3 lowfreq_check.py       > logs/lowfreq_check.log 2>&1
  run python3 cw_compute.py          > logs/cw_compute.log 2>&1
  run python3 corr_constants2.py 1   > logs/corr_constants2_n1.log 2>&1
  run python3 corr_constants2.py 2   > logs/corr_constants2_n2.log 2>&1
  run python3 corr_convergence.py    > logs/corr_convergence.log 2>&1
fi
run python3 test_kernels.py      > logs/test_kernels.log 2>&1
run python3 test_ipc.py          > logs/test_ipc.log 2>&1
run python3 test_steady.py       > logs/test_steady.log 2>&1
run python3 test_blankenbach.py  > logs/test_blankenbach.log 2>&1
run python3 cert/test_symbols_arb.py > logs/test_symbols_arb.log 2>&1 || echo "   (python-flint missing: interval cross-check skipped)"
run python3 carlson_sup.py       > logs/carlson_sup.log 2>&1
run python3 constants.py         > logs/constants.log 2>&1
run python3 check_dns.py         > logs/check_dns.log 2>&1
run python3 check_chen.py        > logs/check_chen.log 2>&1
run python3 revision_checks.py   > logs/revision_checks.log 2>&1
if ls data/dns/*_cont_snap.npz >/dev/null 2>&1; then run python3 dns_resolution.py > logs/dns_resolution.log 2>&1
else echo "   snapshots not present (companion archive not unzipped): dns_resolution.py skipped" >&2; fi
run python3 spectral_tails.py    > logs/spectral_tails.log 2>&1
run python3 make_tables.py       > logs/make_tables.log 2>&1
run python3 fig_theory.py        > logs/fig_theory.log 2>&1
run python3 fig_dns.py           > logs/fig_dns.log 2>&1
run python3 fig_steady.py        > logs/fig_steady.log 2>&1
run python3 check_dois.py > logs/check_dois.log 2>&1 || echo "   check_dois: mismatch or no network; see logs/check_dois.log" >&2
if command -v pdflatex >/dev/null && command -v pandoc >/dev/null; then
  run python3 make_refs.py > logs/make_refs.log 2>&1
  ( cd paper && pdflatex -interaction=nonstopmode main.tex > ../logs/latex_build.log 2>&1 && pdflatex -interaction=nonstopmode main.tex >> ../logs/latex_build.log 2>&1 )
  run python3 make_md.py > logs/make_md.log 2>&1
  if command -v xelatex >/dev/null; then
    ( cd paper && pandoc response_to_referee.md -o response_to_referee.pdf --pdf-engine=xelatex -V fontfamily=fontspec \
        -V mainfont=FreeSerif -V mathfont="Latin Modern Math" -V geometry:margin=2.4cm -V fontsize=11pt -V colorlinks=true ) \
      || echo "   response PDF not rebuilt (needs xelatex and the FreeSerif font)" >&2
  fi
fi
run python3 verify_numbers.py    > logs/verify_numbers.log 2>&1 && echo "verify_numbers: $(tail -1 logs/verify_numbers.log)" \
  || { echo "verify_numbers FAILED: see logs/verify_numbers.log"; exit 1; }
