#!/usr/bin/env bash
# Regenerate everything. Run from the repository root.
#
# Exit status: 0 only if every claim in 00_claims.tsv reconciles. 02_reconcile.R
# exits 1 on any MISMATCH / NOT COMPUTED; that status is carried through to the
# caller (it used to be swallowed by `|| true`, so CI saw success on a mismatch).
set -euo pipefail
Rscript analyses/repro/01_recompute.R
reconcile_status=0
Rscript analyses/repro/02_reconcile.R || reconcile_status=$?
Rscript analyses/repro/03_which_run.R
if [ "$reconcile_status" -ne 0 ]; then
  echo "run_all.sh: claims did NOT reconcile (02_reconcile.R exit $reconcile_status)" >&2
  exit "$reconcile_status"
fi
