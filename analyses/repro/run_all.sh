#!/usr/bin/env bash
# Regenerate everything. Run from the repository root.
set -euo pipefail
Rscript analyses/repro/01_recompute.R
Rscript analyses/repro/02_reconcile.R || true   # non-zero exit means claims disagree
Rscript analyses/repro/03_which_run.R
