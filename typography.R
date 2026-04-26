# VarViz typography tier system — single base_font drives Big/Medium/Small.
# Sourced by server.R and tests/test_typography.R. Multipliers per design doc
# docs/plans/2026-04-26-sk-revisions-design.md (Comment 3(ii), set 3b).

vv_base_font <- 12       # live-app default (pt); cowplot export overrides via base_font argument
vv_big       <- vv_base_font * 1.00
vv_medium    <- vv_base_font * 0.75
vv_small     <- vv_base_font * 0.55

# Backward-compat alias — existing call sites resolve to medium tier until migrated (Task 14).
vv_axis_size <- vv_medium
