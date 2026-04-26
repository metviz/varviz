# VarViz typography tier system — single base_font drives Big/Medium/Small.
# Sourced by server.R and tests/test_typography.R. Multipliers per design doc
# docs/plans/2026-04-26-sk-revisions-design.md (Comment 3(ii), set 3b).

# Unified anchor for both live-app rendering and AppNote static export — the
# protein-landscape plot is built once via reactive() and the same theme
# objects are used for both pathways, so a single base_font applies. SK's
# Comment 3(ii) readability concern motivated bumping from 12 → 14.
vv_base_font <- 14       # pt
vv_big       <- vv_base_font * 1.00
vv_medium    <- vv_base_font * 0.75
vv_small     <- vv_base_font * 0.55
