# Self-check for the headless runner's clinician-evidence inputs and MANE reporting.
# Run from project root: Rscript analyses/tests/test_run_single_variant_inputs.R
#
# Why this exists. run_single_variant.R documents itself as driving
# build_variant_table() "so the output is what the UI would render". That was not
# true for two parameters: denovo_status was pinned to the literal "not_denovo" at
# the call site and segregation was never passed at all, so PS2, PM6 and the PP1
# tiers were unreachable from the headless path no matter what the caller wanted.
# The MANE Select transcript had the same shape of problem in reverse: the runner
# loads gene_data and dbNSFP reads are already MANE-indexed, but the transcript the
# predictors were read from never reached the printed block or the TSV, so an
# export could not say which isoform produced its scores.
#
# BMPR2 p.Glu376Lys is the anchor case. Published 2026 (PMID 42767215) as a
# recurrent de novo gain-of-function allele, it is exactly the variant whose call
# turns on evidence the headless path could not express.
#
# Everything here is offline. The engine functions are extracted from the shipped
# server.R rather than restated, so the test fails if the real points map drifts.

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
RUNNER <- resolve("analyses/run_single_variant.R")
SERVER <- resolve("server.R")

fails <- 0L; checks <- 0L
check <- function(label, ok) {
  ok <- isTRUE(ok)
  checks <<- checks + 1L
  if (!ok) fails <<- fails + 1L
  cat(if (ok) "  ok   " else "  FAIL ", label, "\n", sep = "")
}

runner <- readLines(RUNNER, warn = FALSE)
src    <- readLines(SERVER, warn = FALSE)

# ── Arg parsing: the two new flags must exist and round-trip ───────────────
# DEFAULTS decides which flags are known and which get numeric coercion, so the
# parser is driven, not mirrored.
eval(parse(text = paste(runner[grep("^DEFAULTS <- list\\(", runner):
                               (grep("^DEFAULTS <- list\\(", runner) - 1 +
                                which(grepl("\\)\\s*$", runner[grep("^DEFAULTS <- list\\(", runner):length(runner)]))[1])],
                        collapse = "\n")))
# parse_args validates against these, so they are part of the parser's contract.
eval(parse(text = runner[grep("^DENOVO_VALUES <- ", runner)][1]))
eval(parse(text = runner[grep("^SEG_VALUES <- ", runner)][1]))
pb <- grep("^parse_args <- function", runner)
pe <- pb - 1 + which(runner[pb:length(runner)] == "}")[1]
eval(parse(text = paste(runner[pb:pe], collapse = "\n")))

check("--denovo and --seg are known flags",
      all(c("denovo", "seg") %in% names(DEFAULTS)))
check("both default to the engine's own defaults",
      identical(DEFAULTS$denovo, "not_denovo") && identical(DEFAULTS$seg, "none"))

a <- parse_args(c("BMPR2", "p.Glu376Lys", "--denovo=denovo_confirmed", "--seg=pp1_strong"))
check("--denovo round-trips as a string", identical(a$denovo, "denovo_confirmed"))
check("--seg round-trips as a string",    identical(a$seg, "pp1_strong"))
check("omitted flags keep their defaults",
      identical(parse_args(c("BMPR2", "p.Glu376Lys"))$denovo, "not_denovo"))

# A typo in a value is worse than a typo in a flag name: the flag name errors
# today, but a bad value would silently fall through the engine's switch() and
# score nothing, which reads as "the evidence did not help".
bad_dn <- try(parse_args(c("BMPR2", "p.Glu376Lys", "--denovo=de_novo")), silent = TRUE)
check("an unrecognised --denovo value is rejected", inherits(bad_dn, "try-error"))
bad_sg <- try(parse_args(c("BMPR2", "p.Glu376Lys", "--seg=PP1")), silent = TRUE)
check("an unrecognised --seg value is rejected", inherits(bad_sg, "try-error"))
check("an unknown flag name still errors",
      inherits(try(parse_args(c("BMPR2", "p.Glu376Lys", "--nope=1")), silent = TRUE), "try-error"))

# The accepted sets must be the engine's, not a second list that can drift.
seg_sw <- src[grep("pp1_tag <- switch\\(segregation,", src):
              (grep("pp1_tag <- switch\\(segregation,", src) + 5)]
# Require the `=` so the `pp1_tag <- switch(...)` line itself is not read as an arm.
engine_seg <- c("none", trimws(sub("\\s*=.*$", "", grep("^\\s+pp1\\w*\\s*=", seg_sw, value = TRUE))))
check("--seg accepts exactly the engine's switch() arms",
      setequal(engine_seg, eval(parse(text = runner[grep("^SEG_VALUES <- ", runner)][1]))))
engine_dn <- c("not_denovo", "denovo_assumed", "denovo_confirmed")
check("--denovo accepts exactly the engine's three states",
      setequal(engine_dn, eval(parse(text = runner[grep("^DENOVO_VALUES <- ", runner)][1]))))

# ── Wiring: the call site must pass the parsed values, not literals ────────
# A flag that parses but is dropped at the call site is the exact bug this test
# exists for, so assert the call site, not just the parser.
call_beg <- grep("^vtbl <- build_variant_table\\(", runner)
call_end <- call_beg - 1 + which(grepl("^\\s*consurf_file_name = ", runner[call_beg:length(runner)]))[1]
call <- paste(runner[call_beg:call_end], collapse = " ")
check("denovo_status is wired to the parsed flag",
      grepl("denovo_status\\s*=\\s*p\\$denovo\\b", call))
check("denovo_status is no longer a hardcoded literal",
      !grepl('denovo_status\\s*=\\s*"', call))
check("segregation is passed at all",
      grepl("segregation\\s*=\\s*p\\$seg\\b", call))

# ── MANE: resolved from gene_data, printed, and written to the TSV ─────────
mb <- grep("^mane_tx_for_gene <- function", src)
me <- mb - 1 + which(src[mb:length(src)] == "}")[1]
eval(parse(text = paste(src[mb:me], collapse = "\n")))
gd <- local({ e <- new.env(); load(resolve("data/VarViz.RData"), envir = e); e$gene_data })

check("BMPR2 resolves to its MANE Select transcript",
      identical(mane_tx_for_gene("BMPR2", gd), "ENST00000374580"))
check("the version suffix is dropped, as dbNSFP reports it",
      !grepl("\\.", mane_tx_for_gene("BMPR2", gd)))
check("an unknown symbol yields NA, not an error",
      is.na(mane_tx_for_gene("NOT_A_GENE", gd)))
check("MANE_Transcript reaches the printed block",
      any(grepl('"MANE_Transcript"', runner)))
check("MANE_Transcript reaches the TSV header",
      any(grepl('"mane_transcript"\\s*=', runner)))

# The app's export is the file a clinician actually keeps, so it needs the same
# provenance. Its header already records every other parameter that moves a call;
# the transcript the predictors were read from belongs with them.
hb <- grep('^\\s*"varviz_version"\\s*=\\s*VARVIZ_VERSION,', src)
he <- hb - 1 + which(grepl('^\\s*\\)\\s*$', src[hb:length(src)]))[1]
app_hdr <- paste(src[hb:he], collapse = " ")
check("the app export header carries mane_transcript too",
      grepl('"mane_transcript"\\s*=', app_hdr))
check("it resolves through the same shipped helper, not a second lookup",
      grepl("mane_tx_for_gene\\(", app_hdr))
# Dropped-when-empty is the app's convention for absent fields, but an
# unresolved MANE entry is not absence: the reads fell back to the first array
# element. Both paths must say so rather than omit the line.
check("both paths write \"unresolved\" rather than dropping the field",
      grepl('"unresolved"', app_hdr) && any(grepl('"unresolved"', runner)))

# ── The anchor case, scored by the shipped classifier ──────────────────────
# Observed for BMPR2 p.Glu376Lys on 2026-10-11, VarViz 1.2.3, app defaults.
cb <- grep("^classify_acmg <- function", src)
ce <- cb - 1 + which(src[cb:length(src)] == "}")[1]
eval(parse(text = paste(src[cb:ce], collapse = "\n")))
eval(parse(text = src[grep("^PP1_PP4_CAP <- ", src)]))
tb <- grep("^ACMG_TAG_PTS <- c\\($", src); te <- tb - 1 + which(src[tb:length(src)] == ")")[1]
eval(parse(text = paste(src[tb:te], collapse = "\n")))

BMPR2 <- c("PM1_strong", "PM2", "PP2", "PP3")
check("public evidence alone reproduces the observed 7 points",
      classify_acmg(BMPR2)$pts == 7)
check("and lands Likely Pathogenic",
      classify_acmg(BMPR2)$class == "Likely Pathogenic")
# PS2 is the evidence the paper actually turns on: de novo in at least five trios.
check("confirmed de novo (PS2) carries it to Pathogenic",
      classify_acmg(c(BMPR2, "PS2"))$class == "Pathogenic")
check("PS2 is worth Strong, so 7 + 4 = 11",
      classify_acmg(c(BMPR2, "PS2"))$pts == 11)
# Assumed de novo is Moderate and must NOT be enough on its own — if this ever
# passes, the runner would let unconfirmed paternity reach a Pathogenic call.
check("assumed de novo (PM6) stops short at Likely Pathogenic",
      classify_acmg(c(BMPR2, "PM6"))$class == "Likely Pathogenic")
check("PM6 is worth Moderate, so 7 + 2 = 9",
      classify_acmg(c(BMPR2, "PM6"))$pts == 9)

# ── acmg_blind() must be called with its full argument list ───────────────
# Anchor: CASR p.Gly509Arg on the live app, 2026-10-11. Pass-Full is
# PM1 2 + PM2 1 + PP2 1 + PM5 2 + PP3_moderate 2 = 8. PP3's pre-cap tier is
# strong; the outcome-aware cap lowers it because a strong PP3 would reach 10
# and tip the full pass into Pathogenic. Under blinding PM1 (clinvar_hotspot)
# and PM5 go, the total can no longer reach Pathogenic, and PP3 must therefore
# return to strong: 1 + 1 + 4 = 6, Likely Pathogenic.
#
# acmg_blind() does that only when it is given am_pp3_precap. Called with two
# arguments it silently inherits the full pass's capped PP3 and reports 4,
# VUS-High. The app (server.R:7821, :8707) and the published harness
# (ps_final_harness.R:479) all pass the full set; run_single_variant.R did not,
# so a headless export disagreed with the app on exactly these variants.
eval(parse(text = readLines(resolve("analyses/lib/clinvar_blind.R"), warn = FALSE)))
for (fn in c("cap_am_pathogenic", "acmg_blind")) {
  b <- grep(paste0("^", fn, " <- function"), src)
  e <- b - 1 + which(src[b:length(src)] == "}")[1]
  eval(parse(text = paste(src[b:e], collapse = "\n")))
}

G509 <- c("PM1", "PM2", "PP2", "PM5", "PP3_moderate")
two <- acmg_blind(G509, "clinvar_hotspot")
four <- acmg_blind(G509, "clinvar_hotspot", "PP3_strong", 2L)
check("two-arg acmg_blind inherits the capped PP3 (the defect)",
      two$pts == 4 && two$classification == "VUS-High")
check("four-arg acmg_blind restores the pre-cap tier",
      four$pts == 6 && four$classification == "Likely Pathogenic")
check("and that is what the live app reports for CASR p.Gly509Arg",
      "PP3_strong" %in% four$tags && !any(c("PM1", "PM5") %in% four$tags))

# The regression guard: every call site passes four arguments. A parser that
# works while an argument is dropped at the call site is the failure mode here.
call_sites <- c(
  grep("acmg_blind\\(", runner, value = TRUE),
  grep("acmg_blind\\(", src, value = TRUE)
)
call_sites <- grep("^\\s*#", call_sites, value = TRUE, invert = TRUE)
check("at least three acmg_blind call sites found", length(call_sites) >= 3)
check("no call site passes only two arguments",
      !any(grepl("acmg_blind\\([^,)]+,[^,)]+\\)", call_sites)))
check("the runner reads the pre-cap column",
      any(grepl("ACMG_AM_PP3_Precap", runner)))
check("the runner reads the raised-from column",
      any(grepl("ACMG_AM_Raised_From", runner)))

# ── AlphaMissense source must be deterministic and declared ──────────────
# AlphaFold EBI answers 403 to download.file()'s default libcurl user-agent but
# 200 once HTTPUserAgent is set. Without it the per-gene substitutions CSV never
# arrives, the engine falls back to MyVariant/dbNSFP, and the two sources differ
# enough to change a call: CASR p.Gly509Arg scored AM 0.996 (MyVariant) with no
# PP3 raise recorded, against 0.9962 (AlphaFold) with a raise to strong on the
# live app -- Pass-Blind VUS-High 4 here, Likely Pathogenic 6 there. The fallback
# was silent, so an export could not say which source produced its scores.
check("the runner sets a user-agent before fetching",
      any(grepl("HTTPUserAgent", runner)))
check("the user-agent is set before the pipeline runs",
      min(grep("HTTPUserAgent", runner)) < min(grep("^mane_tx <- |^pfam_d", runner)))
check("AM_Source reaches the printed block",
      any(grepl('"AM_Source"', runner)))
check("AM_Source reaches the TSV header",
      any(grepl('"am_source"\\s*=', runner)))
check("a fallback to a non-AlphaFold source is warned about, not silent",
      any(grepl("AM source fallback|falling back", runner)))

cat(sprintf("\nrun_single_variant inputs: %d check(s), %d failure(s)\n",
            checks, fails))
quit(status = if (fails) 1L else 0L)
