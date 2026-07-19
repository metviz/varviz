# Gene-specific PP3 Calibration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compute a queried gene's own REVEL/AM/CADD Evidence Strength (LR+) from its ClinVar P/LP vs LB/B variants, show it against the genome-wide Pejaver default on the variant card, and allow a gated opt-in override of the PP3 tier.

**Architecture:** A pure, network-free stats core (`gene_calib_stats.R`, sourced by server.R like the existing `typography.R`) holds all LR+ / CI / tier / sweep math and is unit-tested with `Rscript` + `stopifnot`. Thin network functions in server.R fetch the benign ClinVar arm and batch-fetch dbNSFP scores. An assembly function takes the two arms *as data* (so it is unit-testable with fakes), applies leave-one-out for the queried variant, and returns a structured result the card renders. Override is advisory-first, session-only, gated by N.

**Tech Stack:** R, Shiny (`server.R`/`ui.R`), `httr` + `jsonlite` (already used), NCBI E-utilities, MyVariant.info. Tests are standalone `Rscript` files using `stopifnot` — the repo convention (`test_acmg_posterior.R`), no test framework.

## Global Constraints

- Evidence Strength (LR+ / OddsPath) is **prior-free**; prior_P enters only at §4.1 posterior, never here. Verbatim invariant from the spec.
- OddsPath tier boundaries: Supporting ≥ 2.0813, Moderate ≥ 4.33, Strong ≥ 18.7, Very strong ≥ 350. C = 2.0813 = 350^(1/8). Must match `acmg_posterior` (server.R).
- Min N per arm for high confidence = **20**. Below 20: still show numbers + CI, flag "low confidence", override requires explicit confirm. Compute requires both arms ≥ 1.
- Min-sensitivity floor default = **0.90**, user-settable in the parameters panel.
- Default behavior is **advisory only** — never mutate `acmg_res`/PP3 unless the user opts in and gating passes. No auto-override.
- Queried variant that is in ClinVar: **leave-one-out** from its calibration arm; surface its own ClinVar assertion as a **separate** line (no double-count).
- No local dbNSFP file. Bulk scores via MyVariant.info batch POST only.
- Back up `server.R` and `ui.R` to timestamped `.bak` in the same dir before editing them (project rule).
- Predictor cut-points (as encoded in VarViz today): REVEL 0.644/0.773/0.932; CADD 28.1/35; AM 0.564 (supporting only).

---

### Task 1: Pure stats core — LR+, Haldane, Katz CI, tier map, gene-optimal sweep

**Files:**
- Create: `gene_calib_stats.R`
- Test: `test_gene_calib_stats.R`

**Interfaces:**
- Produces:
  - `oddspath_tier(lr) -> character` one of "None"/"Supporting"/"Moderate"/"Strong"/"Very strong"
  - `gene_lr_plus(pos, neg, t) -> list(lr, lo, hi, sens, spec, tp, fp, n_pos, n_neg)` — `pos`/`neg` numeric score vectors, higher = more damaging; `sens`/`spec` are raw (uncorrected), `lr`/`lo`/`hi` use Haldane +0.5.
  - `gene_optimal(pos, neg, min_sens=0.90) -> list from gene_lr_plus + $t` (or `NULL` if none clears the floor)
  - `PREDICTORS` — named list; each has `$col` and `$cuts` (named numeric of tier→threshold)

- [ ] **Step 1: Write the failing test**

```r
# test_gene_calib_stats.R  — run: Rscript test_gene_calib_stats.R
source("gene_calib_stats.R")

# oddspath_tier boundaries (reuse §4.1 C)
stopifnot(
  oddspath_tier(11.83) == "Moderate",   # the Bindra REVEL>=0.7 value
  oddspath_tier(2.5)   == "Supporting",
  oddspath_tier(20)    == "Strong",
  oddspath_tier(400)   == "Very strong",
  oddspath_tier(1.5)   == "None"
)

# perfect separation -> finite LR+ via Haldane, not Inf
r <- gene_lr_plus(pos = c(0.9,0.95,0.99), neg = c(0.1,0.2,0.05), t = 0.5)
stopifnot(is.finite(r$lr), r$lr > 1, r$sens == 1, r$spec == 1)

# hand 2x2: pos>=t: 8/10 (sens .8), neg>=t: 2/10 (spec .8) -> LR+ approx .8/.2 = 4 (Haldane shifts slightly)
pos <- c(rep(1,8), rep(0,2)); neg <- c(rep(1,2), rep(0,8))
r2 <- gene_lr_plus(pos, neg, t = 1)
stopifnot(abs(r2$lr - 4) < 0.6, r2$lo < r2$lr, r2$hi > r2$lr, r2$tp == 8, r2$fp == 2)

# gene_optimal respects the sensitivity floor
pos3 <- c(0.6,0.7,0.8,0.9,0.95); neg3 <- c(0.1,0.2,0.3,0.65,0.72)
opt_hi <- gene_optimal(pos3, neg3, min_sens = 0.90)   # must keep >=90% of pos -> low t
opt_lo <- gene_optimal(pos3, neg3, min_sens = 0.50)   # may pick higher t, higher LR+
stopifnot(!is.null(opt_hi), opt_hi$sens >= 0.90 - 1e-9,
          opt_lo$lr >= opt_hi$lr)

# floor impossible -> NULL
stopifnot(is.null(gene_optimal(pos3, neg3, min_sens = 1.01)))

cat("gene_calib_stats: all checks pass\n")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript test_gene_calib_stats.R`
Expected: FAIL — `cannot open file 'gene_calib_stats.R'` / `could not find function "oddspath_tier"`.

- [ ] **Step 3: Write minimal implementation**

```r
# gene_calib_stats.R — pure Evidence-Strength / LR+ helpers.
# No network, no shiny. Sourced by server.R (like typography.R).
# Evidence Strength (LR+/OddsPath) is prior-free; prior enters only at acmg_posterior.

# OddsPath (LR+) -> Tavtigian Evidence-Strength tier. Boundaries = prior-free odds,
# C = 2.0813 = 350^(1/8), consistent with acmg_posterior in server.R.
oddspath_tier <- function(lr) {
  if (!is.finite(lr) || lr <= 1) return("None")
  if (lr >= 350)    return("Very strong")
  if (lr >= 18.7)   return("Strong")
  if (lr >= 4.33)   return("Moderate")
  if (lr >= 2.0813) return("Supporting")
  "None"
}

# LR+ at threshold t. pos/neg = numeric score vectors (P/LP arm, LB/B arm),
# higher score = more damaging. Haldane-Anscombe +0.5 on the 2x2 for lr/CI;
# raw (uncorrected) sens/spec returned for the sensitivity floor.
gene_lr_plus <- function(pos, neg, t) {
  pos <- pos[!is.na(pos)]; neg <- neg[!is.na(neg)]
  n_pos <- length(pos); n_neg <- length(neg)
  tp <- sum(pos >= t); fn <- n_pos - tp
  fp <- sum(neg >= t); tn <- n_neg - fp
  a <- tp + 0.5; b <- fn + 0.5; c <- fp + 0.5; d <- tn + 0.5   # Haldane cells
  sens_c <- a / (a + b); spec_c <- d / (c + d)
  lr <- sens_c / (1 - spec_c)
  var_ln <- (1 - sens_c)/(sens_c*(a + b)) + spec_c/((1 - spec_c)*(c + d))
  se <- sqrt(var_ln)
  list(lr = lr,
       lo = exp(log(lr) - 1.96*se),
       hi = exp(log(lr) + 1.96*se),
       sens = if (n_pos) tp / n_pos else NA_real_,
       spec = if (n_neg) tn / n_neg else NA_real_,
       tp = tp, fp = fp, n_pos = n_pos, n_neg = n_neg)
}

# gene-optimal threshold: max LR+ subject to raw sens >= min_sens.
# candidate thresholds = distinct pos scores (exact grid). NULL if none clears floor.
gene_optimal <- function(pos, neg, min_sens = 0.90) {
  pos <- pos[!is.na(pos)]; neg <- neg[!is.na(neg)]
  if (length(pos) == 0) return(NULL)
  best <- NULL
  for (t in sort(unique(pos))) {
    r <- gene_lr_plus(pos, neg, t)
    if (is.na(r$sens) || r$sens + 1e-9 < min_sens) next
    if (is.null(best) || r$lr > best$lr) { r$t <- t; best <- r }
  }
  best
}

# Predictor config: column name in the scored arm + Pejaver cut-points (tier -> threshold).
PREDICTORS <- list(
  REVEL = list(col = "revel", cuts = c(Supporting = 0.644, Moderate = 0.773, Strong = 0.932)),
  AM    = list(col = "am",    cuts = c(Supporting = 0.564)),
  CADD  = list(col = "cadd",  cuts = c(Supporting = 28.1, Moderate = 35))
)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript test_gene_calib_stats.R`
Expected: PASS — `gene_calib_stats: all checks pass`.

- [ ] **Step 5: Commit**

```bash
git add gene_calib_stats.R test_gene_calib_stats.R
git commit -m "feat(calib): pure LR+/OddsPath stats core for gene-specific PP3"
```

---

### Task 2: Parse ClinVar `name` → dbNSFP-style hgvsp (join key)

**Files:**
- Modify: `gene_calib_stats.R` (append)
- Test: `test_gene_calib_stats.R` (append)

**Interfaces:**
- Produces: `parse_clinvar_hgvsp(name) -> character` like `"p.R130G"`, or `NA` if no simple missense `p.Xxx###Yyy` is present. Vectorized over `name`.

- [ ] **Step 1: Write the failing test** (append before the final `cat(...)`)

```r
# hgvsp parse: ClinVar name string -> dbNSFP hgvsp (1-letter)
stopifnot(
  parse_clinvar_hgvsp("NM_000314.8(PTEN):c.388C>G (p.Arg130Gly)") == "p.R130G",
  parse_clinvar_hgvsp("NM_000257.4(MYH7):c.1988G>A (p.Arg663His)") == "p.R663H",
  is.na(parse_clinvar_hgvsp("NM_000314.8(PTEN):c.209+1G>A")),        # splice, no p.
  identical(
    parse_clinvar_hgvsp(c("x (p.Gly12Asp)", "no-protein", "y (p.Ter100Cys)")),
    c("p.G12D", NA_character_, NA_character_))                        # Ter ref not a missense start
)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript test_gene_calib_stats.R`
Expected: FAIL — `could not find function "parse_clinvar_hgvsp"`.

- [ ] **Step 3: Write minimal implementation** (append to `gene_calib_stats.R`)

```r
# 3-letter -> 1-letter amino acid. Ter included for alt (nonsense), excluded as a
# missense ref (a p.Ter### start is not a missense variant).
.AA3TO1 <- c(Ala="A",Arg="R",Asn="N",Asp="D",Cys="C",Gln="Q",Glu="E",Gly="G",
             His="H",Ile="I",Leu="L",Lys="K",Met="M",Phe="F",Pro="P",Ser="S",
             Thr="T",Trp="W",Tyr="Y",Val="V",Ter="*")

# ClinVar `name` -> dbNSFP hgvsp like "p.R130G"; NA unless a clean missense p.Xxx###Yyy.
parse_clinvar_hgvsp <- function(name) {
  vapply(name, function(nm) {
    if (is.na(nm)) return(NA_character_)
    g <- regmatches(nm, regexec("p\\.([A-Za-z]{3})([0-9]+)([A-Za-z]{3})", nm))[[1]]
    if (length(g) != 4) return(NA_character_)
    ref <- .AA3TO1[[g[2]]]; alt <- .AA3TO1[[g[4]]]
    if (is.null(ref) || is.null(alt) || ref == "*") return(NA_character_)
    paste0("p.", ref, g[3], alt)
  }, character(1), USE.NAMES = FALSE)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript test_gene_calib_stats.R`
Expected: PASS — `gene_calib_stats: all checks pass`.

- [ ] **Step 5: Commit**

```bash
git add gene_calib_stats.R test_gene_calib_stats.R
git commit -m "feat(calib): parse ClinVar name to dbNSFP hgvsp join key"
```

---

### Task 3: Assembly — arms in, structured calibration out (LOO, validation, optimal)

**Files:**
- Modify: `gene_calib_stats.R` (append)
- Test: `test_gene_calib_stats.R` (append)

**Interfaces:**
- Consumes: `gene_lr_plus`, `gene_optimal`, `oddspath_tier`, `PREDICTORS` (Task 1).
- Produces:
  `gene_calibration(path_arm, benign_arm, query_hgvsp=NULL, min_sens=0.90) -> list`
  where `path_arm`/`benign_arm` are data.frames with columns `hgvsp, revel, am, cadd`.
  Returns:
  ```
  list(
    loo_applied = logical,
    predictors = list( REVEL = list(
        n_pos, n_neg, confident,           # confident = min(n_pos,n_neg) >= 20
        validation = data.frame(tier_cut, threshold, lr, lo, hi, gene_tier,
                                pejaver_tier, agree),
        optimal = list(threshold, lr, lo, hi, gene_tier) | NULL
      ), AM = ..., CADD = ... )
  )
  ```

- [ ] **Step 1: Write the failing test** (append before final `cat`)

```r
mk <- function(hg, rv) data.frame(hgvsp = hg, revel = rv,
                                  am = rv, cadd = rv*40, stringsAsFactors = FALSE)
path_arm   <- mk(sprintf("p.A%dV", 1:30), c(rep(0.95, 25), rep(0.60, 5)))
benign_arm <- mk(sprintf("p.G%dS", 1:30), c(rep(0.10, 27), rep(0.80, 3)))

res <- gene_calibration(path_arm, benign_arm, min_sens = 0.80)
stopifnot(
  res$loo_applied == FALSE,
  res$predictors$REVEL$n_pos == 30, res$predictors$REVEL$n_neg == 30,
  res$predictors$REVEL$confident == TRUE,                       # both arms >= 20
  nrow(res$predictors$REVEL$validation) == 3,                   # 3 REVEL cut-points
  !is.null(res$predictors$REVEL$optimal),
  res$predictors$REVEL$optimal$lr > 1
)

# LOO: dropping a P/LP member reduces n_pos by 1
res2 <- gene_calibration(path_arm, benign_arm, query_hgvsp = "p.A1V", min_sens = 0.80)
stopifnot(res2$loo_applied == TRUE, res2$predictors$REVEL$n_pos == 29)

# small benign arm -> not confident, still computes
res3 <- gene_calibration(path_arm, benign_arm[1:8, ], min_sens = 0.80)
stopifnot(res3$predictors$REVEL$confident == FALSE,
          !is.null(res3$predictors$REVEL$validation))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript test_gene_calib_stats.R`
Expected: FAIL — `could not find function "gene_calibration"`.

- [ ] **Step 3: Write minimal implementation** (append to `gene_calib_stats.R`)

```r
# Assemble per-predictor validation + gene-optimal from two scored arms.
# path_arm/benign_arm: data.frame(hgvsp, revel, am, cadd). Pure (no network).
gene_calibration <- function(path_arm, benign_arm, query_hgvsp = NULL, min_sens = 0.90) {
  loo <- FALSE
  if (!is.null(query_hgvsp) && !is.na(query_hgvsp)) {
    if (query_hgvsp %in% path_arm$hgvsp)   { path_arm   <- path_arm[path_arm$hgvsp   != query_hgvsp, , drop = FALSE]; loo <- TRUE }
    if (query_hgvsp %in% benign_arm$hgvsp) { benign_arm <- benign_arm[benign_arm$hgvsp != query_hgvsp, , drop = FALSE]; loo <- TRUE }
  }
  preds <- lapply(names(PREDICTORS), function(pname) {
    cfg <- PREDICTORS[[pname]]
    pos <- path_arm[[cfg$col]];  neg <- benign_arm[[cfg$col]]
    pos <- pos[!is.na(pos)];     neg <- neg[!is.na(neg)]
    val <- do.call(rbind, lapply(names(cfg$cuts), function(tier) {
      t <- cfg$cuts[[tier]]; r <- gene_lr_plus(pos, neg, t)
      gt <- oddspath_tier(r$lr)
      data.frame(tier_cut = tier, threshold = t, lr = r$lr, lo = r$lo, hi = r$hi,
                 gene_tier = gt, pejaver_tier = tier,
                 agree = identical(gt, tier), stringsAsFactors = FALSE)
    }))
    opt <- gene_optimal(pos, neg, min_sens = min_sens)
    opt_out <- if (is.null(opt)) NULL else
      list(threshold = opt$t, lr = opt$lr, lo = opt$lo, hi = opt$hi,
           gene_tier = oddspath_tier(opt$lr))
    list(n_pos = length(pos), n_neg = length(neg),
         confident = (length(pos) >= 20 && length(neg) >= 20),
         validation = val, optimal = opt_out)
  })
  names(preds) <- names(PREDICTORS)
  list(loo_applied = loo, predictors = preds)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript test_gene_calib_stats.R`
Expected: PASS — `gene_calib_stats: all checks pass`.

- [ ] **Step 5: Commit**

```bash
git add gene_calib_stats.R test_gene_calib_stats.R
git commit -m "feat(calib): assemble per-predictor validation + gene-optimal with LOO"
```

---

### Task 4: Benign ClinVar arm fetch (server.R)

**Files:**
- Modify: `server.R` (add `extract_clinvar_benign` near `extract_clinvar` at 824; add `api_cache$gene_calib` near the namespace block ~72–82)
- Backup: `server.R.bak.<ts>.pre-4.2benign` before editing

**Interfaces:**
- Consumes: existing `cache_get`/`cache_set`, `httr`, `jsonlite`, and the esummary parse pattern at server.R:867–1148.
- Produces: `extract_clinvar_benign(gene_name) -> data.frame` with the same columns `extract_clinvar` returns (incl. `name`, `ClinicalSignificance`), but for LB/B variants.

- [ ] **Step 1: Back up + add cache namespace**

```bash
cp server.R "server.R.bak.$(date +%Y%m%d_%H%M%S).pre-4.2benign"
```

Add to the `api_cache$...` block (server.R ~72–82):

```r
api_cache$gene_calib  <- list()   # gene_name -> assembled gene_calibration() result
```

- [ ] **Step 2: Add the benign fetch function**

Immediately after `extract_clinvar` closes (find its closing brace after 1163), add a sibling that reuses the *same* esummary parse by delegating to a small shared helper. Minimal-risk approach: `extract_clinvar_benign` is `extract_clinvar` with the search queries swapped to benign clinsig. Copy `extract_clinvar`, rename, and change **only** the `search_queries` vector to:

```r
    search_queries <- c(
      paste0(gene_name, '[gene] AND ("benign"[clinsig] OR "likely benign"[clinsig]) AND "homo sapiens"[orgn]'),
      paste0(gene_name, '[gene] AND (benign[clinsig]) AND human[orgn]')
    )
```

and change the cache namespace it reads/writes from `"clinvar"` to `"clinvar_benign"` (add `api_cache$clinvar_benign <- list()` to the namespace block). Leave all esummary parsing identical so `ClinicalSignificance`, `name`, `clinvar_goldstar` populate the same way.

- [ ] **Step 3: Parse-check**

Run: `Rscript -e 'invisible(parse("server.R")); cat("parses OK\n")'`
Expected: `parses OK`.

- [ ] **Step 4: Live smoke test** (network — a gene with benign ClinVar entries)

Run:
```bash
Rscript -e 'source("typography.R"); source("gene_calib_stats.R"); source("server.R", local=new.env())' 2>/dev/null || true
Rscript - <<'RS'
suppressMessages({library(httr); library(jsonlite)})
# minimal harness: source just the function by extracting is impractical; instead
# call the running app path in dev. For CI, assert the search returns >0 IDs:
q <- URLencode('BRCA1[gene] AND ("benign"[clinsig] OR "likely benign"[clinsig]) AND "homo sapiens"[orgn]')
u <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=", q, "&retmax=50&retmode=json")
ids <- fromJSON(content(GET(u), "text", encoding="UTF-8"))$esearchresult$idlist
cat("BRCA1 benign IDs:", length(ids), "\n"); stopifnot(length(ids) > 0)
RS
```
Expected: `BRCA1 benign IDs: <N>` with N > 0. (Confirms the benign query returns results; full parse is exercised in-app.)

- [ ] **Step 5: Commit**

```bash
git add server.R
git commit -m "feat(calib): fetch benign ClinVar arm + gene_calib cache namespace"
```

---

### Task 5: Batch dbNSFP score fetch (server.R)

**Files:**
- Modify: `server.R` (add `fetch_dbnsfp_batch` near `fetch_dbnsfp` ~1923)
- Backup: `server.R.bak.<ts>.pre-4.2batch` before editing

**Interfaces:**
- Consumes: `httr`, `jsonlite`, `parse_clinvar_hgvsp` (Task 2, via sourced `gene_calib_stats.R`).
- Produces: `fetch_dbnsfp_batch(gene_name, hgvsps) -> data.frame(hgvsp, revel, am, cadd)` — one row per matched hgvsp; unmatched hgvsps absent (report count).

- [ ] **Step 1: Back up**

```bash
cp server.R "server.R.bak.$(date +%Y%m%d_%H%M%S).pre-4.2batch"
```

- [ ] **Step 2: Add the batch function**

```r
# Bulk REVEL/AM/CADD for a set of dbNSFP hgvsps via one MyVariant.info POST per <=1000.
# Returns data.frame(hgvsp, revel, am, cadd); unmatched hgvsps are omitted.
fetch_dbnsfp_batch <- function(gene_name, hgvsps) {
  hgvsps <- unique(hgvsps[!is.na(hgvsps)])
  if (length(hgvsps) == 0) return(data.frame(hgvsp=character(), revel=numeric(),
                                             am=numeric(), cadd=numeric(),
                                             stringsAsFactors = FALSE))
  chunks <- split(hgvsps, ceiling(seq_along(hgvsps) / 1000))
  rows <- list()
  for (ch in chunks) {
    body <- list(q = paste(ch, collapse = ","),
                 scopes = "dbnsfp.hgvsp",
                 fields = "dbnsfp.revel,dbnsfp.alphamissense,cadd.phred")
    resp <- tryCatch(httr::POST("https://myvariant.info/v1/query",
                                body = body, encode = "form", httr::timeout(30)),
                     error = function(e) NULL)
    if (is.null(resp) || httr::status_code(resp) != 200) next
    hits <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
                               simplifyVector = FALSE)
    for (h in hits) {
      if (isTRUE(h$notfound)) next
      rows[[length(rows) + 1]] <- data.frame(
        hgvsp = h$query %||% NA_character_,
        revel = as.numeric(h$dbnsfp$revel        %||% NA),
        am    = as.numeric(h$dbnsfp$alphamissense %||% NA),
        cadd  = as.numeric(h$cadd$phred           %||% NA),
        stringsAsFactors = FALSE)
    }
    Sys.sleep(0.2)  # rate-limit courtesy, matches existing dbNSFP path
  }
  if (length(rows) == 0) return(data.frame(hgvsp=character(), revel=numeric(),
                                           am=numeric(), cadd=numeric(),
                                           stringsAsFactors = FALSE))
  do.call(rbind, rows)
}
```

Note: `%||%` — confirm it exists in server.R; per session memory it does **not**. If absent, add near the top of `gene_calib_stats.R`:

```r
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a)==1 && is.na(a))) b else a
```

- [ ] **Step 3: Parse-check**

Run: `Rscript -e 'invisible(parse("server.R")); cat("parses OK\n")'`
Expected: `parses OK`.

- [ ] **Step 4: Live smoke test** (network)

Run:
```bash
Rscript - <<'RS'
suppressMessages({library(httr); library(jsonlite)})
`%||%` <- function(a,b) if (is.null(a)||length(a)==0||(length(a)==1&&is.na(a))) b else a
body <- list(q="p.R130G,p.R663H", scopes="dbnsfp.hgvsp",
             fields="dbnsfp.revel,dbnsfp.alphamissense,cadd.phred")
r <- POST("https://myvariant.info/v1/query", body=body, encode="form", timeout(30))
cat("status:", status_code(r), "\n")
hits <- fromJSON(content(r,"text",encoding="UTF-8"), simplifyVector=FALSE)
cat("hits:", length(hits), " first revel:", (hits[[1]]$dbnsfp$revel %||% NA), "\n")
stopifnot(status_code(r) == 200, length(hits) >= 1)
RS
```
Expected: `status: 200`, `hits: >=1` with a numeric REVEL. (Confirms batch POST shape + field paths.)

- [ ] **Step 5: Commit**

```bash
git add server.R gene_calib_stats.R
git commit -m "feat(calib): batch dbNSFP REVEL/AM/CADD fetch via MyVariant.info POST"
```

---

### Task 6: Live wrapper — fetch both arms, score, assemble, cache (server.R)

**Files:**
- Modify: `server.R` (add `gene_calibration_live` after Task 5's function)
- Backup: `server.R.bak.<ts>.pre-4.2live`

**Interfaces:**
- Consumes: `extract_clinvar` (824), `extract_clinvar_benign` (Task 4), `fetch_dbnsfp_batch` (Task 5), `parse_clinvar_hgvsp` + `gene_calibration` (stats core), `cache_get`/`cache_set`.
- Produces: `gene_calibration_live(gene_name, query_hgvsp=NULL, min_sens=0.90) -> list` — the `gene_calibration()` result plus `$match_rate` (fraction of ClinVar variants scored) and `$n_raw = list(path, benign)` (pre-match counts). Cached in `api_cache$gene_calib` keyed by `gene_name` (arms only; LR+/optimal recomputed per `min_sens`/`query_hgvsp`).

- [ ] **Step 1: Back up**

```bash
cp server.R "server.R.bak.$(date +%Y%m%d_%H%M%S).pre-4.2live"
```

- [ ] **Step 2: Add the wrapper**

```r
# Live gene-specific calibration: fetch both ClinVar arms, score via dbNSFP batch,
# join on hgvsp, then assemble. Caches the scored arms per gene; recomputes stats
# per (query_hgvsp, min_sens) which are cheap.
gene_calibration_live <- function(gene_name, query_hgvsp = NULL, min_sens = 0.90) {
  arms <- cache_get("gene_calib", gene_name)
  if (is.null(arms)) {
    pv <- extract_clinvar(gene_name)          # P/LP arm
    bn <- extract_clinvar_benign(gene_name)   # LB/B arm
    add_scores <- function(df) {
      if (is.null(df) || !nrow(df)) return(data.frame(hgvsp=character(), revel=numeric(),
                                                      am=numeric(), cadd=numeric()))
      hg <- parse_clinvar_hgvsp(df$name)
      sc <- fetch_dbnsfp_batch(gene_name, hg)
      merge(data.frame(hgvsp = hg, stringsAsFactors = FALSE), sc, by = "hgvsp")
    }
    arms <- list(path = add_scores(pv), benign = add_scores(bn),
                 n_raw = list(path = if (is.null(pv)) 0 else nrow(pv),
                              benign = if (is.null(bn)) 0 else nrow(bn)))
    cache_set("gene_calib", gene_name, arms)
  }
  res <- gene_calibration(arms$path, arms$benign,
                          query_hgvsp = query_hgvsp, min_sens = min_sens)
  scored <- nrow(arms$path) + nrow(arms$benign)
  raw    <- arms$n_raw$path + arms$n_raw$benign
  res$match_rate <- if (raw > 0) scored / raw else 0
  res$n_raw <- arms$n_raw
  res
}
```

- [ ] **Step 3: Parse-check**

Run: `Rscript -e 'invisible(parse("server.R")); cat("parses OK\n")'`
Expected: `parses OK`.

- [ ] **Step 4: Live smoke test** — launch the app on a well-populated gene

Run: `R -e 'shiny::runApp(".", port=7788, launch.browser=FALSE)'` in the background, query a gene with ≥20 P/LP and ≥20 LB/B (e.g. `BRCA1`, `MYH7`) via the app, and confirm the console logs a non-empty `gene_calibration_live` result: both arms scored, `match_rate` printed, `confident == TRUE` for REVEL.
Expected: REVEL validation table has 3 rows with finite LR+ and CIs; `match_rate` > 0.5.
(If launching headless is impractical here, defer this observation to Task 7's card smoke test — the function is pure-composed from already-tested parts.)

- [ ] **Step 5: Commit**

```bash
git add server.R
git commit -m "feat(calib): live wrapper — fetch/score/join both arms, cache per gene"
```

---

### Task 7: Variant-card section + params panel + gated override (server.R, ui.R)

**Files:**
- Modify: `server.R` (card render — add a "Gene-specific calibration" section near the title bar block ~6340–6354; read `input$calib_min_sens`; override toggle handling)
- Modify: `ui.R` (min-sensitivity numeric input in the parameters panel; per-card override checkbox is rendered in server.R HTML)
- Backup: both files to `.bak.<ts>.pre-4.2ui`

**Interfaces:**
- Consumes: `gene_calibration_live` (Task 6), `acmg_res` + `acmg_res$pts` in the card scope, `esc()`, `sec_hdr()`, `clinvar_goldstar` (server.R:863/972), `input$calib_min_sens`.
- Produces: rendered HTML section per variant card; session-only PP3 override applied to `acmg_res` **only** when the user opts in and gating passes.

- [ ] **Step 1: Back up**

```bash
cp server.R "server.R.bak.$(date +%Y%m%d_%H%M%S).pre-4.2ui"
cp ui.R     "ui.R.bak.$(date +%Y%m%d_%H%M%S).pre-4.2ui"
```

- [ ] **Step 2: Add min-sensitivity input to `ui.R` parameters panel**

Find the parameters panel (search `ui.R` for the prevalence/penetrance inputs) and add alongside them:

```r
numericInput("calib_min_sens",
             label = "Gene calibration — min sensitivity",
             value = 0.90, min = 0.50, max = 1.00, step = 0.05)
```

- [ ] **Step 3: Render the calibration section in the card (server.R)**

After the `title_bar <- paste0(...)` block (ends ~6354), build and append a section. Insert into the card's row assembly (where `row1`, etc. are concatenated) a new `calib_section`:

```r
        calib_section <- tryCatch({
          q_hgvsp <- parse_clinvar_hgvsp(as.character(r$Variant))  # if r$Variant carries p.change; else NA
          cal <- gene_calibration_live(gene_name = as.character(input$gene),
                                       query_hgvsp = q_hgvsp,
                                       min_sens = input$calib_min_sens %||% 0.90)
          rows <- lapply(names(cal$predictors), function(pn) {
            p <- cal$predictors[[pn]]
            conf <- if (p$confident) "" else
              paste0(' <span style="color:#b45309;">low confidence (n=', p$n_pos, '/', p$n_neg, ')</span>')
            vrows <- apply(p$validation, 1, function(v) {
              agree <- identical(v[["gene_tier"]], v[["pejaver_tier"]])
              paste0('<tr><td>', pn, ' &ge; ', v[["threshold"]], '</td>',
                     '<td>LR+ ', formatC(as.numeric(v[["lr"]]), format="f", digits=1),
                     ' (', formatC(as.numeric(v[["lo"]]),format="f",digits=1), '&ndash;',
                     formatC(as.numeric(v[["hi"]]),format="f",digits=1), ')</td>',
                     '<td>', v[["gene_tier"]], '</td>',
                     '<td>', if (agree) '=' else '&ne;', ' Pejaver ', v[["pejaver_tier"]], '</td></tr>')
            })
            opt <- p$optimal
            optrow <- if (is.null(opt)) '' else
              paste0('<tr><td>', pn, ' optimal &ge; ', formatC(opt$threshold,format="f",digits=3), '</td>',
                     '<td>LR+ ', formatC(opt$lr,format="f",digits=1), '</td><td>', opt$gene_tier, '</td><td>max LR+ @ sens floor</td></tr>')
            paste0(vrows, collapse = "", optrow, conf)
          })
          loo_note <- if (isTRUE(cal$loo_applied))
            '<div style="font-size:10px;color:#64748b;">Queried variant removed from calibration (leave-one-out).</div>' else ''
          clinvar_line <- paste0(
            '<div style="font-size:11px;">This variant in ClinVar: ',
            esc(as.character(r$ClinicalSignificance %||% "—")),
            if (!is.null(r$clinvar_goldstar)) paste0(' (', esc(r$clinvar_goldstar), '&#9733;)') else '',
            ' &mdash; independent of the calibration above (no double-count).</div>')
          paste0(sec_hdr("Gene-specific calibration", "#334155", "&#9878;"),
                 '<tr><td colspan="99"><table style="width:100%;font-size:11px;">',
                 paste0(rows, collapse=""), '</table>', loo_note, clinvar_line,
                 '<div style="font-size:10px;color:#94a3b8;">match rate ',
                 formatC(100*cal$match_rate, format="f", digits=0), '% of ClinVar variants scored.</div>',
                 '</td></tr>')
        }, error = function(e) "")
```

Then include `calib_section` in the final card `paste0(...)` that assembles the rows.

- [ ] **Step 4: Parse-check both files**

Run:
```bash
Rscript -e 'invisible(parse("server.R")); invisible(parse("ui.R")); cat("both parse OK\n")'
```
Expected: `both parse OK`.

- [ ] **Step 5: Live smoke test** — render a card

Launch the app, query a confident gene and a sparse gene:
- Confident gene (BRCA1/MYH7): the card shows the "Gene-specific calibration" section with 3 REVEL rows + optimal, AM/CADD rows, `=`/`≠` vs Pejaver, and the ClinVar line separate from calibration.
- Sparse gene: shows "low confidence (n=X/Y)" and still renders LR+/CI.
Expected: both render without error; low-confidence flag visible on the sparse gene.

- [ ] **Step 6: Commit**

```bash
git add server.R ui.R
git commit -m "feat(calib): variant-card gene-specific calibration section + min-sens param"
```

---

### Task 8: Gated opt-in override (server.R)

**Files:**
- Modify: `server.R` (per-card override checkbox in the calib section HTML + observer that applies the gene-specific tier to the PP3 tag for the session when gating passes)
- Backup: `server.R.bak.<ts>.pre-4.2override`

**Interfaces:**
- Consumes: the rendered `cal` result (Task 6/7), `acmg_res` PP3 tagging path, a `reactiveVal` `calib_override` (gene+predictor+tier), Shiny `observeEvent`.
- Produces: session-only reclassification: when the user enables override for a confident result, the chosen predictor's gene tier replaces the Pejaver PP3 contribution and the ACMG engine re-runs; low-confidence requires a confirm modal first.

- [ ] **Step 1: Back up**

```bash
cp server.R "server.R.bak.$(date +%Y%m%d_%H%M%S).pre-4.2override"
```

- [ ] **Step 2: Add override control + observer**

In `calib_section`, add a checkbox per confident predictor (or a single "apply gene-specific PP3" toggle keyed to the gene). Because the card is HTML (not native inputs), use a `actionButton`/`checkboxInput` rendered via `renderUI` in the parameters panel keyed to the current gene, OR a session `reactiveVal`:

```r
calib_override <- reactiveVal(NULL)   # list(gene, predictor, tier) or NULL

observeEvent(input$calib_apply, {
  cal <- gene_calibration_live(as.character(input$gene),
                               min_sens = input$calib_min_sens %||% 0.90)
  pred <- input$calib_pred                       # selectInput of REVEL/AM/CADD
  p <- cal$predictors[[pred]]
  apply_it <- function() calib_override(list(gene = input$gene, predictor = pred,
                                             tier = p$optimal$gene_tier))
  if (isTRUE(p$confident)) apply_it() else {
    showModal(modalDialog(
      title = "Low-confidence gene calibration",
      paste0("This gene has n=", p$n_pos, "/", p$n_neg,
             " (< 20 per arm). Apply the gene-specific PP3 tier anyway?"),
      footer = tagList(modalButton("Cancel"),
                       actionButton("calib_confirm", "Apply anyway"))))
  }
})
observeEvent(input$calib_confirm, {
  removeModal()
  cal <- gene_calibration_live(as.character(input$gene),
                               min_sens = input$calib_min_sens %||% 0.90)
  p <- cal$predictors[[input$calib_pred]]
  calib_override(list(gene = input$gene, predictor = input$calib_pred,
                      tier = p$optimal$gene_tier))
})
```

In the ACMG PP3 tagging path (the `add_pp3` block ~4882), after computing the Pejaver PP3 level, if `calib_override()` is set for the current gene, replace the PP3 tier with the override tier before the points accumulate:

```r
    ovr <- calib_override()
    if (!is.null(ovr) && identical(ovr$gene, input$gene)) {
      # map ovr$tier ("Supporting"/"Moderate"/"Strong") -> add_pp3 level, replacing Pejaver
      lvl <- c(Supporting="1", Moderate="2", Strong="3")[ovr$tier]
      if (!is.na(lvl)) { acmg_tags <- drop_pp3(acmg_tags); add_pp3(lvl) }
    }
```

Add the `calib_pred` selectInput and `calib_apply` actionButton to `ui.R` beside `calib_min_sens`. (`drop_pp3` — a one-line helper removing existing PP3 tags before re-adding; define near `add_pp3`.)

- [ ] **Step 3: Parse-check**

Run: `Rscript -e 'invisible(parse("server.R")); invisible(parse("ui.R")); cat("both parse OK\n")'`
Expected: `both parse OK`.

- [ ] **Step 4: Live smoke test**

- Confident gene: pick predictor, click Apply → PP3 tier changes to the gene-specific tier, ACMG score/badge update, no modal.
- Sparse gene: click Apply → confirm modal appears; "Apply anyway" applies, Cancel leaves the Pejaver call intact.
- Toggle off / re-query different gene → override clears (keyed to gene).
Expected: all three behave as described; default (no click) leaves the ACMG call identical to pre-feature.

- [ ] **Step 5: Commit**

```bash
git add server.R ui.R
git commit -m "feat(calib): gated opt-in PP3 override (session-only, confirm on low-N)"
```

---

## Self-Review

**Spec coverage:**
- §3 data flow → Tasks 4 (benign arm), 5 (batch dbNSFP), 6 (join + cache). ✓
- §4 statistic (LR+, Haldane, Katz CI, tier map) → Task 1. ✓
- §5 two outputs (validation + gene-optimal) → Tasks 1 (`gene_optimal`), 3 (assembly), 7 (render). ✓
- §6 circularity LOO + separate ClinVar assertion line → Task 3 (LOO), Task 7 (ClinVar line, no double-count). ✓
- §7 gating/confidence → Task 3 (`confident` flag), Task 7 (low-confidence render), Task 8 (confirm modal). ✓
- §8 authority/UI (advisory default, opt-in override, min-sens param) → Tasks 7 (section + param), 8 (override). ✓
- §9 files (server.R, ui.R, backups) → Tasks 4–8. ✓
- §11 testing (unit + live smoke) → unit in Tasks 1–3; live smoke in 4–8. ✓

**Placeholder scan:** No TBD/TODO. Two flagged conditionals resolved inline: `%||%` existence (added to `gene_calib_stats.R` if absent — Task 5 Step 2), and whether `r$Variant` carries the protein change for `query_hgvsp` (Task 7 uses `parse_clinvar_hgvsp` which returns NA safely → LOO simply doesn't fire, no crash).

**Type consistency:** `gene_calibration(path_arm, benign_arm, query_hgvsp, min_sens)` signature identical across Tasks 3, 6. Arm columns `hgvsp/revel/am/cadd` consistent across Tasks 3, 5, 6. `gene_lr_plus` return keys (`lr/lo/hi/sens/spec/tp/fp/n_pos/n_neg`) consistent Tasks 1, 3. Tier strings ("Supporting"/"Moderate"/"Strong"/"Very strong"/"None") consistent across `oddspath_tier`, validation, override map.

**Known risk carried forward:** Task 8 override wiring depends on the exact structure of the `add_pp3`/`acmg_tags` mechanism at ~4882; the implementer must read that block before Step 2 and adapt `drop_pp3`/`add_pp3` to the actual tag representation. Flagged in the task.
