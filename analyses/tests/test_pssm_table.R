# Self-check for the genome-wide PSSM table (analyses/derived/pfam_pssm_human.rds).
# Run: Rscript test_pssm_table.R   (after analyses/15p + 20)
#
# Guards the two defects the first build hit: quantise() flattening the column x
# 20 matrix to a vector (breaks the runtime M[column, mut] lookup), and a
# mega-family (PF00005) silently dropping. Skips cleanly if the table is absent.
RDS <- "analyses/derived/pfam_pssm_human.rds"
if (!file.exists(RDS)) { cat("skip: build the table first (", RDS, " missing)\n", sep=""); quit(status=0) }
source("analyses/lib/pfam_pssm.R")
o <- readRDS(RDS)
Q_MIN <- -12; Q_MAX <- 6
deq <- function(q) (as.numeric(q) + 127) / 254 * (Q_MAX - Q_MIN) + Q_MIN

# ── Every PSSM is a column x 20 matrix with AA names (not a flat vector) ────
stopifnot(length(o$pssm) > 6000, all(c("pssm","map","quant","aa") %in% names(o)))
samp <- o$pssm[unique(round(seq(1, length(o$pssm), length.out = 200)))]
for (m in samp) {
  stopifnot(!is.null(dim(m)), ncol(m) == 20L,
            identical(colnames(m), AA20), nrow(m) >= 1L)
}

# ── The dropped mega-family and the biggest clinical family are present ────
stopifnot("PF00005" %in% names(o$pssm),   # ABC transporter (OOM'd, recovered)
          "PF00069" %in% names(o$pssm),   # protein kinase (largest clinical fam)
          nrow(o$pssm[["PF00005"]]) == 2369L)

# ── Map schema + scale ─────────────────────────────────────────────────────
stopifnot(all(c("id","residue","column","family") %in% colnames(o$map)),
          nrow(o$map) > 1e7,
          all(o$map$family %in% names(o$pssm)))   # no orphan map rows

# ── The stored, quantised table reproduces a known DOLPHIN delta ───────────
# SNCA A53T, PF01387/SYUA_HUMAN residue 53: DOLPHIN reported +1.165. int8
# quantisation costs <=0.07/step, so anything within 0.2 confirms the pipeline.
p  <- o$pssm[["PF01387"]]
mr <- o$map[o$map$family == "PF01387" & grepl("SYUA_HUMAN", o$map$id) & o$map$residue == 53, ]
stopifnot(nrow(mr) >= 1)
Md <- deq(p); dim(Md) <- dim(p); dimnames(Md) <- dimnames(p)
d <- Md[mr$column[1], "T"] - Md[mr$column[1], "A"]
stopifnot(abs(d - 1.165) < 0.2)

cat(sprintf("pssm table: all checks pass (%s families, %s residues, A53T=%.3f vs 1.165)\n",
            format(length(o$pssm), big.mark=","), format(nrow(o$map), big.mark=","), d))
