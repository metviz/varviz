#!/usr/bin/env Rscript
# 26_build_varviz_sqlite.R — build the cloud-runtime data artifacts.
# Moves the two per-gene "elephant" tables out of RAM. Two deployment profiles:
#
#   FULL STORE (paid tier, >1 GB bundle) -> data/varviz_store.sqlite
#     * pssm     : one row per Pfam family, raw-packed int8 score matrix as BLOB
#     * pssm_map : (id, residue, family, column), indexed on (id, residue)
#     * ccrs     : the CCRS table, indexed on ensembl_gene_name and uniprotAcc
#     Both tables lazy -> ~47 MB RAM, but the 1.5 GB bundle needs a paid plan.
#
#   FREE-TIER DEPLOY (1 GB bundle cap) -> two files copied into varvizR/data/:
#     * data/ccrs.sqlite              — ccrs lazy (SQLite), ~746 MB
#     * data/pfam_pssm_human_deploy.rds — PSSM stays in-memory but PRE-BAKED to
#       raw+factor + gzip so it loads in ~1 s (vs ~15 s for the integer rds) at
#       cold start; copy it to varvizR/data/pfam_pssm_human.rds.
#     Per-worker RAM ~600 MB, bundle ~800 MB. This is what's live on shinyapps.
#
# gene_data stays in RData (9 MB). The canonical data/pfam_pssm_human.rds keeps
# its integer form (this script + tests read it); only the deploy copy is baked.
suppressMessages({library(DBI); library(RSQLite)})

RDS  <- "data/pfam_pssm_human.rds"
RDAT <- "data/VarViz.RData"
OUT  <- "data/varviz_store.sqlite"
if (file.exists(OUT)) unlink(OUT)

t <- readRDS(RDS)                       # integer matrices, char map
e <- new.env(); load(RDAT, envir = e)   # ccrs_data, gene_data
con <- dbConnect(SQLite(), OUT)
dbExecute(con, "PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;")

# ---- meta (quant bounds, aa order) ----
dbExecute(con, "CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT)")
meta <- data.frame(key=c("quant_min","quant_max","aa","built","source"),
                   value=c(t$quant[["min"]], t$quant[["max"]], paste(t$aa,collapse=""),
                           as.character(Sys.Date()), "26_build_varviz_sqlite.R"),
                   stringsAsFactors=FALSE)
dbWriteTable(con, "meta", meta, append=TRUE)

# ---- pssm: raw-packed matrices (q+127) as BLOB, column-major, nr rows x 20 aa ----
dbExecute(con, "CREATE TABLE pssm(family TEXT PRIMARY KEY, nr INTEGER, nc INTEGER, dat BLOB)")
ins <- dbSendStatement(con, "INSERT INTO pssm VALUES(:family,:nr,:nc,:dat)")
fams <- names(t$pssm)
for (i in seq_along(fams)) {
  m <- t$pssm[[i]]
  raw <- as.raw(as.integer(m) + 127L)   # [-127,127] -> [0,254]
  dbBind(ins, list(family=fams[i], nr=nrow(m), nc=ncol(m), dat=list(raw)))
}
dbClearResult(ins)

# ---- pssm_map: (id, residue, family, column) indexed on (id,residue) ----
dbWriteTable(con, "pssm_map",
             data.frame(id=t$map$id, residue=as.integer(t$map$residue),
                        family=t$map$family, column=as.integer(t$map$column),
                        stringsAsFactors=FALSE))
dbExecute(con, "CREATE INDEX idx_map ON pssm_map(id, residue)")

# ---- ccrs: whole table, indexed on the two query keys ----
dbWriteTable(con, "ccrs", as.data.frame(e$ccrs_data))
dbExecute(con, "CREATE INDEX idx_ccrs_gene ON ccrs(ensembl_gene_name)")
dbExecute(con, "CREATE INDEX idx_ccrs_uacc ON ccrs(uniprotAcc)")

dbExecute(con, "VACUUM")
cat(sprintf("pssm families: %d | pssm_map rows: %d | ccrs rows: %d\n",
            length(fams), nrow(t$map), nrow(e$ccrs_data)))
dbDisconnect(con)
cat(sprintf("wrote %s (%.0f MB)\n", OUT, file.info(OUT)$size/1e6))

# ── Free-tier deploy artifacts (what actually ships to shinyapps) ─────────────
source("analyses/lib/pssm_lookup.R")

# (a) ccrs-only SQLite — ccrs stays lazy, pssm does not.
CCRS_OUT <- "data/ccrs.sqlite"; if (file.exists(CCRS_OUT)) unlink(CCRS_OUT)
cc <- dbConnect(SQLite(), CCRS_OUT); dbExecute(cc, "PRAGMA journal_mode=OFF")
dbWriteTable(cc, "ccrs", as.data.frame(e$ccrs_data))
dbExecute(cc, "CREATE INDEX idx_ccrs_gene ON ccrs(ensembl_gene_name)")
dbExecute(cc, "CREATE INDEX idx_ccrs_uacc ON ccrs(uniprotAcc)")
dbExecute(cc, "VACUUM"); dbDisconnect(cc)
cat(sprintf("wrote %s (%.0f MB)\n", CCRS_OUT, file.info(CCRS_OUT)$size/1e6))

# (b) pre-baked PSSM rds — pssm_table_load() re-packs the int8 matrices to raw
# and factorizes the map; saving that keeps the conversion OUT of the cold-start
# path (readRDS of a raw+factor gzip ~1 s vs ~15 s for the integer+convert).
# gzip (not xz): larger file but ~1 s load vs ~4 s, and cold-start latency is the
# thing we are optimising. Copy this to varvizR/data/pfam_pssm_human.rds.
PSSM_OUT <- "data/pfam_pssm_human_deploy.rds"
saveRDS(pssm_table_load(RDS), PSSM_OUT, compress = "gzip")
cat(sprintf("wrote %s (%.0f MB) — copy to varvizR/data/pfam_pssm_human.rds\n",
            PSSM_OUT, file.info(PSSM_OUT)$size / 1e6))
