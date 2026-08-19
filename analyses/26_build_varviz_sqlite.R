#!/usr/bin/env Rscript
# 26_build_varviz_sqlite.R — build the lazy-load store for the cloud runtime.
# Moves the two per-gene "elephant" tables out of RAM and onto disk, keyed for
# random access so a request loads only the queried gene's slice:
#   * pssm     : one row per Pfam family, raw-packed int8 score matrix as BLOB
#   * pssm_map : (id, residue, family, column), indexed on (id, residue)
#   * ccrs     : the CCRS table, indexed on ensembl_gene_name and uniprotAcc
# gene_data stays in RData (9 MB, negligible). Output: data/varviz_store.sqlite
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
