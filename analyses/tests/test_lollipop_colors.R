# Guards the lollipop colour decision in pfamplot().
#
# Why this exists. Lollipop colour reads like a category key and is not one.
# pfamplot() colours the stem by the REFERENCE amino acid's chemical class and the
# dot by the VARIANT amino acid's class, and sets legend.position = "none", so no
# key is drawn and nothing in that code path ever receives a mechanism label. A
# variant such as p.A364E renders red because glutamate is acidic, whatever
# functional group it belongs to.
#
# The fix adds opt-in category colouring. A caller may attach mech_color (and
# optionally mech_label) to the highlight frame; when absent, behaviour is exactly
# the shipped amino-acid-class palette. The decision is extracted into
# lollipop_colors() so it can be asserted against the shipped function rather than
# against a copy of the rule in this file, which is the failure mode of commit
# 4c93a92 (a test re-implemented the Pathogenic cap and passed while the shipped
# rule was dead code).
#
# The first block below is the regression guard that matters: the default path
# must keep emitting the amino-acid-class colours, because every existing figure
# and every live app session depends on it.
#
# Run from project root: Rscript analyses/tests/test_lollipop_colors.R

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
src_path <- resolve("server.R")

# Pull just lollipop_colors() out of server.R rather than sourcing the whole app,
# which would load shiny and the full data stack for one pure function. This is
# the same readLines/parse/eval idiom as analyses/tests/test_acmg_posterior.R.
# The eval input is a fixed line range of this repository's own tracked source,
# not user or network data, and the point of reading it from server.R is to assert
# against the shipped function rather than a copy of it in this file.
src <- readLines(src_path, warn = FALSE)
beg <- grep("^lollipop_colors <- function", src)
if (length(beg) != 1L) {
  cat("FAIL: lollipop_colors() not found in server.R (expected exactly one definition, found ",
      length(beg), ")\n", sep = "")
  quit(status = 1L)
}
end <- beg - 1L + which(src[beg:length(src)] == "}")[1]
eval(parse(text = paste(src[beg:end], collapse = "\n")))

fails <- 0L
check <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

# Chemical classes as the shipped palette defines them: acidic #e04040,
# basic #4472c4, polar #2ca02c, aromatic #e67e22, aliphatic #8e44ad, proline #8B4513.
# One-letter notation, because the shipped rule reads the first and LAST character
# of the name. That is the form the app's own exports use; see the three-letter
# block below for what the rule does with "p.Ala364Glu".
h_plain <- data.frame(
  Mutation = c("p.A364E", "p.L87P", "p.S247F", "p.V104I"),
  prot_pos = c(364L, 87L, 247L, 104L),
  stringsAsFactors = FALSE
)

cat("-- default path: amino-acid chemical class, no mechanism, no key --\n")
d <- lollipop_colors(h_plain)

check(isFALSE(d$mech), "mech is FALSE when the highlight frame carries no mech_color")
check(is.null(d$pal), "no palette is returned on the default path")
# A364E: ref Ala aliphatic -> #8e44ad ; var Glu acidic -> #e04040
check(identical(d$stem[1], "#8e44ad"), "p.A364E stem takes the reference Ala aliphatic colour")
check(identical(d$dot[1],  "#e04040"), "p.A364E dot takes the variant Glu acidic colour")
# L87P: ref Leu aliphatic -> #8e44ad ; var Pro -> #8B4513
check(identical(d$stem[2], "#8e44ad"), "p.L87P stem takes the reference Leu aliphatic colour")
check(identical(d$dot[2],  "#8B4513"), "p.L87P dot takes the variant Pro colour")
# S247F: ref Ser polar -> #2ca02c ; var Phe aromatic -> #e67e22
check(identical(d$stem[3], "#2ca02c"), "p.S247F stem takes the reference Ser polar colour")
check(identical(d$dot[3],  "#e67e22"), "p.S247F dot takes the variant Phe aromatic colour")
# The colour carries no group information: an externally-assigned category does
# not determine the dot colour on the default path.
check(!identical(d$dot[1], "#e67e22"),
      "dot colour is set by residue chemistry, not by any assigned category")

cat("-- known limitation: three-letter notation is not parsed --\n")
# Recorded rather than fixed. The rule takes the last character of the name, so
# "p.Ala364Glu" yields "U" and a grey dot, while the stem is right only because
# "A" is both the first letter of "Ala" and alanine's one-letter code. Callers
# passing three-letter notation silently get grey dots, so the behaviour is pinned
# here. The opt-in path does not parse residues and is unaffected.
d3 <- lollipop_colors(data.frame(Mutation = "p.Ala364Glu", prot_pos = 364L,
                                 stringsAsFactors = FALSE))
check(identical(d3$dot[1], "#555555"),
      "three-letter notation gives a grey dot (parses 'u' of 'Glu', not 'E')")

cat("-- unknown residue falls back to grey, as shipped --\n")
d_unk <- lollipop_colors(data.frame(Mutation = "p.Xaa1Zzz", prot_pos = 1L,
                                    stringsAsFactors = FALSE))
check(identical(d_unk$stem[1], "#555555"), "unrecognised reference residue gives #555555")
check(identical(d_unk$dot[1],  "#555555"), "unrecognised variant residue gives #555555")

cat("-- opt-in path: mechanism colours and a key --\n")
h_mech <- h_plain
h_mech$mech_color <- c("#d62728", "#d62728", "#ff7f0e", "#7f7f7f")
h_mech$mech_label <- c("LOF", "LOF", "GOF", "Uncertain")
m <- lollipop_colors(h_mech)

check(isTRUE(m$mech), "mech is TRUE when mech_color is supplied")
check(identical(m$stem, h_mech$mech_color), "stems follow the supplied mechanism colours")
check(identical(m$dot,  h_mech$mech_color), "dots follow the supplied mechanism colours")
check(identical(m$label, h_mech$mech_label), "labels are carried through for the key")
check(is.character(m$pal) && !is.null(names(m$pal)),
      "a named label-to-colour palette is returned so a key can be drawn")
check(identical(sort(names(m$pal)), c("GOF", "LOF", "Uncertain")),
      "the palette has one entry per distinct mechanism label")
check(identical(unname(m$pal[["GOF"]]), "#ff7f0e"), "the GOF label maps to its supplied colour")

cat("-- mech_color present but all NA is treated as absent --\n")
h_na <- h_plain
h_na$mech_color <- NA_character_
check(isFALSE(lollipop_colors(h_na)$mech),
      "an all-NA mech_color column does not switch on mechanism colouring")

cat("-- mech_label may be omitted; colours then label themselves --\n")
h_nolab <- h_plain
h_nolab$mech_color <- c("#d62728", "#d62728", "#ff7f0e", "#7f7f7f")
n <- lollipop_colors(h_nolab)
check(isTRUE(n$mech), "mech is TRUE with mech_color and no mech_label")
check(identical(sort(unique(n$label)), sort(unique(h_nolab$mech_color))),
      "labels fall back to the colour values when mech_label is absent")

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
