library(testthat)

# testthat::test_file() changes cwd to tests/ during the test body, so resolve
# the source path relative to whichever cwd we're running under.
src_path <- if (file.exists("analyses/lib/clinvar_blind.R")) {
  "analyses/lib/clinvar_blind.R"
} else {
  "../analyses/lib/clinvar_blind.R"
}
source(src_path, local = TRUE)

test_that("strip_clinvar_tags removes the five direct circular criteria", {
  tags <- c("PVS1", "PS1", "PS1_moderate", "PM5", "PP5", "BP6",
            "PM2", "PP3", "BP4", "BA1")
  pm1_pathway <- character(0)
  out <- strip_clinvar_tags(tags, pm1_pathway = pm1_pathway)

  expect_false("PS1"          %in% out)
  expect_false("PS1_moderate" %in% out)
  expect_false("PM5"          %in% out)
  expect_false("PP5"          %in% out)
  expect_false("BP6"          %in% out)
  expect_true(all(c("PVS1", "PM2", "PP3", "BP4", "BA1") %in% out))
})

test_that("PM1 stripped only when pathway is clinvar_hotspot", {
  out_ccrs <- strip_clinvar_tags(c("PM1", "PM2"), pm1_pathway = "ccrs")
  expect_true("PM1" %in% out_ccrs)

  out_dom <- strip_clinvar_tags(c("PM1", "PM2"), pm1_pathway = "uniprot_domain")
  expect_true("PM1" %in% out_dom)

  out_hot <- strip_clinvar_tags(c("PM1", "PM2"), pm1_pathway = "clinvar_hotspot")
  expect_false("PM1" %in% out_hot)
  expect_true("PM2" %in% out_hot)
})

test_that("PP5/BP6 strength variants all stripped", {
  tags <- c("PP5", "PP5_moderate", "PP5_supporting",
            "BP6", "BP6_moderate", "BP6_supporting", "PP3")
  out <- strip_clinvar_tags(tags, pm1_pathway = character(0))
  expect_equal(out, "PP3")
})

test_that("empty input returns empty vector", {
  expect_equal(
    strip_clinvar_tags(character(0), pm1_pathway = character(0)),
    character(0)
  )
})
