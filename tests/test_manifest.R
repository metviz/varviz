library(testthat)

manifest_path <- if (file.exists("analyses/lib/manifest.R")) "analyses/lib/manifest.R" else "../analyses/lib/manifest.R"
source(manifest_path, local = TRUE)

test_that("write_manifest produces a JSON file with required fields", {
  tmpdir <- tempfile("manifest_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))

  raw_bytes <- charToRaw("hello world")
  out_path <- write_manifest(
    out_dir       = tmpdir,
    label         = "unit_test",
    url           = "https://example.org/api",
    params        = list(query = "foo", limit = 10L),
    response_raw  = raw_bytes,
    http_status   = 200L
  )

  expect_true(file.exists(out_path))
  m <- jsonlite::fromJSON(out_path, simplifyVector = FALSE)

  expect_identical(m$url, "https://example.org/api")
  expect_identical(m$params$query, "foo")
  expect_equal(m$params$limit, 10L)
  expect_identical(m$http_status, 200L)
  expect_equal(m$response_bytes, length(raw_bytes))
  # SHA256 of "hello world" is well-known
  expect_identical(m$response_sha256,
                   "b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9")
  expect_match(m$timestamp_utc, "^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}Z$")
  expect_true(nchar(m$r_version) > 0)
  expect_true("jsonlite" %in% names(m$package_versions))
})

test_that("write_manifest creates the output directory if missing", {
  tmpdir <- file.path(tempfile("manifest_test_"), "nested", "dir")
  on.exit(unlink(dirname(dirname(tmpdir)), recursive = TRUE))

  out_path <- write_manifest(
    out_dir       = tmpdir,
    label         = "create_dir_test",
    url           = "https://example.org",
    params        = list(),
    response_raw  = raw(),
    http_status   = 200L
  )
  expect_true(file.exists(out_path))
})
