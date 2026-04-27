# Manifest helper — writes a JSON record per API fetch for reproducibility.
# Required fields: label, url, params, timestamp_utc, response_sha256, http_status,
# response_bytes, r_version, platform, package_versions.

suppressMessages({
  library(jsonlite)
  library(digest)
})

write_manifest <- function(out_dir, label, url, params, response_raw, http_status) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  fname <- sprintf("%s__%s.json",
                   format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"),
                   label)
  out_path <- file.path(out_dir, fname)

  pkgs <- c("jsonlite", "digest", "httr2", "dplyr", "ggplot2")
  pkg_versions <- vapply(pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
  }, character(1))

  manifest <- list(
    label            = label,
    url              = url,
    params           = params,
    timestamp_utc    = ts,
    http_status      = as.integer(http_status),
    response_bytes   = length(response_raw),
    response_sha256  = digest::digest(response_raw, algo = "sha256", serialize = FALSE),
    r_version        = R.version.string,
    platform         = R.version$platform,
    package_versions = as.list(pkg_versions[!is.na(pkg_versions)])
  )

  jsonlite::write_json(manifest, out_path,
                       auto_unbox = TRUE, pretty = TRUE, na = "null")
  out_path
}
