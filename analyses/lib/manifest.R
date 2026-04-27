# Manifest helper — writes a JSON record per API fetch for reproducibility.
# Required fields: label, url, params, timestamp_utc, response_sha256, http_status,
# response_bytes, r_version, platform, package_versions.

suppressMessages({
  library(jsonlite)
  library(digest)
})

write_manifest <- function(out_dir, label, url, params, response_raw, http_status) {
  if (!dir.exists(out_dir)) {
    ok <- dir.create(out_dir, recursive = TRUE)
    if (!ok && !dir.exists(out_dir)) stop("write_manifest: cannot create out_dir: ", out_dir)
  }

  now <- Sys.time()
  ts    <- format(now, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  safe_label <- gsub("[^A-Za-z0-9._-]", "_", label)
  fname <- sprintf("%s__%s.json",
                   format(now, "%Y%m%dT%H%M%SZ", tz = "UTC"),
                   safe_label)
  out_path <- file.path(out_dir, fname)

  loaded_pkgs <- setdiff(loadedNamespaces(),
                         c("base", "stats", "graphics", "grDevices",
                           "utils", "datasets", "methods"))
  pkg_versions <- setNames(
    vapply(loaded_pkgs, function(p) as.character(utils::packageVersion(p)), character(1)),
    loaded_pkgs
  )

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
    package_versions = as.list(pkg_versions)
  )

  jsonlite::write_json(manifest, out_path,
                       auto_unbox = TRUE, pretty = TRUE, na = "null")
  out_path
}
