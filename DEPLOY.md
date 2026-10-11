<!--
Reference copy. The operational copy lives beside the deploy staging directory
at ~/tools/varvizR/DEPLOY.md and is what a deploy is run against; that tree is
not a git repository, which is why this copy exists. Keep the two in step.

Sanitised for a public repository: home-directory paths are ~-relative and the
hosting plan's cost is omitted. Nothing operational is left out.
-->

# Deploy staging for varviz.shinyapps.io

This is a trimmed runtime subset, not a checkout. The development repository is
`~/tools/varviz/`; only the files the running app needs live here.

## Sync, then deploy

    cd ~/tools/varviz
    for f in server.R ui.R www/help.html typography.R gene_calib_stats.R \
             analyses/lib/pssm_lookup.R analyses/lib/dolphin.R \
             analyses/lib/clinvar_blind.R VERSION; do
      cp "$f" "~/tools/varvizR/$f"
    done

    Rscript -e 'rsconnect::deployApp("~/tools/varvizR",
                  appName="varviz", account="varviz", forceUpdate=TRUE)'

Confirm from rsconnect's own output — "Deployment complete" — and then from the
live app, not from a process check:

    curl -s https://varviz.shinyapps.io/varviz/ | grep -oE 'v[0-9]+\.[0-9]+\.[0-9]+'

## Every file server.R sources unguarded

`server.R` has five `source()` calls and none is wrapped in `tryCatch` or
guarded by `file.exists`. A missing one does not degrade the app; it crashes
every worker at startup, and the failure surfaces only after the upload.

    typography.R
    gene_calib_stats.R
    analyses/lib/dolphin.R
    analyses/lib/pssm_lookup.R
    analyses/lib/clinvar_blind.R

`clinvar_blind.R` was absent here through the 1.1.2 era because `server.R` did
not yet source it. The 1.2.2 sync added it. Re-check this list whenever
`server.R` changes:

    grep -n 'source(' server.R

## Keep backups out of this directory

`.rscignore` matches exact names only — it has no glob support in rsconnect
1.10.1, so the `*.bak.*` line it used to carry never excluded anything and the
1.2.2 deploy shipped five backup files. Backups now go to
`~/tools/varvizR_staging_backups/`, outside the bundle.

This directory is **not** a git repository, so nothing here is recoverable from
a checkout. Its current `.rscignore` is reproduced below; if the staging tree is
ever rebuilt, recreate it by hand or the excluded files start shipping again
with no error to say so.

    VarViz.CASR.pdf
    DEPLOY.md

`DEPLOY.md` excludes itself: it is deploy procedure, plan pricing and instance
sizing, which belong beside the bundle and not inside it. The bundle is 17 files
with this in place and 18 without, so the file count in the deploy log is the
quickest check that the exclusion is still working.

## Keep www/ small

Everything under `www/` is served at the URL root and ships in every bundle,
including files no page references. `www/Agasthya.jpg` was 1381x1548 and 1.4 MB
against 20 KB and 36 KB for the two photographs beside it (2026-10-11, b4fdf4d);
it is now 800 px on the long edge. None of the three headshots is referenced from
`ui.R`, `server.R` or `help.html`, but an unreferenced file under `www/` is still
a live URL somebody may have linked, so resize rather than delete.

## Plan and instance size

Basic plan, instance `xxlarge` = 4 GB. Ongoing, not a stopgap.
A single-variant render has peaked above 2 GB. The lazy `data/ccrs.sqlite`
store and the repacked PSSM table are what brought a worker from ~1 GB to
~600 MB at rest; do not revert either. Bundle is ~806 MB on disk, ~250 MB
compressed.
