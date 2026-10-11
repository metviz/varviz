<!--
Canonical copy: the one in the repository root. The deploy staging tree carries
a byte-identical copy, synced by the loop below, because that tree is not a git
repository and would otherwise lose this file with no way to restore it. If the
two ever differ, the repository wins.

    cmp ./DEPLOY.md "$HOME/tools/varvizR/DEPLOY.md"

Written for a public repository: home-directory paths go through $HOME and the
hosting plan's cost is omitted. Nothing operational is left out.
-->

# Deploy staging for varviz.shinyapps.io

This is a trimmed runtime subset, not a checkout. The development repository is
`~/tools/varviz/`; only the files the running app needs live here.

## Sync, then deploy

    cd "$HOME/tools/varviz"
    for f in server.R ui.R www/help.html typography.R gene_calib_stats.R \
             analyses/lib/pssm_lookup.R analyses/lib/dolphin.R \
             analyses/lib/clinvar_blind.R VERSION DEPLOY.md; do
      cp "$f" "$HOME/tools/varvizR/$f"
    done

`$HOME` rather than `~`: bash does not expand a tilde inside double quotes, so
`"~/tools/varvizR/$f"` silently creates a directory called `~`. R expands `~` in
a path itself, which is why the `deployApp` call below can keep it.

`DEPLOY.md` is in that list so this file cannot drift from the repository copy.
It is `.rscignore`d, so syncing it does not put it in the bundle.

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

The staging directory is **not** a git repository, so nothing in it is
recoverable from a checkout. Its `.rscignore` is reproduced below; if that tree
is ever rebuilt, recreate the file by hand or the excluded entries start
shipping again with no error to say so. This document is covered by the sync
loop above, but `.rscignore` is not, because it is the thing that decides what
reaches the bundle and should not be rewritten by a routine copy.

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
