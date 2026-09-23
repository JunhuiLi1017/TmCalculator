#!/usr/bin/env bash
# ===========================================================================
# setup_cluster_env.sh -- conda environment for the benchmark scripts
#
#   bash inst/scripts/setup_cluster_env.sh [env_name] [git_ref]
#   bash inst/scripts/setup_cluster_env.sh tmcalc dev
#
# Four steps:
#   1. conda environment with R 4.4.1
#   2. remotes, ps and BiocManager from CRAN
#   3. TmCalculator from GitHub (brings inst/scripts with it)
#   4. BSgenome.Hsapiens.UCSC.hg38 from Bioconductor
#
# Everything except r-base comes from R rather than from conda, so the library
# is a single ordinary R library and the versions are the ones BiocManager
# pins for this release. The Bioconductor dependencies are compiled from
# source, which is why the environment also carries a toolchain: conda's R is
# built against conda's libraries, and the module system's gcc links against
# the system ones, which produces packages that load and then segfault.
#
# WHAT THIS INSTALLS IS WHAT WAS PUSHED. Unpushed commits are not in the
# tarball GitHub serves, and a benchmark run against a stale build looks
# exactly like a valid one until the numbers are compared with the
# manuscript. The verification section checks for the specific fixes this
# work depends on rather than trusting the version string, which does not
# change when a commit is added.
# ===========================================================================
set -euo pipefail

ENV_NAME="${1:-tmcalc}"
GIT_REF="${2:-dev}"
REPO="JunhuiLi1017/TmCalculator"
CRAN="https://cloud.r-project.org"

command -v conda >/dev/null || { echo "conda not on PATH"; exit 1; }
# `conda activate` needs the shell hook in a non-interactive shell.
source "$(conda info --base)/etc/profile.d/conda.sh"

# --- 1. environment -------------------------------------------------------
echo "=== 1/4  creating environment '${ENV_NAME}' with R 4.4.1 ==="
conda create -y -n "${ENV_NAME}" -c conda-forge \
  r-base=4.4.1 compilers make pkg-config \
  zlib bzip2 xz libcurl openssl libxml2

conda activate "${ENV_NAME}"
R --vanilla -q -e 'cat("R", as.character(getRversion()), "at", R.home(), "\n")'

# A stray ~/.R/Makevars written for the module toolchain overrides the
# compiler conda just installed and is the usual cause of Bioconductor
# packages failing to link here.
if [ -f "${HOME}/.R/Makevars" ]; then
  echo
  echo "WARNING: ~/.R/Makevars exists and will override the conda toolchain."
  echo "         Move it aside before continuing if the builds below fail."
fi

# --- 2. CRAN packages -----------------------------------------------------
echo
echo "=== 2/4  installing remotes, ps, BiocManager ==="
R --vanilla -q -e "install.packages(c('remotes','ps','BiocManager'), \
                                    repos='${CRAN}')"

# --- 3. TmCalculator ------------------------------------------------------
echo
echo "=== 3/4  installing TmCalculator from ${REPO}@${GIT_REF} ==="
# The dependency solver needs both halves of the repository list. CRAN has to
# be set FIRST: BiocManager::repositories() adds the Bioconductor repositories
# to whatever getOption("repos") already holds, and under --vanilla that is
# the unresolved placeholder "@CRAN@". Left that way the Bioconductor half
# resolves and the CRAN half does not, so the install fails on Rcpp, BH, snow
# and the rest of the ordinary CRAN dependencies while never complaining about
# a Bioconductor package.
#
# R_REMOTES_USE_PAK=false keeps remotes on its own installer. pak reports the
# same failure as an opaque "error in pak subprocess", and it resolves the
# whole dependency graph up front, so one unreachable repository stops
# everything rather than just the packages that need it.
R_REMOTES_USE_PAK=false R --vanilla -q -e "
  options(repos = c(CRAN = '${CRAN}'))
  options(repos = BiocManager::repositories())
  stopifnot(!any(getOption('repos') == '@CRAN@'))
  print(getOption('repos'))
  remotes::install_github('${REPO}', ref = '${GIT_REF}',
                          upgrade = 'never', build_vignettes = FALSE)"

# --- 4. genome ------------------------------------------------------------
echo
echo "=== 4/4  installing BSgenome.Hsapiens.UCSC.hg38 (about 850 MB) ==="
R --vanilla -q -e "options(repos = c(CRAN = '${CRAN}')); \
   BiocManager::install('BSgenome.Hsapiens.UCSC.hg38', \
                        ask = FALSE, update = FALSE)"

# --- verification ---------------------------------------------------------
echo
echo "=== verification ==="
R --vanilla -q -e '
  ok <- TRUE
  say <- function(label, pass, detail = "") {
    cat(sprintf("  [%s] %-42s %s\n", if (pass) "ok" else "FAIL", label, detail))
    if (!pass) assign("ok", FALSE, envir = .GlobalEnv)
  }
  say("R version", getRversion() >= "4.4.0", as.character(getRversion()))

  suppressPackageStartupMessages(library(TmCalculator))
  say("TmCalculator version", TRUE,
      as.character(utils::packageVersion("TmCalculator")))

  ## The compiled core must be present. An install that silently fell back to
  ## the R path, or a load_all build compiled at -O0, makes every timing in
  ## the sweep wrong by roughly a factor of six.
  say("compiled nearest-neighbor core",
      exists("cpp_tm_nn_dhds", envir = asNamespace("TmCalculator")))

  ## Reverse-complement completion of the nearest-neighbor table. A key
  ## XY/WZ and its character reversal are the same duplex read from the
  ## opposite strand and must carry identical parameters. Four of these rows
  ## were transposed in earlier builds, so this distinguishes a current
  ## install from a stale one more reliably than the version number.
  tbl <- TmCalculator:::get_table("DNA_NN_SantaLucia_2004")
  rev_key <- function(k) paste(rev(strsplit(k, "", fixed = TRUE)[[1]]),
                               collapse = "")
  pairs <- c("TT/AA", "AC/TG", "AG/TC", "TC/AG", "TG/AC", "CC/GG")
  say("RC table completion",
      all(vapply(pairs, function(k)
        isTRUE(all.equal(unname(tbl[k, ]), unname(tbl[rev_key(k), ]))),
        logical(1))))

  say("gc_content() exported (gc() retired)",
      "gc_content" %in% getNamespaceExports("TmCalculator") &&
      !("gc" %in% getNamespaceExports("TmCalculator")))

  say("BPPARAM defaults to NULL",
      is.null(formals(TmCalculator::tm_calculate)$BPPARAM))

  say("hg38 genome package",
      requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE))
  say("ps (per-task memory)", requireNamespace("ps", quietly = TRUE))

  scripts <- c("bench_tm_calculate.R", "bench_tm_calculate.lsf",
               "bench_tm_calculate_local.sh")
  present <- vapply(scripts, function(f)
    nzchar(system.file("scripts", f, package = "TmCalculator")), logical(1))
  say("benchmark scripts installed", all(present),
      paste(scripts[!present], collapse = " "))

  cat(sprintf("\n  library:  %s\n  scripts:  %s\n",
              .libPaths()[1], system.file("scripts", package = "TmCalculator")))
  if (!ok) quit(status = 1)
'

cat <<EOF

Environment ready:  conda activate ${ENV_NAME}

Smoke test on an interactive allocation before submitting the real sweep,
which runs for hours:

  bsub -Is -n 4 -R "span[hosts=1]" -W 0:30 -q interactive bash
  conda activate ${ENV_NAME}
  Rscript "\$(R --vanilla -q -s -e 'cat(system.file("scripts/bench_tm_calculate.R", package="TmCalculator"))')" \\
    --workers 1,2 --reps 1 --outdir smoke_cluster

Two things to check in that output before going further:
  * df -T on the scratch path must report a LOCAL filesystem. nfs, gpfs or
    lustre means staging did nothing, and several workers reading the .2bit
    at once measure storage contention instead of the software.
  * n_windows must be 14,687,330 at every worker count. If the counts differ
    the segment boundaries are wrong and every timing below it is meaningless.
EOF
