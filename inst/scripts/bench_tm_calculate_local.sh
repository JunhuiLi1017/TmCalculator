#!/bin/bash
# ===========================================================================
# Run the bench_tm_calculate.R sweep on a laptop, with memory sampling.
#
#   bash inst/scripts/bench_tm_calculate_local.sh
#
# This is the macOS and Linux counterpart of bench_tm_calculate.lsf. The R
# driver does not measure memory itself, for the reason given in its header:
# tm_calculate() has no hook inside its tasks and is not going to grow one
# for a benchmark. Memory is therefore sampled from outside, and the LSF
# script's sampler reads /proc, which macOS does not have. This one uses ps,
# which both systems do, and writes the same three columns the driver joins
# on, so the two environments produce comparable CSVs.
#
# WHAT IS DIFFERENT FROM THE CLUSTER RUN, and it belongs in any figure
# caption that puts the two side by side:
#
#   * pss_mb is a copy of rss_mb rather than a measurement. macOS has no
#     proportional set size. On the cluster the two agreed to 1.002x, which
#     is what says PSOCK workers share essentially nothing, so copying it
#     here is defensible; it is still not measured, and the column is kept
#     only so the file format matches.
#   * The machine is not quiet. A cluster job at least has a slot
#     reservation; a laptop has a browser. Close what you can, put the
#     machine on mains power, and disable App Nap and sleep for the duration
#     (caffeinate does that below on macOS). Thermal throttling on a laptop
#     is real and will show up as a wide range at the higher worker counts,
#     which is a result about laptops rather than noise to be hidden.
#   * There is no memory quota to compare against. The cluster run's point
#     was that 16 GB was enough; here 16 GB is simply all there is, and the
#     interesting number is how close the sweep comes to it.
#
# Everything the run writes goes to results_laptop/ under the current
# directory. Expect 60 to 90 minutes: six worker counts by three
# repetitions, the one-worker runs being about eight minutes each.
# ===========================================================================
set -euo pipefail

OUTDIR="$(pwd)/results_laptop"
mkdir -p "$OUTDIR"

# --- locate the R driver ---------------------------------------------------
# Unlike the LSF script this one keeps its $0, so the sibling file is the
# first place to look; the other two cover being run from the package root
# or against an installed copy.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_R=""
for p in "$HERE/bench_tm_calculate.R" \
         "inst/scripts/bench_tm_calculate.R" \
         "scripts/bench_tm_calculate.R"; do
  [ -f "$p" ] && { BENCH_R="$p"; break; }
done
if [ -z "$BENCH_R" ]; then
  BENCH_R="$(Rscript -e 'cat(system.file("scripts", "bench_tm_calculate.R",
                                         package = "TmCalculator"))' \
             2>/dev/null || true)"
fi
[ -f "$BENCH_R" ] || { echo "ERROR: cannot locate bench_tm_calculate.R" >&2; exit 1; }

# --- fail now, not in an hour ----------------------------------------------
echo "R         : $(command -v R)"
R --version | head -1
Rscript -e 'for (p in c("TmCalculator", "BSgenome.Hsapiens.UCSC.hg38",
                        "BiocParallel"))
              if (!requireNamespace(p, quietly = TRUE)) stop("missing: ", p)
            if (!"regions" %in% names(formals(TmCalculator::tm_calculate)))
              stop("this TmCalculator predates the merged tm_calculate()")
            cat("deps OK; TmCalculator ",
                as.character(packageVersion("TmCalculator")), "\n")'

# --- machine state, recorded rather than assumed ---------------------------
echo "=== machine ($(date -Iseconds 2>/dev/null || date)) ==="
uname -srm
if [ "$(uname -s)" = "Darwin" ]; then
  sysctl -n machdep.cpu.brand_string 2>/dev/null || true
  echo "physical cores : $(sysctl -n hw.physicalcpu)"
  echo "memory (GB)    : $(( $(sysctl -n hw.memsize) / 1024 / 1024 / 1024 ))"
  pmset -g ps | head -1          # on mains power, or on battery and throttled?
else
  lscpu | grep -E "^(Model name|Core|Socket|Thread)" || true
  free -g | head -2 || true
fi
echo "output dir     : $OUTDIR"
echo

# --- memory sampling -------------------------------------------------------
# ps reports RSS in kilobytes on both systems. Sampled every 10 s to match
# the cluster run, so both series are floors on the peak by the same margin.
MEMLOG="$OUTDIR/mem_sampling.tsv"
mem_sample_once() {
  local rss=0 max1=0 n=0 r
  for pid in $(pgrep -x R 2>/dev/null || true); do
    r=$(ps -o rss= -p "$pid" 2>/dev/null | tr -d ' ')
    [ -z "$r" ] && continue
    rss=$((rss + r)); n=$((n + 1))
    [ "$r" -gt "$max1" ] && max1=$r
  done
  # pss_mb repeats rss_mb: see the header.
  printf "%s\t%d\t%d\t%d\t%d\n" "$(date +%s)" "$n" \
         "$((rss / 1024))" "$((rss / 1024))" "$((max1 / 1024))"
}

printf "unix_time\tn_proc\trss_mb\tpss_mb\tmax1_mb\n" > "$MEMLOG"
( while :; do mem_sample_once >> "$MEMLOG"; sleep 10; done ) &
SAMPLER_PID=$!
trap 'kill "$SAMPLER_PID" 2>/dev/null || true' EXIT

# --- run -------------------------------------------------------------------
# Every flag matches bench_tm_calculate.lsf. If one of them drifts, the two
# environments stop being comparable and the figure that puts them on the
# same axes stops meaning anything, so change them in both files or neither.
RUN=(Rscript "$BENCH_R"
     --outdir  "$OUTDIR"
     --workers 6,5,4,3,2,1
     --reps    3
     --window  200
     --slide   200
     --unit    segment
     --segsize 50e6
     --genome  BSgenome.Hsapiens.UCSC.hg38
     --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY
     --nn-table DNA_NN_Breslauer_1986
     --na      50
     --memlog  "$MEMLOG")

# caffeinate keeps the machine awake for the duration. A laptop that sleeps
# mid-sweep does not fail, it produces one absurd wall time and no warning.
if command -v caffeinate >/dev/null 2>&1; then
  caffeinate -dimsu "${RUN[@]}"
else
  "${RUN[@]}"
fi

kill "$SAMPLER_PID" 2>/dev/null || true

# --- memory envelope over the whole sweep ----------------------------------
echo
echo "=== memory footprint over the whole sweep ==="
awk -F'\t' 'NR>1 {
      if ($3 > mr) { mr = $3; mrn = $2 }
      if ($5 > m1) { m1 = $5 }
      n++
    }
    END {
      printf "  samples                 : %d (10 s apart)\n", n
      printf "  peak sum of RSS         : %d MB\n", mr
      printf "  peak single R process   : %d MB\n", m1
      printf "  R processes at RSS peak : %d\n", mrn
    }' "$MEMLOG"
echo "  per-sample log          : $MEMLOG"
echo
echo "Per-run memory is joined onto bench_tm_calculate.csv by timestamp;"
echo "this block is the envelope over every run together."
echo
echo "done: $(date)"
