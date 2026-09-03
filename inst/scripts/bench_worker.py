#!/usr/bin/env python3
"""bench_worker.py -- Biopython Bio.SeqUtils.MeltingTemp on one input size.

Invoked by benchmark_tools.R; the same interface as bench_worker.R so the
driver can time both the same way.

    python3 bench_worker.py <seqfile> <n> <outcsv>

Conditions match the R workers: SantaLucia 2004 nearest-neighbor parameters
(Biopython's DNA_NN4), 50 mM Na+, 25 nM of each strand, and the
Schildkraut/Lifson salt correction (Biopython saltcorr=1, i.e.
16.6*log10[Na+], the same formula as TmCalculator's "Schildkraut2010").
"""
import csv
import sys
import time

from Bio.SeqUtils import MeltingTemp as mt

seqfile, n, outcsv = sys.argv[1], int(sys.argv[2]), sys.argv[3]

with open(seqfile) as fh:
    seqs = [next(fh).strip() for _ in range(n)]

t0 = time.perf_counter()
tm = []
for s in seqs:
    try:
        tm.append(mt.Tm_NN(s, nn_table=mt.DNA_NN4,
                           Na=50, K=0, Tris=0, Mg=0, dNTPs=0,
                           dnac1=25, dnac2=25, selfcomp=False,
                           saltcorr=1))
    except Exception:
        tm.append(float("nan"))
elapsed = time.perf_counter() - t0

with open(outcsv, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["i", "Tm"])
    for i, v in enumerate(tm, start=1):
        w.writerow([i, v])

n_na = sum(1 for v in tm if v != v)
print(f"BENCH_ELAPSED {elapsed:.3f}")
print("BENCH_PREP 0.000")
print(f"BENCH_N {len(tm)}")
print(f"BENCH_NA {n_na}")
