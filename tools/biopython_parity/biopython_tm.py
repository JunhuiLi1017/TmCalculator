#!/usr/bin/env python3
"""Compute Biopython melting temperatures for a case file.

Worker half of the TmCalculator/Biopython parity harness. Reads the case grid
written by compare_biopython.R, calls Bio.SeqUtils.MeltingTemp.Tm_NN once per
case, and writes one row per case with either a temperature or the error text.

    python3 biopython_tm.py cases.csv biopython.csv

Parameter sets are named with TmCalculator's names in the case file and
translated here, so the two sides can never silently pair up different tables.
Only the sets Biopython actually ships are mappable; a case naming anything
else is reported as an error rather than quietly skipped.

strict=True is deliberate: where a stack is missing from both tables Biopython
raises and TmCalculator returns NA, and the comparison should see that as the
agreement it is rather than compare against an imputed number.
"""

import argparse
import csv
import sys

try:
    from Bio.SeqUtils import MeltingTemp as mt
except ImportError:  # pragma: no cover
    sys.exit("Biopython is not installed. pip install biopython")

NN_TABLES = {
    "DNA_NN_Breslauer_1986": mt.DNA_NN1,
    "DNA_NN_Sugimoto_1996": mt.DNA_NN2,
    "DNA_NN_Allawi_1998": mt.DNA_NN3,
    "DNA_NN_SantaLucia_2004": mt.DNA_NN4,
    "RNA_NN_Freier_1986": mt.RNA_NN1,
    "RNA_NN_Xia_1998": mt.RNA_NN2,
    "RNA_NN_Chen_2012": mt.RNA_NN3,
    "RNA_DNA_NN_Sugimoto_1995": mt.R_DNA_NN1,
}
TMM_TABLES = {"DNA_TMM_Bommarito_2000": mt.DNA_TMM1}
IMM_TABLES = {"DNA_IMM_Peyret_1999": mt.DNA_IMM1}
DE_TABLES = {
    "DNA_DE_Bommarito_2000": mt.DNA_DE1,
    "RNA_DE_Turner_2010": mt.RNA_DE1,
}

# TmCalculator salt_method -> Biopython saltcorr method number
SALT = {
    "none": 0,
    "Schildkraut2010": 1,
    "Wetmur1991": 2,
    "SantaLucia1996": 3,
    "SantaLucia1998-1": 4,
    "SantaLucia1998-2": 5,
    "Owczarzy2004": 6,
    "Owczarzy2008": 7,
}


def as_bool(x):
    return str(x).strip().upper() in ("TRUE", "T", "1", "YES")


def run_case(row):
    nn = NN_TABLES.get(row["nn_table"])
    if nn is None:
        raise KeyError("no Biopython equivalent of " + row["nn_table"])
    c_seq = row["c_seq"] or None
    return mt.Tm_NN(
        row["seq"],
        check=False,
        strict=True,
        c_seq=c_seq,
        shift=int(row["shift"]),
        nn_table=nn,
        tmm_table=TMM_TABLES[row["tmm_table"]],
        imm_table=IMM_TABLES[row["imm_table"]],
        de_table=DE_TABLES[row["de_table"]],
        dnac1=float(row["dnac1"]),
        dnac2=float(row["dnac2"]),
        selfcomp=as_bool(row["selfcomp"]),
        Na=float(row["Na"]),
        K=float(row["K"]),
        Tris=float(row["Tris"]),
        Mg=float(row["Mg"]),
        dNTPs=float(row["dNTPs"]),
        saltcorr=SALT[row["salt_method"]],
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("cases")
    ap.add_argument("out")
    args = ap.parse_args()

    with open(args.cases, newline="") as fh:
        rows = list(csv.DictReader(fh))

    n_ok = n_err = 0
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "bio_tm", "bio_error"])
        for row in rows:
            try:
                w.writerow([row["id"], repr(float(run_case(row))), ""])
                n_ok += 1
            except Exception as exc:                      # noqa: BLE001
                w.writerow([row["id"], "", type(exc).__name__ + ": " + str(exc)])
                n_err += 1

    print(
        "biopython: {} cases, {} computed, {} raised".format(
            len(rows), n_ok, n_err
        ),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
