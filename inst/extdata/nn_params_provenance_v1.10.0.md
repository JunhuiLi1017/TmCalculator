# Nearest-neighbor parameters added in v1.10.0: provenance

Provenance record for the parameter sets added in this release: the two
raised in review (Banerjee 2020, Zuber 2022) and the two molecular-crowding
sets from the same group (Ghosh 2020, 2023). Both were
transcribed from the published tables; the "package key" columns give the
mapping onto this package's row-name convention
(`5'XY3' / 3'WZ5'`, RNA U written as T).

Units: ΔH in kcal/mol, ΔS in cal/(mol·K) (eu). The package stores
`(ΔH, ΔS)` in the two columns `left`, `right` and derives ΔG internally,
so only ΔH and ΔS are transcribed here.

---

## 1. Banerjee et al. (2020) — RNA/DNA hybrid, 100 mM NaCl

Banerjee D, Tateishi-Karimata H, Ohyama T, Ghosh S, Endoh T, Takahashi S,
Sugimoto N. *Improved nearest-neighbor parameters for the stability of
RNA/DNA hybrids under a physiological condition.* Nucleic Acids Research
2020;48(21):12042–12054. doi:10.1093/nar/gkaa572 — Table 2.

**Status: implemented** as `RNA_DNA_NN_Banerjee_2020` (v1.10.0), with
`attr(, "salt_mM") = 100` so salt correction is skipped automatically at
the fitted condition.

Paper keys are written `rXY/dWZ` with the DNA strand 5'→3'; the package
writes the bottom strand 3'→5', so `rAC/dGT` → `AC/TG`.

| Paper key | Package key | ΔH | ΔS |
|---|---|---|---|
| rAA/dTT | AA/TT | −7.8 | −22.9 |
| rAC/dGT | AC/TG | −10.1 | −27.3 |
| rAG/dCT | AG/TC | −9.4 | −26.2 |
| rAU/dAT | AT/TA | −5.8 | −17.5 |
| rCA/dTG | CA/GT | −9.8 | −27.4 |
| rCC/dGG | CC/GG | −9.5 | −24.8 |
| rCG/dCG | CG/GC | −9.0 | −24.3 |
| rCU/dAG | CT/GA | −6.1 | −17.9 |
| rGA/dTC | GA/CT | −8.6 | −22.7 |
| rGC/dGC | GC/CG | −10.6 | −27.7 |
| rGG/dCC | GG/CC | −13.3 | −35.7 |
| rGU/dAC | GT/CA | −9.3 | −25.5 |
| rUA/dTA | TA/AT | −6.6 | −19.7 |
| rUC/dGA | TC/AG | −6.5 | −16.3 |
| rUG/dCA | TG/AC | −8.9 | −23.3 |
| rUU/dAA | TT/AA | −7.4 | −24.3 |
| init. rG−dC / rC−dG | init_G/C | 0 | −4.9 |
| init. rA−dT / rU−dA | init_A/T | 0 | −7.0 |

Model notes. Initiation is expressed **per duplex end** only: there is no
global initiation term, so `init`, `init_oneG/C`, `init_allA/T`,
`init_5T/A` and `sym` are all zero. The package core already applies
`init_A/T × (number of AT ends)` and `init_G/C × (number of GC ends)`,
which reproduces the published model exactly.

Consistency check. For the sixteen stacks, ΔH − 310.15·ΔS/1000 reproduces
the paper's tabulated ΔG°37 to the printed precision (e.g. rGG/dCC:
−13.3 + 11.07 = −2.23 ≈ −2.3). The two initiation rows do **not**:
the paper lists ΔG°37 = +2.0 and +2.6, whereas ΔH = 0 with
ΔS = −4.9 and −7.0 gives +1.52 and +2.17. Since the package computes Tm
from ΔH and ΔS, the tabulated ΔH/ΔS values are used as printed. This
discrepancy appears to be in the source table.

---

## 2. Zuber et al. (2022) — RNA/RNA, improved end effects

Zuber J, Schroeder SJ, Sun H, Turner DH, Mathews DH. *Nearest neighbor
rules for RNA helix folding thermodynamics: improved end effects.*
Nucleic Acids Research 2022;50(9):5251–5262. doi:10.1093/nar/gkac261 —
Tables 1A and 1B, "New Model" columns.

**Status: not yet implemented** (see the end-effect note below).

Despite the title, these are duplex parameters: they are fitted from
optical melting of short RNA duplexes using the two-state relation
1/Tm = R·ln(Ct/4)/ΔH° + ΔS°/ΔH°, and the comparison column of Table 1A is
Xia et al. (1998), i.e. this set is a direct successor to
`RNA_NN_Xia_1998`. The stacks are therefore usable for duplex Tm without
modification.

### 2A. Watson–Crick–Franklin stacks (Table 1A)

The paper's convention (top 5'→3' / bottom 3'→5') is the package's
convention, so mapping is U→T plus, for three stacks, the equivalent
representative already used by the package (a stack read from the
opposite strand is the same stack).

| Paper key | Package key | ΔH | ΔS |
|---|---|---|---|
| AA/UU | AA/TT | −7.44 | −20.98 |
| AU/UA | AT/TA | −8.91 | −25.22 |
| UA/AU | TA/AT | −9.16 | −25.40 |
| CA/GU | CA/GT | −10.47 | −27.08 |
| AC/UG | GT/CA | −11.98 | −31.37 |
| AG/UC | CT/GA | −9.34 | −23.66 |
| GA/CU | GA/CT | −13.75 | −36.53 |
| CG/GC | CG/GC | −9.61 | −23.46 |
| GC/CG | GC/CG | −16.52 | −42.13 |
| CC/GG | GG/CC | −13.94 | −34.41 |
| Initiation | init | +4.66 | +1.78 |
| Symmetry | sym | 0 | −1.38 |
| AU End on AU | (see note) | +4.36 | +13.35 |
| AU End on CG | (see note) | +3.17 | +8.79 |

Verified: ΔH − 310.15·ΔS/1000 reproduces every tabulated ΔG°37 (e.g.
GC/CG: −16.52 + 13.07 = −3.45 ≈ −3.46; Initiation: 4.66 − 0.55 = 4.11 ≈
4.10). The three re-mapped rows are confirmed against the paper's own
"1998 Model" column, which matches the values already stored in
`RNA_NN_Xia_1998` (AC/UG ↔ our `GT/CA` = −11.40; AG/UC ↔ our `CT/GA` =
−10.48; CC/GG ↔ our `GG/CC` = −13.39).

### The end-effect problem

Xia 1998 applies one terminal-AU penalty per AU end, which the package
stores in `init_A/T` (3.72, 10.5) and applies as
`init_A/T × (number of AT ends)`. Zuber 2022 splits this into two terms
that depend on the **penultimate** base pair: +4.36/+13.35 when the AU
end sits on an AU pair, +3.17/+8.79 when it sits on a CG pair. The
current core cannot express this, which is why the set is not yet
shipped: entering either value alone would silently misrepresent the
model.

The fix is small and well defined. The condition "terminal AU pair with
penultimate X pair" is exactly a property of the **terminal dinucleotide
stack**, and the core already computes those keys (`keys_fr[1]` and
`.right_key()`) for the dangling-end and terminal-mismatch lookups. What
is needed is an additional optional table looked up on those same keys
whose contribution is **added without trimming** the terminal base pair
(unlike `de_tbl`/`tmm_tbl`, which consume it). Scheduled in `ROADMAP.md`.

### 2B. GU-containing stacks (Table 1B)

| Paper key | Package key | ΔH | ΔS |
|---|---|---|---|
| GC/UG | GT/CG | −14.73 | −40.32 |
| CU/GG | CT/GG | −9.26 | −23.64 |
| GG/CU | GG/CT | −12.41 | −34.23 |
| CG/GU | CG/GT | −5.64 | −14.83 |
| AU/UG | AT/TG | −9.23 | −27.32 |
| GA/UU | TT/AG | −10.58 | −32.19 |
| UG/GU | TG/GT | −8.76 | −27.04 |
| UA/GU | TG/AT | −2.72 | −8.08 |
| GG/UU | GG/TT | −9.06 | −28.57 |
| GU/UG | GT/TG | −7.66 | −24.11 |
| AG/UU | AG/TT | −5.10 | −16.53 |
| AU End on GU | (end table) | +5.16 | +18.96 |
| GU End on CG | (end table) | +3.91 | +12.17 |
| GU End on AU | (end table) | +3.65 | +12.78 |
| GU End on GU | (end table) | +6.23 | +22.47 |

These eleven stacks are the same eleven GU rows carried by
`RNA_NN_Chen_2012`, under the strand-symmetry equivalences noted above
(e.g. paper `GC/UG` ≡ package `GT/CG`), so Table 1B slots into the
existing row set once the end table exists.

Out of scope: Table 1B also lists a `GGUC/CUGG` term (−32.49, −92.57),
a four-nucleotide motif correction rather than a nearest-neighbor stack.
It cannot be represented in a dinucleotide model and is not planned.

---

## 3. Ready-to-paste table for Zuber 2022 (WCF)

For use in `R/zzz.R` once the end table exists. Row order follows
`rownames(DNA_NN_Breslauer_1986)`, as the other RNA sets do. The
`init_A/T` value below is the AU-end-on-CG term; it is a placeholder and
must be replaced by the two-term end table before release.

```r
  RNA_NN_Zuber_2022 <- matrix(c(
      4.66,   1.78,     # init
      3.17,   8.79,     # init_A/T   PLACEHOLDER: AU end on CG only
      0.00,   0.00,     # init_G/C
      0.00,   0.00,     # init_oneG/C
      0.00,   0.00,     # init_allA/T
      0.00,   0.00,     # init_5T/A
      0.00,  -1.38,     # sym
     -7.44, -20.98,     # AA/TT
     -8.91, -25.22,     # AT/TA
     -9.16, -25.40,     # TA/AT
    -10.47, -27.08,     # CA/GT
    -11.98, -31.37,     # GT/CA
     -9.34, -23.66,     # CT/GA
    -13.75, -36.53,     # GA/CT
     -9.61, -23.46,     # CG/GC
    -16.52, -42.13,     # GC/CG
    -13.94, -34.41      # GG/CC
  ), ncol = 2, byrow = TRUE,
     dimnames = list(rownames(DNA_NN_Breslauer_1986), nn_col))
```


---

## 4. Ghosh et al. (2023) — RNA/RNA under molecular crowding

Ghosh S, Takahashi S, Banerjee D, Ohyama T, Endoh T, Tateishi-Karimata H,
Sugimoto N. *Nearest-neighbor parameters for the prediction of RNA duplex
stability in diverse in vitro and cellular-like crowding conditions.*
Nucleic Acids Research 2023;51(9):4101-4111. doi:10.1093/nar/gkad020 —
Table 1 (40 wt% PEG200, 100 mM NaCl).

**Status: implemented** as `RNA_NN_Ghosh_2023_PEG200`, `salt_mM = 100`.

Paper keys match the package convention directly (U written as T), in the
same order as `nn_row_std`.

| Package key | ΔH | ΔS |
|---|---|---|
| AA/TT | −10.0 | −30.4 |
| AT/TA | −10.1 | −30.8 |
| TA/AT | −11.1 | −31.5 |
| CA/GT | −12.1 | −32.1 |
| GT/CA | −10.7 | −28.7 |
| CT/GA | −11.2 | −30.4 |
| GA/CT | −11.7 | −30.5 |
| CG/GC | −11.1 | −28.8 |
| GC/CG | −13.8 | −34.6 |
| GG/CC | −14.8 | −38.4 |
| init | +4.6 | −2.9 |
| init_A/T (per terminal AU) | +6.5 | +18.2 |
| sym (self-complementary) | 0 | −1.4 |

All twelve tabulated ΔG°37 values are reproduced by ΔH − 310.15·ΔS/1000
to the printed precision.

---

## 5. Ghosh et al. (2020) — DNA/DNA under molecular crowding

Ghosh S, Takahashi S, Ohyama T, Endoh T, Tateishi-Karimata H, Sugimoto N.
*Nearest-neighbor parameters for predicting DNA duplex stability in diverse
molecular crowding conditions.* PNAS 2020;117(25):14194-14201.
doi:10.1073/pnas.1920886117 — Supporting Information Tables S7 (enthalpy)
and S8 (entropy), 40 wt% PEG200 with 100 mM NaCl.

**Status: implemented** as `DNA_NN_Ghosh_2020_PEG200`, `salt_mM = 100`.

**Important transcription note.** Tables S7 and S8 have two parameter
columns. The first, ΔH°/ΔS°\[cation\], is the no-cosolute reference at
100 mM NaCl (identical to Table S3). The second,
ΔH°/ΔS°\[40 wt% PEG 200\], is the **crowding increment, not an absolute
value** — several entries are positive, which is impossible for a stacking
term. The values below are the sum of the two columns.

| Package key | ΔH\[cation\] | ΔH increment | **ΔH used** | ΔS\[cation\] | ΔS increment | **ΔS used** |
|---|---|---|---|---|---|---|
| AA/TT | −7.9 | +1.4 | **−6.5** | −23.3 | +4.1 | **−19.2** |
| AT/TA | −7.2 | −2.2 | **−9.4** | −21.3 | −8.1 | **−29.4** |
| TA/AT | −7.2 | +2.9 | **−4.3** | −22.0 | +8.7 | **−13.3** |
| CA/GT | −8.5 | −4.6 | **−13.1** | −23.4 | −15.4 | **−38.8** |
| GT/CA | −8.4 | −0.8 | **−9.2** | −23.2 | −3.6 | **−26.8** |
| CT/GA | −7.8 | +4.4 | **−3.4** | −21.5 | +13.6 | **−7.9** |
| GA/CT | −8.2 | +3.3 | **−4.9** | −23.4 | +10.4 | **−13.0** |
| CG/GC | −10.6 | +4.2 | **−6.4** | −28.2 | +12.1 | **−16.1** |
| GC/CG | −9.8 | +5.6 | **−4.2** | −25.0 | +15.7 | **−9.3** |
| GG/CC | −8.0 | +4.0 | **−4.0** | −20.4 | +11.5 | **−8.9** |
| init_G/C (per GC end) | +0.1 | −10.2 | **−10.1** | −2.8 | −32.3 | **−35.1** |
| init_A/T (per AT end) | +2.3 | −5.2 | **−2.9** | +4.1 | −16.8 | **−12.7** |
| sym | — | — | **0** | — | — | **−1.4** |

Validation: re-predicting duplexes from Table S4 with the derived values
reproduces the paper's own predicted column, e.g.

| Duplex | ΔH derived | ΔH published | ΔS derived | ΔS published |
|---|---|---|---|---|
| d(GGCAGTTC) | −65.5 | −65.7 | −194.1 | −195.0 |
| d(CGCTGTAG) | −64.2 | −64.3 | −190.3 | −191.4 |
| d(GGACGTCC) (self-comp.) | −62.8 | −63.0 | −185.1 | −185.7 |
| d(CGTCGACG) (self-comp.) | −67.6 | −67.9 | −199.5 | −200.6 |
| d(ATGAGCTCAT) (self-comp.) | −71.6 | −71.6 | −214.3 | −215.1 |

Residual differences come from rounding the one-decimal tabulated inputs
across nine or more summed terms.

Initiation in this model is per duplex end (per G-C end and per A-T end),
matching the package's `init_G/C` / `init_A/T` handling; the global `init`
term and the `oneG/C`, `allA/T`, `5T/A` terms are therefore zero.
