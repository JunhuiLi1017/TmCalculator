# Draft replies for issues #8 and #9

Both fixed in 1.1.2, on GitHub now, CRAN to follow.

**Note on the affected range.** 1.0.5–1.1.1 is correct for **#8** only. **#9 is
not a regression** — the missing `R >= 6` branch is in every release in the git
history, back to 1.0.2 (2022). The two replies say different things about this
on purpose.

---

## Reply to #8 — strand asymmetry with an internal mismatch

Confirmed, and fixed in 1.1.2 (on GitHub now; a CRAN release will follow).
Thank you — the report was precise enough to reproduce immediately, and the
regression test you suggested is now in the package.

**Affected releases: 1.0.5 through 1.1.1** — that is 1.0.5, 1.0.6, 1.0.7,
1.0.8, 1.0.9, 1.1.0 and 1.1.1. Everything up to and including **1.0.4 is
correct**. Your four calls on 1.0.4:

```
f("GCATCGTAGGCTAGCT", "CGTAGCATCCGCTCGA")   # 59.9731
f("AGCTCGCCTACGATGC", "TCGATCGGATGCTACG")   # 59.9731
```

1.1.2 returns the same 59.9731 for both. That the pre-regression
implementation and the fixed one agree to four decimals, by different code
paths and with different shipped tables, is the strongest check we have that
the fix is right rather than merely symmetric.

### Your diagnosis of the symptom is exactly right; the mechanism is different

There is no swap-based fallback. There is **no fallback at all** — `Tbl::get()`
in `src/tm_nn_core.cpp` was a single `unordered_map::find`, and the R reference
path was a single `%in% rownames()`. A miss contributed zero, exactly as you
inferred from the arithmetic.

Why the retry disappeared is the part worth recording. Up to 1.0.4 the loop was

```r
if      (neighbors     %in% imm_table_name) ...
else if (rev_neighbors %in% imm_table_name) ...
else if (neighbors     %in% nn_table_name)  ...
else if (rev_neighbors %in% nn_table_name)  ...
else stop(...)
```

with `rev_neighbors` built by full reversal — the rule you identified. One
commit (2026-05-26) vectorised that loop and, in the same change, replaced the
retry with a completed parameter table: `.complete_nn_rc()`, which fills in the
six reverse-orientation rows the published nearest-neighbor tables omit.
**The completion was applied to the nearest-neighbor tables only.** The
internal-mismatch, terminal-mismatch and dangling-end tables were left in one
orientation with nothing left to look up the other.

That is why perfectly paired duplexes were untouched — they only use the NN
table, where the completion covered the gap — and why the fault survived seven
releases: every test in the package used perfectly paired duplexes. Your
invariant is the one thing that would have caught it, which is exactly your
point about the 1.1.0 transposition.

### A second fault underneath, which the invariant alone does not catch

Chasing this turned up an independent problem in the terminal-mismatch table.
Its keys carry the **penultimate** pair first and the **terminal** pair second
(`"AA/TA"` is a Watson-Crick pair then a mismatch), while the walk built the
terminal key with the terminal pair first. Two consequences:

- A genuine terminal mismatch never matched — 0 times across 400 random
  duplexes carrying one — and was scored with the internal-mismatch parameters
  instead. On one 16-mer, G·A, G·G and G·T at position 1 all returned the
  identical 62.1434, i.e. the terminal stack contributed nothing whatever the
  mismatch was.
- A duplex whose terminal pair **is** Watson-Crick but whose next pair is not
  produces a terminal-first key of exactly the shape the table stores, so it
  matched **spuriously** and collected a penalty it had not earned.

The second one is symmetric — both readings are wrong in the same way — so the
strand invariant passes on it. It needs its own test, which 1.1.2 has.

### What 1.1.2 does

- The reversed spelling is retried for the NN, internal-mismatch and
  dangling-end tables, as you suggested (`rev_key()` on the whole string).
- The terminal-mismatch key is **built** in the table's orientation instead —
  the reversal of the first stack at the left-hand end, the last stack as-is at
  the right — and that table is excluded from the retry, since for it the
  orientation carries meaning rather than being two spellings of one thing.
- A stack present in both `nn_table` and `imm_table` is taken from the NN table
  only. They used to be added together, which double-counted the G·U wobble
  stacks an RNA set shares with `DNA_IMM_Peyret_1999`. No DNA set overlaps that
  table in either orientation, so no DNA value moves.
- `init_5T/A` is charged once per strand whose 5' end is T, rather than for the
  top strand only. Zero in all 31 shipped sets, so nothing computed moves; a
  user-supplied table with a non-zero value used to break the invariant.
- Dangling ends and the Zuber 2022 end-effect table were checked and were
  already correct; `.right_key()` was already putting those keys in the right
  orientation.

`tests/testthat/test_regressions_1_1_2.R` contains your four calls verbatim,
your invariant over perfect / internal-mismatch / terminal-mismatch /
two-mismatch / dangling-end duplexes, and ten terminal-mismatch cases pinning
which of them may reach that table.

---

## Reply to #9 — Owczarzy2008 returns NA when sqrt([Mg2+])/[Mon] >= 6

Confirmed, and fixed in 1.1.2 (on GitHub now; a CRAN release will follow).
Thank you — the diagnosis was exact, including the Biopython reference.

**This one is not a regression.** Unlike #8, it is in every release in the git
history: v1.0.2 (2022) has the identical

```r
if (R < 0.22) { ... } else if (R >= 0.22 && R < 6.0) { ... }
```

with no `else`, so `corr` is never assigned. The symptom has changed over time
— older versions error out where 1.1.x returns NA — but the gap is the same
code in the same place. Anyone who has ever used `Owczarzy2008` with a
low-monovalent buffer is affected, not only users of 1.0.5–1.1.1.

The missing branch is now the published expression with the constants
unmodified (a = 3.92, b = −0.911, c = 6.26, d = 1.42, e = −48.2, f = 52.5,
g = 8.31), exactly as `Bio.SeqUtils.MeltingTemp.salt_correction()` method 7
does; only the competing regime reparameterises a, d and g on [Mon]. Your
example now gives:

```r
tm_nn(to_genomic_ranges("GCATCGTAGGCTAGCTTGCA"),
      salt_method = "Owczarzy2008", Na = 1, Mg = 5, dNTPs = 0)
#  before: Tm NA, GC NA
#  now:    Tm 63.162, GC 55
```

Two related things came out of fixing it.

**A worse sibling, in the same function.** `Owczarzy2008` with magnesium but
**no** monovalent cation at all — `Na = 0, K = 0, Tris = 0, Mg = 5` — fell into
the guard that returns a zero correction when [Mon] is zero. That guard is
right for the six methods that take log([Mon]) and wrong for this one: in the
divalent-dominated regime [Mon] drops out of the expression entirely, so the
correction is defined. This was worse than the NA you found, because nothing
marked the result as untrustworthy. It now returns the divalent form, and gives
the same answer as `Na = 1, Mg = 5`, which is what the model says it should.

**Your point about the silent NA is taken.** `tm_nn()` now warns whenever it
returns NA, naming the cause — a sequence the model cannot evaluate, or a salt
correction undefined at the requested conditions — and `GC` no longer follows
`Tm` into NA, since base composition is a property of the sequence rather than
of the thermodynamic model.

While in the area we also removed `Owczarzy2004` and `Owczarzy2008` from
`tm_gc()`. Those corrections apply to the reciprocal of the melting temperature
in kelvin, referenced to the same duplex in 1 M Na⁺, and carry a 1/(2(N−1))
duplex-length term of their own; `tm_gc()` was adding them to a Celsius value,
so the Tm came back essentially uncorrected. They remain available in
`tm_nn()`, where they belong.
