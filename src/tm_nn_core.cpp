// Fast nearest-neighbor dH/dS accumulation for tm_nn().
//
// This is a line-by-line port of the string-scanning part of the R function
// .tm_nn_core() (R/tm_nn.R). It handles shift/padding with '.', dangling
// ends, terminal and internal mismatches, and the initiation terms, and
// returns per-sequence delta_H, delta_S and A/C/G/T counts. Everything
// downstream (two-state Tm formula, salt correction, chemical correction)
// stays in R, vectorized, so the published correction functions remain the
// single source of truth for those formulas.
//
// Semantics preserved from the R implementation:
//  * A key present in both the NN and IMM tables is taken from the NN table
//    only; a stack has one delta_H and delta_S.
//  * A stack key is looked up in both orientations (see rev_key below),
//    because a stack and its character reversal are the same physical stack
//    and the published tables store only one of the two. A key neither table
//    holds marks the sequence (nostack = 1) rather than contributing zero.
//  * A missing initiation row ('init', 'init_5T/A', ...) marks the sequence
//    as failed (ok = 0), matching the R behavior where the subscript
//    error/NA was converted to NA by tryCatch.
//  * Base counts are taken on the ORIGINAL input sequence (before padding
//    and trimming), as .GC_fast()/gc() were called on seq_str.

#include <Rcpp.h>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace {

// A key "XY/WZ" denotes 5'-XY-3' paired with 3'-WZ-5'. Reading the same
// physical stack from the other strand reverses the whole key: the bottom
// strand 3'-WZ-5' read 5'->3' is "ZW", and the top strand read 3'->5' is "YX".
// So "XY/WZ" and its character reversal "ZW/YX" name the same stack, which is
// the rule Biopython's Tm_NN applies at lookup time (neighbors[::-1]).
inline std::string rev_key(const std::string& k) {
  return std::string(k.rbegin(), k.rend());
}

struct Tbl {
  std::unordered_map<std::string, std::pair<double, double> > m;
  // Whether a miss may be retried against the reversed key. True for every
  // table that stores each stack in one orientation only, which is all of
  // them except the Zuber end-effect table: that one already lists both
  // orientations of most of its keys, so a reversed lookup there could return
  // a value belonging to a different stack rather than the same one.
  bool rev_ok;
  Tbl() : rev_ok(true) {}
  bool get(const std::string& key, double& dh, double& ds) const {
    std::unordered_map<std::string, std::pair<double, double> >::const_iterator
        it = m.find(key);
    if (it == m.end()) {
      if (!rev_ok || key.size() < 2) return false;
      it = m.find(rev_key(key));
      if (it == m.end()) return false;
    }
    dh = it->second.first;
    ds = it->second.second;
    return true;
  }
};

Tbl make_tbl(const List& L, bool rev_ok = true) {
  Tbl t;
  t.rev_ok = rev_ok;
  // An absent or empty key vector yields an empty table (used for parameter
  // sets that define no end effects); guard against R NULL reaching here.
  if (!L.containsElementNamed("keys") ||
      Rf_isNull(L["keys"])) return t;
  CharacterVector keys = L["keys"];
  NumericVector dh = L["dh"];
  NumericVector ds = L["ds"];
  const R_xlen_t n = keys.size();
  t.m.reserve(static_cast<size_t>(n) * 2);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (CharacterVector::is_na(keys[i])) continue;
    t.m[std::string(keys[i])] = std::make_pair(dh[i], ds[i]);
  }
  return t;
}

// .right_key(): cseq[len], cseq[len-1], '/', seq[len], seq[len-1] (1-based)
inline std::string right_key(const std::string& s, const std::string& c) {
  const size_t n = s.size();
  std::string k;
  if (n < 2 || c.size() < 2) return k;  // empty: matches nothing
  k.reserve(5);
  k += c[n - 1];
  k += c[n - 2];
  k += '/';
  k += s[n - 1];
  k += s[n - 2];
  return k;
}

inline bool starts_dotdot(const std::string& x) {
  return x.size() >= 2 && x[0] == '.' && x[1] == '.';
}
inline bool ends_dotdot(const std::string& x) {
  const size_t n = x.size();
  return n >= 2 && x[n - 1] == '.' && x[n - 2] == '.';
}

}  // namespace

// [[Rcpp::export]]
NumericMatrix cpp_tm_nn_dhds(CharacterVector seqs, CharacterVector cseqs,
                             int shift, List nn, List tmm, List imm, List de,
                             bool self_comp_eff, List end) {
  const R_xlen_t nseq = seqs.size();
  if (cseqs.size() != nseq)
    stop("'seqs' and 'cseqs' must have the same length");

  const Tbl nn_t  = make_tbl(nn);
  // Orientation carries meaning for terminal mismatches: see the TMM block.
  const Tbl tmm_t = make_tbl(tmm, false);
  const Tbl imm_t = make_tbl(imm);
  const Tbl de_t  = make_tbl(de);
  // Empty for all reference-salt sets. No reversed retry: see Tbl::rev_ok.
  const Tbl end_t = make_tbl(end, false);

  // columns: dh, ds, nA, nC, nG, nT, len, ok, nostack
  NumericMatrix out(nseq, 9);
  colnames(out) = CharacterVector::create("dh", "ds", "nA", "nC", "nG", "nT",
                                          "len", "ok", "nostack");

  for (R_xlen_t si = 0; si < nseq; ++si) {
    double dh = 0.0, ds = 0.0;
    bool ok = true;
    // A stack that no table defines, as opposed to a missing initiation row
    // or a sequence too short to stack. Reported separately because it means
    // something different: the duplex is outside the model rather than
    // outside the parameter set.
    bool nostack = false;

    if (CharacterVector::is_na(seqs[si]) || CharacterVector::is_na(cseqs[si])) {
      out(si, 7) = 0.0;
      continue;
    }
    const std::string raw(seqs[si]);
    const std::string craw(cseqs[si]);

    // Clean: uppercase, keep only A/C/G/T/I. This ports the per-sequence
    // part of check_filter_seq(method = "tm_nn"); doing it here means the
    // regex pass over every window runs on the workers instead of serially
    // in the main process.
    std::string orig;
    orig.reserve(raw.size());
    for (size_t i = 0; i < raw.size(); ++i) {
      char ch = raw[i];
      if (ch >= 'a' && ch <= 'z') ch = static_cast<char>(ch - 32);
      switch (ch) {
        case 'A': case 'C': case 'G': case 'T': case 'I': orig += ch; break;
        default: break;
      }
    }
    std::string corig;
    corig.reserve(craw.size());
    for (size_t i = 0; i < craw.size(); ++i) {
      char ch = craw[i];
      if (ch >= 'a' && ch <= 'z') ch = static_cast<char>(ch - 32);
      switch (ch) {
        case 'A': case 'C': case 'G': case 'T': case 'I': corig += ch; break;
        default: break;
      }
    }
    std::string s(orig);
    std::string c(corig);

    // Base counts on the original sequence (for GC in R)
    long nA = 0, nC = 0, nG = 0, nT = 0;
    for (size_t i = 0; i < orig.size(); ++i) {
      switch (orig[i]) {
        case 'A': ++nA; break;
        case 'C': ++nC; break;
        case 'G': ++nG; break;
        case 'T': ++nT; break;
        default: break;
      }
    }
    out(si, 2) = static_cast<double>(nA);
    out(si, 3) = static_cast<double>(nC);
    out(si, 4) = static_cast<double>(nG);
    out(si, 5) = static_cast<double>(nT);
    out(si, 6) = static_cast<double>(orig.size());

    if (orig.size() < 2) {
      out(si, 7) = 0.0;
      continue;
    }

    // -- Shift / length padding with '.' (port of the R block) --------------
    if (shift != 0 || s.size() != c.size()) {
      if (shift > 0) {
        s.insert(0, static_cast<size_t>(shift), '.');
      } else if (shift < 0) {
        c.insert(0, static_cast<size_t>(-shift), '.');
      }
      if (c.size() > s.size()) s.append(c.size() - s.size(), '.');
      if (c.size() < s.size()) c.append(s.size() - c.size(), '.');
      while (starts_dotdot(s) || starts_dotdot(c)) {
        s.erase(0, 1);
        c.erase(0, 1);
      }
      while (ends_dotdot(s) || ends_dotdot(c)) {
        s.erase(s.size() - 1, 1);
        c.erase(c.size() - 1, 1);
      }
    }

    const size_t n0 = s.size();
    if (n0 < 2 || c.size() != n0) {
      out(si, 7) = 0.0;
      continue;
    }

    // -- Dinucleotide keys "s_i s_{i+1} / c_i c_{i+1}" ----------------------
    std::vector<std::string> keys;
    keys.reserve(n0 - 1);
    for (size_t i = 0; i + 1 < n0; ++i) {
      std::string k;
      k.reserve(5);
      k += s[i];
      k += s[i + 1];
      k += '/';
      k += c[i];
      k += c[i + 1];
      keys.push_back(k);
    }
    size_t lo = 0, hi = keys.size();  // active range [lo, hi)

    std::string key_left = keys[lo];
    std::string key_right = right_key(s, c);
    double th, ts;

    // -- Dangling ends ------------------------------------------------------
    if (de_t.get(key_left, th, ts)) {
      dh += th; ds += ts;
      ++lo;
      s.erase(0, 1);
      c.erase(0, 1);
    }
    if (de_t.get(key_right, th, ts)) {
      dh += th; ds += ts;
      if (hi > lo) --hi;
      s.erase(s.size() - 1, 1);
      c.erase(c.size() - 1, 1);
    }

    // -- Terminal mismatches ------------------------------------------------
    // The TMM tables are keyed with the PENULTIMATE pair first and the
    // terminal pair second: "AA/TA" is a Watson-Crick pair followed by a
    // mismatch. That is the orientation you get by reading along whichever
    // strand runs 5'->3' towards the duplex end -- the top strand at the
    // right-hand end, so keys[hi-1] is already in it, and the bottom strand
    // at the left-hand end, so the key there is the reversal of keys[lo].
    //
    // Unlike the stacking tables, the orientation carries meaning here and a
    // reversed retry would be wrong, not merely redundant: a duplex whose
    // terminal pair is Watson-Crick but whose penultimate pair is not has a
    // terminal-first key of exactly the shape the table stores, so probing
    // that spelling collects a terminal-mismatch penalty the molecule has not
    // earned. Hence rev_ok = false on tmm_t.
    if (lo < hi && tmm_t.get(rev_key(keys[lo]), th, ts)) {
      dh += th; ds += ts;
      ++lo;
      if (!s.empty()) s.erase(0, 1);
      if (!c.empty()) c.erase(0, 1);
    }
    if (lo < hi && tmm_t.get(keys[hi - 1], th, ts)) {
      dh += th; ds += ts;
      --hi;
      if (!s.empty()) s.erase(s.size() - 1, 1);
      if (!c.empty()) c.erase(c.size() - 1, 1);
    }

    // -- End effects (Zuber 2022 style) -------------------------------------
    // Terminal terms that depend on the penultimate base pair, i.e. on the
    // terminal dinucleotide stack. Unlike dangling ends and terminal
    // mismatches, these are ADDED without consuming the terminal pair, so
    // the stack itself is still counted below. Empty table for parameter
    // sets that use a single per-end penalty (stored in init_A/T instead).
    if (!end_t.m.empty()) {
      if (lo < hi && end_t.get(keys[lo], th, ts)) { dh += th; ds += ts; }
      // right end: .right_key() rewrites the terminal stack in the same
      // (terminal pair first) orientation used for the left end
      const std::string ekey_r = right_key(s, c);
      if (!ekey_r.empty() && end_t.get(ekey_r, th, ts)) { dh += th; ds += ts; }
    }

    // -- Initiation terms ---------------------------------------------------
    if (nn_t.get("init", th, ts)) { dh += th; ds += ts; } else ok = false;

    // The 5'-T penalty is due once for each strand whose 5' end is T. The top
    // strand's 5' end is s[0]; the bottom strand's is c[n-1], the last base of
    // the complement. Charging only the first made Tm depend on which strand
    // was handed over as the sequence. Every shipped parameter set carries
    // this row as zero, so no shipped result moves; a user table with a
    // non-zero value used to break the strand symmetry and now does not.
    //
    // Biopython's Tm_NN tests seq.endswith("A") for the second term, which is
    // the same thing only when the terminal pair is Watson-Crick. Reading the
    // complement directly is exact and symmetric by construction, since the
    // other strand's first base IS this one.
    int t5 = 0;
    if (!s.empty() && s[0] == 'T') ++t5;
    if (!c.empty() && c[c.size() - 1] == 'T') ++t5;
    if (t5 > 0) {
      if (nn_t.get("init_5T/A", th, ts)) { dh += th * t5; ds += ts * t5; }
      else ok = false;
    }

    const char first_base = s.empty() ? '\0' : s[0];
    const char last_base  = s.empty() ? '\0' : s[s.size() - 1];
    int gc_ends = 0;
    if (first_base == 'G' || first_base == 'C') ++gc_ends;
    if (last_base == 'G' || last_base == 'C') ++gc_ends;
    const int at_ends = 2 - gc_ends;

    if (gc_ends == 0) {
      if (nn_t.get("init_allA/T", th, ts)) { dh += th; ds += ts; }
      else ok = false;
    } else {
      if (nn_t.get("init_oneG/C", th, ts)) { dh += th; ds += ts; }
      else ok = false;
    }
    if (nn_t.get("init_A/T", th, ts)) { dh += th * at_ends; ds += ts * at_ends; }
    else ok = false;
    if (nn_t.get("init_G/C", th, ts)) { dh += th * gc_ends; ds += ts * gc_ends; }
    else ok = false;

    // -- Stacking / internal mismatch lookups -------------------------------
    // One stack, one parameter. The two tables are consulted in order rather
    // than both added: a stack has a single delta_H and delta_S, and where the
    // two tables overlap they describe it in different chemistries. That
    // happens for the G.U wobble stacks of the RNA sets, whose spellings also
    // occur in the DNA internal-mismatch table (Peyret 1999), where they mean
    // a DNA G.T mismatch. The nearest-neighbor set wins because it is the one
    // chosen for the molecule. No DNA set overlaps the mismatch table in
    // either orientation, so this changes nothing for DNA.
    // A stack neither table defines is not treated as contributing zero. Of
    // the 256 dinucleotide stacks over A/C/G/T, 116 have parameters and the
    // remaining 140 all carry two adjacent mismatches, which the two-state
    // nearest-neighbor model does not describe: the published sets measure a
    // mismatch flanked by Watson-Crick pairs. Only the three tandem G.T
    // stacks (GG/TT, GT/TG, TG/GT, from Allawi and SantaLucia 1997) have
    // measured values, which is the same coverage MELTING 5 reports for DNA.
    // Scoring the rest as zero overstates stability silently; the sequence is
    // marked instead and reported as NA with a warning, which is what MELTING
    // and Biopython do at whole-call granularity.
    for (size_t i = lo; i < hi; ++i) {
      if (nn_t.get(keys[i], th, ts))       { dh += th; ds += ts; }
      else if (imm_t.get(keys[i], th, ts)) { dh += th; ds += ts; }
      else { ok = false; nostack = true; }
    }

    // -- Symmetry correction (flag precomputed in R) ------------------------
    if (self_comp_eff) {
      if (nn_t.get("sym", th, ts)) { dh += th; ds += ts; }
      else ok = false;
    }

    out(si, 0) = dh;
    out(si, 1) = ds;
    out(si, 7) = ok ? 1.0 : 0.0;
    out(si, 8) = nostack ? 1.0 : 0.0;
  }

  return out;
}
