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
//  * A key present in both the NN and IMM tables contributes from both
//    (independent lookups, as in R).
//  * Keys not found in any table are silently skipped.
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

struct Tbl {
  std::unordered_map<std::string, std::pair<double, double> > m;
  bool get(const std::string& key, double& dh, double& ds) const {
    std::unordered_map<std::string, std::pair<double, double> >::const_iterator
        it = m.find(key);
    if (it == m.end()) return false;
    dh = it->second.first;
    ds = it->second.second;
    return true;
  }
  bool has(const std::string& key) const { return m.find(key) != m.end(); }
};

Tbl make_tbl(const List& L) {
  Tbl t;
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
  const Tbl tmm_t = make_tbl(tmm);
  const Tbl imm_t = make_tbl(imm);
  const Tbl de_t  = make_tbl(de);
  const Tbl end_t = make_tbl(end);   // empty for all reference-salt sets

  // columns: dh, ds, nA, nC, nG, nT, len, ok
  NumericMatrix out(nseq, 8);
  colnames(out) = CharacterVector::create("dh", "ds", "nA", "nC", "nG", "nT",
                                          "len", "ok");

  for (R_xlen_t si = 0; si < nseq; ++si) {
    double dh = 0.0, ds = 0.0;
    bool ok = true;

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
      key_left = (lo < hi) ? keys[lo] : std::string();
      s.erase(0, 1);
      c.erase(0, 1);
    }
    if (de_t.get(key_right, th, ts)) {
      dh += th; ds += ts;
      if (hi > lo) --hi;
      s.erase(s.size() - 1, 1);
      c.erase(c.size() - 1, 1);
      key_right = right_key(s, c);
    }

    // -- Terminal mismatches ------------------------------------------------
    if (!key_left.empty() && tmm_t.get(key_left, th, ts)) {
      dh += th; ds += ts;
      ++lo;
      if (!s.empty()) s.erase(0, 1);
      if (!c.empty()) c.erase(0, 1);
    }
    if (!key_right.empty() && tmm_t.get(key_right, th, ts)) {
      dh += th; ds += ts;
      if (hi > lo) --hi;
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

    if (!s.empty() && s[0] == 'T') {
      if (nn_t.get("init_5T/A", th, ts)) { dh += th; ds += ts; }
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
    for (size_t i = lo; i < hi; ++i) {
      if (nn_t.get(keys[i], th, ts))  { dh += th; ds += ts; }
      if (imm_t.get(keys[i], th, ts)) { dh += th; ds += ts; }
    }

    // -- Symmetry correction (flag precomputed in R) ------------------------
    if (self_comp_eff) {
      if (nn_t.get("sym", th, ts)) { dh += th; ds += ts; }
      else ok = false;
    }

    out(si, 0) = dh;
    out(si, 1) = ds;
    out(si, 7) = ok ? 1.0 : 0.0;
  }

  return out;
}
