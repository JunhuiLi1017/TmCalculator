// Per-sequence base counts for gc(), tm_gc() and tm_wallace().
//
// The composition-based Tm methods used to split every sequence into a
// character vector with seqinr::s2c() and then scan it five times with %in%,
// once inside gc() and again inside salt_correct(). On the E. coli case study
// that made tm_gc() 77x slower than the compiled nearest-neighbor path.
//
// This counts every IUPAC code in one pass per sequence, reading the string
// bytes directly through CHAR(STRING_ELT(...)) so that no std::string is
// constructed and no R character object is allocated per base. The ambiguity
// apportioning stays in R (.gc_vec), operating on these counts, so that the
// published formula remains readable and is what the equivalence tests in
// tests/testthat/test_gc_vec.R check.
//
// Semantics preserved from gc():
//  * Input is uppercased.
//  * The accepted alphabet is A B C D G H I K M N R S T V W Y. Anything else
//    is tallied in the "other" column so that R can raise the same warning.
//  * N and I are accepted but contribute to neither the GC nor the AT side.
//    They are reported in their own columns so that the denominator actually
//    used by gc() -- A+C+G+T -- is inspectable rather than implicit: a window
//    of 199 N and one G is 100% GC by that definition, and the caller has to
//    be able to see why.
//  * NA input yields an all-NA row.

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix cpp_base_counts(CharacterVector seqs) {
  const R_xlen_t n = seqs.size();
  const int NC = 18;                 // len A C G T S W M K R Y V H D B N I other
  IntegerMatrix out(n, NC);

  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP s = STRING_ELT(seqs, i);
    if (s == NA_STRING) {
      for (int j = 0; j < NC; ++j) out(i, j) = NA_INTEGER;
      continue;
    }

    int cnt[NC] = {0};
    const char* p = CHAR(s);
    for (; *p != '\0'; ++p) {
      char c = *p;
      if (c >= 'a' && c <= 'z') c = (char)(c - 32);
      ++cnt[0];                      // len
      switch (c) {
        case 'A': ++cnt[1];  break;
        case 'C': ++cnt[2];  break;
        case 'G': ++cnt[3];  break;
        case 'T': ++cnt[4];  break;
        case 'S': ++cnt[5];  break;
        case 'W': ++cnt[6];  break;
        case 'M': ++cnt[7];  break;
        case 'K': ++cnt[8];  break;
        case 'R': ++cnt[9];  break;
        case 'Y': ++cnt[10]; break;
        case 'V': ++cnt[11]; break;
        case 'H': ++cnt[12]; break;
        case 'D': ++cnt[13]; break;
        case 'B': ++cnt[14]; break;
        case 'N': ++cnt[15]; break;  // accepted; contributes to neither side
        case 'I': ++cnt[16]; break;  // accepted (inosine); likewise
        default:  ++cnt[17]; break;  // outside gc()'s alphabet
      }
    }
    for (int j = 0; j < NC; ++j) out(i, j) = cnt[j];
  }

  colnames(out) = CharacterVector::create(
    "len", "A", "C", "G", "T", "S", "W", "M", "K", "R", "Y",
    "V", "H", "D", "B", "N", "I", "other");
  return out;
}
