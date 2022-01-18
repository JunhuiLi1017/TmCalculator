% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Tm_GC.R
\name{Tm_GC}
\alias{Tm_GC}
\title{Calculate the melting temperature using empirical formulas based on GC content}
\usage{
Tm_GC(
  ntseq,
  ambiguous = FALSE,
  userset = NULL,
  variant = c("Primer3Plus", "Chester1993", "QuikChange", "Schildkraut1965",
    "Wetmur1991_MELTING", "Wetmur1991_RNA", "Wetmur1991_RNA/DNA", "vonAhsen2001"),
  Na = 0,
  K = 0,
  Tris = 0,
  Mg = 0,
  dNTPs = 0,
  saltcorr = c("Schildkraut2010", "Wetmur1991", "SantaLucia1996", "SantaLucia1998-1",
    "Owczarzy2004", "Owczarzy2008"),
  mismatch = TRUE,
  DMSO = 0,
  fmd = 0,
  DMSOfactor = 0.75,
  fmdfactor = 0.65,
  fmdmethod = c("concentration", "molar")
)
}
\arguments{
\item{ntseq}{Sequence (5' to 3') of one strand of the nucleic acid duplex as string or vector of characters.}

\item{ambiguous}{Ambiguous bases are taken into account to compute the G and C content when ambiguous is TRUE.}

\item{userset}{A vector of four coefficient values. Usersets override value sets.}

\item{variant}{Empirical constants coefficient with 8 variant: Chester1993, QuikChange, Schildkraut1965, Wetmur1991_MELTING, Wetmur1991_RNA, Wetmur1991_RNA/DNA, Primer3Plus and vonAhsen2001}

\item{Na}{Millimolar concentration of Na, default is 0}

\item{K}{Millimolar concentration of K, default is 0}

\item{Tris}{Millimolar concentration of Tris, default is 0}

\item{Mg}{Millimolar concentration of Mg, default is 0}

\item{dNTPs}{Millimolar concentration of dNTPs, default is 0}

\item{saltcorr}{Salt correction method should be chosen when provide 'userset'. Options are "Schildkraut2010", "Wetmur1991","SantaLucia1996","SantaLucia1998-1","Owczarzy2004","Owczarzy2008". Note that "SantaLucia1998-2" is not available for this function.}

\item{mismatch}{If 'True' (default) every 'X' in the sequence is counted as mismatch}

\item{DMSO}{Percent DMSO}

\item{fmd}{Formamide concentration in percentage (fmdmethod="concentration") or molar (fmdmethod="molar").}

\item{DMSOfactor}{Coeffecient of Tm decreases per percent DMSO. Default=0.75 von Ahsen N (2001) <PMID:11673362>. Other published values are 0.5, 0.6 and 0.675.}

\item{fmdfactor}{Coeffecient of Tm decrease per percent formamide. Default=0.65. Several papers report factors between 0.6 and 0.72.}

\item{fmdmethod}{"concentration" method for formamide concentration in percentage and "molar" for formamide concentration in molar}
}
\description{
Calculate the melting temperature using empirical formulas based on GC content with different options
}
\details{
Empirical constants coefficient with 8 variant:

Chester1993: Tm = 69.3 + 0.41(Percentage_GC) - 650/N
  
QuikChange: Tm = 81.5 + 0.41(Percentage_GC) - 675/N - Percentage_mismatch
  
Schildkraut1965: Tm = 81.5 + 0.41(Percentage_GC) - 675/N + 16.6 x log[Na+]
  
Wetmur1991_MELTING: Tm = 81.5 + 0.41(Percentage_GC) - 500/N + 16.6 x log([Na+]/(1.0 + 0.7 x [Na+])) - Percentage_mismatch
  
Wetmur1991_RNA: Tm = 78 + 0.7(Percentage_GC) - 500/N + 16.6 x log([Na+]/(1.0 + 0.7 x [Na+])) - Percentage_mismatch
  
Wetmur1991_RNA/DNA: Tm = 67 + 0.8(Percentage_GC) - 500/N + 16.6 x log([Na+]/(1.0 + 0.7 x [Na+])) - Percentage_mismatch
  
Primer3Plus: Tm = 81.5 + 0.41(Percentage_GC) - 600/N + 16.6 x log[Na+]
  
vonAhsen2001: Tm = 77.1 + 0.41(Percentage_GC) - 528/N + 11.7 x log[Na+]
}
\examples{

ntseq <- c("ATCGTGCGTAGCAGTACGATCAGTAG")
out <- Tm_GC(ntseq,ambiguous=TRUE,variant="Primer3Plus",Na=50,mismatch=TRUE)
out
out$Tm
out$Options

}
\references{
Marmur J , Doty P . Determination of the base composition of deoxyribonucleic acid from its thermal denaturation temperature.[J]. Journal of Molecular Biology, 1962, 5(1):109-118.

Schildkraut C . Dependence of the melting temperature of DNA on salt concentration[J]. Biopolymers, 2010, 3(2):195-208.

Wetmur J G . DNA Probes: Applications of the Principles of Nucleic Acid Hybridization[J]. CRC Critical Reviews in Biochemistry, 1991, 26(3-4):33.

Untergasser A , Cutcutache I , Koressaar T , et al. Primer3--new capabilities and interfaces[J]. Nucleic Acids Research, 2012, 40(15):e115-e115.

von Ahsen N, Wittwer CT, Schutz E , et al. Oligonucleotide melting temperatures under PCR conditions: deoxynucleotide Triphosphate and Dimethyl sulfoxide concentrations with comparison to alternative empirical formulas. Clin Chem 2001, 47:1956-1961.
}
\author{
Junhui Li
}