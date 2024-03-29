% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Tm_NN.R
\name{Tm_NN}
\alias{Tm_NN}
\title{Calculate melting temperature using nearest neighbor thermodynamics}
\usage{
Tm_NN(
  ntseq,
  ambiguous = FALSE,
  comSeq = NULL,
  shift = 0,
  nn_table = c("DNA_NN4", "DNA_NN1", "DNA_NN2", "DNA_NN3", "RNA_NN1", "RNA_NN2",
    "RNA_NN3", "R_DNA_NN1"),
  tmm_table = "DNA_TMM1",
  imm_table = "DNA_IMM1",
  de_table = c("DNA_DE1", "RNA_DE1"),
  dnac1 = 25,
  dnac2 = 25,
  selfcomp = FALSE,
  Na = 0,
  K = 0,
  Tris = 0,
  Mg = 0,
  dNTPs = 0,
  saltcorr = c("Schildkraut2010", "Wetmur1991", "SantaLucia1996", "SantaLucia1998-1",
    "SantaLucia1998-2", "Owczarzy2004", "Owczarzy2008"),
  DMSO = 0,
  fmd = 0,
  DMSOfactor = 0.75,
  fmdfactor = 0.65,
  fmdmethod = c("concentration", "molar")
)
}
\arguments{
\item{ntseq}{Sequence (5' to 3') of one strand of the nucleic acid duplex
as string or vector of characters.}

\item{ambiguous}{Ambiguous bases are taken into account to compute the G
and C content when ambiguous is TRUE.Default is FALSE.}

\item{comSeq}{Complementary sequence. The sequence of the template/target in 3'->5' direction}

\item{shift}{Shift of the primer/probe sequence on the template/target sequence, default=0. 
for example: when shift=0, the first nucleotide base at 5` end of primer align to first one at 3` 
end of template. When shift=-1, the second nucleotide base at 5` end of primer align to first one at 3` 
end of template.

When shift=1, the first nucleotide base at 5` end of primer align to second one at 3` end 
of template. The shift parameter is necessary to align  primer/probe and template/target 
if they have different lengths or if they should have dangling ends.}

\item{nn_table}{Thermodynamic NN values, eight tables are implemented.

For DNA/DNA hybridizations:
  DNA_NN1,DNA_NN2,DNA_NN3,DNA_NN4

For RNA/RNA hybridizations:
  RNA_NN1,RNA_NN2,RNA_NN3

For RNA/DNA hybridizations:
  R_DNA_NN1}

\item{tmm_table}{Thermodynamic values for terminal mismatches. Default: DNA_TMM1}

\item{imm_table}{Thermodynamic values for internal mismatches, may include insosine mismatches. Default: DNA_IMM1}

\item{de_table}{Thermodynamic values for dangling ends. DNA_DE1(default) and RNA_DE1}

\item{dnac1}{Concentration of the higher concentrated strand [nM]. Typically this will
be the primer (for PCR) or the probe. Default=25.}

\item{dnac2}{Concentration of the lower concentrated strand [nM].}

\item{selfcomp}{Sequence self-complementary, default=False. If 'True' 
the primer is thought binding to itself, thus dnac2 is not considered.}

\item{Na}{Millimolar concentration of Na, default is 0}

\item{K}{Millimolar concentration of K, default is 0}

\item{Tris}{Millimolar concentration of Tris, default is 0}

\item{Mg}{Millimolar concentration of Mg, default is 0}

\item{dNTPs}{Millimolar concentration of dNTPs, default is 0}

\item{saltcorr}{Salt correction method should be chosen when provide 'userset' 
Options are "Schildkraut2010", "Wetmur1991","SantaLucia1996","SantaLucia1998-1",
"SantaLucia1998-2","Owczarzy2004","Owczarzy2008". Note that NA means no salt correction.}

\item{DMSO}{Percent DMSO}

\item{fmd}{Formamide concentration in percentage (fmdmethod="concentration") or molar (fmdmethod="molar").}

\item{DMSOfactor}{Coeffecient of Tm decreases per percent DMSO. Default=0.75 von Ahsen N (2001) <PMID:11673362>. Other published values are 0.5, 0.6 and 0.675.}

\item{fmdfactor}{Coeffecient of Tm decrease per percent formamide. Default=0.65. Several papers report factors between 0.6 and 0.72.}

\item{fmdmethod}{"concentration" method for formamide concentration in percentage and "molar" for formamide concentration in molar.}
}
\description{
Calculate melting temperature using nearest neighbor thermodynamics
}
\details{
DNA_NN1: Breslauer K J (1986) <doi:10.1073/pnas.83.11.3746>
 
 DNA_NN2: Sugimoto N (1996) <doi:10.1093/nar/24.22.4501>
 
 DNA_NN3: Allawi H (1998) <doi:10.1093/nar/26.11.2694>
 
 DNA_NN4: SantaLucia J (2004) <doi:10.1146/annurev.biophys.32.110601.141800>
 
 RNA_NN1: Freier S (1986) <doi:10.1073/pnas.83.24.9373>
 
 RNA_NN2: Xia T (1998) <doi:10.1021/bi9809425>
 
 RNA_NN3: Chen JL (2012) <doi:10.1021/bi3002709>
 
 R_DNA_NN1: Sugimoto N (1995)<doi:10.1016/S0048-9697(98)00088-6>
 
 DNA_TMM1: Bommarito S (2000)  <doi:10.1093/nar/28.9.1929>
 
 DNA_IMM1: Peyret N (1999) <doi:10.1021/bi9825091> & Allawi H T (1997) <doi:10.1021/bi962590c> & Santalucia N (2005) <doi:10.1093/nar/gki918>
 
 DNA_DE1: Bommarito S (2000) <doi:10.1093/nar/28.9.1929>
 
 RNA_DE1: Turner D H (2010) <doi:10.1093/nar/gkp892>
}
\examples{

ntseq <- c("AAAATTTTTTTCCCCCCCCCCCCCCGGGGGGGGGGGGTGTGCGCTGC")
out <- Tm_NN(ntseq,Na=50)
out
out$Options

}
\references{
Breslauer K J , Frank R , Blocker H , et al. Predicting DNA duplex stability from the base sequence.[J]. Proceedings of the National Academy of Sciences, 1986, 83(11):3746-3750.

Sugimoto N , Nakano S , Yoneyama M , et al. Improved Thermodynamic Parameters and Helix Initiation Factor to Predict Stability of DNA Duplexes[J]. Nucleic Acids Research, 1996, 24(22):4501-5.

Allawi, H. Thermodynamics of internal C.T mismatches in DNA[J]. Nucleic Acids Research, 1998, 26(11):2694-2701.

Hicks L D , Santalucia J . The thermodynamics of DNA structural motifs.[J]. Annual Review of Biophysics & Biomolecular Structure, 2004, 33(1):415-440.

Freier S M , Kierzek R , Jaeger J A , et al. Improved free-energy parameters for predictions of RNA duplex stability.[J]. Proceedings of the National Academy of Sciences, 1986, 83(24):9373-9377.

Xia T , Santalucia , J , Burkard M E , et al. Thermodynamic Parameters for an Expanded Nearest-Neighbor Model for Formation of RNA Duplexes with Watson-Crick Base Pairs,[J]. Biochemistry, 1998, 37(42):14719-14735.

Chen J L , Dishler A L , Kennedy S D , et al. Testing the Nearest Neighbor Model for Canonical RNA Base Pairs: Revision of GU Parameters[J]. Biochemistry, 2012, 51(16):3508-3522.

Bommarito S, Peyret N, Jr S L. Thermodynamic parameters for DNA sequences with dangling ends[J]. Nucleic Acids Research, 2000, 28(9):1929-1934.

Turner D H , Mathews D H . NNDB: the nearest neighbor parameter database for predicting stability of nucleic acid secondary structure[J]. Nucleic Acids Research, 2010, 38(Database issue):D280-D282.

Sugimoto N , Nakano S I , Katoh M , et al. Thermodynamic Parameters To Predict Stability of RNA/DNA Hybrid Duplexes[J]. Biochemistry, 1995, 34(35):11211-11216.

Allawi H, SantaLucia J: Thermodynamics and NMR of internal G-T mismatches in DNA. Biochemistry 1997, 36:10581-10594.

Santalucia N E W J . Nearest-neighbor thermodynamics of deoxyinosine pairs in DNA duplexes[J]. Nucleic Acids Research, 2005, 33(19):6258-67.

Peyret N , Seneviratne P A , Allawi H T , et al. Nearest-Neighbor Thermodynamics and NMR of DNA Sequences with Internal A-A, C-C, G-G, and T-T Mismatches, [J]. Biochemistry, 1999, 38(12):3468-3477.
}
\author{
Junhui Li
}
