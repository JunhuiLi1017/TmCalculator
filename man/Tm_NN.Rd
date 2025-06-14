% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tm_nn.R
\name{tm_nn}
\alias{tm_nn}
\title{Calculate melting temperature using nearest neighbor thermodynamics}
\usage{
tm_nn(
  gr_seq,
  ambiguous = FALSE,
  shift = 0,
  nn_table = c("DNA_NN_SantaLucia_2004", "DNA_NN_Breslauer_1986", "DNA_NN_Sugimoto_1996",
    "DNA_NN_Allawi_1998", "RNA_NN_Freier_1986", "RNA_NN_Xia_1998", "RNA_NN_Chen_2012",
    "RNA_DNA_NN_Sugimoto_1995"),
  tmm_table = "DNA_TMM_Bommarito_2000",
  imm_table = "DNA_IMM_Peyret_1999",
  de_table = c("DNA_DE_Bommarito_2000", "RNA_DE_Turner_2010"),
  dnac_high = 25,
  dnac_low = 25,
  self_comp = FALSE,
  Na = 50,
  K = 0,
  Tris = 0,
  Mg = 0,
  dNTPs = 0,
  salt_corr_method = c("Schildkraut2010", "Wetmur1991", "SantaLucia1996",
    "SantaLucia1998-1", "Owczarzy2004", "Owczarzy2008"),
  DMSO = 0,
  formamide_value_unit = list(value = 0, unit = "percent"),
  dmso_factor = 0.75,
  formamide_factor = 0.65
)
}
\arguments{
\item{gr_seq}{Pre-processed sequence(s) in 5' to 3' direction. This should be the output from
to_genomic_ranges() function.}

\item{ambiguous}{Logical value controlling how ambiguous bases are handled:
- TRUE: Ambiguous bases (e.g., N, R, Y) are included in calculations
- FALSE (default): Ambiguous bases are excluded from calculations}

\item{shift}{Integer value controlling the alignment offset between primer and template sequences.
  Visual representation of different shift values:

  shift = 0 (default):
  Primer:    5' ATGCG 3'
  Template:  3' TACGC 5'

  shift = -1:
  Primer:    5' ATGCG 3'
  Template:  3'  TACGC 5'
              ^

  shift = 1:
  Primer:    5'  ATGCG 3'
  Template:  3' TACGC 5'
              ^

  The shift parameter is necessary when:
  - Sequences have different lengths
  - Dangling ends are required
  - Specific alignment positions are needed}

\item{nn_table}{Thermodynamic nearest-neighbor parameters for different nucleic acid hybridizations.
  Eight parameter sets are available, organized by hybridization type:

  DNA/DNA hybridizations:
  - "DNA_NN_Breslauer_1986": Original DNA/DNA parameters
  - "DNA_NN_Sugimoto_1996": Improved DNA/DNA parameters
  - "DNA_NN_Allawi_1998": DNA/DNA parameters with internal mismatch corrections
  - "DNA_NN_SantaLucia_2004": Updated DNA/DNA parameters

  RNA/RNA hybridizations:
  - "RNA_NN_Freier_1986": Original RNA/RNA parameters
  - "RNA_NN_Xia_1998": Improved RNA/RNA parameters
  - "RNA_NN_Chen_2012": Updated RNA/RNA parameters with GU pair corrections

  RNA/DNA hybridizations:
  - "RNA_DNA_NN_Sugimoto_1995": RNA/DNA hybridization parameters}

\item{tmm_table}{Thermodynamic parameters for terminal mismatches. Default: "DNA_TMM_Bommarito_2000"
These parameters account for mismatches at the ends of the duplex.}

\item{imm_table}{Thermodynamic parameters for internal mismatches. Default: "DNA_IMM_Peyret_1999"
These parameters account for mismatches within the duplex, including inosine mismatches.}

\item{de_table}{Thermodynamic parameters for dangling ends. Default: "DNA_DE_Bommarito_2000"
Available options:
- "DNA_DE_Bommarito_2000": Parameters for DNA dangling ends
- "RNA_DE_Turner_2010": Parameters for RNA dangling ends}

\item{dnac_high}{Concentration of the higher concentrated strand in nM. Default: 25
Typically this is the primer (for PCR) or the probe concentration.}

\item{dnac_low}{Concentration of the lower concentrated strand in nM. Default: 25
This is typically the template concentration.}

\item{self_comp}{Logical value indicating if the sequence is self-complementary:
- TRUE: Sequence can bind to itself, dnac_low is ignored
- FALSE (default): Sequence binds to a different complementary sequence}

\item{Na}{Millimolar concentration of sodium ions. Default: 50}

\item{K}{Millimolar concentration of potassium ions. Default: 0}

\item{Tris}{Millimolar concentration of Tris buffer. Default: 0}

\item{Mg}{Millimolar concentration of magnesium ions. Default: 0}

\item{dNTPs}{Millimolar concentration of deoxynucleotide triphosphates. Default: 0}

\item{salt_corr_method}{Method for calculating salt concentration corrections to the melting temperature.
Available options:
- "Schildkraut2010": Updated salt correction method
- "Wetmur1991": Classic salt correction method
- "SantaLucia1996": DNA-specific salt correction
- "SantaLucia1998-1": Improved DNA salt correction
- "SantaLucia1998-2": Alternative DNA salt correction
- "Owczarzy2004": Comprehensive salt correction
- "Owczarzy2008": Updated comprehensive salt correction
Note: Setting to NA disables salt correction}

\item{DMSO}{Percent DMSO concentration in the reaction mixture. Default: 0
DMSO can lower the melting temperature of nucleic acid duplexes.}

\item{formamide_value_unit}{A list containing formamide concentration information:
- value: numeric value of formamide concentration
- unit: character string specifying the unit ("percent" or "molar")
Default: list(value=0, unit="percent")}

\item{dmso_factor}{Coefficient of melting temperature (Tm) decrease per percent DMSO.
Default: 0.75 (von Ahsen N, 2001, PMID:11673362)
Other published values: 0.5, 0.6, 0.675}

\item{formamide_factor}{Coefficient of melting temperature (Tm) decrease per percent formamide.
Default: 0.65
Literature reports values ranging from 0.6 to 0.72}
}
\description{
Calculate melting temperature using nearest neighbor thermodynamics. The function checks if all
sequence combinations in the input sequence are present in the thermodynamic parameter tables
before performing calculations.
}
\details{
DNA_NN_Breslauer_1986: Breslauer K J (1986) <doi:10.1073/pnas.83.11.3746>
 
 DNA_NN_Sugimoto_1996: Sugimoto N (1996) <doi:10.1093/nar/24.22.4501>
 
 DNA_NN_Allawi_1998: Allawi H (1998) <doi:10.1093/nar/26.11.2694>
 
 DNA_NN_SantaLucia_2004: SantaLucia J (2004) <doi:10.1146/annurev.biophys.32.110601.141800>
 
 RNA_NN_Freier_1986: Freier S (1986) <doi:10.1073/pnas.83.24.9373>
 
 RNA_NN_Xia_1998: Xia T (1998) <doi:10.1021/bi9809425>
 
 RNA_NN_Chen_2012: Chen JL (2012) <doi:10.1021/bi3002709>
 
 RNA_DNA_NN_Sugimoto_1995: Sugimoto N (1995)<doi:10.1016/S0048-9697(98)00088-6>
 
 DNA_TMM_Bommarito_2000: Bommarito S (2000)  <doi:10.1093/nar/28.9.1929>
 
 DNA_IMM_Peyret_1999: Peyret N (1999) <doi:10.1021/bi9825091> & Allawi H T (1997) <doi:10.1021/bi962590c> & Santalucia N (2005) <doi:10.1093/nar/gki918>
 
 DNA_DE_Bommarito_2000: Bommarito S (2000) <doi:10.1093/nar/28.9.1929>
 
 RNA_DE_Turner_2010: Turner D H (2010) <doi:10.1093/nar/gkp892>
}
\examples{

input_seq <- c("AAAATTTTTTTCCCCCCCCCCCCCCGGGGGGGGGGGGTGTGCGCTGC",
"AAAATTTTTTTCCCCCCCCCCCCCCGGGGGGGGGGGGTGTGCGCTGC")
gr_seq <- to_genomic_ranges(input_seq)
out <- tm_nn(gr_seq, Na=50)
out
out$tm_nn$Options

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
