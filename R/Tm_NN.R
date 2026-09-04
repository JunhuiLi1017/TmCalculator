#' Calculate melting temperature using nearest neighbor thermodynamics
#'
#' Calculates melting temperature (Tm) from nearest-neighbor (NN) thermodynamics,
#' summing the stacking enthalpies and entropies of adjacent base-pair steps and
#' applying initiation, symmetry, salt and chemical corrections. Terminal
#' mismatches, internal mismatches and dangling ends are supported through
#' dedicated parameter tables. The function verifies that every dinucleotide step
#' in the input sequence is present in the selected tables before calculating.
#'
#' @section Choosing a parameter set:
#' Parameter sets fall into two families that are handled differently.
#'
#' The two families are distinguished by one mechanical criterion: a set is
#' condition-specific if and only if it carries a \code{salt_mM} attribute.
#'
#' \strong{Reference-salt sets} (no \code{salt_mM}; Breslauer 1986,
#' Sugimoto 1996, Allawi 1998, SantaLucia 2004, Freier 1986, Xia 1998,
#' Chen 2012, Zuber 2022, Sugimoto 1995) were fitted at a single reference
#' sodium concentration, and other conditions are reached by applying one of
#' the \code{salt_method} correction formulas.
#'
#' \strong{Condition-specific sets} (the Weber/VarGibbs series, Banerjee 2020,
#' and the molecular-crowding sets of Ghosh 2020 and Ghosh 2023) were instead
#' fitted directly at a stated sodium concentration and are intended to
#' \emph{replace} salt correction rather than be corrected. When the requested
#' \code{Na} matches the set's \code{salt_mM} value,
#' salt correction is skipped automatically; when it does not, the correction is
#' still applied but a warning is issued, since correcting an already
#' condition-specific set double-counts the ionic effect. Whether a correction
#' was applied is recorded in the returned \code{options}.
#'
#' As a rule of thumb, pick the set whose fitted salt is closest to your
#' experimental condition rather than correcting a distant one.
#' 
#' @param gr_seq Pre-processed sequence(s) in 5' to 3' direction. This should be the output from
#'   to_genomic_ranges() function.
#' 
#' @param ambiguous Logical value controlling how ambiguous bases are handled:
#'   - TRUE: Ambiguous bases (e.g., N, R, Y) are included in calculations
#'   - FALSE (default): Ambiguous bases are excluded from calculations
#' 
#' @param shift Integer value controlling the alignment offset between primer and template sequences.
#'   Visual representation of different shift values:
#' 
#'   shift = 0 (default):
#'   Primer:    5' ATGCG 3'
#'   Template:  3' TACGC 5'
#' 
#'   shift = -1:
#'   Primer:    5' ATGCG 3'
#'   Template:  3'  TACGC 5'
#'               ^
#' 
#'   shift = 1:
#'   Primer:    5'  ATGCG 3'
#'   Template:  3' TACGC 5'
#'               ^
#' 
#'   The shift parameter is necessary when:
#'   - Sequences have different lengths
#'   - Dangling ends are required
#'   - Specific alignment positions are needed
#' 
#' @param nn_table Thermodynamic nearest-neighbor parameters for different nucleic acid hybridizations.
#'   Parameter sets are listed below by hybridization type.
#'   Sets marked with a sodium concentration were fitted at that condition and
#'   are not salt-corrected further (see the "Choosing a parameter set" section).
#'
#'   DNA/DNA hybridizations, reference salt:
#'   - "DNA_NN_Breslauer_1986": Original DNA/DNA parameters
#'   - "DNA_NN_Sugimoto_1996": Improved DNA/DNA parameters
#'   - "DNA_NN_Allawi_1998": DNA/DNA parameters with internal mismatch corrections
#'   - "DNA_NN_SantaLucia_2004": Unified DNA/DNA parameters (default)
#'
#'   DNA/DNA hybridizations, melting-temperature optimized (Weber 2015):
#'   - "DNA_NN_Weber_2015": Combined dataset, 1020 mM. Recommended when a
#'     salt-optimized DNA set is wanted at high salt
#'   - "DNA_NN_Weber_OW04_69", "...119", "...220", "...621", "...1020":
#'     fitted independently at 69, 119, 220, 621 and 1020 mM sodium
#'
#'   DNA/DNA under molecular crowding (cell-like rather than dilute solution):
#'   - "DNA_NN_Ghosh_2020_PEG200": fitted in 40 wt% PEG200 with 100 mM NaCl
#'
#'   RNA/RNA hybridizations, reference salt:
#'   - "RNA_NN_Freier_1986": Original RNA/RNA parameters
#'   - "RNA_NN_Xia_1998": Improved RNA/RNA parameters
#'   - "RNA_NN_Chen_2012": Updated RNA/RNA parameters with GU pair corrections
#'   - "RNA_NN_Zuber_2022": Successor to Xia 1998 with improved end effects.
#'     The terminal-AU penalty is replaced by end terms that depend on the
#'     penultimate base pair, applied automatically from a companion table
#'
#'   RNA/RNA under molecular crowding (cell-like rather than dilute solution):
#'   - "RNA_NN_Ghosh_2023_PEG200": fitted in 40 wt% PEG200 with 100 mM NaCl,
#'     and shown to describe duplexes in an intracellular cation composition
#'
#'   RNA/RNA hybridizations, salt-optimized (Ferreira 2019). VIF (variable
#'   initiation factors) gave better cross-validation than FIF (fixed):
#'   - "RNA_NN_Weber_VIF_71", "...121", "...221", "...621", "...1021"
#'   - "RNA_NN_Weber_FIF_71", "...121", "...221", "...621", "...1021"
#'
#'   RNA/DNA hybridizations:
#'   - "RNA_DNA_NN_Sugimoto_1995": RNA/DNA hybridization parameters
#'   - "RNA_DNA_NN_Weber_2019_FT": curve-fitting derived, 1000 mM. Best
#'     performing high-salt hybrid set in Basilio Barbosa (2019)
#'   - "RNA_DNA_NN_Weber_2019_VH": van't Hoff derived, 1000 mM
#'   - "RNA_DNA_NN_Weber_2019_LS": low salt, 100 mM
#'   - "RNA_DNA_NN_Banerjee_2020": improved hybrid parameters fitted at a
#'     physiological condition (100 mM NaCl), Banerjee et al. (2020)
#' 
#' @param tmm_table Thermodynamic parameters for terminal mismatches. Default: "DNA_TMM_Bommarito_2000"
#'   These parameters account for mismatches at the ends of the duplex.
#' 
#' @param imm_table Thermodynamic parameters for internal mismatches. Default: "DNA_IMM_Peyret_1999"
#'   These parameters account for mismatches within the duplex, including inosine mismatches.
#' 
#' @param de_table Thermodynamic parameters for dangling ends. Default: "DNA_DE_Bommarito_2000"
#'   Available options:
#'   - "DNA_DE_Bommarito_2000": Parameters for DNA dangling ends
#'   - "RNA_DE_Turner_2010": Parameters for RNA dangling ends
#' 
#' @param dnac_high Concentration of the higher concentrated strand in nM. Default: 25
#'   Typically this is the primer (for PCR) or the probe concentration.
#' 
#' @param dnac_low Concentration of the lower concentrated strand in nM. Default: 25
#'   This is typically the template concentration.
#' 
#' @param self_comp Logical value indicating if the sequence is self-complementary:
#'   - TRUE: Sequence can bind to itself, dnac_low is ignored
#'   - FALSE (default): Sequence binds to a different complementary sequence
#' 
#' @param Na Millimolar concentration of sodium ions. Default: 50
#' 
#' @param K Millimolar concentration of potassium ions. Default: 0
#' 
#' @param Tris Millimolar concentration of Tris buffer. Default: 0
#' 
#' @param Mg Millimolar concentration of magnesium ions. Default: 0
#' 
#' @param dNTPs Millimolar concentration of deoxynucleotide triphosphates. Default: 0
#' 
#' @param salt_method Salt correction method. Options are:
#'   Available options:
#'   - "Schildkraut2010": Updated salt correction method
#'   - "Wetmur1991": Classic salt correction method
#'   - "SantaLucia1996": DNA-specific salt correction
#'   - "SantaLucia1998-1": Improved DNA salt correction
#'   - "Owczarzy2004": Comprehensive salt correction
#'   - "Owczarzy2008": Updated comprehensive salt correction
#'   - "none": Disables salt correction entirely
#'   Note: Parameter sets fitted at a specific sodium concentration (those
#'   carrying a "salt_mM" attribute, e.g. the Weber/VarGibbs series and
#'   Banerjee 2020) already
#'   account for salt. When the requested \code{Na} matches the concentration
#'   such a set was fitted at, salt correction is skipped automatically; when
#'   it does not, a warning is issued.
#' 
#' @param DMSO Percent DMSO concentration in the reaction mixture. Default: 0
#'   DMSO can lower the melting temperature of nucleic acid duplexes.
#' 
#' @param formamide_unit Formamide concentration as `list(value, unit)`. Default: list(value = 0, unit = "percent")
#'   - value: numeric value of formamide concentration
#'   - unit: character string specifying the unit ("percent" or "molar")
#'   Default: list(value=0, unit="percent")
#' 
#' @param dmso_factor Coefficient of melting temperature (Tm) decrease per percent DMSO.
#'   Default: 0.75 (von Ahsen N, 2001, PMID:11673362)
#'   Other published values: 0.5, 0.6, 0.675
#' 
#' @param formamide_factor Coefficient of melting temperature (Tm) decrease per percent formamide.
#'   Default: 0.65
#'   Literature reports values ranging from 0.6 to 0.72
#'
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object
#'   specifying the parallel backend. \code{BiocParallel::SnowParam(n)} is
#'   recommended; \code{BiocParallel::MulticoreParam} can be slower than
#'   serial on large inputs because copy-on-write interacts badly with R's
#'   garbage collector (see \code{\link{tm_calculate}}). The default,
#'   \code{BiocParallel::SerialParam()}, runs serially. Sequences are split
#'   into one chunk per worker, so parallelization pays off for large inputs
#'   (e.g. genome-wide windows) rather than a handful of primers.
#'
#' @details
#' 
#'  DNA_NN_Breslauer_1986: Breslauer K J (1986) <doi:10.1073/pnas.83.11.3746>
#'  
#'  DNA_NN_Sugimoto_1996: Sugimoto N (1996) <doi:10.1093/nar/24.22.4501>
#'  
#'  DNA_NN_Allawi_1998: Allawi H (1998) <doi:10.1093/nar/26.11.2694>
#'  
#'  DNA_NN_SantaLucia_2004: SantaLucia J (2004) <doi:10.1146/annurev.biophys.32.110601.141800>
#'  
#'  RNA_NN_Freier_1986: Freier S (1986) <doi:10.1073/pnas.83.24.9373>
#'  
#'  RNA_NN_Xia_1998: Xia T (1998) <doi:10.1021/bi9809425>
#'  
#'  RNA_NN_Chen_2012: Chen JL (2012) <doi:10.1021/bi3002709>
#'  
#'  RNA_DNA_NN_Sugimoto_1995: Sugimoto N (1995)<doi:10.1016/S0048-9697(98)00088-6>
#'
#'  The following sets were derived by melting-temperature optimization and are
#'  fitted at the sodium concentration given in parentheses. They are not
#'  salt-corrected further; see the \dQuote{Choosing a parameter set} section.
#'
#'  DNA_NN_Weber_2015 (1020 mM): Weber G (2015) <doi:10.1093/bioinformatics/btu751>
#'
#'  DNA_NN_Weber_OW04_69 (69 mM), DNA_NN_Weber_OW04_119 (119 mM),
#'  DNA_NN_Weber_OW04_220 (220 mM), DNA_NN_Weber_OW04_621 (621 mM),
#'  DNA_NN_Weber_OW04_1020 (1020 mM): Weber G (2015) <doi:10.1093/bioinformatics/btu751>
#'
#'  RNA_NN_Weber_VIF_71 (71 mM), RNA_NN_Weber_VIF_121 (121 mM),
#'  RNA_NN_Weber_VIF_221 (221 mM), RNA_NN_Weber_VIF_621 (621 mM),
#'  RNA_NN_Weber_VIF_1021 (1021 mM): Ferreira I (2019) <doi:10.1016/j.chemphys.2019.01.016>,
#'  variable initiation factors
#'
#'  RNA_NN_Weber_FIF_71 (71 mM), RNA_NN_Weber_FIF_121 (121 mM),
#'  RNA_NN_Weber_FIF_221 (221 mM), RNA_NN_Weber_FIF_621 (621 mM),
#'  RNA_NN_Weber_FIF_1021 (1021 mM): Ferreira I (2019) <doi:10.1016/j.chemphys.2019.01.016>,
#'  fixed initiation factors
#'
#'  RNA_DNA_NN_Weber_2019_FT (1000 mM), RNA_DNA_NN_Weber_2019_VH (1000 mM),
#'  RNA_DNA_NN_Weber_2019_LS (100 mM): Basilio Barbosa V (2019) <doi:10.1016/j.bpc.2019.106189>
#'
#'  RNA_DNA_NN_Banerjee_2020 (100 mM): Banerjee D (2020) <doi:10.1093/nar/gkaa572>
#'
#'  RNA_NN_Zuber_2022: Zuber J (2022) <doi:10.1093/nar/gkac261>
#'
#'  RNA_NN_Ghosh_2023_PEG200 (100 mM, 40 wt% PEG200): Ghosh S (2023) <doi:10.1093/nar/gkad020>
#'
#'  DNA_NN_Ghosh_2020_PEG200 (100 mM, 40 wt% PEG200): Ghosh S (2020) <doi:10.1073/pnas.1920886117>
#'
#'  DNA_TMM_Bommarito_2000: Bommarito S (2000)  <doi:10.1093/nar/28.9.1929>
#'  
#'  DNA_IMM_Peyret_1999: Peyret N (1999) <doi:10.1021/bi9825091> & Allawi H T (1997) <doi:10.1021/bi962590c> & Santalucia N (2005) <doi:10.1093/nar/gki918>
#'  
#'  DNA_DE_Bommarito_2000: Bommarito S (2000) <doi:10.1093/nar/28.9.1929>
#'  
#'  RNA_DE_Turner_2010: Turner D H (2010) <doi:10.1093/nar/gkp892>
#' 
#' @references 
#' 
#' Breslauer K J , Frank R , Blocker H , et al. Predicting DNA duplex stability from the base sequence.[J]. Proceedings of the National Academy of Sciences, 1986, 83(11):3746-3750.
#' 
#' Sugimoto N , Nakano S , Yoneyama M , et al. Improved Thermodynamic Parameters and Helix Initiation Factor to Predict Stability of DNA Duplexes[J]. Nucleic Acids Research, 1996, 24(22):4501-5.
#' 
#' Allawi, H. Thermodynamics of internal C.T mismatches in DNA[J]. Nucleic Acids Research, 1998, 26(11):2694-2701.
#' 
#' Hicks L D , Santalucia J . The thermodynamics of DNA structural motifs.[J]. Annual Review of Biophysics & Biomolecular Structure, 2004, 33(1):415-440.
#' 
#' Freier S M , Kierzek R , Jaeger J A , et al. Improved free-energy parameters for predictions of RNA duplex stability.[J]. Proceedings of the National Academy of Sciences, 1986, 83(24):9373-9377.
#' 
#' Xia T , Santalucia , J , Burkard M E , et al. Thermodynamic Parameters for an Expanded Nearest-Neighbor Model for Formation of RNA Duplexes with Watson-Crick Base Pairs,[J]. Biochemistry, 1998, 37(42):14719-14735.
#' 
#' Chen J L , Dishler A L , Kennedy S D , et al. Testing the Nearest Neighbor Model for Canonical RNA Base Pairs: Revision of GU Parameters[J]. Biochemistry, 2012, 51(16):3508-3522.
#' 
#' Bommarito S, Peyret N, Jr S L. Thermodynamic parameters for DNA sequences with dangling ends[J]. Nucleic Acids Research, 2000, 28(9):1929-1934.
#' 
#' Turner D H , Mathews D H . NNDB: the nearest neighbor parameter database for predicting stability of nucleic acid secondary structure[J]. Nucleic Acids Research, 2010, 38(Database issue):D280-D282.
#' 
#' Sugimoto N , Nakano S I , Katoh M , et al. Thermodynamic Parameters To Predict Stability of RNA/DNA Hybrid Duplexes[J]. Biochemistry, 1995, 34(35):11211-11216.
#' 
#' Allawi H, SantaLucia J: Thermodynamics and NMR of internal G-T mismatches in DNA. Biochemistry 1997, 36:10581-10594.
#' 
#' Weber G. Optimization method for obtaining nearest-neighbour DNA entropies and enthalpies directly from melting temperatures. Bioinformatics, 2015, 31(6):871-877.
#' 
#' Ferreira I, Jolley E A, Znosko B M, Weber G. Replacing salt correction factors with optimized RNA nearest-neighbour enthalpy and entropy parameters. Chemical Physics, 2019, 521:69-76.
#' 
#' Basilio Barbosa V, de Oliveira Martins E, Weber G. Nearest-neighbour parameters optimized for melting temperature prediction of DNA/RNA hybrids at high and low salt concentrations. Biophysical Chemistry, 2019, 251:106189.
#' 
#' Banerjee D, Tateishi-Karimata H, Ohyama T, et al. Improved nearest-neighbor parameters for the stability of RNA/DNA hybrids under a physiological condition. Nucleic Acids Research, 2020, 48(21):12042-12054.
#' 
#' Zuber J, Schroeder S J, Sun H, Turner D H, Mathews D H. Nearest neighbor rules for RNA helix folding thermodynamics: improved end effects. Nucleic Acids Research, 2022, 50(9):5251-5262.
#' 
#' Ghosh S, Takahashi S, Ohyama T, et al. Nearest-neighbor parameters for predicting DNA duplex stability in diverse molecular crowding conditions. Proceedings of the National Academy of Sciences, 2020, 117(25):14194-14201.
#' 
#' Ghosh S, Takahashi S, Banerjee D, et al. Nearest-neighbor parameters for the prediction of RNA duplex stability in diverse in vitro and cellular-like crowding conditions. Nucleic Acids Research, 2023, 51(9):4101-4111.
#' 
#' Santalucia N E W J . Nearest-neighbor thermodynamics of deoxyinosine pairs in DNA duplexes[J]. Nucleic Acids Research, 2005, 33(19):6258-67.
#' 
#' Peyret N , Seneviratne P A , Allawi H T , et al. Nearest-Neighbor Thermodynamics and NMR of DNA Sequences with Internal A-A, C-C, G-G, and T-T Mismatches, [J]. Biochemistry, 1999, 38(12):3468-3477.
#'
#' Weber G. Optimization method for obtaining nearest-neighbour DNA entropies and enthalpies directly from melting temperatures[J]. Bioinformatics, 2015, 31(6):871-877.
#'
#' Ferreira I, Jolley E A, Znosko B M, et al. Replacing salt correction factors with optimized RNA nearest-neighbour enthalpy and entropy parameters[J]. Chemical Physics, 2019, 521:69-76.
#'
#' Basilio Barbosa V, de Oliveira Martins E, Weber G. Nearest-neighbour parameters optimized for melting temperature prediction of DNA/RNA hybrids at high and low salt concentrations[J]. Biophysical Chemistry, 2019, 251:106189.
#'
#' Banerjee D, Tateishi-Karimata H, Ohyama T, Ghosh S, Endoh T, Takahashi S, Sugimoto N. Improved nearest-neighbor parameters for the stability of RNA/DNA hybrids under a physiological condition[J]. Nucleic Acids Research, 2020, 48(21):12042-12054.
#'
#' Zuber J, Schroeder S J, Sun H, Turner D H, Mathews D H. Nearest neighbor rules for RNA helix folding thermodynamics: improved end effects[J]. Nucleic Acids Research, 2022, 50(9):5251-5262.
#'
#' Ghosh S, Takahashi S, Banerjee D, Ohyama T, Endoh T, Tateishi-Karimata H, Sugimoto N. Nearest-neighbor parameters for the prediction of RNA duplex stability in diverse in vitro and cellular-like crowding conditions[J]. Nucleic Acids Research, 2023, 51(9):4101-4111.
#'
#' Ghosh S, Takahashi S, Ohyama T, Endoh T, Tateishi-Karimata H, Sugimoto N. Nearest-neighbor parameters for predicting DNA duplex stability in diverse molecular crowding conditions[J]. Proceedings of the National Academy of Sciences, 2020, 117(25):14194-14201.
#'
#' @return A \code{TmCalculator} list with:
#'   \item{\code{gr}}{The input \code{GRanges} with metadata columns \code{Tm}
#'     and \code{GC} (melting temperature in \eqn{^{\circ}}C and GC percent).
#'     Sequences that could not be evaluated receive \code{NA}.}
#'   \item{\code{options}}{The calculation parameters actually used, including
#'     the parameter tables and their citations, ion and additive
#'     concentrations, \code{Salt correction method}, and two fields recording
#'     how salt was handled: \code{Salt correction applied} (logical) and
#'     \code{Parameter set fitted at [Na+] (mM)}, which is \code{NA} for
#'     reference-salt sets. Regions skipped for containing \code{N} are listed
#'     in \code{Skipped regions containing N}.}
#'
#' @seealso \code{\link{tm_calculate}} for a single entry point to the
#'   nearest-neighbor, GC-content and Wallace methods.
#'
#' @author Junhui Li
#'
#' @importFrom BiocParallel bplapply bpnworkers SerialParam
#'
#' @examples
#'
#' input_seq <- c("AAAATTTTTTTCCCCCCCCCCCCCCGGGGGGGGGGGGTGTGCGCTGC",
#' "AAAATTTTTTTCCCCCCCCCCCCCCGGGGGGGGGGGGTGTGCGCTGC")
#' seqs <- to_genomic_ranges(input_seq)
#' out <- tm_nn(seqs, Na=50)
#' out
#'
#' # A parameter set fitted at a stated sodium concentration. Because Na
#' # matches the concentration the set was fitted at, salt correction is
#' # skipped automatically rather than applied on top of it.
#' out_ls <- tm_nn(seqs, nn_table = "RNA_DNA_NN_Weber_2019_LS", Na = 100)
#' out_ls$options[["Salt correction applied"]]
#' out_ls$options[["Parameter set fitted at [Na+] (mM)"]]
#'
#' @export tm_nn

tm_nn <- function(gr_seq,
                  ambiguous     = FALSE,
                  shift         = 0,
                  nn_table      = c("DNA_NN_SantaLucia_2004",
                                     "DNA_NN_Ghosh_2020_PEG200",
                                     "DNA_NN_Breslauer_1986",
                                     "DNA_NN_Sugimoto_1996",
                                     "DNA_NN_Allawi_1998",
                                     "RNA_NN_Freier_1986",
                                     "RNA_NN_Xia_1998",
                                     "RNA_NN_Chen_2012",
                                     "RNA_NN_Zuber_2022",
                                     "RNA_NN_Ghosh_2023_PEG200",
                                     "RNA_DNA_NN_Sugimoto_1995",
                                     "DNA_NN_Weber_2015",
                                     "DNA_NN_Weber_OW04_69",
                                     "DNA_NN_Weber_OW04_119",
                                     "DNA_NN_Weber_OW04_220",
                                     "DNA_NN_Weber_OW04_621",
                                     "DNA_NN_Weber_OW04_1020",
                                     "RNA_NN_Weber_VIF_71",
                                     "RNA_NN_Weber_VIF_121",
                                     "RNA_NN_Weber_VIF_221",
                                     "RNA_NN_Weber_VIF_621",
                                     "RNA_NN_Weber_VIF_1021",
                                     "RNA_NN_Weber_FIF_71",
                                     "RNA_NN_Weber_FIF_121",
                                     "RNA_NN_Weber_FIF_221",
                                     "RNA_NN_Weber_FIF_621",
                                     "RNA_NN_Weber_FIF_1021",
                                     "RNA_DNA_NN_Weber_2019_FT",
                                     "RNA_DNA_NN_Weber_2019_VH",
                                     "RNA_DNA_NN_Weber_2019_LS",
                                     "RNA_DNA_NN_Banerjee_2020"),
                  tmm_table      = "DNA_TMM_Bommarito_2000",
                  imm_table      = "DNA_IMM_Peyret_1999",
                  de_table       = c("DNA_DE_Bommarito_2000",
                                     "RNA_DE_Turner_2010"),
                  dnac_high      = 25,
                  dnac_low       = 25,
                  self_comp      = FALSE,
                  Na             = 50,
                  K              = 0,
                  Tris           = 0,
                  Mg             = 0,
                  dNTPs          = 0,
                  salt_method    = c("Schildkraut2010",
                                        "Wetmur1991",
                                        "SantaLucia1996",
                                        "SantaLucia1998-1",
                                        "Owczarzy2004",
                                        "Owczarzy2008",
                                        "none"),
                  DMSO           = 0,
                  formamide_unit = list(value = 0, unit = "percent"),
                  dmso_factor    = 0.75,
                  formamide_factor     = 0.65,
                  BPPARAM        = BiocParallel::SerialParam()) {

  # -- Validate args once ----------------------------------------------------
  nn_table <- match.arg(nn_table)
  tmm_table <- match.arg(tmm_table)
  imm_table <- match.arg(imm_table)
  de_table <- match.arg(de_table)
  salt_method <- match.arg(salt_method)

  # -- Load tables once (from package sysdata or .TM_CONSTANTS) -------------
  # In the final package, replace with: tbl <- .TM_CONSTANTS$NN[[nn_table]]
  # For now, build once in this call (still 100x faster than per-sequence):
  nn_tbl <- get_table(nn_table)   # internal helper (see below)
  tmm_tbl <- get_table(tmm_table)
  imm_tbl <- get_table(imm_table)
  de_tbl <- get_table(de_table)
  # Companion end-effect table, if the selected parameter set ships one
  # (currently only Zuber 2022). Empty matrix otherwise, which leaves the
  # calculation identical to previous releases.
  end_tbl <- get_end_table(nn_table)

  # -- Salt-correction guard -------------------------------------------------
  # Some parameter sets (the Weber/VarGibbs series) are fitted AT a specific
  # sodium concentration and are intended to replace salt correction rather
  # than be corrected. Applying a correction on top of them double-counts the
  # ionic effect and yields silently wrong Tm values. Sets that carry the
  # attribute "salt_mM" are handled here.
  salt_fn_eff <- if (identical(salt_method, "none")) NULL else salt_method
  tbl_salt <- attr(nn_tbl, "salt_mM")
  if (!is.null(tbl_salt) && !is.na(tbl_salt) && !is.null(salt_fn_eff)) {
    if (isTRUE(all.equal(as.numeric(tbl_salt), as.numeric(Na),
                         tolerance = 0.02, scale = 1))) {
      # Already fitted at this salt: skip correction silently.
      salt_fn_eff <- NULL
    } else {
      warning("Parameter set '", nn_table, "' was fitted at ", tbl_salt,
              " mM [Na+], but Na = ", Na, " mM was requested. The '",
              salt_method, "' correction is being applied on top of a ",
              "condition-specific parameter set, which is approximate. ",
              "Prefer the set fitted nearest your condition, or set ",
              "salt_method = \"none\".", call. = FALSE)
    }
  }

  # -- Pairwise N filtering (fast path; the per-base cleaning that used to
  # -- happen here now runs inside the C++ core, on the workers) -------------
  has_n <- .col_has_n(gr_seq$sequence) | .col_has_n(gr_seq$complement)
  if (any(has_n)) {
    warning(
      sprintf(
        "Skipped %d region(s) because sequence or complement contains 'N'. See your_output$options$'Skipped regions containing N' for details",
        sum(has_n)
      ),
      call. = FALSE
    )
  }
  gr_seq_dropoff <- gr_seq[has_n]
  gr_seq <- gr_seq[!has_n]
  if (length(gr_seq) == 0) {
    stop("No valid regions left for tm_nn calculation after filtering sequences with 'N'.")
  }

  # -- Process all sequences (chunked, optionally in parallel) ---------------
  n   <- length(gr_seq)

  # Keep the raw columns (possibly DNAStringSet); slices are coerced to
  # character per chunk so PSOCK workers never receive an XVector whose
  # serialization would drag the whole shared pool along.
  all_seqs  <- mcols(gr_seq)$sequence
  all_cseqs <- mcols(gr_seq)$complement

  chunk_res <- .bp_map_chunks(
    n = n,
    make_chunk = function(idx) list(sequence = as.character(all_seqs[idx]),
                                    complement = as.character(all_cseqs[idx])),
    worker = .tm_nn_chunk,
    BPPARAM = BPPARAM,
    ambiguous = ambiguous, shift = shift,
    nn_tbl = nn_tbl, tmm_tbl = tmm_tbl, imm_tbl = imm_tbl, de_tbl = de_tbl,
    end_tbl = end_tbl,
    dnac_high = dnac_high, dnac_low = dnac_low, self_comp = self_comp,
    Na = Na, K = K, Tris = Tris, Mg = Mg, dNTPs = dNTPs,
    salt_fn = salt_fn_eff, DMSO = DMSO, dmso_factor = dmso_factor,
    formamide_factor = formamide_factor, formamide_unit = formamide_unit
  )
  tm <- chunk_res$Tm
  gc <- chunk_res$GC
  if (!"GC" %in% names(GenomicRanges::mcols(gr_seq))) {
    gr_seq$GC <- gc
  }
  gr_seq$Tm <- tm
  
  nn_table_list <- list("DNA_NN_Breslauer_1986" = "Breslauer K J (1986) <doi:10.1073/pnas.83.11.3746>",
                        "DNA_NN_Sugimoto_1996" = "Sugimoto N (1996) <doi:10.1093/nar/24.22.4501>",
                        "DNA_NN_Allawi_1998" = "Allawi H (1998) <doi:10.1093/nar/26.11.2694>",
                        "DNA_NN_SantaLucia_2004" = "SantaLucia J (2004) <doi:10.1146/annurev.biophys.32.110601.141800>",
                        "RNA_NN_Freier_1986" = "Freier S (1986) <doi:10.1073/pnas.83.24.9373>",
                        "RNA_NN_Xia_1998" = "Xia T (1998) <doi:10.1021/bi9809425>",
                        "RNA_NN_Chen_2012" = "Chen JL (2012) <doi:10.1021/bi3002709>",
                        "RNA_DNA_NN_Sugimoto_1995" = "Sugimoto N (1995)<doi:10.1016/S0048-9697(98)00088-6>",
                        "DNA_TMM_Bommarito_2000" = "Bommarito S (2000)  <doi:10.1093/nar/28.9.1929>",
                        "DNA_IMM_Peyret_1999" = "Peyret N (1999) <doi:10.1021/bi9825091> & Allawi H T (1997) <doi:10.1021/bi962590c> & Santalucia N (2005) <doi:10.1093/nar/gki918>",
                        "DNA_DE_Bommarito_2000" = "Bommarito S (2000) <doi:10.1093/nar/28.9.1929>",
                        "RNA_DE_Turner_2010" = "Turner D H (2010) <doi:10.1093/nar/gkp892>",
                        "DNA_NN_Weber_2015" = "Weber G (2015) <doi:10.1093/bioinformatics/btu751>",
                        "DNA_NN_Weber_OW04_69" = "Weber G (2015) <doi:10.1093/bioinformatics/btu751>, 69 mM [Na+]",
                        "DNA_NN_Weber_OW04_119" = "Weber G (2015) <doi:10.1093/bioinformatics/btu751>, 119 mM [Na+]",
                        "DNA_NN_Weber_OW04_220" = "Weber G (2015) <doi:10.1093/bioinformatics/btu751>, 220 mM [Na+]",
                        "DNA_NN_Weber_OW04_621" = "Weber G (2015) <doi:10.1093/bioinformatics/btu751>, 621 mM [Na+]",
                        "DNA_NN_Weber_OW04_1020" = "Weber G (2015) <doi:10.1093/bioinformatics/btu751>, 1020 mM [Na+]",
                        "RNA_NN_Weber_VIF_71" = "Ferreira I (2019) <doi:10.1016/j.chemphys.2019.01.016>, VIF, 71 mM [Na+]",
                        "RNA_NN_Weber_VIF_121" = "Ferreira I (2019) <doi:10.1016/j.chemphys.2019.01.016>, VIF, 121 mM [Na+]",
                        "RNA_NN_Weber_VIF_221" = "Ferreira I (2019) <doi:10.1016/j.chemphys.2019.01.016>, VIF, 221 mM [Na+]",
                        "RNA_NN_Weber_VIF_621" = "Ferreira I (2019) <doi:10.1016/j.chemphys.2019.01.016>, VIF, 621 mM [Na+]",
                        "RNA_NN_Weber_VIF_1021" = "Ferreira I (2019) <doi:10.1016/j.chemphys.2019.01.016>, VIF, 1021 mM [Na+]",
                        "RNA_NN_Weber_FIF_71" = "Ferreira I (2019) <doi:10.1016/j.chemphys.2019.01.016>, FIF, 71 mM [Na+]",
                        "RNA_NN_Weber_FIF_121" = "Ferreira I (2019) <doi:10.1016/j.chemphys.2019.01.016>, FIF, 121 mM [Na+]",
                        "RNA_NN_Weber_FIF_221" = "Ferreira I (2019) <doi:10.1016/j.chemphys.2019.01.016>, FIF, 221 mM [Na+]",
                        "RNA_NN_Weber_FIF_621" = "Ferreira I (2019) <doi:10.1016/j.chemphys.2019.01.016>, FIF, 621 mM [Na+]",
                        "RNA_NN_Weber_FIF_1021" = "Ferreira I (2019) <doi:10.1016/j.chemphys.2019.01.016>, FIF, 1021 mM [Na+]",
                        "RNA_DNA_NN_Weber_2019_FT" = "Basilio Barbosa V (2019) <doi:10.1016/j.bpc.2019.106189>, curve fitting, 1000 mM [Na+]",
                        "RNA_DNA_NN_Weber_2019_VH" = "Basilio Barbosa V (2019) <doi:10.1016/j.bpc.2019.106189>, van't Hoff, 1000 mM [Na+]",
                        "RNA_DNA_NN_Weber_2019_LS" = "Basilio Barbosa V (2019) <doi:10.1016/j.bpc.2019.106189>, low salt, 100 mM [Na+]",
                        "RNA_DNA_NN_Banerjee_2020" = "Banerjee D (2020) <doi:10.1093/nar/gkaa572>, physiological condition, 100 mM [Na+]",
                        "RNA_NN_Zuber_2022" = "Zuber J (2022) <doi:10.1093/nar/gkac261>, improved end effects",
                        "RNA_NN_Ghosh_2023_PEG200" = "Ghosh S (2023) <doi:10.1093/nar/gkad020>, 40 wt% PEG200 crowding, 100 mM [Na+]",
                        "DNA_NN_Ghosh_2020_PEG200" = "Ghosh S (2020) <doi:10.1073/pnas.1920886117>, 40 wt% PEG200 crowding, 100 mM [Na+]")
  
  # Create result list with proper structure
  
  gr_seq <- .normalize_tm_gc_metadata(gr_seq)

  result_list <- list(
    gr = gr_seq,
    options = list("Ambiguous" = ambiguous,
                   "Shift" = shift,
                   "Thermodynamic NN values" = paste0(nn_table, ": ", nn_table_list[[nn_table]]), 
                   "Thermodynamic values for terminal mismatches" = paste0(tmm_table,": ",nn_table_list[[tmm_table]]), 
                   "Thermodynamic values for internal mismatches" = paste0(imm_table,": ",nn_table_list[[imm_table]]),
                   "Thermodynamic values for dangling ends" = paste0(de_table,": ",nn_table_list[[de_table]]), 
                   "Concentration of the higher concentrated strand" = dnac_high,
                   "Concentration of the lower concentrated strand" = dnac_low, 
                   "Sequence self-complementary" = self_comp, 
                   "Na" = Na,
                   "K" = K,
                   "Tris" = Tris,
                   "Mg" = Mg,
                   "dNTPs" = dNTPs,
                   "Salt correction method" = salt_method,
                   "Salt correction applied" = !is.null(salt_fn_eff),
                   "Parameter set fitted at [Na+] (mM)" =
                     if (is.null(attr(nn_tbl, "salt_mM"))) NA_real_
                     else as.numeric(attr(nn_tbl, "salt_mM")),
                   "Percent of DMSO" = DMSO,
                   "Formamide concentration" = formamide_unit$value,
                   "DMSO factor" = dmso_factor,
                   "Formamide concentration unit" = formamide_unit$unit,
                   "Formamide factor" = formamide_factor,
                   "Skipped regions containing N" = gr_seq_dropoff)
  )
  
  class(result_list) <- c("TmCalculator", "list")
  attr(result_list, "nonhidden") <- "gr"
  return(result_list)
}

# -- Fast per-column N detection ----------------------------------------------
# vcountPattern works directly on XStringSet without character coercion.
#' @keywords internal
.col_has_n <- function(x) {
  if (inherits(x, "XStringSet")) {
    Biostrings::vcountPattern("N", x) > 0
  } else {
    grepl("N", toupper(as.character(x)), fixed = FALSE)
  }
}

# -- Table -> C++ handoff -----------------------------------------------------
#' @keywords internal
.tbl_to_cpp <- function(tbl) {
  keys <- rownames(tbl)
  # A zero-row table (e.g. the end table of a set without end effects) has
  # NULL rownames; hand C++ a zero-length character vector instead of NULL.
  if (is.null(keys)) keys <- character(0)
  list(keys = keys,
       dh = as.numeric(tbl[, 1]),
       ds = as.numeric(tbl[, 2]))
}

# -- Chunk worker: NN Tm/GC over a block of sequences (Rcpp-backed) -----------
# Called by .bp_map_chunks(), either directly (serial) or on a BiocParallel
# worker. `chunk` is list(sequence=, complement=) for this worker's block.
#
# The C++ core (src/tm_nn_core.cpp) uppercases each sequence, strips
# characters outside A/C/G/T/I (the former check_filter_seq step, now on the
# workers), accumulates delta_H/delta_S and base counts; the two-state Tm
# formula plus salt and chemical corrections are applied here, vectorized.
# After cleaning, sequences contain only A/C/G/T/I, so the `ambiguous` flag
# cannot change GC values on this path.
#' @keywords internal
.tm_nn_chunk <- function(chunk, ambiguous, shift, nn_tbl, tmm_tbl, imm_tbl,
                         de_tbl, end_tbl, dnac_high, dnac_low, self_comp,
                         Na, K, Tris, Mg, dNTPs, salt_fn,
                         DMSO, dmso_factor, formamide_factor, formamide_unit) {
  self_comp_eff <- isTRUE(self_comp) && ("sym" %in% rownames(nn_tbl))

  res <- cpp_tm_nn_dhds(
    as.character(chunk$sequence), as.character(chunk$complement),
    as.integer(shift),
    .tbl_to_cpp(nn_tbl), .tbl_to_cpp(tmm_tbl),
    .tbl_to_cpp(imm_tbl), .tbl_to_cpp(de_tbl),
    self_comp_eff, .tbl_to_cpp(end_tbl)
  )

  dh  <- res[, "dh"]
  ds  <- res[, "ds"]
  nGC <- res[, "nG"] + res[, "nC"]
  acgt <- res[, "nA"] + res[, "nC"] + res[, "nG"] + res[, "nT"]
  len <- res[, "len"]
  ok  <- res[, "ok"] > 0

  # One definition of GC throughout the package: (G+C)/(A+C+G+T), i.e. gc()
  # semantics. Previously this path reported (G+C)/length and salt-corrected
  # with (G+C)/(A+C+G+T); the two diverge whenever inosine is present, since
  # I counts in the length but is not a determinable base.
  gc_salt <- ifelse(acgt > 0, 100 * nGC / acgt, NA_real_)

  k <- if (self_comp_eff) dnac_high * 1e-9 else
    (dnac_high - (dnac_low / 2.0)) * 1e-9
  R <- 1.987

  corr_salt <- NULL
  if (!is.null(salt_fn) && !identical(salt_fn, "none")) {
    corr_salt <- .salt_correct_vec(Na = Na, K = K, Tris = Tris, Mg = Mg,
                                   dNTPs = dNTPs, method = salt_fn,
                                   gc_pct = gc_salt, seq_len = len)
    if (identical(salt_fn, "SantaLucia1998-2")) {
      ds <- ds + corr_salt
    }
    tm <- (1000 * dh) / (ds + (R * log(k))) - 273.15
    if (salt_fn %in% c("Schildkraut2010", "Wetmur1991",
                       "SantaLucia1996", "SantaLucia1998-1")) {
      tm <- tm + corr_salt
    }
    if (salt_fn %in% c("Owczarzy2004", "Owczarzy2008")) {
      tm <- (1 / (1 / (tm + 273.15) + corr_salt) - 273.15)
    }
  } else {
    tm <- (1000 * dh) / (ds + (R * log(k))) - 273.15
  }

  tm <- tm + .chem_correct_vec(DMSO = DMSO,
                               formamide_unit = formamide_unit,
                               dmso_factor = dmso_factor,
                               formamide_factor = formamide_factor,
                               pt_gc = gc_salt)

  # Reproduce the R core's failure semantics: any sequence whose evaluation
  # errored (short/missing init rows) or whose salt correction is undefined
  # (e.g. Owczarzy2008 with sqrt(Mg)/mon >= 6) gets NA for BOTH Tm and GC,
  # exactly as tryCatch() around .tm_nn_core() did.
  bad <- !ok
  if (!is.null(corr_salt)) {
    bad <- bad | is.na(corr_salt)
  }
  tm[bad] <- NA_real_
  gc_out <- gc_salt
  gc_out[bad] <- NA_real_
  list(Tm = unname(tm), GC = unname(gc_out))
}

# -- Pure-R chunk worker, kept as reference implementation --------------------
# Used by unit tests to verify the Rcpp path reproduces the original R
# results exactly; not called in normal operation.
#' @keywords internal
.tm_nn_chunk_r <- function(chunk, ambiguous, shift, nn_tbl, tmm_tbl, imm_tbl,
                           de_tbl, end_tbl, dnac_high, dnac_low, self_comp,
                           Na, K, Tris, Mg, dNTPs, salt_fn,
                           DMSO, dmso_factor, formamide_factor, formamide_unit) {
  m  <- length(chunk$sequence)
  tm <- rep(NA_real_, m)
  gc <- rep(NA_real_, m)
  for (j in seq_len(m)) {
    seq_str  <- chunk$sequence[j]
    cseq_str <- chunk$complement[j]
    if (is.na(seq_str) || is.na(cseq_str)) {
      next
    }
    # Same cleaning as the C++ core: uppercase, keep only A/C/G/T/I
    seq_str  <- gsub("[^ACGTI]", "", toupper(seq_str), perl = TRUE)
    cseq_str <- gsub("[^ACGTI]", "", toupper(cseq_str), perl = TRUE)
    if (nchar(seq_str) < 2L) {
      next
    }
    result <- tryCatch(
      .tm_nn_core(seq_str, cseq_str, ambiguous, shift, nn_tbl = nn_tbl,
                  tmm_tbl = tmm_tbl, imm_tbl = imm_tbl, de_tbl = de_tbl,
                  end_tbl = end_tbl,
                  dnac_high, dnac_low, self_comp,
                  Na, K, Tris, Mg, dNTPs, salt_fn = salt_fn,
                  DMSO, dmso_factor, formamide_factor, formamide_unit),
      error = function(e) list(Tm = NA_real_, GC = NA_real_)
    )
    tm[j] <- result$Tm
    gc[j] <- result$GC
  }
  list(Tm = unname(tm), GC = unname(gc))
}

# -- Core single-sequence NN computation --------------------------------------
# All table lookups are vectorized (no per-base loop)

.tm_nn_core <- function(seq_str, cseq_str, ambiguous, shift, nn_tbl, tmm_tbl, imm_tbl, de_tbl, dnac_high, dnac_low, self_comp,
                        end_tbl = matrix(numeric(0), nrow = 0, ncol = 2),
                        Na, K, Tris, Mg, dNTPs, salt_fn, DMSO, dmso_factor, formamide_factor, formamide_unit) {


  tmp_seq <- seq_str
  tmp_cseq <- cseq_str
  delta_h <- 0
  delta_s <- 0
  
  # for de end
  if(shift != 0 || nchar(tmp_seq) != nchar(tmp_cseq)) {
    if(shift > 0) {
      tmp_seq <- paste0(paste(rep('.', shift), collapse = ""), tmp_seq)
    }else{
      tmp_cseq <- paste0(paste(rep('.', abs(shift)), collapse = ""), tmp_cseq)
    }
    if(nchar(tmp_cseq)>nchar(tmp_seq)){
      tmp_seq <- paste0(tmp_seq, paste(rep('.', nchar(tmp_cseq) - nchar(tmp_seq)), collapse = ""))
    }
    if(nchar(tmp_cseq)<nchar(tmp_seq)){
      tmp_cseq <- paste0(tmp_cseq, paste(rep('.', nchar(tmp_seq) - nchar(tmp_cseq)), collapse = ""))
    }
    
    while(substring(tmp_seq, 1, 2) == ".." || substring(tmp_cseq, 1, 2) == "..") {
      tmp_seq <- substring(tmp_seq, 2, nchar(tmp_seq))
      tmp_cseq <- substring(tmp_cseq, 2, nchar(tmp_cseq))
    }
    
    while(substring(tmp_seq, nchar(tmp_seq)-1, nchar(tmp_seq)) == ".." || substring(tmp_cseq, nchar(tmp_cseq)-1, nchar(tmp_cseq)) == "..") {
      tmp_seq <- substring(tmp_seq, 1, nchar(tmp_seq)-1)
      tmp_cseq <- substring(tmp_cseq, 1, nchar(tmp_cseq)-1)
    }
  }
  
  # -- Build all dinucleotide keys at once ----------------------------------
  # substring() is faster than strsplit -> paste for large n
  n     <- nchar(tmp_seq)
  n_int <- n - 1L
  
  fwd  <- substring(tmp_seq, 1:n_int, 2:(n_int+1))          # forward strand
  #bwd <- substring(paste(rev(strsplit(tmp_seq, "", fixed=TRUE)[[1]]), collapse=""), 1:n_int, 2:(n_int+1))
  cfwd <- substring(tmp_cseq, 1:n_int, 2:(n_int+1))
  #cbwd <- substring(paste(rev(strsplit(tmp_cseq, "", fixed=TRUE)[[1]]), collapse=""), 1:n_int, 2:(n_int+1))

  keys_fr <- paste0(fwd, "/", cfwd)
  #keys_rf <- paste0(cbwd, "/", bwd)
  

  keys_t_left <- keys_fr[1]
  keys_t_right <- .right_key(tmp_seq, tmp_cseq, n)
  
  #for dang end
  if(keys_t_left %in% rownames(de_tbl)) {
    delta_h <- de_tbl[keys_t_left,1] + delta_h
    delta_s <- de_tbl[keys_t_left,2] + delta_s
    keys_fr <- keys_fr[-1]
    keys_t_left <- keys_fr[1]
    tmp_seq  <- substring(tmp_seq,  2, n)
    tmp_cseq <- substring(tmp_cseq, 2, n)
  }
  
  if (keys_t_right %in% rownames(de_tbl)) {
    delta_h <- de_tbl[keys_t_right, 1] + delta_h
    delta_s <- de_tbl[keys_t_right, 2] + delta_s
    keys_fr <- keys_fr[-length(keys_fr)]
    n <- nchar(tmp_seq) - 1L
    tmp_seq <- substring(tmp_seq, 1, n)
    tmp_cseq <- substring(tmp_cseq, 1, n)
    keys_t_right <- .right_key(tmp_seq, tmp_cseq, n)
  }
  
  # for terminal mismatch
  if(keys_t_left %in% rownames(tmm_tbl)) {
    delta_h <- tmm_tbl[keys_t_left, 1] + delta_h
    delta_s <- tmm_tbl[keys_t_left, 2] + delta_s
    keys_fr <- keys_fr[-1]
    n <- nchar(tmp_seq)
    tmp_seq  <- substring(tmp_seq,  2, n)
    tmp_cseq <- substring(tmp_cseq, 2, n)
  }
  
  if(keys_t_right %in% rownames(tmm_tbl)) {
    delta_h <- tmm_tbl[keys_t_right, 1] + delta_h
    delta_s <- tmm_tbl[keys_t_right, 2] + delta_s
    #keys_rf <- keys_rf[-1]
    keys_fr <- keys_fr[-length(keys_fr)]
    n <- nchar(tmp_seq)-1
    tmp_seq <- substring(tmp_seq, 1, n)
    tmp_cseq <- substring(tmp_cseq, 1, n)
  }
  
  # for end effects (Zuber 2022): added, terminal pair NOT consumed
  if (nrow(end_tbl) > 0L) {
    if (length(keys_fr) > 0L && keys_fr[1] %in% rownames(end_tbl)) {
      delta_h <- end_tbl[keys_fr[1], 1] + delta_h
      delta_s <- end_tbl[keys_fr[1], 2] + delta_s
    }
    key_e_right <- .right_key(tmp_seq, tmp_cseq, nchar(tmp_seq))
    if (key_e_right %in% rownames(end_tbl)) {
      delta_h <- end_tbl[key_e_right, 1] + delta_h
      delta_s <- end_tbl[key_e_right, 2] + delta_s
    }
  }

  # for initial of nearest neighbor
  delta_h <- nn_tbl['init', 1] + delta_h
  delta_s <- nn_tbl['init', 2] + delta_s
  
  if(substring(tmp_seq, 1, 1) == 'T'){
    delta_h <- nn_tbl['init_5T/A', 1] + delta_h
    delta_s <- nn_tbl['init_5T/A', 2] + delta_s
  }
  #if(substring(tmp_seq, 1, 1) == 'A'){
  #  delta_h <- nn_tbl['init_5T/A', 1] + delta_h
  #  delta_s <- nn_tbl['init_5T/A', 2] + delta_s
  #}
  
  # -- Initiation parameters -------------------------------------------------
  first_base <- substr(tmp_seq, 1, 1)
  last_base  <- substr(tmp_seq, nchar(tmp_seq), nchar(tmp_seq))
  gc_ends    <- sum(c(first_base, last_base) %in% c("G","C"))
  at_ends    <- 2L - gc_ends
  
  if(gc_ends == 0){
    delta_h <- nn_tbl['init_allA/T', 1] + delta_h
    delta_s <- nn_tbl['init_allA/T', 2] + delta_s
  }else{
    delta_h <- nn_tbl['init_oneG/C', 1] + delta_h
    delta_s <- nn_tbl['init_oneG/C', 2] + delta_s
  }
  
  delta_h <- nn_tbl['init_A/T', 1] * at_ends + delta_h
  delta_s <- nn_tbl['init_A/T', 2] * at_ends + delta_s
  delta_h <- nn_tbl['init_G/C', 1] * gc_ends + delta_h
  delta_s <- nn_tbl['init_G/C', 2] * gc_ends + delta_s
  
  # -- Vectorized table lookup -----------------------------------------------
  # for nn table
  matched_nn_fr <- keys_fr %in% rownames(nn_tbl)
  #matched_nn_rf <- keys_rf %in% rownames(nn_tbl)
  #which_nn_fr <- which(matched_nn_fr)
  #which_nn_rf <- which(matched_nn_rf)
  #n_int <- length(keys_fr)
  #pos_nn_rf <- which_nn_rf[!which_nn_rf %in% (n_int - which_nn_fr + 1)]

  delta_h <- sum(nn_tbl[keys_fr[matched_nn_fr], 1]) + delta_h
  delta_s <- sum(nn_tbl[keys_fr[matched_nn_fr], 2]) + delta_s
  #delta_h <- sum(nn_tbl[keys_rf[pos_nn_rf], 1]) + delta_h
  #delta_s <- sum(nn_tbl[keys_rf[pos_nn_rf], 2]) + delta_s
  
  # for imm table
  matched_imm_fr <- keys_fr %in% rownames(imm_tbl)
  #matched_imm_rf <- keys_rf %in% rownames(imm_tbl)
  
  #if(any(c(matched_imm_fr,matched_imm_rf))){
  if(any(c(matched_imm_fr))){
    #which_imm_fr <- which(matched_imm_fr)
    #which_imm_rf <- which(matched_imm_rf)
    #pos_imm_rf <- which_imm_rf[!which_imm_rf %in% (n_int - which_imm_fr + 1)]
    
    delta_h <- sum(imm_tbl[keys_fr[matched_imm_fr], 1]) + delta_h
    delta_s <- sum(imm_tbl[keys_fr[matched_imm_fr], 2]) + delta_s
    
    #delta_h <- sum(imm_tbl[keys_rf[pos_imm_rf], 1]) + delta_h
    #delta_s <- sum(imm_tbl[keys_rf[pos_imm_rf], 2]) + delta_s
  }
  
  # -- Symmetry correction ---------------------------------------------------
  k <- (dnac_high-(dnac_low/2.0))*1e-9
  
  if (self_comp && "sym" %in% rownames(nn_tbl)) {
    k <- dnac_high*1e-9
    delta_h <- delta_h + nn_tbl["sym", 1]
    delta_s <- delta_s + nn_tbl["sym", 2]
  }

  R <- 1.987
  if(!is.null(salt_fn) && !identical(salt_fn, "none")){
    corr_salt <- salt_correct(Na = Na, 
                                 K = K,
                                 Tris = Tris,
                                 Mg = Mg,
                                 dNTPs = dNTPs,
                                 method = salt_fn,
                                 input_seq = seq_str,
                                 ambiguous = ambiguous)
    if(salt_fn == "SantaLucia1998-2"){
      delta_s <- corr_salt+delta_s
    }
    tm <- (1000 * delta_h) / (delta_s + (R * (log(k)))) - 273.15
    if (salt_fn %in% c("Schildkraut2010", "Wetmur1991",
                           "SantaLucia1996", "SantaLucia1998-1")) {
      tm <- tm + corr_salt
    }
    if(salt_fn %in% c("Owczarzy2004","Owczarzy2008")){
      tm <- (1 / (1 / (tm + 273.15) + corr_salt) - 273.15)
    }
  } else {
    tm <- (1000 * delta_h) / (delta_s + (R * (log(k)))) - 273.15
  }
  pt_gc <- .gc_vec(seq_str, ambiguous = ambiguous)
  corr_chem <- chem_correct(
    DMSO = DMSO,
    formamide_unit = formamide_unit,
    dmso_factor = dmso_factor,
    formamide_factor = formamide_factor,
    pt_gc = pt_gc
  )
  tm <- tm + corr_chem
  return(list(Tm = tm, GC = pt_gc))
}

# Package-private cache for thermodynamic tables (not .GlobalEnv)
.TM_TABLE_CACHE <- new.env(parent = emptyenv())

# -- Helper: companion end-effect table for a parameter set ------------------
# Returns the "<nn_table>_END" entry when the set defines penultimate-pair
# dependent end terms (Zuber 2022), otherwise a zero-row matrix.
#' @keywords internal
get_end_table <- function(table_name) {
  end_name <- paste0(table_name, "_END")
  if (!is.null(.TM_CONSTANTS[[end_name]])) {
    return(get_table(end_name))
  }
  matrix(numeric(0), nrow = 0, ncol = 2,
         dimnames = list(NULL, c("left", "right")))
}

# -- Helper: get table (package data or build) ------------------------------
get_table <- function(table_name) {
  if (!exists(table_name, envir = .TM_TABLE_CACHE, inherits = FALSE)) {
    assign(table_name, .TM_CONSTANTS[[table_name]], envir = .TM_TABLE_CACHE)
  }
  get(table_name, envir = .TM_TABLE_CACHE, inherits = FALSE)
}


.rev_str <- function(s) paste(rev(strsplit(s, "", fixed=TRUE)[[1]]), collapse="")

# -----------------------------------------------------------------------------
# FIX 12: Fast complement
# Drop-in replacement for complement()
# -----------------------------------------------------------------------------

#' Fast complement and reverse complement
#'
#' Uses chartr() instead of character-by-character translation.
#' ~10x faster than the seqinr-style implementation for long sequences.
#'
#' @param seq_str Character string or vector.
#' @param rev Logical. Return reverse complement? Default FALSE.
#' @return Character string(s) with complement.
#' @export
complement_fast <- function(seq_str, rev = FALSE) {
  # chartr handles both upper and lower case simultaneously
  comp <- chartr("ACGTacgtWSRYKMBVDH", "TGCAtgcaWSYRMKVBHD", seq_str)
  if (rev) {
    vapply(comp, function(s) {
      paste(rev(strsplit(s, "", fixed = TRUE)[[1]]), collapse = "")
    }, character(1L), USE.NAMES = FALSE)
  } else {
    comp
  }
}


.right_key <- function(seq, cseq, len) {
  paste0(
    substr(cseq, len,    len),       # cseq last base
    substr(cseq, len-1L, len-1L),   # cseq second-to-last
    "/",
    substr(seq,  len,    len),       # seq last base
    substr(seq,  len-1L, len-1L)    # seq second-to-last
  )
}



# -----------------------------------------------------------------------------
# QUICK BENCHMARK to verify gains
# -----------------------------------------------------------------------------

#' @keywords internal
benchmark_tm_nn <- function(n_seqs = 1000L, seq_len = 200L) {
  set.seed(42)
  bases <- c("A", "C", "G", "T")
  seqs <- vapply(seq_len(n_seqs), function(i) {
    paste(sample(bases, seq_len, replace = TRUE), collapse = "")
  }, character(1))

  cat(sprintf(
    "Benchmarking per-row tm_nn() vs one tm_nn() on a multi-row GRanges (%d seq x %d bp)\n\n",
    n_seqs, seq_len
  ))

  t_loop <- system.time({
    tm_loop <- vapply(seq_along(seqs), function(i) {
      gr_i <- to_genomic_ranges(seqs[i])
      out <- tm_nn(gr_i, Na = 50, salt_method = "Owczarzy2004")
      as.numeric(S4Vectors::mcols(out$gr)$Tm[1])
    }, numeric(1))
  })

  gr_all <- to_genomic_ranges(seqs)
  t_batch <- system.time({
    out_b <- tm_nn(gr_all, Na = 50, salt_method = "Owczarzy2004")
    tm_batch <- as.numeric(S4Vectors::mcols(out_b$gr)$Tm)
  })

  cat(sprintf("  Per-sequence (loop):       %.2f sec\n", t_loop["elapsed"]))
  cat(sprintf("  One tm_nn() on GRanges:    %.2f sec\n", t_batch["elapsed"]))
  if (t_batch["elapsed"] > 0) {
    cat(sprintf(
      "  Speedup:                   %.1fx\n\n",
      unname(t_loop["elapsed"] / t_batch["elapsed"])
    ))
  }

  valid <- !is.na(tm_loop) & !is.na(tm_batch)
  if (sum(valid) > 0) {
    max_diff <- max(abs(tm_loop[valid] - tm_batch[valid]))
    cat(sprintf(
      "  Max Tm difference: %.6f deg C %s\n",
      max_diff,
      if (max_diff < 0.01) "(PASS)" else "(CHECK!)"
    ))
  }

  invisible(list(
    loop = tm_loop, batch = tm_batch,
    speedup = if (t_batch["elapsed"] > 0) {
      unname(t_loop["elapsed"] / t_batch["elapsed"])
    } else {
      NA_real_
    }
  ))
}
