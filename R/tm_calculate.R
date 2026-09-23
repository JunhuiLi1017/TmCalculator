#' Calculate melting temperature using multiple methods
#'
#' Calculates nucleic acid melting temperature (Tm) by one of three methods, and
#' returns the result as a \code{GRanges} object so that Tm can be used directly
#' as a quantitative genomic feature alongside other assays:
#' \itemize{
#'   \item \strong{Nearest neighbor} (\code{tm_nn}) sums the stacking free
#'     energies of adjacent base-pair steps using experimentally derived
#'     enthalpy and entropy parameters, and so resolves sequences of identical
#'     base composition but different order. It supports salt and chemical
#'     corrections, mismatches and dangling ends, and is the default.
#'   \item \strong{GC content} (\code{tm_gc}) computes Tm as an empirical
#'     function of GC percentage with corrections for length and ionic strength.
#'     Cheaper than NN, but blind to sequence order.
#'   \item \strong{Wallace rule} (\code{tm_wallace}) assigns a fixed
#'     contribution per base. It is calibrated for short oligonucleotides,
#'     typically 14 to 20 bp, and ignores sequence context and reaction
#'     conditions, so it is not appropriate for long sequences or for
#'     genome-wide windows.
#' }
#'
#' The input sequence is processed once and passed to the selected method, which
#' is faster than calling the individual functions separately.
#'
#' @section Salt handling:
#' Most nearest-neighbor parameter sets were fitted at a single reference sodium
#' concentration, and other conditions are reached through the
#' \code{salt_method} correction formulas. The Weber/VarGibbs sets were instead
#' fitted directly at a stated sodium concentration and are meant to replace
#' salt correction. When such a set is selected and \code{Na} matches the
#' concentration it was fitted at, correction is skipped automatically; when it
#' does not, correction is applied with a warning. See \code{\link{tm_nn}} for
#' details.
#' 
#' @section Available Options:
#' 
#' \strong{Method Selection:}
#' \itemize{
#'   \item \code{method}: c("tm_nn", "tm_gc", "tm_wallace")
#' }
#' 
#' \strong{Nearest Neighbor (NN) Method Options:}
#' \itemize{
#'   \item \code{nn_table}:
#'     \itemize{
#'       \item DNA/DNA: "DNA_NN_Breslauer_1986", "DNA_NN_Sugimoto_1996",
#'         "DNA_NN_Allawi_1998", "DNA_NN_SantaLucia_2004" (default)
#'       \item DNA/DNA, molecular crowding: "DNA_NN_Ghosh_2020_PEG200"
#'         (40 wt% PEG200, 100 mM NaCl)
#'       \item DNA/DNA, salt-optimized: "DNA_NN_Weber_2015" (1020 mM),
#'         "DNA_NN_Weber_OW04_69", "..._119", "..._220", "..._621",
#'         "..._1020" (fitted at 69 to 1020 mM sodium)
#'       \item RNA/RNA: "RNA_NN_Freier_1986", "RNA_NN_Xia_1998",
#'         "RNA_NN_Chen_2012", "RNA_NN_Zuber_2022" (improved end effects),
#'         "RNA_NN_Ghosh_2023_PEG200" (molecular crowding, cell-like)
#'       \item RNA/RNA, salt-optimized: "RNA_NN_Weber_VIF_71", "..._121",
#'         "..._221", "..._621", "..._1021" and the corresponding
#'         "RNA_NN_Weber_FIF_*" sets
#'       \item RNA/DNA: "RNA_DNA_NN_Sugimoto_1995",
#'         "RNA_DNA_NN_Weber_2019_FT", "RNA_DNA_NN_Weber_2019_VH" (1000 mM),
#'         "RNA_DNA_NN_Weber_2019_LS" (100 mM),
#'         "RNA_DNA_NN_Banerjee_2020" (100 mM)
#'     }
#'   \item \code{tmm_table} (Terminal Mismatches):
#'     \itemize{
#'       \item "DNA_TMM_Bommarito_2000" (default)
#'     }
#'   \item \code{imm_table} (Internal Mismatches):
#'     \itemize{
#'       \item "DNA_IMM_Peyret_1999" (default)
#'     }
#'   \item \code{de_table} (Dangling Ends):
#'     \itemize{
#'       \item "DNA_DE_Bommarito_2000" (default)
#'       \item "RNA_DE_Turner_2010"
#'     }
#' }
#' 
#' \strong{GC Method Options:}
#' \itemize{
#'   \item \code{variant}:
#'     \itemize{
#'       \item "Primer3Plus" (default)
#'       \item "Chester1993"
#'       \item "QuikChange"
#'       \item "Schildkraut1965"
#'       \item "Wetmur1991_MELTING"
#'       \item "Wetmur1991_RNA"
#'       \item "Wetmur1991_RNA/DNA"
#'       \item "vonAhsen2001"
#'     }
#' }
#' 
#' \strong{Salt Correction Options:}
#' \itemize{
#'   \item \code{salt_method}:
#'     \itemize{
#'       \item "Schildkraut2010" (default)
#'       \item "Wetmur1991"
#'       \item "SantaLucia1996"
#'       \item "SantaLucia1998-1"
#'       \item "Owczarzy2004" (\code{method = "tm_nn"} only)
#'       \item "Owczarzy2008" (\code{method = "tm_nn"} only)
#'       \item "none" (also selected automatically when \code{nn_table} was
#'         fitted at the requested \code{Na})
#'     }
#' }
#'
#' With \code{method = "tm_gc"} the salt term belongs to the published formula
#' selected by \code{variant}, so naming a different one is ignored there, with
#' a warning, unless \code{userset} is supplied. \code{"none"} and \code{NA}
#' drop the correction and are honoured on either path.
#' 
#' \strong{Formamide Unit Options:}
#' \itemize{
#'   \item \code{formamide_unit$unit}:
#'     \itemize{
#'       \item "percent" (default)
#'       \item "molar"
#'     }
#' }
#' 
#' \strong{Other Parameters:}
#' \itemize{
#'   \item \code{ambiguous}: TRUE/FALSE (default: FALSE)
#'   \item \code{shift}: Integer value (default: 0)
#'   \item \code{dnac_high}: Numeric value in nM (default: 25)
#'   \item \code{dnac_low}: Numeric value in nM (default: 25)
#'   \item \code{self_comp}: TRUE/FALSE (default: FALSE)
#'   \item \code{Na}: Millimolar concentration (default: 50)
#'   \item \code{K}: Millimolar concentration (default: 0)
#'   \item \code{Tris}: Millimolar concentration (default: 0)
#'   \item \code{Mg}: Millimolar concentration (default: 0)
#'   \item \code{dNTPs}: Millimolar concentration (default: 0)
#'   \item \code{DMSO}: Percent concentration (default: 0)
#'   \item \code{dmso_factor}: Numeric value (default: 0.75)
#'   \item \code{formamide_factor}: Numeric value (default: 0.65)
#'   \item \code{mismatch}: TRUE/FALSE (default: TRUE)
#' }
#' 
#' @param input_seq Where the sequence comes from. One of:
#'   \itemize{
#'     \item \strong{Sequences}, as a character vector in 5' to 3' direction,
#'       e.g. \code{c("ATGCG", "GGCCA")}. Names, when present, become the
#'       \code{seqnames} of the result; unnamed sequences are keyed by their
#'       position in the vector.
#'     \item \strong{An installed BSgenome package}, by name, e.g.
#'       \code{"BSgenome.Hsapiens.UCSC.hg38"}. It is named rather than passed
#'       as a loaded object because each worker opens the genome for itself,
#'       so no sequence crosses between processes. See
#'       \code{BSgenome::available.genomes()} for what exists, and install the
#'       package before calling.
#'     \item \strong{A FASTA file}, by path; gzipped files are read directly.
#'     \item \strong{A \code{GRanges}} carrying a \code{sequence} metadata
#'       column. The complement is derived from it when absent.
#'   }
#'   \code{regions} selects from any of the four, and means the same thing in
#'   each: the identifier before the colon is resolved against whatever names
#'   the source itself offers, and falls back to position. A BSgenome or a
#'   FASTA file with no \code{regions} is taken whole, which for a genome
#'   means its standard chromosomes.
#'
#'   Also accepted, and unchanged from earlier versions: a character vector of
#'   coordinate strings \code{"chr:start-end:strand:species"}, for example
#'   \code{"chr1:100-200:+:BSgenome.Hsapiens.UCSC.hg38"}, where strand
#'   defaults to \code{"+"}. This form carries its own genome in every
#'   element, so it is read directly and ignores \code{regions}; to profile
#'   coordinates against a genome, prefer passing the genome here and the
#'   coordinates as \code{regions}.
#' 
#' @param complement_seq Complementary sequence(s) in 3' to 5' direction. If not provided,
#'   the function will automatically generate it from input_seq. This is the template/target
#'   sequence that the input sequence will hybridize with. Can be provided as input_seq format besides A NULL value(default)
#' 
#' @param method Method(s) to use for Tm calculation. Can be one or more of:
#'   - "tm_nn": Nearest Neighbor thermodynamics (default)
#'   - "tm_gc": GC content-based method
#'   - "tm_wallace": Wallace rule
#'   Default: c("tm_nn", "tm_gc", "tm_wallace")
#' 
#' @param ambiguous Logical. If TRUE, ambiguous bases are taken into account when computing
#'   the G and C content. The function handles various ambiguous bases (S, W, M, K, R, Y, V, H, D, B)
#'   by proportionally distributing their contribution to GC content based on their possible
#'   nucleotide compositions. Default: FALSE
#' 
#' @param shift Integer value controlling the alignment offset between primer and template sequences.
#'   Only applicable for the NN method. Default: 0
#' 
#' @param nn_table Thermodynamic nearest-neighbor parameters for different nucleic acid hybridizations.
#'   Only applicable for the NN method. Sets whose name encodes a sodium
#'   concentration were fitted at that condition and are not salt-corrected
#'   again. See \code{\link{tm_nn}} for the full list and guidance on choosing
#'   between them. Default: "DNA_NN_SantaLucia_2004"
#' 
#'
#'   Alternatively, supply a matrix or data.frame of parameters directly. This
#'   is the route for parameter sets the package does not ship, in particular
#'   sets covering modified bases such as 5-methylcytosine. Requirements:
#'   \itemize{
#'     \item numeric, with columns 1 and 2 read as delta H (kcal/mol) and
#'       delta S (cal/mol/K); further columns are ignored;
#'     \item row names giving the parameter keys, e.g. \code{"AA/TT"},
#'       \code{"init"}, \code{"init_A/T"}, \code{"sym"};
#'     \item every key of a built-in reference set must be present. The
#'       reference is named by \code{attr(x, "reference")}, or defaults to
#'       the first built-in listed for the argument, which is a DNA/DNA set;
#'       RNA and hybrid tables should therefore set the attribute. Extra keys
#'       beyond the reference are kept, which is how a modified-base set adds
#'       stacks rather than replacing them.
#'   }
#'   The supplied table is reordered to the reference key order before use, so
#'   that two tables differing only in row order give identical results. A
#'   missing key would otherwise contribute zero to the calculation instead of
#'   raising an error, which is why the full key set is required. Keys that
#'   disagree with their reverse complement produce a warning: expected for
#'   modified bases, a transposition error otherwise.
#'
#'   Two optional attributes are honoured. \code{attr(x, "salt_mM")} marks a
#'   set as fitted at a stated sodium concentration, which suppresses the salt
#'   correction at that concentration exactly as for the built-in sets fitted
#'   this way; without it the table is treated as a reference-condition set and
#'   \code{salt_method} is applied. \code{attr(x, "end_table")} supplies a
#'   companion penultimate-pair end-effect table.
#' @param tmm_table Thermodynamic parameters for terminal mismatches. Only applicable for the NN method.
#'   Default: "DNA_TMM_Bommarito_2000"
#' 
#' @param imm_table Thermodynamic parameters for internal mismatches. Only applicable for the NN method.
#'   Default: "DNA_IMM_Peyret_1999"
#' 
#' @param de_table Thermodynamic parameters for dangling ends. Only applicable for the NN method.
#'   Default: "DNA_DE_Bommarito_2000"
#' 
#' @param dnac_high Concentration of the higher concentrated strand in nM. Only applicable for the NN method.
#'   Default: 25
#' 
#' @param dnac_low Concentration of the lower concentrated strand in nM. Only applicable for the NN method.
#'   Default: 25
#' 
#' @param self_comp Logical value indicating if the sequence is self-complementary. Only applicable
#'   for the NN method. Default: FALSE
#' 
#' @param variant Empirical constants coefficient for GC method. Only applicable for the GC method.
#'   Default: "Primer3Plus"
#' 
#' @param userset A vector of four coefficient values for GC method. Only applicable for the GC method.
#'   Usersets override value sets. Default: NULL
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
#' @param salt_method Salt correction method for Tm. Default: "Schildkraut2010"
#'   Available options:
#'   - "none": Disables salt correction. Also selected automatically when the
#'     chosen \code{nn_table} was fitted at the requested \code{Na}.
#'   - "Schildkraut2010": Updated salt correction method
#'   - "Wetmur1991": Classic salt correction method
#'   - "SantaLucia1996": DNA-specific salt correction
#'   - "SantaLucia1998-1": Improved DNA salt correction
#'   - "Owczarzy2004": Comprehensive salt correction
#'   - "Owczarzy2008": Updated comprehensive salt correction
#'   Default: "Schildkraut2010"
#'
#'   With \code{method = "tm_gc"} the salt term is part of the published
#'   formula selected by \code{variant}, so naming a \emph{different} one is
#'   ignored there, with a warning, unless \code{userset} is supplied.
#'   \code{"none"} (or \code{NA}) is not a substitution but a request to drop
#'   the correction, and is honoured on either path without a warning. The two
#'   Owczarzy corrections are not available to \code{"tm_gc"} at all: they
#'   apply to the reciprocal of the melting temperature in kelvin, referenced
#'   to the same duplex in 1 M Na+, and carry a duplex-length term the
#'   GC-content formulas already have.
#'
#' @param DMSO Percent DMSO concentration in the reaction mixture. Default: 0
#' 
#' @param formamide_unit Formamide concentration as `list(value, unit)`. Default: list(value = 0, unit = "percent")
#'   - value: Numeric value of formamide concentration
#'   - unit: Either "percent" or "molar"
#' 
#' @param dmso_factor Coefficient of Tm decreases per percent DMSO. Default: 0.75
#'   Other published values are 0.5, 0.6 and 0.675.
#' 
#' @param formamide_factor Tm decrease per percent formamide. Default: 0.65
#'   Several papers report factors between 0.6 and 0.72.
#' 
#' @param mismatch Logical. If TRUE, every '.' in the sequence is counted as a mismatch.
#'   Only applicable for the GC method. Default: TRUE
#'
#' @param regions What to take from \code{input_seq}. \code{NULL}, the
#'   default, means all of it, except for a \pkg{BSgenome}, where it means
#'   \code{GenomeInfoDb::standardChromosomes()} of that genome, which for
#'   GRCh38 includes chrM. Otherwise: names or numbers
#'   (\code{c("chr1", "chr2")}, \code{1:2}), coordinate strings
#'   \code{"name:start-end"} with commas and scientific notation accepted
#'   (\code{"chr1:5,000,000-6e6"}), a mixture of the two, or a
#'   \code{GRanges}.
#'
#'   The identifier before the colon is resolved against whatever names the
#'   source itself offers, and falls back to position. So \code{"chr1"} is a
#'   chromosome in a \pkg{BSgenome}, a record in a FASTA file and a
#'   \code{seqname} in a \code{GRanges}, and \code{"1:1-200"} is the first
#'   200 bases of the first sequence in an unnamed character vector. On a
#'   \pkg{BSgenome} the \code{chr} prefix is added or removed as the genome
#'   requires, since that is a convention rather than information; elsewhere
#'   names are matched exactly.
#' @param window Window width in base pairs. \code{NULL}, the default, gives
#'   one melting temperature per region, which is what short records such as
#'   probes, primers and oligonucleotides call for; a region longer than 1 Mb
#'   with \code{window = NULL} is an error rather than one meaningless Tm.
#'   Regions shorter than \code{window} are returned whole.
#' @param slide Step between window starts, defaulting to \code{window} for a
#'   non-overlapping tiling. Ignored when \code{window} is \code{NULL}.
#' @param unit How regions become tasks: \code{"segment"} cuts them into
#'   pieces of about \code{segment_size} bp, \code{"region"} makes one task
#'   per region. Segments are faster and need less memory per worker, because
#'   no worker then holds a whole large chromosome.
#' @param segment_size Task size in base pairs when \code{unit = "segment"},
#'   rounded down to a multiple of \code{slide} so that the window grid is
#'   the one an unsegmented run would produce. Default 50 Mb.
#' @param BPPARAM A \code{BiocParallelParam} from \pkg{BiocParallel}, e.g.
#'   \code{SnowParam(workers = 5)}, to spread the tasks over processes.
#'   \code{NULL}, the default, runs them here, with no dependency on
#'   \pkg{BiocParallel}. Parallelism divides the work by region, never the
#'   sequences of one region: each task opens the source itself, so only a
#'   name and a coordinate pair cross between processes.
#' @param keep_sequence Keep the \code{sequence} and \code{complement}
#'   columns. The default keeps them for sequences the caller supplied and
#'   drops them for a genome or a file, where they run to roughly 500 MB per
#'   large chromosome.
#' @param tmpdir Directory for the temporary FASTA file written when
#'   sequences are supplied directly and there is tiling or parallelism to
#'   do. Worth setting on a cluster, where \code{tempdir()} is often a small
#'   partition.
#' @param verbose Report the task and window counts.
#'
#' @details
#' The three methods differ in resolution and in the range of sequence lengths
#' over which they are calibrated, so they are not interchangeable.
#'
#' \code{tm_nn} is the appropriate default. Because it sums sequence-dependent
#' stacking terms, it distinguishes sequences of identical GC content but
#' different base order, which the other two cannot. Its parameters were derived
#' from short duplexes under a two-state assumption; when applied to long
#' sequences or to fixed-width genomic windows the resulting value is best read
#' as a relative measure of local thermodynamic stability for comparison across
#' windows, rather than as an absolute experimental melting temperature.
#'
#' \code{tm_gc} computes Tm from GC percentage using one of several published
#' empirical formulas selected by \code{variant}, with corrections for length
#' and ionic strength. It extends to longer sequences at low computational cost
#' but cannot resolve base order.
#'
#' \code{tm_wallace} applies the 2 + 4 rule, adding 2 \eqn{^{\circ}}C per A or T
#' and 4 \eqn{^{\circ}}C per G or C. It is calibrated for short oligonucleotides,
#' typically 14 to 20 bp, and takes no account of sequence context, salt or
#' chemical additives. Accuracy degrades quickly with length, so it is retained
#' for compatibility rather than recommended for genome-scale work.
#'
#' Salt and chemical corrections apply to \code{tm_nn} and \code{tm_gc} only.
#' The input sequence is parsed and validated once and reused by the selected
#' method, which is faster than calling the individual functions directly.
#' 
#' @return A \code{TmCalculator} list with:
#'   \item{\code{gr}}{The input \code{GRanges} with metadata columns \code{Tm}
#'     and \code{GC} (melting temperature in \eqn{^{\circ}}C and GC percent).}
#'   \item{\code{options}}{Calculation parameters and method information. For
#'     the nearest-neighbor method this includes \code{Salt correction applied}
#'     (logical) and \code{Parameter set fitted at [Na+] (mM)}, which record
#'     whether a salt correction was actually performed.}
#' 
#' @encoding UTF-8
#' @author Junhui Li
#' 
#' @export
#' 
#' @importFrom GenomeInfoDb genome
#' 
#' @examples
#' \dontrun{
#' input_seq <- c("chr1:1000100-1000150:+:BSgenome.Hsapiens.UCSC.hg38")
#' result <- tm_calculate(
#'   input_seq,
#'   method = "tm_nn",
#'   nn_table = "DNA_NN_SantaLucia_2004",
#'   salt_method = "Owczarzy2008"
#' )
#'
#' # A hybrid parameter set fitted at 100 mM sodium. Salt correction is
#' # skipped automatically because Na matches the fitted condition.
#' result_ls <- tm_calculate(
#'   input_seq,
#'   method = "tm_nn",
#'   nn_table = "RNA_DNA_NN_Weber_2019_LS",
#'   Na = 100
#' )
#'
#' # Genome scale. The source is named rather than loaded, because each
#' # worker opens it for itself; `regions` says what to cover and `window`
#' # says at what resolution.
#' hg38 <- "BSgenome.Hsapiens.UCSC.hg38"
#' tm_calculate(hg38, regions = "chr21:10e6-20e6", window = 200, slide = 200)
#'
#' # Five processes. Tasks are divided by region, never by splitting the
#' # sequences of one region, so only coordinates cross between them.
#' library(BiocParallel)
#' whole <- tm_calculate(hg38, window = 200, slide = 200,
#'                       BPPARAM = SnowParam(workers = 5))
#' whole$gr
#'
#' # The same, for a FASTA file and for sequences already in R. With
#' # window = NULL, the default, each record returns a single Tm, which is
#' # what a file of probes or primers calls for.
#' tm_calculate("probes.fa", BPPARAM = SnowParam(workers = 4))
#' tm_calculate(c("ACGTGCTAGCTAGCTAGC", "GGCCATATATGCGC"))
#' }
#'
#' @seealso \code{\link{tm_nn}} for the nearest-neighbor method and the full
#'   list of thermodynamic parameter sets.
#'
#' @export tm_calculate
tm_calculate <- function(input_seq,
                        method = c("tm_nn", "tm_gc", "tm_wallace"),
                        complement_seq = NULL,
                        ambiguous = FALSE,
                        shift = 0,
                        nn_table = c("DNA_NN_SantaLucia_2004",
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
                        tmm_table = "DNA_TMM_Bommarito_2000",
                        imm_table = "DNA_IMM_Peyret_1999",
                        de_table = c("DNA_DE_Bommarito_2000",
                                    "RNA_DE_Turner_2010"),
                        dnac_high = 25,
                        dnac_low = 25,
                        self_comp = FALSE,
                        variant = c("Primer3Plus",
                                    "Chester1993",
                                    "QuikChange",
                                    "Schildkraut1965",
                                    "Wetmur1991_MELTING",
                                    "Wetmur1991_RNA",
                                    "Wetmur1991_RNA/DNA",
                                    "vonAhsen2001"),
                        userset = NULL,
                        Na = 50,
                        K = 0,
                        Tris = 0,
                        Mg = 0,
                        dNTPs = 0,
                        salt_method = c("Schildkraut2010",
                                            "Wetmur1991",
                                            "SantaLucia1996",
                                            "SantaLucia1998-1",
                                            "Owczarzy2004",
                                            "Owczarzy2008",
                                            "none"),
                        DMSO = 0,
                        formamide_unit = list(value = 0, unit = "percent"),
                        dmso_factor = 0.75,
                        formamide_factor = 0.65,
                        mismatch = TRUE,
                        regions = NULL,
                        window = NULL,
                        slide = window,
                        unit = c("segment", "region"),
                        segment_size = 50e6,
                        BPPARAM = NULL,
                        keep_sequence = NULL,
                        tmpdir = tempdir(),
                        verbose = FALSE) {
  # Read before anything reassigns it: missing() is only reliable while the
  # argument is untouched, and `salt_method <- match.arg(salt_method)` below
  # would make it FALSE for every call. NULL counts as not named, because it
  # is the spelling tm_gc() documents for "use the formula's own"; match.arg()
  # would otherwise turn it into the first candidate and warn about a method
  # the caller never asked for.
  salt_named  <- !missing(salt_method) && !is.null(salt_method)
  method      <- match.arg(method, several.ok = FALSE)
  unit        <- match.arg(unit)
  # tm_gc() spells "no salt correction" as NA as well as "none"; match.arg()
  # cannot express NA, so the two are made the same thing here rather than
  # left as a disagreement between the two entry points.
  if (length(salt_method) == 1L && is.na(salt_method)) salt_method <- "none"
  # Validated once and passed down as a scalar. Without this the full default
  # candidate vector reached tm_gc(), whose own match.arg() has no "none"
  # choice and failed with "'arg' must be of length 1" for any tm_gc call
  # relying on defaults. Same reasoning for variant.
  salt_method <- match.arg(salt_method)
  variant     <- match.arg(variant)

  # A GC-content formula carries its own salt term, so salt_method applies to
  # tm_gc() only when userset is supplied. Saying so here, once, is the only
  # place it can be said: this function runs in the calling process, whereas
  # tm_gc() runs once per task under a BPPARAM.
  if (identical(method, "tm_gc") && salt_named) {
    if (salt_method %in% c("Owczarzy2004", "Owczarzy2008")) {
      # Rejected here rather than inside tm_gc() so that the message reaches
      # the caller instead of surfacing as one worker's error.
      stop("`salt_method = \"", salt_method, "\"` is not available for ",
           "method = \"tm_gc\". The Owczarzy corrections apply to the ",
           "reciprocal of the melting temperature in kelvin, referenced to ",
           "the same duplex in 1 M Na+, and carry a duplex-length term of ",
           "their own, which the GC-content formulas already have. Use ",
           "method = \"tm_nn\" for them.", call. = FALSE)
    }
    # "none" is not a substitution but a request to drop the correction, and
    # tm_gc() honours it on either path, so it is not warned about.
    if (is.null(userset) && !identical(salt_method, "none")) {
      own <- get_table("GC_VARTAB")[variant, "salt_correct"]
      if (!identical(salt_method, own)) {
        carries <- if (is.na(own)) "no salt term of its own"
                   else paste0("the '", own, "' salt term")
        warning("variant '", variant, "' carries ", carries, ", so ",
                "`salt_method = \"", salt_method, "\"` is ignored by tm_gc(). ",
                "Supply `userset` to choose the correction yourself.",
                call. = FALSE)
      }
    }
  }

  # Everything that describes the thermodynamic model and nothing that
  # describes where the sequence comes from. Each task hands this back to
  # tm_calculate(), so it must not carry regions, window or BPPARAM: those
  # would make the call recurse instead of compute.
  model <- list(method = method, ambiguous = ambiguous, shift = shift,
                nn_table = nn_table, tmm_table = tmm_table,
                imm_table = imm_table, de_table = de_table,
                dnac_high = dnac_high, dnac_low = dnac_low,
                self_comp = self_comp, Na = Na, K = K, Tris = Tris, Mg = Mg,
                dNTPs = dNTPs, userset = userset, variant = variant,
                salt_method = salt_method, DMSO = DMSO,
                formamide_unit = formamide_unit, dmso_factor = dmso_factor,
                formamide_factor = formamide_factor, mismatch = mismatch)

  # -- the direct route -------------------------------------------------------
  # Sequences already in hand, nothing to tile, one process. This is what the
  # function has always done for a character vector or a GRanges, and it stays
  # the shortest path through it: no temporary file, no task machinery.
  a_source <- is.character(input_seq) && length(input_seq) == 1L &&
    ((file.exists(input_seq) && !dir.exists(input_seq)) ||
       requireNamespace(input_seq, quietly = TRUE))
  if (!a_source && is.null(regions) && is.null(window) && is.null(BPPARAM)) {
    gr <- if (methods::is(input_seq, "GRanges")) .tm_complete_gr(input_seq)
          else to_genomic_ranges(input_seq = input_seq,
                                 complement_seq = complement_seq)
    return(.tm_model(gr, model))
  }

  # -- the profiling route ----------------------------------------------------
  src <- .tm_source(input_seq, complement_seq)
  if (is.null(keep_sequence))
    # Worth keeping for sequences the caller already had; ruinous for a
    # genome, where the columns run to about 500 MB per large chromosome.
    keep_sequence <- src$kind %in% c("sequences", "granges")

  # Resolved against the source the caller gave, not against the staged copy
  # below. `regions` names the caller's own sequences, and those names do not
  # always survive staging: a GRanges with two ranges on chr1 cannot put
  # "chr1" on two FASTA records, where a record name has to identify one
  # record. Resolving first also means a GRanges source is selected by
  # overlap, which is what having coordinates makes possible.
  req <- .tm_regions(regions, src)

  if (src$kind %in% c("sequences", "granges")) {
    # Staged to a file so that the workers read the sequences rather than
    # receive them, which is what moves window construction and result
    # assembly into the worker as well. Serially this costs one write and one
    # read; in parallel it replaces a serialisation that costs more.
    seqs <- if (src$kind == "granges")
      stats::setNames(as.character(GenomicRanges::mcols(src$gr)$sequence),
                      as.character(GenomeInfoDb::seqnames(src$gr)))
      else src$seqs
    path <- .spill_fasta(seqs, tmpdir)
    on.exit(unlink(path), add = TRUE)
    # Records are written in input order, so req$idx still addresses the
    # right one; req$name keeps the caller's name for the output.
    src  <- .tm_source(path)
  }

  if (is.null(window)) {
    # One window per region is right for a probe and absurd for a chromosome:
    # the nearest-neighbour model is not calibrated at that length, and the
    # single Tm it would return means nothing.
    longest <- max(req$end - req$start + 1)
    if (longest > 1e6)
      stop("'window' is NULL, which asks for one melting temperature per ",
           "region, but the longest region is ", format(longest, big.mark = ","),
           " bp.\n  Set 'window' (and 'slide') to tile it, for example ",
           "window = 200, slide = 200.")
  }
  step  <- if (is.null(slide)) 1L else as.integer(slide)
  tasks <- .tm_tasks(req, src, unit, segment_size, step)
  gr    <- .tm_run(tasks, src, window, slide, model, BPPARAM,
                   keep_sequence, verbose)

  # What the options report is what was applied, not what was asked for. On
  # the direct route that is already true, because the returned object is the
  # one tm_gc() built; here the model is reported instead, and for a built-in
  # GC variant the salt term it carries is the variant's own rather than
  # whatever salt_method resolved to. `model` itself is left alone: the
  # workers have already run against it.
  opts <- model
  if (identical(method, "tm_gc"))
    opts$salt_method <- if (identical(salt_method, "none")) NA_character_
                        else if (is.null(userset))
                          get_table("GC_VARTAB")[variant, "salt_correct"]
                        else salt_method

  result <- list(gr = gr,
                 options = c(opts, list(window = window, slide = slide,
                                        unit = unit, n_tasks = length(tasks))))
  class(result) <- c("TmCalculator", "list")
  attr(result, "nonhidden") <- "gr"

  # A data.frame representation is available lazily via result$df
  # (see `$.TmCalculator` in print.TmCalculator.R); converting
  # genome-scale GRanges eagerly here cost seconds per call.
  result
}


# ---------------------------------------------------------------------------
#' Fill in a GRanges that carries sequences but not their complements
#'
#' Every other input form reaches \code{\link{to_genomic_ranges}}, which
#' derives the complement. A \code{GRanges} handed straight to
#' \code{\link{tm_calculate}} skipped that step, so an object built with a
#' \code{sequence} column and nothing else (which is exactly what
#' \code{regions} treats as a source) reached the compiled core with some
#' sequences and no complements and failed there, on a message about
#' \code{cseqs} that says nothing about what the caller did.
#'
#' The complement is the plain base-for-base one, not the reverse
#' complement: \code{tm_nn()} reads the two strands in register, so
#' reversing one would pair every position with the wrong partner.
#'
#' @param gr A \code{GRanges}.
#' @return The same object, with a \code{complement} column.
#' @keywords internal
.tm_complete_gr <- function(gr) {
  mc <- GenomicRanges::mcols(gr)
  if (is.null(mc$sequence))
    stop("A GRanges given as 'input_seq' must carry a 'sequence' column.\n",
         "  To profile coordinates against a genome, pass the genome as ",
         "'input_seq' and the coordinates as 'regions'.")
  if (is.null(mc$complement))
    gr$complement <- generate_complement(as.character(mc$sequence))
  gr
}

# ---------------------------------------------------------------------------
#' Apply the selected model to windows that already carry their sequence
#'
#' The method dispatch that used to be the whole of \code{tm_calculate()}.
#' It is separate now because both routes through the function end here: the
#' direct one, and every task of the profiling one.
#' @param gr Windows with \code{sequence} and \code{complement} columns.
#' @param model Model arguments, as assembled by \code{\link{tm_calculate}}.
#' @return A \code{TmCalculator} object.
#' @keywords internal
.tm_model <- function(gr, model) {
  with(model, {
    # Calculate Tm using each selected method
    if ("tm_nn" %in% method) {
      result <- tm_nn(
        gr_seq = gr,
        ambiguous = ambiguous,
        shift = shift,
        nn_table = nn_table,
        tmm_table = tmm_table,
        imm_table = imm_table,
        de_table = de_table,
        dnac_high = dnac_high,
        dnac_low = dnac_low,
        self_comp = self_comp,
        Na = Na,
        K = K,
        Tris = Tris,
        Mg = Mg,
        dNTPs = dNTPs,
        salt_method = salt_method,
        DMSO = DMSO,
        formamide_unit = formamide_unit,
        dmso_factor = dmso_factor,
        formamide_factor = formamide_factor
      )
    }
  
    if ("tm_gc" %in% method) {
      result <- tm_gc(
        gr_seq = gr,
        ambiguous = ambiguous,
        userset = userset,
        variant = variant,
        Na = Na,
        K = K,
        Tris = Tris,
        Mg = Mg,
        dNTPs = dNTPs,
        # NULL lets tm_gc() take the correction published with the variant.
        # tm_calculate() has already warned if the caller named a different
        # one; repeating that warning here would repeat it once per task.
        # "none" is forwarded, because dropping the correction is a request
        # tm_gc() honours rather than a substitution it refuses.
        salt_method = if (is.null(userset) && !identical(salt_method, "none"))
                        NULL else salt_method,
        mismatch = mismatch,
        DMSO = DMSO,
        formamide_unit = formamide_unit,
        dmso_factor = dmso_factor,
        formamide_factor = formamide_factor
      )
    }
  
    if ("tm_wallace" %in% method) {
      result <- tm_wallace(
        gr_seq = gr,
        ambiguous = ambiguous
      )
    }
    result
  })
}
