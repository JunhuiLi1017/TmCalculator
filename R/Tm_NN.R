#' Calculate melting temperature using nearest neighbor thermodynamics
#' 
#' Calculate melting temperature using nearest neighbor thermodynamics. The function checks if all
#' sequence combinations in the input sequence are present in the thermodynamic parameter tables
#' before performing calculations.
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
#'   Eight parameter sets are available, organized by hybridization type:
#' 
#'   DNA/DNA hybridizations:
#'   - "DNA_NN_Breslauer_1986": Original DNA/DNA parameters
#'   - "DNA_NN_Sugimoto_1996": Improved DNA/DNA parameters
#'   - "DNA_NN_Allawi_1998": DNA/DNA parameters with internal mismatch corrections
#'   - "DNA_NN_SantaLucia_2004": Updated DNA/DNA parameters
#' 
#'   RNA/RNA hybridizations:
#'   - "RNA_NN_Freier_1986": Original RNA/RNA parameters
#'   - "RNA_NN_Xia_1998": Improved RNA/RNA parameters
#'   - "RNA_NN_Chen_2012": Updated RNA/RNA parameters with GU pair corrections
#' 
#'   RNA/DNA hybridizations:
#'   - "RNA_DNA_NN_Sugimoto_1995": RNA/DNA hybridization parameters
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
#'   - "SantaLucia1998-2": Alternative DNA salt correction
#'   - "Owczarzy2004": Comprehensive salt correction
#'   - "Owczarzy2008": Updated comprehensive salt correction
#'   Note: Setting to NA disables salt correction
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
#' Santalucia N E W J . Nearest-neighbor thermodynamics of deoxyinosine pairs in DNA duplexes[J]. Nucleic Acids Research, 2005, 33(19):6258-67.
#' 
#' Peyret N , Seneviratne P A , Allawi H T , et al. Nearest-Neighbor Thermodynamics and NMR of DNA Sequences with Internal A-A, C-C, G-G, and T-T Mismatches, [J]. Biochemistry, 1999, 38(12):3468-3477.
#' 
#' @author Junhui Li
#' 
#' @examples
#' 
#' input_seq <- c("AAAATTTTTTTCCCCCCCCCCCCCCGGGGGGGGGGGGTGTGCGCTGC",
#' "AAAATTTTTTTCCCCCCCCCCCCCCGGGGGGGGGGGGTGTGCGCTGC")
#' seqs <- to_genomic_ranges(input_seq)
#' out <- tm_nn(seqs, Na=50)
#' out
#' 
#' @export tm_nn

tm_nn <- function(gr_seq,
                  ambiguous     = FALSE,
                  shift         = 0,
                  nn_table      = c("DNA_NN_SantaLucia_2004",
                                     "DNA_NN_Breslauer_1986",
                                     "DNA_NN_Sugimoto_1996",
                                     "DNA_NN_Allawi_1998",
                                     "RNA_NN_Freier_1986",
                                     "RNA_NN_Xia_1998",
                                     "RNA_NN_Chen_2012",
                                     "RNA_DNA_NN_Sugimoto_1995"),
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
                                        "Owczarzy2008"),
                  DMSO           = 0,
                  formamide_unit = list(value = 0, unit = "percent"),
                  dmso_factor    = 0.75,
                  formamide_factor     = 0.65) {

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
  
  # Process sequence with pairwise N filtering
  region_ids <- names(gr_seq)
  if (is.null(region_ids)) {
    region_ids <- rep("", length(gr_seq))
  }
  empty_id <- is.na(region_ids) | region_ids == ""
  if (any(empty_id)) {
    region_ids[empty_id] <- paste0(
      as.character(GenomicRanges::seqnames(gr_seq))[empty_id], ":",
      GenomicRanges::start(gr_seq)[empty_id], "-",
      GenomicRanges::end(gr_seq)[empty_id]
    )
  }
  
  filtered_seq <- check_filter_seq(
    list(
      sequence = gr_seq$sequence,
      complement = gr_seq$complement,
      region_ids = region_ids
    ),
    method = "tm_nn"
  )
  
  gr_seq_dropoff <- gr_seq[!filtered_seq$kept]
  gr_seq <- gr_seq[filtered_seq$kept]
  gr_seq$sequence <- filtered_seq$sequence
  gr_seq$complement <- filtered_seq$complement
  
  # -- Process all sequences -------------------------------------------------
  n   <- length(gr_seq)
  tm  <- numeric(n)
  gc <- numeric(n)

  all_seqs  <- as.character(mcols(gr_seq)$sequence)
  all_cseqs <- as.character(mcols(gr_seq)$complement)
  
  for (i in seq_len(n)) {
    seq_str  <- all_seqs[i]
    cseq_str <- all_cseqs[i]
    
    if (nchar(seq_str) < 2L ) {
      tm[i] <- NA_real_
      gc[i] <- NA_real_
      next
    }
    result <- tryCatch(
      .tm_nn_core(seq_str, cseq_str, ambiguous, shift, nn_tbl=nn_tbl, tmm_tbl=tmm_tbl, imm_tbl=imm_tbl, de_tbl=de_tbl, dnac_high, dnac_low, self_comp,
                  Na, K, Tris, Mg, dNTPs, salt_fn=salt_method, DMSO, dmso_factor, formamide_factor, formamide_unit),
      error = function(e) NA_real_
    )
    tm[i] <- result$Tm
    gc[i] <- result$GC
  }
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
                        "RNA_DE_Turner_2010" = "Turner D H (2010) <doi:10.1093/nar/gkp892>")
  
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

# -- Core single-sequence NN computation --------------------------------------
# All table lookups are vectorized (no per-base loop)

.tm_nn_core <- function(seq_str, cseq_str, ambiguous, shift, nn_tbl, tmm_tbl, imm_tbl, de_tbl, dnac_high, dnac_low, self_comp,
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
  if(!is.null(salt_fn)){
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
  pt_gc <- .GC_fast(seq_str, ambiguous = ambiguous)
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

# -- Helper: get table (package data or build) ------------------------------
get_table <- function(table_name) {
  if (!exists(table_name, envir = .TM_TABLE_CACHE, inherits = FALSE)) {
    assign(table_name, .TM_CONSTANTS[[table_name]], envir = .TM_TABLE_CACHE)
  }
  get(table_name, envir = .TM_TABLE_CACHE, inherits = FALSE)
}


# -----------------------------------------------------------------------------
# FIX 5: Fast GC calculation
# Add to R/GC.R as an internal helper .GC_fast()
# -----------------------------------------------------------------------------

#' @keywords internal
.GC_fast <- function(seq_upper, ambiguous = FALSE) {
  # seq_upper: already uppercased character string
  n <- nchar(seq_upper)
  if (n == 0L) return(NA_real_)

  if (!ambiguous) {
    nGC <- n - nchar(gsub("[GC]", "", seq_upper, perl = TRUE))
    return(100 * nGC / n)
  }

  # Ambiguous: count each IUPAC code's GC contribution
  # G=1, C=1, S(G+C)=1, Y(C+T)=0.5, R(A+G)=0.5,
  # K(G+T)=0.5, M(A+C)=0.5, B(CGT)=2/3, D(AGT)=1/3, H(ACT)=1/3, V(ACG)=2/3
  nG <- n - nchar(gsub("G", "", seq_upper, fixed = TRUE))
  nC <- n - nchar(gsub("C", "", seq_upper, fixed = TRUE))
  nS <- n - nchar(gsub("S", "", seq_upper, fixed = TRUE))   # G+C = 1.0
  nY <- n - nchar(gsub("Y", "", seq_upper, fixed = TRUE))   # C+T = 0.5
  nR <- n - nchar(gsub("R", "", seq_upper, fixed = TRUE))   # A+G = 0.5
  nK <- n - nchar(gsub("K", "", seq_upper, fixed = TRUE))   # G+T = 0.5
  nM <- n - nchar(gsub("M", "", seq_upper, fixed = TRUE))   # A+C = 0.5
  nB <- n - nchar(gsub("B", "", seq_upper, fixed = TRUE))   # C+G+T = 2/3
  nV <- n - nchar(gsub("V", "", seq_upper, fixed = TRUE))   # A+C+G = 2/3

  gc_count <- nG + nC + nS + 0.5*(nY + nR + nK + nM) + (2/3)*(nB + nV)
  100 * gc_count / n
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
