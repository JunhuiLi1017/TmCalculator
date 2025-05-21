#' Calculate melting temperature using nearest neighbor thermodynamics
#' 
#' Calculate melting temperature using nearest neighbor thermodynamics
#' 
#' @param input_seq Input sequence in 5' to 3' direction. Can be provided as either:
#'   - A character string (e.g., "ATGCG")
#'   - A vector of characters (e.g., c("A","T","G","C","G"))
#' 
#' @param ambiguous Logical value controlling how ambiguous bases are handled:
#'   - TRUE: Ambiguous bases (e.g., N, R, Y) are included in calculations
#'   - FALSE (default): Ambiguous bases are excluded from calculations
#' 
#' @param complement_seq Complementary sequence in 3' to 5' direction. If not provided,
#'   the function will automatically generate it from input_seq. This is the template/target
#'   sequence that the input sequence will hybridize with.
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
#'              ^
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
#' @param Na Millimolar concentration of sodium ions. Default: 0
#' 
#' @param K Millimolar concentration of potassium ions. Default: 0
#' 
#' @param Tris Millimolar concentration of Tris buffer. Default: 0
#' 
#' @param Mg Millimolar concentration of magnesium ions. Default: 0
#' 
#' @param dNTPs Millimolar concentration of deoxynucleotide triphosphates. Default: 0
#' 
#' @param salt_correct Method for calculating salt concentration corrections to the melting temperature.
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
#' @param formamide_value_unit A list containing formamide concentration information:
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
#' input_seq <- c("AAAATTTTTTTCCCCCCCCCCCCCCGGGGGGGGGGGGTGTGCGCTGC")
#' out <- Tm_NN(input_seq, Na=50)
#' out
#' out$Options
#' 
#' @export Tm_NN
Tm_NN <- function(input_seq,
                  ambiguous=FALSE,
                  complement_seq=NULL,
                  shift=0,
                  nn_table=c("DNA_NN_SantaLucia_2004",
                             "DNA_NN_Breslauer_1986",
                             "DNA_NN_Sugimoto_1996",
                             "DNA_NN_Allawi_1998",
                             "RNA_NN_Freier_1986",
                             "RNA_NN_Xia_1998",
                             "RNA_NN_Chen_2012",
                             "RNA_DNA_NN_Sugimoto_1995"),
                  tmm_table="DNA_TMM_Bommarito_2000",
                  imm_table="DNA_IMM_Peyret_1999",
                  de_table=c("DNA_DE_Bommarito_2000",
                             "RNA_DE_Turner_2010"),
                  dnac_high=25,
                  dnac_low=25,
                  self_comp=FALSE, 
                  Na=0,
                  K=0, 
                  Tris=0, 
                  Mg=0, 
                  dNTPs=0, 
                  salt_correct=c("Schildkraut2010",
                             "Wetmur1991",
                             "SantaLucia1996",
                             "SantaLucia1998-1",
                             "SantaLucia1998-2",
                             "Owczarzy2004",
                             "Owczarzy2008"),
                  DMSO=0,
                  formamide_value_unit=list(value=0, unit="percent"),
                  dmso_factor=0.75,
                  formamide_factor=0.65){
  nn_table <- match.arg(nn_table)
  tmm_table <- match.arg(tmm_table)
  imm_table <- match.arg(imm_table)
  de_table <- match.arg(de_table)
  salt_correct <- match.arg(salt_correct)
  
  # Get tables from the thermodynamic_tables data
  imm_table_list <- thermodynamic_tables[[imm_table]]
  nn_table_list <- thermodynamic_tables[[nn_table]]
  tmm_table_list <- thermodynamic_tables[[tmm_table]]
  de_table_list <- thermodynamic_tables[[de_table]]
  
  imm_table_name <- rownames(imm_table_list)
  nn_table_name <- rownames(nn_table_list)
  tmm_table_name <- rownames(tmm_table_list)
  de_table_name <- rownames(de_table_list)
  
  my_seq <- check_filter(input_seq,method = "Tm_NN")
  my_seq_c2s <- c2s(my_seq)
  pt_gc <- GC(my_seq,ambiguous = ambiguous)
  if(is.null(complement_seq)){
    complement_seq <- complement(my_seq_c2s,FALSE)
  }
  my_cseq <- check_filter(complement_seq,method = "Tm_NN")
  
  tmp_seq <- my_seq
  tmp_cseq <- my_cseq
  delta_h <- 0
  delta_s <- 0
  d_h <- 1
  d_s <- 2
  if(shift!=0 | length(my_seq)!=length(my_cseq)){
    if(shift>0){
      tmp_seq <- append(rep('.',shift),my_seq)
    }else{
      tmp_cseq <- append(rep('.',abs(shift)),my_cseq)
    }
    if(length(tmp_cseq)>length(tmp_seq)){
      tmp_seq <- append(tmp_seq,rep('.',length(tmp_cseq)-length(tmp_seq)))
    }
    if(length(tmp_cseq)<length(tmp_seq)){
      tmp_cseq <- append(tmp_cseq,rep('.',length(tmp_seq)-length(tmp_cseq)))
    }
    while(all(tmp_seq[1:2]==".") | all(tmp_cseq[1:2]==".")){
      tmp_seq <- tmp_seq[-1]
      tmp_cseq <- tmp_cseq[-1]
    }
    while(all(tmp_seq[(length(tmp_seq)-1):length(tmp_seq)]==".") | all(tmp_cseq[(length(tmp_cseq)-1):length(tmp_cseq)]==".")){
      tmp_seq <- tmp_seq[-length(tmp_seq)]
      tmp_cseq <- tmp_cseq[-length(tmp_cseq)]
    } 
    if(tmp_seq[1]=="." | tmp_cseq[1]=="."){
      left_de <- paste0(c2s(tmp_seq[1:2]),'/',c2s(tmp_cseq[1:2]),collapse='')
      if(left_de %in% de_table_name){
        delta_h <- de_table_list[left_de,d_h]+delta_h
        delta_s <- de_table_list[left_de,d_s]+delta_s
      }else{
        stop("No such combination in de_table_list")
      }
      tmp_seq <- tmp_seq[-1]
      tmp_cseq <- tmp_cseq[-1]
    }
    if(tmp_seq[length(tmp_seq)] == '.' | tmp_cseq[length(tmp_cseq)] == '.'){
      right_de <- paste0(c2s(tmp_cseq[c(length(tmp_cseq),(length(tmp_cseq)-1))]),'/',
                         c2s(tmp_seq[c(length(tmp_seq),(length(tmp_seq)-1))]),collapse='')
      if(right_de %in% de_table_name){
        delta_h <- de_table_list[right_de,d_h]+delta_h
        delta_s <- de_table_list[right_de,d_s]+delta_s
      }else{
        stop("No such combination in de_table_list")
      }
      tmp_seq <- tmp_seq[-length(tmp_seq)]
      tmp_cseq <- tmp_cseq[-length(tmp_cseq)]
    }
  }
  left_tmm <- paste0(c2s(tmp_cseq[c(2,1)]),'/',c2s(tmp_seq[c(2,1)]),collapse='')
  if(left_tmm %in% tmm_table_name){
    delta_h <- tmm_table_list[left_tmm,d_h]+delta_h
    delta_s <- tmm_table_list[left_tmm,d_s]+delta_s
    tmp_seq <- tmp_seq[-1]
    tmp_cseq <- tmp_cseq[-1]
  }
  right_tmm <- paste0(c2s(tmp_seq[(length(tmp_seq)-1):length(tmp_seq)]),'/',
                      c2s(tmp_cseq[(length(tmp_cseq)-1):length(tmp_cseq)]),collapse='')
  if(right_tmm %in% tmm_table_name){
    delta_h <- tmm_table_list[right_tmm,d_h]+delta_h
    delta_s <- tmm_table_list[right_tmm,d_s]+delta_s
    tmp_seq <- tmp_seq[-length(tmp_seq)]
    tmp_cseq <- tmp_cseq[-length(tmp_cseq)]
  }
  delta_h <- nn_table_list['init',d_h]+delta_h
  delta_s <- nn_table_list['init',d_s]+delta_s
  
  if(GC(my_seq) == 0){
    delta_h <- nn_table_list['init_allA/T',d_h]+delta_h
    delta_s <- nn_table_list['init_allA/T',d_s]+delta_s
  }else{
    delta_h <- nn_table_list['init_oneG/C',d_h]+delta_h
    delta_s <- nn_table_list['init_oneG/C',d_s]+delta_s
  }

  if(my_seq[1] == 'T'){
    delta_h <- nn_table_list['init_5T/A',d_h]+delta_h
    delta_s <- nn_table_list['init_5T/A',d_s]+delta_s
  }
  if(my_seq[1] == 'A'){
    delta_h <- nn_table_list['init_5T/A',d_h]+delta_h
    delta_s <- nn_table_list['init_5T/A',d_s]+delta_s
  }
  
  ends <- c(my_seq[1],my_seq[length(my_seq)])
  at <- sum(ends %in% 'A')+sum(ends %in% 'T')
  gc <- sum(ends %in% 'G')+sum(ends %in% 'C')
  delta_h <- nn_table_list['init_A/T',d_h]*at+delta_h
  delta_s <- nn_table_list['init_A/T',d_s]*at+delta_s
  delta_h <- nn_table_list['init_G/C',d_h]*gc+delta_h
  delta_s <- nn_table_list['init_G/C',d_s]*gc+delta_s

  for(bn in 1:(length(tmp_seq)-1)){
    neighbors <- paste0(c2s(tmp_seq[bn:(bn+1)]),'/',c2s(tmp_cseq[bn:(bn+1)]),collapse='')
    rev_neighbors <- paste0(c2s(tmp_cseq[(bn+1):bn]),'/',c2s(tmp_seq[(bn+1):bn]),collapse='')
    if(neighbors %in% imm_table_name){
      delta_h <- imm_table_list[neighbors,d_h]+delta_h
      delta_s <- imm_table_list[neighbors,d_s]+delta_s
    }else if(rev_neighbors %in% imm_table_name){
      delta_h <- imm_table_list[rev_neighbors,d_h]+delta_h
      delta_s <- imm_table_list[rev_neighbors,d_s]+delta_s
    }else if(neighbors %in% nn_table_name){
      delta_h <- nn_table_list[neighbors,d_h]+delta_h
      delta_s <- nn_table_list[neighbors,d_s]+delta_s
    }else if(rev_neighbors %in% nn_table_name){
      delta_h <- nn_table_list[rev_neighbors,d_h]+delta_h
      delta_s <- nn_table_list[rev_neighbors,d_s]+delta_s
    }else{
      stop("No such combination in de_table_list")
    }
  }
  k <- (dnac_high-(dnac_low/2.0))*1e-9
  if(self_comp==TRUE){
    k <- dnac_high*1e-9
    delta_h <- nn_table_list['sym',d_h]
    delta_s <- nn_table_list['sym',d_s]
  }
  R <- 1.987
  if(!is.null(salt_correct)){
    corr_salt = salt_correction(Na=Na,K=K,Tris=Tris,Mg=Mg,dNTPs=dNTPs,method=salt_correct,input_seq=my_seq,ambiguous = ambiguous)
    if(salt_correct == "SantaLucia1998-2"){
      delta_s <- corr_salt+delta_s
    }
    tm <- (1000*delta_h)/(delta_s+(R*(log(k))))-273.15
    if(salt_correct %in% c("Schildkraut2010","Wetmur1991","SantaLucia1996","SantaLucia1998-1")){
      tm <- tm+corr_salt
    }
    if(salt_correct %in% c("Owczarzy2004","Owczarzy2008")){
      tm <- (1/(1/(tm+273.15)+corr_salt)-273.15)
    }
  }else{
    tm <- (1000*delta_h)/(delta_s+(R*(log(k))))-273.15
  }
  
  corr_chem <- chem_correction(DMSO=DMSO,formamide_value_unit=formamide_value_unit,dmso_factor=dmso_factor,formamide_factor=formamide_factor,pt_gc=pt_gc)
  tm <- tm + corr_chem
  
  nn_table_list <- list("DNA_NN_Breslauer_1986"="Breslauer K J (1986) <doi:10.1073/pnas.83.11.3746>",
                      "DNA_NN_Sugimoto_1996"="Sugimoto N (1996) <doi:10.1093/nar/24.22.4501>",
                      "DNA_NN_Allawi_1998"="Allawi H (1998) <doi:10.1093/nar/26.11.2694>",
                      "DNA_NN_SantaLucia_2004"="SantaLucia J (2004) <doi:10.1146/annurev.biophys.32.110601.141800>",
                      "RNA_NN_Freier_1986"="Freier S (1986) <doi:10.1073/pnas.83.24.9373>",
                      "RNA_NN_Xia_1998"="Xia T (1998) <doi:10.1021/bi9809425>",
                      "RNA_NN_Chen_2012"="Chen JL (2012) <doi:10.1021/bi3002709>",
                      "RNA_DNA_NN_Sugimoto_1995"="Sugimoto N (1995)<doi:10.1016/S0048-9697(98)00088-6>",
                      "DNA_TMM_Bommarito_2000"="Bommarito S (2000)  <doi:10.1093/nar/28.9.1929>",
                      "DNA_IMM_Peyret_1999"="Peyret N (1999) <doi:10.1021/bi9825091> & Allawi H T (1997) <doi:10.1021/bi962590c> & Santalucia N (2005) <doi:10.1093/nar/gki918>",
                      "DNA_DE_Bommarito_2000"="Bommarito S (2000) <doi:10.1093/nar/28.9.1929>",
                      "RNA_DE_Turner_2010"="Turner D H (2010) <doi:10.1093/nar/gkp892>")

  result_list <- vector('list',2L)
  names(result_list) <- c("Tm","Options")
  result_list$Tm <- as.numeric(tm)
  result_list$Options <- list("Sequence"=input_seq,
                            "Check filter"=c2s(my_seq),
                            "Ambiguous"=ambiguous,
                            "Complement Sequence" = complement_seq,  
                            "Shift"= shift, 
                            "Thermodynamic NN values" = paste0(nn_table,": ",nn_table_list[[nn_table]]), 
                            "Thermodynamic values for terminal mismatches" = paste0(tmm_table,": ",nn_table_list[[tmm_table]]), 
                            "Thermodynamic values for internal mismatches" = paste0(imm_table,": ",nn_table_list[[imm_table]]),
                            "Thermodynamic values for dangling ends" = paste0(de_table,": ",nn_table_list[[de_table]]), 
                            "Concentration of the higher concentrated strand" = dnac_high,
                            "Concentration of the lower concentrated strand" = dnac_low, 
                            "Sequence self-complementary" = self_comp, 
                            "Na"=Na,
                            "K"=K,
                            "Tris"=Tris,
                            "Mg"=Mg,
                            "dNTPs"=dNTPs,
                            "Salt correction method"=salt_correct,
                            "Percent of DMSO"=DMSO,
                            "Formamide concentration"=formamide_value_unit$value,
                            "Coeffecient of Tm decreases per percent DMSO"=dmso_factor,
                            "Method for formamide concentration"=formamide_value_unit$unit,
                            "Coefficient of Tm decrease per percent formamide"=formamide_factor,
                            "Percent of GC"=pt_gc)
  class(result_list) <- c("TmCalculator","list")
  attr(result_list, "nonhidden") <- "Tm"
  return(result_list)
}
