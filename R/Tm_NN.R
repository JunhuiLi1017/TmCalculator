#' Calculate melting temperature using nearest neighbor thermodynamics
#' 
#' Calculate melting temperature using nearest neighbor thermodynamics
#' 
#' @param ntseq Sequence (5' to 3') of one strand of the nucleic acid duplex
#'   as string or vector of characters.
#' 
#' @param ambiguous Ambiguous bases are taken into account to compute the G
#'  and C content when ambiguous is TRUE.Default is FALSE.
#' 
#' @param comSeq Complementary sequence. The sequence of the template/target in 3'->5' direction
#' 
#' @param shift Shift of the primer/probe sequence on the template/target sequence, default=0. 
#' for example: when shift=0, the first nucleotide base at 5` end of primer align to first one at 3` 
#' end of template. When shift=-1, the second nucleotide base at 5` end of primer align to first one at 3` 
#' end of template.
#' 
#' When shift=1, the first nucleotide base at 5` end of primer align to second one at 3` end 
#' of template. The shift parameter is necessary to align  primer/probe and template/target 
#' if they have different lengths or if they should have dangling ends.
#' 
#' @param nn_table Thermodynamic NN values, eight tables are implemented.
#' 
#' For DNA/DNA hybridizations:
#'   DNA_NN1,DNA_NN2,DNA_NN3,DNA_NN4
#' 
#' For RNA/RNA hybridizations:
#'   RNA_NN1,RNA_NN2,RNA_NN3
#' 
#' For RNA/DNA hybridizations:
#'   R_DNA_NN1
#' 
#' @param tmm_table Thermodynamic values for terminal mismatches. Default: DNA_TMM1
#' 
#' @param imm_table Thermodynamic values for internal mismatches, may include insosine mismatches. Default: DNA_IMM1
#' 
#' @param de_table Thermodynamic values for dangling ends. DNA_DE1(default) and RNA_DE1
#' 
#' @param dnac1 Concentration of the higher concentrated strand [nM]. Typically this will
#'  be the primer (for PCR) or the probe. Default=25.
#' 
#' @param dnac2 Concentration of the lower concentrated strand [nM].
#' 
#' @param selfcomp Sequence self-complementary, default=False. If 'True' 
#' the primer is thought binding to itself, thus dnac2 is not considered.
#' 
#' @param Na Millimolar concentration of Na, default is 0
#' 
#' @param K Millimolar concentration of K, default is 0
#' 
#' @param Tris Millimolar concentration of Tris, default is 0
#' 
#' @param Mg Millimolar concentration of Mg, default is 0
#' 
#' @param dNTPs Millimolar concentration of dNTPs, default is 0
#' 
#' @param saltcorr Salt correction method should be chosen when provide 'userset' 
#' Options are "Schildkraut2010", "Wetmur1991","SantaLucia1996","SantaLucia1998-1",
#' "SantaLucia1998-2","Owczarzy2004","Owczarzy2008". Note that NA means no salt correction. 
#' 
#' @param DMSO Percent DMSO
#' 
#' @param fmd Formamide concentration in percentage (fmdmethod="concentration") or molar (fmdmethod="molar").
#' 
#' @param DMSOfactor Coeffecient of Tm decreases per percent DMSO. Default=0.75 von Ahsen N (2001) <PMID:11673362>. Other published values are 0.5, 0.6 and 0.675.
#' 
#' @param fmdfactor Coeffecient of Tm decrease per percent formamide. Default=0.65. Several papers report factors between 0.6 and 0.72.
#' 
#' @param fmdmethod "concentration" method for formamide concentration in percentage and "molar" for formamide concentration in molar.
#' 
#' @details 
#' 
#'  DNA_NN1: Breslauer K J (1986) <doi:10.1073/pnas.83.11.3746>
#'  
#'  DNA_NN2: Sugimoto N (1996) <doi:10.1093/nar/24.22.4501>
#'  
#'  DNA_NN3: Allawi H (1998) <doi:10.1093/nar/26.11.2694>
#'  
#'  DNA_NN4: SantaLucia J (2004) <doi:10.1146/annurev.biophys.32.110601.141800>
#'  
#'  RNA_NN1: Freier S (1986) <doi:10.1073/pnas.83.24.9373>
#'  
#'  RNA_NN2: Xia T (1998) <doi:10.1021/bi9809425>
#'  
#'  RNA_NN3: Chen JL (2012) <doi:10.1021/bi3002709>
#'  
#'  R_DNA_NN1: Sugimoto N (1995)<doi:10.1016/S0048-9697(98)00088-6>
#'  
#'  DNA_TMM1: Bommarito S (2000)  <doi:10.1093/nar/28.9.1929>
#'  
#'  DNA_IMM1: Peyret N (1999) <doi:10.1021/bi9825091> & Allawi H T (1997) <doi:10.1021/bi962590c> & Santalucia N (2005) <doi:10.1093/nar/gki918>
#'  
#'  DNA_DE1: Bommarito S (2000) <doi:10.1093/nar/28.9.1929>
#'  
#'  RNA_DE1: Turner D H (2010) <doi:10.1093/nar/gkp892>
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
#' ntseq <- c("AAAATTTTTTTCCCCCCCCCCCCCCGGGGGGGGGGGGTGTGCGCTGC")
#' out <- Tm_NN(ntseq,Na=50)
#' out
#' out$Options
#' 
#' @export Tm_NN
Tm_NN <- function(ntseq,
                  ambiguous=FALSE,
                  comSeq=NULL,
                  shift=0,
                  nn_table=c("DNA_NN4",
                             "DNA_NN1",
                             "DNA_NN2",
                             "DNA_NN3",
                             "RNA_NN1",
                             "RNA_NN2",
                             "RNA_NN3",
                             "R_DNA_NN1"),
                  tmm_table="DNA_TMM1",
                  imm_table="DNA_IMM1",
                  de_table=c("DNA_DE1",
                             "RNA_DE1"),
                  dnac1=25,
                  dnac2=25,
                  selfcomp=FALSE, 
                  Na=0,
                  K=0, 
                  Tris=0, 
                  Mg=0, 
                  dNTPs=0, 
                  saltcorr=c("Schildkraut2010",
                             "Wetmur1991",
                             "SantaLucia1996",
                             "SantaLucia1998-1",
                             "SantaLucia1998-2",
                             "Owczarzy2004",
                             "Owczarzy2008"),
                  DMSO=0,
                  fmd=0, 
                  DMSOfactor=0.75,
                  fmdfactor=0.65,
                  fmdmethod=c("concentration","molar")){
  nn_table <- match.arg(nn_table)
  tmm_table <- match.arg(tmm_table)
  imm_table <- match.arg(imm_table)
  de_table <- match.arg(de_table)
  saltcorr <- match.arg(saltcorr)
  DNA_NN1 <- matrix(c(0,0,0,0,0,0,0,-16.8,0,-20.1,0,0,0,-1.3,-9.1,-24,-8.6,-23.9,-6,-16.9,-5.8,-12.9,
                      -6.5,-17.3,-7.8,-20.8,-5.6,-13.5,-11.9,-27.8,-11.1,-26.7,-11,-26.6),ncol=2,byrow = TRUE)
  rownames(DNA_NN1) <- c("init","init_A/T","init_G/C","init_oneG/C","init_allA/T","init_5T/A","sym","AA/TT",
                         "AT/TA","TA/AT","CA/GT","GT/CA","CT/GA","GA/CT","CG/GC","GC/CG","GG/CC")
  colnames(DNA_NN1) <- c("left","right")
  DNA_NN2 <- matrix(c(0.6,-9,0,0,0,0,0,0,0,0,0,0,0,-1.4,-8,-21.9,-5.6,-15.2,-6.6,-18.4,-8.2,-21,-9.4,
                      -25.5,-6.6,-16.4,-8.8,-23.5,-11.8,-29,-10.5,-26.4,-10.9,-28.4),ncol=2,byrow = TRUE)
  rownames(DNA_NN2) <- rownames(DNA_NN1)
  colnames(DNA_NN2) <- c("left","right")
  DNA_NN3 <- matrix(c(0,0,2.3,4.1,0.1,-2.8,0,0,0,0,0,0,0,-1.4,-7.9,-22.2,-7.2,-20.4,-7.2,-21.3,-8.5,
                      -22.7,-8.4,-22.4,-7.8,-21,-8.2,-22.2,-10.6,-27.2,-9.8,-24.4,-8,-19.9),ncol=2,byrow = TRUE)
  rownames(DNA_NN3) <- rownames(DNA_NN1)
  colnames(DNA_NN3) <- c("left","right")
  DNA_NN4 <- matrix(c(0.2,-5.7,2.2,6.9,0,0,0,0,0,0,0,0,0,-1.4,-7.6,-21.3,-7.2,-20.4,-7.2,-20.4,-8.5,
                      -22.7,-8.4,-22.4,-7.8,-21,-8.2,-22.2,-10.6,-27.2,-9.8,-24.4,-8,-19),ncol=2,byrow = TRUE)
  rownames(DNA_NN4) <- rownames(DNA_NN1)
  colnames(DNA_NN4) <- c("left","right")
  RNA_NN1 <- matrix(c(0,-10.8,0,0,0,0,0,0,0,0,0,0,0,-1.4,-6.6,-18.4,-5.7,-15.5,-8.1,-22.6,-10.5,-27.8,
                      -10.2,-26.2,-7.6,-19.2,-13.3,-35.5,-8,-19.4,-14.2,-34.9,-12.2,-29.7),ncol=2,byrow = TRUE)
  rownames(RNA_NN1) <- rownames(DNA_NN1)
  colnames(RNA_NN1) <- c("left","right")
  RNA_NN2 <- matrix(c(3.61,-1.5,3.72,10.5,0,0,0,0,0,0,0,0,0,-1.4,-6.82,-19,-9.38,-26.7,-7.69,-20.5,-10.44,-26.9,
                      -11.4,-29.5,-10.48,-27.1,-12.44,-32.5,-10.64,-26.7,-14.88,-36.9,-13.39,-32.7),ncol=2,byrow = TRUE)
  rownames(RNA_NN2) <- rownames(DNA_NN1)
  colnames(RNA_NN2) <- c("left","right")
  RNA_NN3 <- matrix(c(6.4,6.99,3.85,11.04,0,0,0,0,0,0,0,0,0,-1.4,-7.09,-19.8,-9.11,-25.8,-8.5,-22.9,-11.03,-28.8,
                      -11.98,-31.3,-10.9,-28.5,-13.21,-34.9,-10.88,-27.4,-16.04,-40.6,-14.18,-35,-13.83,-46.9,
                      -17.82,-56.7,-3.96,-11.6,-0.96,-1.8,-10.38,-31.8,-12.64,-38.9,-7.39,-21,-5.56,-13.9,-9.44,
                      -24.7,-7.03,-16.8,-11.09,-28.8),ncol=2,byrow = TRUE)
  rownames(RNA_NN3) <- c("init","init_A/T","init_G/C","init_oneG/C","init_allA/T","init_5T/A","sym","AA/TT",
                         "AT/TA","TA/AT","CA/GT","GT/CA","CT/GA","GA/CT","CG/GC","GC/CG","GG/CC","GT/TG",
                         "GG/TT","AG/TT","TG/AT","TT/AG","TG/GT","AT/TG","CG/GT","CT/GG","GG/CT","GT/CG")
  colnames(RNA_NN3) <- c("left","right")
  R_DNA_NN1 <- matrix(c(1.9,-3.9,0,0,0,0,0,0,0,0,0,0,0,0,-11.5,-36.4,-7.8,-21.6,-7,-19.7,-8.3,-23.9,-10.4,-28.4,
                        -12.8,-31.9,-16.3,-47.1,-9.1,-23.5,-8.6,-22.9,-8,-17.1,-9.3,-23.2,-5.9,-12.3,-7.8,-23.2,
                        -5.5,-13.5,-9,-26.1,-7.8,-21.9),ncol=2,byrow = TRUE)
  rownames(R_DNA_NN1) <- c("init","init_A/T","init_G/C","init_oneG/C","init_allA/T","init_5T/A","sym","AA/TT",
                           "AC/TG","AG/TC","AT/TA","CA/GT","CC/GG","CG/GC","CT/GA","GA/CT","GC/CG","GG/CC",
                           "GT/CA","TA/AT","TC/AG","TG/AC","TT/AA")
  colnames(R_DNA_NN1) <- c("left","right")
  DNA_IMM1 <- matrix(c(1,0.9,-2.5,-8.3,-4.1,-11.7,-2.8,-8,3.3,10.4,5.8,16.3,-4.4,-12.3,4.1,9.5,-0.1,-1.7,-1.4,
                       -6.2,-1.3,-5.3,-0.6,-2.3,-0.7,-2.3,-0.7,-2.3,-4,-13.2,-0.6,-1,0.5,3.2,0.7,0.7,3,7.4,0.7,
                       0.2,-1.2,-6.2,-0.8,-4.5,-1.5,-6.1,2.3,5.4,5.2,13.5,1.2,0.7,1,0.7,2.3,4.6,5.3,14.6,1.9,
                       3.7,0.6,-0.6,5.2,14.2,-0.7,-3.8,3.4,8,7.6,20.2,1.2,1.7,-0.9,-4.2,-2.9,-9.8,4.7,12.9,0,
                       -4.4,-1.5,-7.2,3.6,8.9,6.1,16.4,-3.1,-9.5,-4.9,-15.3,-6,-15.8,1.6,3.6,-2.7,-10.8,-5,-15.8,
                       -2.2,-8.4,0.2,-1.5,-8.9,-25.5,-5.9,-17.4,-8.8,-25.4,-4.9,-13.9,-5.4,-13.7,-6.8,-19.1,-8.3,
                       -23.8,-5,-12.6,-8.3,-25,-3.4,-11.2,-0.7,-2.6,-1.3,-4.6,2.6,8.9,-7.8,-21.1,-7,-20,-7.6,
                       -20.2,0.49,-0.7,-6.5,-22,-5.6,-18.7,-0.8,-4.3,-1,-2.4,-3.5,-10.6,0.1,-1,-4.3,-12.1,-4.9,
                       -15.8,-1.9,-8.5,0.1,-1.8,1,1,7.1,21.3,-1.1,-3.2,5.8,16.9,-7.6,-22,-3.3,-11.9,0.1,-2.3,
                       1.3,3,-0.5,-1.3),ncol=2,byrow = TRUE)
  rownames(DNA_IMM1) <- c("AG/TT","AT/TG","CG/GT","CT/GG","GG/CT","GG/TT","GT/CG","GT/TG","TG/AT","TG/GT",
                          "TT/AG","AA/TG","AG/TA","CA/GG","CG/GA","GA/CG","GG/CA","TA/AG","TG/AA","AC/TT",
                          "AT/TC","CC/GT","CT/GC","GC/CT","GT/CC","TC/AT","TT/AC","AA/TC","AC/TA","CA/GC",
                          "CC/GA","GA/CC","GC/CA","TA/AC","TC/AA","AA/TA","CA/GA","GA/CA","TA/AA","AC/TC",
                          "CC/GC","GC/CC","TC/AC","AG/TG","CG/GG","GG/CG","TG/AG","AT/TT","CT/GT","GT/CT",
                          "TT/AT","AI/TC","TI/AC","AC/TI","TC/AI","CI/GC","GI/CC","CC/GI","GC/CI","AI/TA",
                          "TI/AA","AA/TI","TA/AI","CI/GA","GI/CA","CA/GI","GA/CI","AI/TT","TI/AT","AT/TI",
                          "TT/AI","CI/GT","GI/CT","CT/GI","GT/CI","AI/TG","TI/AG","AG/TI","TG/AI","CI/GG",
                          "GI/CG","CG/GI","GG/CI","AI/TI","TI/AI","CI/GI","GI/CI")
  colnames(DNA_IMM1) <- c("left","right")
  DNA_TMM1 <- matrix(c(-3.1,-7.8,-2.5,-6.3,-4.3,-10.7,-8,-22.5,-0.1,0.5,-0.7,-1.3,-2.1,-5.1,-3.9,-10.6,-1.1,
                       -2.1,-1.1,-2.7,-3.8,-9.5,-0.7,-19.2,-2.4,-6.5,-3.2,-8.9,-6.1,-16.9,-7.4,-21.2,-1.6,-4,
                       -1.8,-3.8,-2.6,-5.9,-2.7,-6,-5,-13.8,-3.2,-7.1,-2.3,-5.9,-2.7,-7,-0.9,-1.7,-2.3,-6.3,
                       -3.2,-8,-3.9,-10.6,-4.9,-13.5,-3,-7.8,-2.5,-6.3,-0.7,-1.2,-1.9,-4.4,-2.5,-5.9,-3.9,
                       -9.6,-6,-15.5,-4.3,-11.1,-4.6,-11.4,-2,-4.7,-2.4,-5.8,-3.2,-8.7,-3.5,-9.4,-3.8,-9,-6.6,
                       -18.7,-5.7,-15.9,-5.9,-16.1,-3.9,-10.5,-3.6,-9.8),ncol=2,byrow = TRUE)
  rownames(DNA_TMM1) <- c("AA/TA","TA/AA","CA/GA","GA/CA","AC/TC","TC/AC","CC/GC","GC/CC","AG/TG","TG/AG",
                          "CG/GG","GG/CG","AT/TT","TT/AT","CT/GT","GT/CT","AA/TC","AC/TA","CA/GC","CC/GA",
                          "GA/CC","GC/CA","TA/AC","TC/AA","AC/TT","AT/TC","CC/GT","CT/GC","GC/CT","GT/CC",
                          "TC/AT","TT/AC","AA/TG","AG/TA","CA/GG","CG/GA","GA/CG","GG/CA","TA/AG","TG/AA",
                          "AG/TT","AT/TG","CG/GT","CT/GG","GG/CT","GT/CG","TG/AT","TT/AG")
  colnames(DNA_TMM1) <- c("left","right")
  DNA_DE1 <- matrix(c(0.2,2.3,-6.3,-17.1,-3.7,-10,-2.9,-7.6,0.6,3.3,-4.4,-12.6,-4,-11.9,-4.1,-13,-1.1,-1.6,
                      -5.1,-14,-3.9,-10.9,-4.2,-15,-6.9,-20,-4,-10.9,-4.9,-13.8,-0.2,-0.5,-0.7,-0.8,-2.1,
                      -3.9,-5.9,-16.5,-0.5,-1.1,4.4,14.9,-0.2,-0.1,-2.6,-7.4,4.7,14.2,-1.6,-3.6,-3.9,-11.2,
                      -3.2,-10.4,-4.1,-13.1,2.9,10.4,-4.4,-13.1,-5.2,-15,-3.8,-12.6),ncol=2,byrow = TRUE)
  rownames(DNA_DE1) <- c("AA/.T","AC/.G","AG/.C","AT/.A","CA/.T","CC/.G","CG/.C","CT/.A","GA/.T","GC/.G",
                         "GG/.C","GT/.A","TA/.T","TC/.G","TG/.C","TT/.A",".A/AT",".C/AG",".G/AC",".T/AA",
                         ".A/CT",".C/CG",".G/CC",".T/CA",".A/GT",".C/GG",".G/GC",".T/GA",".A/TT",".C/TG",
                         ".G/TC",".T/TA")
  colnames(DNA_DE1) <- c("left","right")
  RNA_DE1 <- matrix(c(-4.9,-13.2,-0.9,-1.3,-5.5,-15.1,-2.3,-5.5,-9,-23.5,-4.1,-10.6,-8.6,-22.2,-7.5,-20.31,
                      -7.4,-20.3,-2.8,-7.7,-6.4,-16.4,-3.6,-9.7,-4.9,-13.2,-0.9,-1.3,-5.5,-15.1,-2.3,-5.5,
                      -5.7,-16.1,-0.7,-1.9,-5.8,-16.4,-2.2,-6.8,-5.7,-16.1,-0.7,-1.9,-5.8,-16.4,-2.2,-6.8,
                      -0.5,-0.6,6.9,22.6,0.6,2.6,0.6,2.6,-1.6,-4.5,0.7,3.2,-4.6,-14.8,-0.4,-1.3,-2.4,-6.1,
                      3.3,11.6,0.8,3.2,-1.4,-4.2,-0.5,-0.6,6.9,22.6,0.6,2.6,0.6,2.6,1.6,6.1,2.2,8.1,0.7,3.5,
                      3.1,10.6,1.6,6.1,2.2,8.1,0.7,3.5,3.1,10.6),ncol=2,byrow = TRUE)
  rownames(RNA_DE1) <- c(".T/AA",".T/CA",".T/GA",".T/TA",".G/AC",".G/CC",".G/GC",".G/TC",".C/AG",".C/CG",
                         ".C/GG",".C/TG",".T/AG",".T/CG",".T/GG",".T/TG",".A/AT",".A/CT",".A/GT",".A/TT",
                         ".G/AT",".G/CT",".G/GT",".G/TT","AT/.A","CT/.A","GT/.A","TT/.A","AG/.C","CG/.C",
                         "GG/.C","TG/.C","AC/.G","CC/.G","GC/.G","TC/.G","AT/.G","CT/.G","GT/.G","TT/.G",
                         "AA/.T","CA/.T","GA/.T","TA/.T","AG/.T","CG/.T","GG/.T","TG/.T")
  colnames(RNA_DE1) <- c("left","right")
  #creat list for all table data
  TableList <- list(DNA_NN1,DNA_NN2,DNA_NN3,DNA_NN4,RNA_NN1,RNA_NN2,RNA_NN3,R_DNA_NN1,DNA_IMM1,DNA_TMM1,DNA_DE1,RNA_DE1)
  names(TableList) <- c("DNA_NN1","DNA_NN2","DNA_NN3","DNA_NN4","RNA_NN1","RNA_NN2","RNA_NN3","R_DNA_NN1","DNA_IMM1","DNA_TMM1","DNA_DE1","RNA_DE1")
  
  imm_table_list <- TableList[[imm_table]]
  nn_table_list <- TableList[[nn_table]]
  tmm_table_list <- TableList[[tmm_table]]
  de_table_list <- TableList[[de_table]]
  
  imm_table_name <- rownames(imm_table_list)
  nn_table_name <- rownames(nn_table_list)
  tmm_table_name <-rownames(tmm_table_list)
  de_table_name <- rownames(de_table_list)
  
  mySeq <- check_filter(ntseq,method = "Tm_NN")
  mySeq_c2s <- c2s(mySeq)
  ptGC <- GC(mySeq,ambiguous = ambiguous)
  if(is.null(comSeq)){
    comSeq <- complement(mySeq_c2s,FALSE)
  }
  mycSeq <- check_filter(comSeq,method = "Tm_NN")
  
  tmp_seq <- mySeq
  tmp_cseq <- mycSeq
  delta_h <- 0
  delta_s <- 0
  d_h <- 1
  d_s <- 2
  if(shift!=0 | length(mySeq)!=length(mycSeq)){
    if(shift>0){
      tmp_seq <- append(rep('.',shift),mySeq)
    }else{
      tmp_cseq <- append(rep('.',abs(shift)),mycSeq)
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
  
  if(GC(mySeq) == 0){
    delta_h <- nn_table_list['init_allA/T',d_h]+delta_h
    delta_s <- nn_table_list['init_allA/T',d_s]+delta_s
  }else{
    delta_h <- nn_table_list['init_oneG/C',d_h]+delta_h
    delta_s <- nn_table_list['init_oneG/C',d_s]+delta_s
  }

  if(mySeq[1] == 'T'){
    delta_h <- nn_table_list['init_5T/A',d_h]+delta_h
    delta_s <- nn_table_list['init_5T/A',d_s]+delta_s
  }
  if(mySeq[1] == 'A'){
    delta_h <- nn_table_list['init_5T/A',d_h]+delta_h
    delta_s <- nn_table_list['init_5T/A',d_s]+delta_s
  }
  
  ends <- c(mySeq[1],mySeq[length(mySeq)])
  AT <- sum(ends %in% 'A')+sum(ends %in% 'T')
  GC <- sum(ends %in% 'G')+sum(ends %in% 'C')
  delta_h <- nn_table_list['init_A/T',d_h]*AT+delta_h
  delta_s <- nn_table_list['init_A/T',d_s]*AT+delta_s
  delta_h <- nn_table_list['init_G/C',d_h]*GC+delta_h
  delta_s <- nn_table_list['init_G/C',d_s]*GC+delta_s

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
  k <- (dnac1-(dnac2/2.0))*1e-9
  if(selfcomp==TRUE){
    k <- dnac1*1e-9
    delta_h <- nn_table_list['sym',d_h]
    delta_s <- nn_table_list['sym',d_s]
  }
  R <- 1.987
  if(!is.null(saltcorr)){
    corrSalt = salt_correction(Na=Na,K=K,Tris=Tris,Mg=Mg,dNTPs=dNTPs,method=saltcorr,ntseq=mySeq,ambiguous = ambiguous)
    if(saltcorr == "SantaLucia1998-2"){
      delta_s <- corrSalt+delta_s
    }
    Tm <- (1000*delta_h)/(delta_s+(R*(log(k))))-273.15
    if(saltcorr %in% c("Schildkraut2010","Wetmur1991","SantaLucia1996","SantaLucia1998-1")){
      Tm <- Tm+corrSalt
    }
    if(saltcorr %in% c("Owczarzy2004","Owczarzy2008")){
      Tm <- (1/(1/(Tm+273.15)+corrSalt)-273.15)
    }
  }else{
    Tm <- (1000*delta_h)/(delta_s+(R*(log(k))))-273.15
  }
  
  corrChem <- chem_correction(DMSO=DMSO,fmd=fmd,DMSOfactor=DMSOfactor,fmdmethod=fmdmethod,fmdfactor=fmdfactor,ptGC=ptGC)
  Tm <- Tm + corrChem
  
  NNTableList <- list("DNA_NN1"="Breslauer K J (1986) <doi:10.1073/pnas.83.11.3746>",
                      "DNA_NN2"="Sugimoto N (1996) <doi:10.1093/nar/24.22.4501>",
                      "DNA_NN3"="Allawi H (1998) <doi:10.1093/nar/26.11.2694>",
                      "DNA_NN4"="SantaLucia J (2004) <doi:10.1146/annurev.biophys.32.110601.141800>",
                      "RNA_NN1"="Freier S (1986) <doi:10.1073/pnas.83.24.9373>",
                      "RNA_NN2"="Xia T (1998) <doi:10.1021/bi9809425>",
                      "RNA_NN3"="Chen JL (2012) <doi:10.1021/bi3002709>",
                      "R_DNA_NN1"="Sugimoto N (1995)<doi:10.1016/S0048-9697(98)00088-6>",
                      "DNA_TMM1"="Bommarito S (2000)  <doi:10.1093/nar/28.9.1929>",
                      "DNA_IMM1"="Peyret N (1999) <doi:10.1021/bi9825091> & Allawi H T (1997) <doi:10.1021/bi962590c> & Santalucia N (2005) <doi:10.1093/nar/gki918>",
                      "DNA_DE1"="Bommarito S (2000) <doi:10.1093/nar/28.9.1929>",
                      "RNA_DE1"="Turner D H (2010) <doi:10.1093/nar/gkp892>")

  resultList <- vector('list',2L)
  names(resultList) <- c("Tm","Options")
  resultList$Tm <- as.numeric(Tm)
  resultList$Options <- list("Sequence"=ntseq,
                            "Check filter"=c2s(mySeq),
                            "Ambiguous"=ambiguous,
                            "Complement Sequence" = comSeq,  
                            "Shift"= shift, 
                            "Thermodynamic NN values" = paste0(nn_table,": ",NNTableList[[nn_table]]), 
                            "Thermodynamic values for terminal mismatches" = paste0(tmm_table,": ",NNTableList[[tmm_table]]), 
                            "Thermodynamic values for internal mismatches" = paste0(imm_table,": ",NNTableList[[imm_table]]),
                            "Thermodynamic values for dangling ends" = paste0(de_table,": ",NNTableList[[de_table]]), 
                            "Concentration of the higher concentrated strand" = dnac1,
                            "Concentration of the lower concentrated strand" = dnac2, 
                            "Sequence self-complementary" = selfcomp, 
                            "Na"=Na,
                            "K"=K,
                            "Tris"=Tris,
                            "Mg"=Mg,
                            "dNTPs"=dNTPs,
                            "Salt correlation method"=saltcorr,
                            "Percent of DMSO"=DMSO,
                            "Formamide concentration"=fmd,
                            "Coeffecient of Tm decreases per percent DMSO"=DMSOfactor,
                            "Method for formamide concentration"=fmdmethod,
                            "Coefficient of Tm decrease per percent formamide"=fmdfactor,
                            "Percent of GC"=ptGC)
  class(resultList) <- c("TmCalculator","list")
  attr(resultList, "nonhidden") <- "Tm"
  return(resultList)
}
