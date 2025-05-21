#' Thermodynamic Tables for Nucleic Acid Hybridization
#' 
#' A comprehensive collection of thermodynamic parameters used for calculating melting temperatures
#' of nucleic acid duplexes. The dataset includes parameters for DNA/DNA, RNA/RNA, and RNA/DNA
#' hybridizations, as well as parameters for mismatches and dangling ends.
#' 
#' @format A list containing 12 matrices of thermodynamic parameters:
#' \describe{
#'   \item{DNA_NN_Breslauer_1986}{DNA/DNA nearest neighbor parameters from Breslauer et al. (1986)}
#'   \item{DNA_NN_Sugimoto_1996}{DNA/DNA nearest neighbor parameters from Sugimoto et al. (1996)}
#'   \item{DNA_NN_Allawi_1998}{DNA/DNA nearest neighbor parameters from Allawi et al. (1998)}
#'   \item{DNA_NN_SantaLucia_2004}{DNA/DNA nearest neighbor parameters from SantaLucia (2004)}
#'   \item{RNA_NN_Freier_1986}{RNA/RNA nearest neighbor parameters from Freier et al. (1986)}
#'   \item{RNA_NN_Xia_1998}{RNA/RNA nearest neighbor parameters from Xia et al. (1998)}
#'   \item{RNA_NN_Chen_2012}{RNA/RNA nearest neighbor parameters from Chen et al. (2012)}
#'   \item{RNA_DNA_NN_Sugimoto_1995}{RNA/DNA nearest neighbor parameters from Sugimoto et al. (1995)}
#'   \item{DNA_IMM_Peyret_1999}{DNA internal mismatch parameters from Peyret et al. (1999)}
#'   \item{DNA_TMM_Bommarito_2000}{DNA terminal mismatch parameters from Bommarito et al. (2000)}
#'   \item{DNA_DE_Bommarito_2000}{DNA dangling end parameters from Bommarito et al. (2000)}
#'   \item{RNA_DE_Turner_2010}{RNA dangling end parameters from Turner et al. (2010)}
#' }
#' 
#' Each matrix contains thermodynamic parameters (enthalpy and entropy) for different
#' nucleic acid interactions. The parameters are used in the nearest neighbor model
#' for calculating melting temperatures of nucleic acid duplexes.
#' 
#' @source Various publications as cited in the references
#' 
#' @references
#' Breslauer K J (1986) <doi:10.1073/pnas.83.11.3746>
#' Sugimoto N (1996) <doi:10.1093/nar/24.22.4501>
#' Allawi H (1998) <doi:10.1093/nar/26.11.2694>
#' SantaLucia J (2004) <doi:10.1146/annurev.biophys.32.110601.141800>
#' Freier S (1986) <doi:10.1073/pnas.83.24.9373>
#' Xia T (1998) <doi:10.1021/bi9809425>
#' Chen JL (2012) <doi:10.1021/bi3002709>
#' Sugimoto N (1995) <doi:10.1016/S0048-9697(98)00088-6>
#' Bommarito S (2000) <doi:10.1093/nar/28.9.1929>
#' Peyret N (1999) <doi:10.1021/bi9825091>
#' Allawi H T (1997) <doi:10.1021/bi962590c>
#' Santalucia N (2005) <doi:10.1093/nar/gki918>
#' Turner D H (2010) <doi:10.1093/nar/gkp892>
#' 
#' @examples
#' # Access DNA/DNA nearest neighbor parameters
#' thermodynamic_tables$DNA_NN_SantaLucia_2004
#' 
#' # Access RNA/RNA nearest neighbor parameters
#' thermodynamic_tables$RNA_NN_Chen_2012
#' 
#' # Access DNA internal mismatch parameters
#' thermodynamic_tables$DNA_IMM_Peyret_1999
"thermodynamic_tables"

# DNA/DNA Nearest Neighbor Parameters
DNA_NN_Breslauer_1986 <- matrix(c(0,0,0,0,0,0,0,-16.8,0,-20.1,0,0,0,-1.3,-9.1,-24,-8.6,-23.9,-6,-16.9,-5.8,-12.9,
                      -6.5,-17.3,-7.8,-20.8,-5.6,-13.5,-11.9,-27.8,-11.1,-26.7,-11,-26.6),ncol=2,byrow = TRUE)
rownames(DNA_NN_Breslauer_1986) <- c("init","init_A/T","init_G/C","init_oneG/C","init_allA/T","init_5T/A","sym","AA/TT",
                         "AT/TA","TA/AT","CA/GT","GT/CA","CT/GA","GA/CT","CG/GC","GC/CG","GG/CC")
colnames(DNA_NN_Breslauer_1986) <- c("left","right")

DNA_NN_Sugimoto_1996 <- matrix(c(0.6,-9,0,0,0,0,0,0,0,0,0,0,0,-1.4,-8,-21.9,-5.6,-15.2,-6.6,-18.4,-8.2,-21,-9.4,
                      -25.5,-6.6,-16.4,-8.8,-23.5,-11.8,-29,-10.5,-26.4,-10.9,-28.4),ncol=2,byrow = TRUE)
rownames(DNA_NN_Sugimoto_1996) <- rownames(DNA_NN_Breslauer_1986)
colnames(DNA_NN_Sugimoto_1996) <- c("left","right")

DNA_NN_Allawi_1998 <- matrix(c(0,0,2.3,4.1,0.1,-2.8,0,0,0,0,0,0,0,-1.4,-7.9,-22.2,-7.2,-20.4,-7.2,-21.3,-8.5,
                      -22.7,-8.4,-22.4,-7.8,-21,-8.2,-22.2,-10.6,-27.2,-9.8,-24.4,-8,-19.9),ncol=2,byrow = TRUE)
rownames(DNA_NN_Allawi_1998) <- rownames(DNA_NN_Breslauer_1986)
colnames(DNA_NN_Allawi_1998) <- c("left","right")

DNA_NN_SantaLucia_2004 <- matrix(c(0.2,-5.7,2.2,6.9,0,0,0,0,0,0,0,0,0,-1.4,-7.6,-21.3,-7.2,-20.4,-7.2,-20.4,-8.5,
                      -22.7,-8.4,-22.4,-7.8,-21,-8.2,-22.2,-10.6,-27.2,-9.8,-24.4,-8,-19),ncol=2,byrow = TRUE)
rownames(DNA_NN_SantaLucia_2004) <- rownames(DNA_NN_Breslauer_1986)
colnames(DNA_NN_SantaLucia_2004) <- c("left","right")

# RNA/RNA Nearest Neighbor Parameters
RNA_NN_Freier_1986 <- matrix(c(0,-10.8,0,0,0,0,0,0,0,0,0,0,0,-1.4,-6.6,-18.4,-5.7,-15.5,-8.1,-22.6,-10.5,-27.8,
                      -10.2,-26.2,-7.6,-19.2,-13.3,-35.5,-8,-19.4,-14.2,-34.9,-12.2,-29.7),ncol=2,byrow = TRUE)
rownames(RNA_NN_Freier_1986) <- rownames(DNA_NN_Breslauer_1986)
colnames(RNA_NN_Freier_1986) <- c("left","right")

RNA_NN_Xia_1998 <- matrix(c(3.61,-1.5,3.72,10.5,0,0,0,0,0,0,0,0,0,-1.4,-6.82,-19,-9.38,-26.7,-7.69,-20.5,-10.44,-26.9,
                      -11.4,-29.5,-10.48,-27.1,-12.44,-32.5,-10.64,-26.7,-14.88,-36.9,-13.39,-32.7),ncol=2,byrow = TRUE)
rownames(RNA_NN_Xia_1998) <- rownames(DNA_NN_Breslauer_1986)
colnames(RNA_NN_Xia_1998) <- c("left","right")

RNA_NN_Chen_2012 <- matrix(c(6.4,6.99,3.85,11.04,0,0,0,0,0,0,0,0,0,-1.4,-7.09,-19.8,-9.11,-25.8,-8.5,-22.9,-11.03,-28.8,
                      -11.98,-31.3,-10.9,-28.5,-13.21,-34.9,-10.88,-27.4,-16.04,-40.6,-14.18,-35,-13.83,-46.9,
                      -17.82,-56.7,-3.96,-11.6,-0.96,-1.8,-10.38,-31.8,-12.64,-38.9,-7.39,-21,-5.56,-13.9,-9.44,
                      -24.7,-7.03,-16.8,-11.09,-28.8),ncol=2,byrow = TRUE)
rownames(RNA_NN_Chen_2012) <- c("init","init_A/T","init_G/C","init_oneG/C","init_allA/T","init_5T/A","sym","AA/TT",
                         "AT/TA","TA/AT","CA/GT","GT/CA","CT/GA","GA/CT","CG/GC","GC/CG","GG/CC","GT/TG",
                         "GG/TT","AG/TT","TG/AT","TT/AG","TG/GT","AT/TG","CG/GT","CT/GG","GG/CT","GT/CG")
colnames(RNA_NN_Chen_2012) <- c("left","right")

# RNA/DNA Nearest Neighbor Parameters
RNA_DNA_NN_Sugimoto_1995 <- matrix(c(1.9,-3.9,0,0,0,0,0,0,0,0,0,0,0,0,-11.5,-36.4,-7.8,-21.6,-7,-19.7,-8.3,-23.9,-10.4,-28.4,
                        -12.8,-31.9,-16.3,-47.1,-9.1,-23.5,-8.6,-22.9,-8,-17.1,-9.3,-23.2,-5.9,-12.3,-7.8,-23.2,
                        -5.5,-13.5,-9,-26.1,-7.8,-21.9),ncol=2,byrow = TRUE)
rownames(RNA_DNA_NN_Sugimoto_1995) <- c("init","init_A/T","init_G/C","init_oneG/C","init_allA/T","init_5T/A","sym","AA/TT",
                           "AC/TG","AG/TC","AT/TA","CA/GT","CC/GG","CG/GC","CT/GA","GA/CT","GC/CG","GG/CC",
                           "GT/CA","TA/AT","TC/AG","TG/AC","TT/AA")
colnames(RNA_DNA_NN_Sugimoto_1995) <- c("left","right")

# Internal Mismatch Parameters
DNA_IMM_Peyret_1999 <- matrix(c(1,0.9,-2.5,-8.3,-4.1,-11.7,-2.8,-8,3.3,10.4,5.8,16.3,-4.4,-12.3,4.1,9.5,-0.1,-1.7,-1.4,
                       -6.2,-1.3,-5.3,-0.6,-2.3,-0.7,-2.3,-0.7,-2.3,-4,-13.2,-0.6,-1,0.5,3.2,0.7,0.7,3,7.4,0.7,
                       0.2,-1.2,-6.2,-0.8,-4.5,-1.5,-6.1,2.3,5.4,5.2,13.5,1.2,0.7,1,0.7,2.3,4.6,5.3,14.6,1.9,
                       3.7,0.6,-0.6,5.2,14.2,-0.7,-3.8,3.4,8,7.6,20.2,1.2,1.7,-0.9,-4.2,-2.9,-9.8,4.7,12.9,0,
                       -4.4,-1.5,-7.2,3.6,8.9,6.1,16.4,-3.1,-9.5,-4.9,-15.3,-6,-15.8,1.6,3.6,-2.7,-10.8,-5,-15.8,
                       -2.2,-8.4,0.2,-1.5,-8.9,-25.5,-5.9,-17.4,-8.8,-25.4,-4.9,-13.9,-5.4,-13.7,-6.8,-19.1,-8.3,
                       -23.8,-5,-12.6,-8.3,-25,-3.4,-11.2,-0.7,-2.6,-1.3,-4.6,2.6,8.9,-7.8,-21.1,-7,-20,-7.6,
                       -20.2,0.49,-0.7,-6.5,-22,-5.6,-18.7,-0.8,-4.3,-1,-2.4,-3.5,-10.6,0.1,-1,-4.3,-12.1,-4.9,
                       -15.8,-1.9,-8.5,0.1,-1.8,1,1,7.1,21.3,-1.1,-3.2,5.8,16.9,-7.6,-22,-3.3,-11.9,0.1,-2.3,
                       1.3,3,-0.5,-1.3),ncol=2,byrow = TRUE)
rownames(DNA_IMM_Peyret_1999) <- c("AG/TT","AT/TG","CG/GT","CT/GG","GG/CT","GG/TT","GT/CG","GT/TG","TG/AT","TG/GT",
                          "TT/AG","AA/TG","AG/TA","CA/GG","CG/GA","GA/CG","GG/CA","TA/AG","TG/AA","AC/TT",
                          "AT/TC","CC/GT","CT/GC","GC/CT","GT/CC","TC/AT","TT/AC","AA/TC","AC/TA","CA/GC",
                          "CC/GA","GA/CC","GC/CA","TA/AC","TC/AA","AA/TA","CA/GA","GA/CA","TA/AA","AC/TC",
                          "CC/GC","GC/CC","TC/AC","AG/TG","CG/GG","GG/CG","TG/AG","AT/TT","CT/GT","GT/CT",
                          "TT/AT","AI/TC","TI/AC","AC/TI","TC/AI","CI/GC","GI/CC","CC/GI","GC/CI","AI/TA",
                          "TI/AA","AA/TI","TA/AI","CI/GA","GI/CA","CA/GI","GA/CI","AI/TT","TI/AT","AT/TI",
                          "TT/AI","CI/GT","GI/CT","CT/GI","GT/CI","AI/TG","TI/AG","AG/TI","TG/AI","CI/GG",
                          "GI/CG","CG/GI","GG/CI","AI/TI","TI/AI","CI/GI","GI/CI")
colnames(DNA_IMM_Peyret_1999) <- c("left","right")

# Terminal Mismatch Parameters
DNA_TMM_Bommarito_2000 <- matrix(c(-3.1,-7.8,-2.5,-6.3,-4.3,-10.7,-8,-22.5,-0.1,0.5,-0.7,-1.3,-2.1,-5.1,-3.9,-10.6,-1.1,
                       -2.1,-1.1,-2.7,-3.8,-9.5,-0.7,-19.2,-2.4,-6.5,-3.2,-8.9,-6.1,-16.9,-7.4,-21.2,-1.6,-4,
                       -1.8,-3.8,-2.6,-5.9,-2.7,-6,-5,-13.8,-3.2,-7.1,-2.3,-5.9,-2.7,-7,-0.9,-1.7,-2.3,-6.3,
                       -3.2,-8,-3.9,-10.6,-4.9,-13.5,-3,-7.8,-2.5,-6.3,-0.7,-1.2,-1.9,-4.4,-2.5,-5.9,-3.9,
                       -9.6,-6,-15.5,-4.3,-11.1,-4.6,-11.4,-2,-4.7,-2.4,-5.8,-3.2,-8.7,-3.5,-9.4,-3.8,-9,-6.6,
                       -18.7,-5.7,-15.9,-5.9,-16.1,-3.9,-10.5,-3.6,-9.8),ncol=2,byrow = TRUE)
rownames(DNA_TMM_Bommarito_2000) <- c("AA/TA","TA/AA","CA/GA","GA/CA","AC/TC","TC/AC","CC/GC","GC/CC","AG/TG","TG/AG",
                          "CG/GG","GG/CG","AT/TT","TT/AT","CT/GT","GT/CT","AA/TC","AC/TA","CA/GC","CC/GA",
                          "GA/CC","GC/CA","TA/AC","TC/AA","AC/TT","AT/TC","CC/GT","CT/GC","GC/CT","GT/CC",
                          "TC/AT","TT/AC","AA/TG","AG/TA","CA/GG","CG/GA","GA/CG","GG/CA","TA/AG","TG/AA",
                          "AG/TT","AT/TG","CG/GT","CT/GG","GG/CT","GT/CG","TG/AT","TT/AG")
colnames(DNA_TMM_Bommarito_2000) <- c("left","right")

# Dangling End Parameters
DNA_DE_Bommarito_2000 <- matrix(c(0.2,2.3,-6.3,-17.1,-3.7,-10,-2.9,-7.6,0.6,3.3,-4.4,-12.6,-4,-11.9,-4.1,-13,-1.1,-1.6,
                      -5.1,-14,-3.9,-10.9,-4.2,-15,-6.9,-20,-4,-10.9,-4.9,-13.8,-0.2,-0.5,-0.7,-0.8,-2.1,
                      -3.9,-5.9,-16.5,-0.5,-1.1,4.4,14.9,-0.2,-0.1,-2.6,-7.4,4.7,14.2,-1.6,-3.6,-3.9,-11.2,
                      -3.2,-10.4,-4.1,-13.1,2.9,10.4,-4.4,-13.1,-5.2,-15,-3.8,-12.6),ncol=2,byrow = TRUE)
rownames(DNA_DE_Bommarito_2000) <- c("AA/.T","AC/.G","AG/.C","AT/.A","CA/.T","CC/.G","CG/.C","CT/.A","GA/.T","GC/.G",
                         "GG/.C","GT/.A","TA/.T","TC/.G","TG/.C","TT/.A",".A/AT",".C/AG",".G/AC",".T/AA",
                         ".A/CT",".C/CG",".G/CC",".T/CA",".A/GT",".C/GG",".G/GC",".T/GA",".A/TT",".C/TG",
                         ".G/TC",".T/TA")
colnames(DNA_DE_Bommarito_2000) <- c("left","right")

RNA_DE_Turner_2010 <- matrix(c(-4.9,-13.2,-0.9,-1.3,-5.5,-15.1,-2.3,-5.5,-9,-23.5,-4.1,-10.6,-8.6,-22.2,-7.5,-20.31,
                      -7.4,-20.3,-2.8,-7.7,-6.4,-16.4,-3.6,-9.7,-4.9,-13.2,-0.9,-1.3,-5.5,-15.1,-2.3,-5.5,
                      -5.7,-16.1,-0.7,-1.9,-5.8,-16.4,-2.2,-6.8,-5.7,-16.1,-0.7,-1.9,-5.8,-16.4,-2.2,-6.8,
                      -0.5,-0.6,6.9,22.6,0.6,2.6,0.6,2.6,-1.6,-4.5,0.7,3.2,-4.6,-14.8,-0.4,-1.3,-2.4,-6.1,
                      3.3,11.6,0.8,3.2,-1.4,-4.2,-0.5,-0.6,6.9,22.6,0.6,2.6,0.6,2.6,1.6,6.1,2.2,8.1,0.7,3.5,
                      3.1,10.6,1.6,6.1,2.2,8.1,0.7,3.5,3.1,10.6),ncol=2,byrow = TRUE)
rownames(RNA_DE_Turner_2010) <- c(".T/AA",".T/CA",".T/GA",".T/TA",".G/AC",".G/CC",".G/GC",".G/TC",".C/AG",".C/CG",
                         ".C/GG",".C/TG",".T/AG",".T/CG",".T/GG",".T/TG",".A/AT",".A/CT",".A/GT",".A/TT",
                         ".G/AT",".G/CT",".G/GT",".G/TT","AT/.A","CT/.A","GT/.A","TT/.A","AG/.C","CG/.C",
                         "GG/.C","TG/.C","AC/.G","CC/.G","GC/.G","TC/.G","AT/.G","CT/.G","GT/.G","TT/.G",
                         "AA/.T","CA/.T","GA/.T","TA/.T","AG/.T","CG/.T","GG/.T","TG/.T")
colnames(RNA_DE_Turner_2010) <- c("left","right")

# Create a list of all tables for easy access
thermodynamic_tables <- list(
  DNA_NN_Breslauer_1986 = DNA_NN_Breslauer_1986,
  DNA_NN_Sugimoto_1996 = DNA_NN_Sugimoto_1996,
  DNA_NN_Allawi_1998 = DNA_NN_Allawi_1998,
  DNA_NN_SantaLucia_2004 = DNA_NN_SantaLucia_2004,
  RNA_NN_Freier_1986 = RNA_NN_Freier_1986,
  RNA_NN_Xia_1998 = RNA_NN_Xia_1998,
  RNA_NN_Chen_2012 = RNA_NN_Chen_2012,
  RNA_DNA_NN_Sugimoto_1995 = RNA_DNA_NN_Sugimoto_1995,
  DNA_IMM_Peyret_1999 = DNA_IMM_Peyret_1999,
  DNA_TMM_Bommarito_2000 = DNA_TMM_Bommarito_2000,
  DNA_DE_Bommarito_2000 = DNA_DE_Bommarito_2000,
  RNA_DE_Turner_2010 = RNA_DE_Turner_2010
)

# Export the tables
usethis::use_data(thermodynamic_tables, overwrite = TRUE) 