#' Thermodynamic Tables for Nucleic Acid Hybridization
#'
#' A comprehensive collection of thermodynamic parameters used for calculating
#' melting temperatures of nucleic acid duplexes. The dataset includes
#' parameters for DNA/DNA, RNA/RNA, RNA/DNA and 2'-O-methylRNA/RNA
#' hybridizations, plus parameters for mismatches, dangling ends, internal /
#' bulge loops, GU wobble, CNG repeats, inosine, hydroxyadenine, azobenzene
#' and locked nucleic acids (LNA).
#'
#' Coverage mirrors the rmelting (Bioconductor) vignette sections 4.2.1 - 4.4.
#' Tables whose numeric values are not yet ported from primary literature are
#' included in the list as \code{NULL} placeholders so that
#' \code{\link{tm_nn}} can still dispatch on the table name; the corresponding
#' contribution is silently skipped after a one-time package-startup notice.
#'
#' @format A named list of two-column matrices with
#' \code{colnames = c("left", "right")} (left = \eqn{\Delta H},
#' right = \eqn{\Delta S}). Rownames are dinucleotide step keys such as
#' \code{"AT/TA"}, mismatch / dangling-end keys, or feature-specific keys
#' (e.g. \code{"CAG_4"} for the (CAG)4 CNG entry).
#'
#' \strong{Populated tables}
#' \describe{
#'   \item{DNA_NN_Breslauer_1986}{DNA/DNA NN, Breslauer et al. (1986)}
#'   \item{DNA_NN_Sugimoto_1996}{DNA/DNA NN, Sugimoto et al. (1996)}
#'   \item{DNA_NN_Allawi_1998}{DNA/DNA NN, Allawi (1998)}
#'   \item{DNA_NN_SantaLucia_2004}{DNA/DNA NN, SantaLucia & Hicks (2004)}
#'   \item{RNA_NN_Freier_1986}{RNA/RNA NN, Freier (1986)}
#'   \item{RNA_NN_Xia_1998}{RNA/RNA NN, Xia (1998)}
#'   \item{RNA_NN_Chen_2012}{RNA/RNA NN with GU, Chen / Serra (2012)}
#'   \item{RNA_DNA_NN_Sugimoto_1995}{RNA/DNA hybrid NN, Sugimoto (1995)}
#'   \item{DNA_IMM_Peyret_1999}{DNA single internal mismatch, Peyret (1999)}
#'   \item{DNA_TMM_Bommarito_2000}{DNA terminal mismatch, Bommarito (2000)}
#'   \item{DNA_DE_Bommarito_2000}{DNA single dangling end, Bommarito (2000)}
#'   \item{RNA_DE_Turner_2010}{RNA single dangling end, Turner (2010)}
#' }
#'
#' \strong{Registered placeholders (values TBD - populate from primary literature)}
#' \describe{
#'   \item{DNA_NN_AllSan_1997}{Allawi & SantaLucia (1997) Biochemistry 36:10581}
#'   \item{DNA_NN_SantaLucia_1996}{SantaLucia et al. (1996) Biochemistry 35:3555}
#'   \item{DNA_NN_Tanaka_2004}{Tanaka et al. (2004) Biochemistry 43:7143}
#'   \item{MeRNA_RNA_NN_Kierzek_2006}{Kierzek et al. (2006) Biochemistry 45:581}
#'   \item{DNA_RNA_IMM_Watkins_2011}{Watkins et al. (2011) Nucleic Acids Res 39:1894}
#'   \item{RNA_IMM_Lu_2006}{Lu et al. (2006) Nucleic Acids Res 34:4912}
#'   \item{RNA_IMM_Davis_Znosko_2007}{Davis & Znosko (2007) Biochemistry 46:13425}
#'   \item{RNA_IMM_Davis_Znosko_2008}{Davis & Znosko (2008) Biochemistry 47:10178}
#'   \item{DNA_TandemMM_AllSanPey}{Allawi & SantaLucia (1997/1998); Peyret (1999)}
#'   \item{RNA_TandemMM_Mathews_1999}{Mathews et al. (1999) JMB 288:911}
#'   \item{DNA_DE_Ohmichi_2002}{Ohmichi et al. (2002) JACS 124:10367}
#'   \item{RNA_DE_Ohmichi_2002}{Ohmichi et al. (2002) JACS 124:10367}
#'   \item{RNA_DE_Serra_2008}{O'Toole / Serra (2006, 2008)}
#'   \item{DNA_DDE_Ohmichi_2002}{Ohmichi et al. (2002) JACS 124:10367}
#'   \item{RNA_DDE_Ohmichi_2002}{Ohmichi et al. (2002) JACS 124:10367}
#'   \item{RNA_DDE_OToole_2005}{O'Toole et al. (2005) Biochemistry 44:14914}
#'   \item{RNA_DDE_OToole_2006}{O'Toole et al. (2006) NAR 34:3338}
#'   \item{DNA_LDE_Ohmichi_2002}{Ohmichi et al. (2002) JACS 124:10367}
#'   \item{RNA_LDE_Ohmichi_2002}{Ohmichi et al. (2002) JACS 124:10367}
#'   \item{DNA_IL_SantaLucia_2004}{SantaLucia & Hicks (2004)}
#'   \item{RNA_IL_Lu_2006}{Lu et al. (2006) NAR 34:4912}
#'   \item{RNA_IL_Badhwar_2007}{Badhwar et al. (2007) Biochemistry 46:14715}
#'   \item{DNA_SBL_Tanaka_2007}{Tan & Chen (2007) Biophys J 92:3615}
#'   \item{DNA_SBL_SantaLucia_2004}{SantaLucia & Hicks (2004)}
#'   \item{RNA_SBL_Lu_2006}{Lu et al. (2006) NAR 34:4912}
#'   \item{RNA_SBL_Blose_2007}{Blose et al. (2007) Biochemistry 46:15123}
#'   \item{DNA_LBL_SantaLucia_2004}{SantaLucia & Hicks (2004)}
#'   \item{RNA_LBL_Lu_2006}{Mathews (1999); Lu et al. (2006)}
#'   \item{RNA_GU_Mathews_1999}{Mathews & Turner (1999) JMB 288:911}
#'   \item{RNA_GU_Chen_2012}{Chen / Serra (2012)}
#'   \item{DNA_CNG_Broda_2005}{Broda et al. (2005) Biochemistry 44:13851
#'     - rownames use \code{paste0("C", N, "G_", n)}, e.g. \code{"CAG_4"}.}
#'   \item{DNA_INO_SantaLucia_2005}{Watkins & SantaLucia (2005) NAR 33:6258}
#'   \item{RNA_INO_Wright_2007}{Wright et al. (2007) Biochemistry 46:4625}
#'   \item{DNA_HA_Kawakami_2001}{Kawakami / Sugimoto (2001) Biochemistry 40:14040}
#'   \item{DNA_AZB_Asanuma_2005}{Asanuma et al. (2005) Angew Chem Int Ed 44:2747}
#'   \item{DNA_LNA_Owczarzy_2011}{Owczarzy et al. (2011) Biochemistry 50:9352}
#'   \item{DNA_LNA_McTigue_2004}{McTigue et al. (2004) Biochemistry 43:5388}
#'   \item{DNA_cLNA_Owczarzy_2011}{Owczarzy et al. (2011) Biochemistry 50:9352}
#'   \item{DNA_cLNA_MM_Owczarzy_2011}{Owczarzy et al. (2011) Biochemistry 50:9352}
#' }
#'
#' To add the numeric values for any placeholder table, construct a matrix
#' with two columns (\eqn{\Delta H} kcal/mol, \eqn{\Delta S} cal/(mol*K)) and
#' the dinucleotide-step rownames described above, then assign it:
#' \preformatted{
#' thermodynamic_nn_params$DNA_CNG_Broda_2005 <- new_table
#' }
#'
#' @source Various publications as cited above.
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
#' Allawi H T & SantaLucia J (1997) <doi:10.1021/bi962590c>
#' SantaLucia J (2005) <doi:10.1093/nar/gki918>
#' Turner D H (2010) <doi:10.1093/nar/gkp892>
#' Tanaka F (2004) Biochemistry 43:7143
#' Kierzek E (2006) Biochemistry 45:581
#' Watkins N E (2011) Nucleic Acids Res 39:1894
#' Lu Z J (2006) Nucleic Acids Res 34:4912
#' Davis A R & Znosko B M (2007, 2008) Biochemistry
#' Mathews D H (1999) <doi:10.1006/jmbi.1999.2700>
#' Ohmichi T (2002) J Am Chem Soc 124:10367
#' O'Toole A S (2005, 2006) Biochemistry / Nucleic Acids Res
#' Badhwar J (2007) Biochemistry 46:14715
#' Tan Z J & Chen S J (2007) Biophys J 92:3615
#' Blose J M (2007) Biochemistry 46:15123
#' Broda M (2005) <doi:10.1021/bi0501447>
#' Wright D J (2007) Biochemistry 46:4625
#' Kawakami J / Sugimoto N (2001) <doi:10.1021/bi010918b>
#' Asanuma H (2005) Angew Chem Int Ed 44:2747
#' McTigue P M (2004) Biochemistry 43:5388
#' Owczarzy R (2011) Biochemistry 50:9352
#'
#' @examples
#' # Access DNA/DNA nearest neighbor parameters
#' thermodynamic_nn_params$DNA_NN_SantaLucia_2004
#'
#' # Access DNA internal mismatch parameters
#' thermodynamic_nn_params$DNA_IMM_Peyret_1999
#'
#' # See which tables are still placeholders (NULL)
#' names(Filter(is.null, thermodynamic_nn_params))
"thermodynamic_nn_params"

# =============================================================================
#  POPULATED TABLES
# =============================================================================

# ---- DNA/DNA Nearest Neighbor Parameters ------------------------------------
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

# ---- RNA/RNA Nearest Neighbor Parameters ------------------------------------
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

# ---- RNA/DNA Nearest Neighbor Parameters ------------------------------------
RNA_DNA_NN_Sugimoto_1995 <- matrix(c(1.9,-3.9,0,0,0,0,0,0,0,0,0,0,0,0,-11.5,-36.4,-7.8,-21.6,-7,-19.7,-8.3,-23.9,-10.4,-28.4,
                        -12.8,-31.9,-16.3,-47.1,-9.1,-23.5,-8.6,-22.9,-8,-17.1,-9.3,-23.2,-5.9,-12.3,-7.8,-23.2,
                        -5.5,-13.5,-9,-26.1,-7.8,-21.9),ncol=2,byrow = TRUE)
rownames(RNA_DNA_NN_Sugimoto_1995) <- c("init","init_A/T","init_G/C","init_oneG/C","init_allA/T","init_5T/A","sym","AA/TT",
                           "AC/TG","AG/TC","AT/TA","CA/GT","CC/GG","CG/GC","CT/GA","GA/CT","GC/CG","GG/CC",
                           "GT/CA","TA/AT","TC/AG","TG/AC","TT/AA")
colnames(RNA_DNA_NN_Sugimoto_1995) <- c("left","right")

# ---- Internal Mismatch Parameters (DNA, Peyret 1999) ------------------------
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

# ---- Terminal Mismatch Parameters (DNA, Bommarito 2000) ---------------------
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

# ---- Dangling End Parameters (DNA, Bommarito 2000) --------------------------
DNA_DE_Bommarito_2000 <- matrix(c(0.2,2.3,-6.3,-17.1,-3.7,-10,-2.9,-7.6,0.6,3.3,-4.4,-12.6,-4,-11.9,-4.1,-13,-1.1,-1.6,
                      -5.1,-14,-3.9,-10.9,-4.2,-15,-6.9,-20,-4,-10.9,-4.9,-13.8,-0.2,-0.5,-0.7,-0.8,-2.1,
                      -3.9,-5.9,-16.5,-0.5,-1.1,4.4,14.9,-0.2,-0.1,-2.6,-7.4,4.7,14.2,-1.6,-3.6,-3.9,-11.2,
                      -3.2,-10.4,-4.1,-13.1,2.9,10.4,-4.4,-13.1,-5.2,-15,-3.8,-12.6),ncol=2,byrow = TRUE)
rownames(DNA_DE_Bommarito_2000) <- c("AA/.T","AC/.G","AG/.C","AT/.A","CA/.T","CC/.G","CG/.C","CT/.A","GA/.T","GC/.G",
                         "GG/.C","GT/.A","TA/.T","TC/.G","TG/.C","TT/.A",".A/AT",".C/AG",".G/AC",".T/AA",
                         ".A/CT",".C/CG",".G/CC",".T/CA",".A/GT",".C/GG",".G/GC",".T/GA",".A/TT",".C/TG",
                         ".G/TC",".T/TA")
colnames(DNA_DE_Bommarito_2000) <- c("left","right")

# ---- Dangling End Parameters (RNA, Turner 2010) -----------------------------
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


# =============================================================================
#  PLACEHOLDER TABLES (numeric values TBD - register from primary literature)
# =============================================================================
#
# Each new table is registered in the master list below as NULL. tm_nn()
# tolerates NULL entries: .get_tm_table() returns NULL for the table, emits a
# one-time package-startup notice, and the corresponding contribution is
# silently skipped. To populate a table, construct a two-column matrix
# (left = dH kcal/mol, right = dS cal/(mol*K)) with dinucleotide-step keys
# as rownames, then assign:
#
#   thermodynamic_nn_params$<TableName> <- new_matrix
#
# Expected key conventions:
#   - NN steps:            "AT/TA", "GC/CG", ... plus init / sym keys
#   - Mismatch steps:      "AG/TT", "CC/GA", ... (same dinuc/dinuc form)
#   - Tandem mismatches:   "AA/TT", "GG/CC" with both bases mismatched
#   - Dangling ends:       ".A/AT", "AA/.T", ".." pads on the dangle side
#   - Internal loops:      either dinucleotide keys or "IL_<size>"
#   - Bulge loops:         "SBL_<base>" (single) or "LBL_<size>" (long)
#   - CNG repeats:         "C<base>G_<n>" (e.g. "CAG_4")
#   - Inosine:             "AI/TC", "TI/AC", ... ("I" = inosine)
#   - Hydroxyadenine:      "AH/TT", "HA/TT", ... ("H" = 8-OH-A; "*" alias)
#   - Azobenzene:          step keys including "X_T" / "X_C"
#   - LNA:                 "LA/T", "LG/C", ... ("L" prefix = locked base)
#   - Consecutive LNA:     "LALA/TT", "LGLA/TC", ...
#   - Consecutive LNA+MM:  same form with one mismatched pair
# =============================================================================


# ---- Master list ------------------------------------------------------------
thermodynamic_nn_params <- list(
  # ---- Populated -----------------------------------------------------------
  DNA_NN_Breslauer_1986     = DNA_NN_Breslauer_1986,
  DNA_NN_Sugimoto_1996      = DNA_NN_Sugimoto_1996,
  DNA_NN_Allawi_1998        = DNA_NN_Allawi_1998,
  DNA_NN_SantaLucia_2004    = DNA_NN_SantaLucia_2004,
  RNA_NN_Freier_1986        = RNA_NN_Freier_1986,
  RNA_NN_Xia_1998           = RNA_NN_Xia_1998,
  RNA_NN_Chen_2012          = RNA_NN_Chen_2012,
  RNA_DNA_NN_Sugimoto_1995  = RNA_DNA_NN_Sugimoto_1995,
  DNA_IMM_Peyret_1999       = DNA_IMM_Peyret_1999,
  DNA_TMM_Bommarito_2000    = DNA_TMM_Bommarito_2000,
  DNA_DE_Bommarito_2000     = DNA_DE_Bommarito_2000,
  RNA_DE_Turner_2010        = RNA_DE_Turner_2010,

  # ---- Additional NN tables (TBD) -----------------------------------------
  DNA_NN_AllSan_1997        = NULL,   # Allawi & SantaLucia (1997) Biochemistry 36:10581
  DNA_NN_SantaLucia_1996    = NULL,   # SantaLucia et al. (1996) Biochemistry 35:3555
  DNA_NN_Tanaka_2004        = NULL,   # Tanaka et al. (2004) Biochemistry 43:7143
  MeRNA_RNA_NN_Kierzek_2006 = NULL,   # Kierzek et al. (2006) Biochemistry 45:581

  # ---- Single (internal) mismatches (4.2.3) -------------------------------
  DNA_RNA_IMM_Watkins_2011  = NULL,   # Watkins et al. (2011) NAR 39:1894
  RNA_IMM_Lu_2006           = NULL,   # Lu et al. (2006) NAR 34:4912
  RNA_IMM_Davis_Znosko_2007 = NULL,   # Davis & Znosko (2007) Biochemistry 46:13425
  RNA_IMM_Davis_Znosko_2008 = NULL,   # Davis & Znosko (2008) Biochemistry 47:10178

  # ---- Tandem mismatches (4.2.4) ------------------------------------------
  DNA_TandemMM_AllSanPey    = NULL,   # Allawi & SantaLucia (1997/1998); Peyret (1999)
  RNA_TandemMM_Mathews_1999 = NULL,   # Mathews et al. (1999) JMB 288:911

  # ---- Single dangling end (4.2.5) - additional sources -------------------
  DNA_DE_Ohmichi_2002       = NULL,   # Ohmichi et al. (2002) JACS 124:10367
  RNA_DE_Ohmichi_2002       = NULL,   # Ohmichi et al. (2002) JACS 124:10367
  RNA_DE_Serra_2008         = NULL,   # O'Toole / Serra (2006, 2008)

  # ---- Double dangling end (4.2.6) ----------------------------------------
  DNA_DDE_Ohmichi_2002      = NULL,   # Ohmichi et al. (2002) JACS 124:10367
  RNA_DDE_Ohmichi_2002      = NULL,   # Ohmichi et al. (2002) JACS 124:10367
  RNA_DDE_OToole_2005       = NULL,   # O'Toole et al. (2005) Biochemistry 44:14914
  RNA_DDE_OToole_2006       = NULL,   # O'Toole et al. (2006) NAR 34:3338

  # ---- Long dangling end (4.2.7) ------------------------------------------
  DNA_LDE_Ohmichi_2002      = NULL,   # Ohmichi et al. (2002) JACS 124:10367
  RNA_LDE_Ohmichi_2002      = NULL,   # Ohmichi et al. (2002) JACS 124:10367

  # ---- Internal loop (4.2.8) ----------------------------------------------
  DNA_IL_SantaLucia_2004    = NULL,   # SantaLucia & Hicks (2004)
  RNA_IL_Lu_2006            = NULL,   # Lu et al. (2006) NAR 34:4912
  RNA_IL_Badhwar_2007       = NULL,   # Badhwar et al. (2007) Biochemistry 46:14715

  # ---- Single bulge loop (4.2.9) ------------------------------------------
  DNA_SBL_Tanaka_2007       = NULL,   # Tan & Chen (2007) Biophys J 92:3615
  DNA_SBL_SantaLucia_2004   = NULL,   # SantaLucia & Hicks (2004)
  RNA_SBL_Lu_2006           = NULL,   # Lu et al. (2006) NAR 34:4912
  RNA_SBL_Blose_2007        = NULL,   # Blose et al. (2007) Biochemistry 46:15123

  # ---- Long bulge loop (4.2.10) -------------------------------------------
  DNA_LBL_SantaLucia_2004   = NULL,   # SantaLucia & Hicks (2004)
  RNA_LBL_Lu_2006           = NULL,   # Mathews (1999); Lu et al. (2006)

  # ---- GU wobble (4.2.2) --------------------------------------------------
  RNA_GU_Mathews_1999       = NULL,   # Mathews & Turner (1999) JMB 288:911
  RNA_GU_Chen_2012          = NULL,   # Chen / Serra (2012) Biochemistry 51:3508

  # ---- CNG triplet repeats (4.2.11) ---------------------------------------
  DNA_CNG_Broda_2005        = NULL,   # Broda et al. (2005) Biochemistry 44:13851

  # ---- Inosine (4.2.12) ---------------------------------------------------
  DNA_INO_SantaLucia_2005   = NULL,   # Watkins & SantaLucia (2005) NAR 33:6258
  RNA_INO_Wright_2007       = NULL,   # Wright et al. (2007) Biochemistry 46:4625

  # ---- Hydroxyadenine (4.2.13) --------------------------------------------
  DNA_HA_Kawakami_2001      = NULL,   # Kawakami / Sugimoto (2001) Biochemistry 40:14040

  # ---- Azobenzene (4.2.14) ------------------------------------------------
  DNA_AZB_Asanuma_2005      = NULL,   # Asanuma et al. (2005) Angew Chem Int Ed 44:2747

  # ---- Single LNA (4.2.15) ------------------------------------------------
  DNA_LNA_Owczarzy_2011     = NULL,   # Owczarzy et al. (2011) Biochemistry 50:9352
  DNA_LNA_McTigue_2004      = NULL,   # McTigue et al. (2004) Biochemistry 43:5388

  # ---- Consecutive LNA (4.3) ----------------------------------------------
  DNA_cLNA_Owczarzy_2011    = NULL,   # Owczarzy et al. (2011) Biochemistry 50:9352

  # ---- Consecutive LNA + single MM (4.4) ----------------------------------
  DNA_cLNA_MM_Owczarzy_2011 = NULL    # Owczarzy et al. (2011) Biochemistry 50:9352
)

# Export the tables
save(thermodynamic_nn_params, file = "data/thermodynamic_nn_params.RData", version = 2)
