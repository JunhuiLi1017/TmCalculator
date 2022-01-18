#' Corrections of melting temperature with salt ions
#' 
#' Corrections coefficient of melting temperature or entropy with different operations
#' 
#' @param Na Millimolar concentration of Na
#' 
#' @param K Millimolar concentration of K
#' 
#' @param Tris Millimolar concentration of Tris
#' 
#' @param Mg Millimolar concentration of Mg
#' 
#' @param dNTPs Millimolar concentration of dNTPs
#' 
#' @param method Method to be applied including "Schildkraut2010", "Wetmur1991","SantaLucia1996", "SantaLucia1998-1", "SantaLucia1998-2","Owczarzy2004","Owczarzy2008". First fourth methods correct Tm, fifth method corrects deltaS, sixth and seventh methods correct 1/Tm. See details for the method description.
#' 
#' @param ntseq Sequence (5' to 3') of one strand of the nucleic acid duplex as string or vector of characters.
#' 
#' @param ambiguous Ambiguous bases are taken into account to compute the G and C content when ambiguous is TRUE.
#' 
#' @details 
#' 
#' The methods are:
#'   
#' 1 Schildkraut C (2010) <doi:10.1002/bip.360030207>
#'   
#' 2 Wetmur J G (1991) <doi:10.3109/10409239109114069>
#'   
#' 3 SantaLucia J (1996) <doi:10.1021/bi951907q>
#'   
#' 4 SantaLucia J (1998) <doi:10.1073/pnas.95.4.1460>
#'   
#' 5 SantaLucia J (1998) <doi:10.1073/pnas.95.4.1460>
#'   
#' 6 Owczarzy R (2004) <doi:10.1021/bi034621r>
#'   
#' 7 Owczarzy R (2008) <doi:10.1021/bi702363u>
#' 
#' methods 1-4: Tm(new) = Tm(old) + correction
#' 
#' method 5: deltaS(new) = deltaS(old) + correction
#' 
#' methods 6+7: Tm(new) = 1/(1/Tm(old) + correction)
#' 
#' @references 
#' 
#' Schildkraut C . Dependence of the melting temperature of DNA on salt concentration[J]. Biopolymers, 2010, 3(2):195-208.
#' 
#' Wetmur J G . DNA Probes: Applications of the Principles of Nucleic Acid Hybridization[J]. CRC Critical Reviews in Biochemistry, 1991, 26(3-4):3
#' 
#' Santalucia , J , Allawi H T , Seneviratne P A . Improved Nearest-Neighbor Parameters for Predicting DNA Duplex Stability, [J]. Biochemistry, 1996, 35(11):3555-3562.
#' 
#' SantaLucia, J. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics[J]. Proceedings of the National Academy of Sciences, 1998, 95(4):1460-1465.
#' 
#' Owczarzy R , You Y , Moreira B G , et al. Effects of Sodium Ions on DNA Duplex Oligomers: Improved Predictions ofMelting Temperatures[J]. Biochemistry, 2004, 43(12):3537-3554.
#' 
#' Owczarzy R , Moreira B G , You Y , et al. Predicting Stability of DNA Duplexes in Solutions Containing Magnesium and Monovalent Cations[J]. Biochemistry, 2008, 47(19):5336-5353.
#' 
#' @author Junhui Li
#' 
#' @examples
#' 
#' ntseq <- c("acgtTGCAATGCCGTAWSDBSYXX")
#' salt_correction(Na=390, K=20, Tris=0, Mg=10, dNTPs=25, method="Owczarzy2008", ntseq)
#' 
#' @export salt_correction

salt_correction <- function(Na=0,
                            K=0,
                            Tris=0,
                            Mg=0,
                            dNTPs=0,
                            method=c("Schildkraut2010",
                                     "Wetmur1991",
                                     "SantaLucia1996",
                                     "SantaLucia1998-1",
                                     "SantaLucia1998-2",
                                     "Owczarzy2004",
                                     "Owczarzy2008"), 
                            ntseq,
                            ambiguous = FALSE){
  method <- match.arg(method)
  if (method %in% c("SantaLucia1998-2","Owczarzy2004","Owczarzy2008")){
    if(is.null(ntseq)){
      stop("'ntseq' should not be NULL when method is one of 'SantaLucia1998-2','Owczarzy2004','Owczarzy2008'")
    }else{
      if(length(ntseq)==1){
        ntseq <- s2c(ntseq)
        mySeq <- toupper(ntseq)
        nSeq <- length(mySeq)
        ptGC <- GC(mySeq,ambiguous=ambiguous)
      }
    }
  }
  if(Na < 0 | K < 0 | Tris < 0 | Mg < 0 | dNTPs < 0){
    stop("all parameters 'Na','K','Tris','Mg','dNTP' should not be less than 0")
  }
  Mon <- Na+K+Tris/2
  mg <- Mg/1000
  if (sum(c(K,Mg,Tris,dNTPs)) > 0 & method != "Owczarzy2008" & dNTPs < Mg){
    Mon <- Mon+120*sqrt(Mg-dNTPs)
  }
  mon <- Mon/1000
  if (mon == 0){
    stop("total ion concentration of zero is not allowed in this method")
  }
  corr <-  0
  if (method == "Schildkraut2010"){
    corr <- 16.6*log10(mon)
  }else if (method == "Wetmur1991"){
    corr <- 16.6*log10((mon)/(1.0+0.7*(mon)))
  }else if (method == "SantaLucia1996"){
    corr <- 12.5*log10(mon)
  }else if (method == "SantaLucia1998-1"){
    corr <- 11.7*log10(mon)
  }else if (method == "SantaLucia1998-2"){
    corr <- 0.368*(nSeq-1)*log(mon)
  }else if (method == "Owczarzy2004"){
    corr <- (4.29*ptGC/100-3.95)*1e-5*log(mon)+9.40e-6*log(mon) ^ 2
  }else if(method == "Owczarzy2008"){
    m7 <- c(3.92, -0.911, 6.26, 1.42, -48.2, 52.5, 8.31)
    dntps <- dNTPs*1e-3
    ka = 3e4
    mg <- (sqrt((ka*dntps-ka*mg+1)**2+4*ka*mg)-(ka*dntps-ka*mg+1))/(2*ka)
    R <- if (Mon > 0) sqrt(mg)/mon
    if (R < 0.22){
      corr <- (4.29*ptGC/100-3.95)*1e-5*log(mon)+9.40e-6*log(mon)**2
    }else if (R >= 0.22 && R < 6.0){
      m7[1] <- 3.92*(0.843-0.352*sqrt(mon)*log(mon))
      m7[4] <- 1.42*(1.279-4.03e-3*log(mon)-8.03e-3*log(mon)**2)
      m7[7] <- 8.31*(0.486-0.258*log(mon)+5.25e-3*log(mon)**3)
      corr <- (m7[1]+m7[2]*log(mg)+(ptGC/100)*(m7[3]+m7[4]*log(mg))+(1/(2.0*(nSeq-1))) *(m7[5]+m7[6]*log(mg)+m7[7]*log(mg)**2))*1e-5
    }
  }
  return(corr)
}
