#' Corrections of melting temperature with salt concentration
#' 
#' Apply corrections to melting temperature calculations based on salt concentrations.
#' Different correction methods are available for various experimental conditions.
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
#' @param method Method for calculating salt concentration corrections to the melting temperature.
#'   Available options:
#'   - "Schildkraut2010": Updated salt correction method
#'   - "Wetmur1991": Classic salt correction method
#'   - "SantaLucia1996": DNA-specific salt correction
#'   - "SantaLucia1998-1": Improved DNA salt correction
#'   - "SantaLucia1998-2": Alternative DNA salt correction (requires input_seq)
#'   - "Owczarzy2004": Comprehensive salt correction (requires input_seq)
#'   - "Owczarzy2008": Updated comprehensive salt correction (requires input_seq)
#'   Note: Setting to NA disables salt correction
#' 
#' @param input_seq Sequence (5' to 3') of one strand of the nucleic acid duplex. Can be provided as either:
#'   - A character string (e.g., "ATGCG")
#'   - A path to a FASTA file containing the sequence(s)
#'   Required for methods: "SantaLucia1998-2", "Owczarzy2004", and "Owczarzy2008"
#' 
#' @param ambiguous Logical. If TRUE, ambiguous bases are taken into account when computing the G and C content.
#'   The function handles various ambiguous bases (S, W, M, K, R, Y, V, H, D, B) by proportionally
#'   distributing their contribution to GC content based on their possible nucleotide compositions.
#' 
#' @details
#' 
#' Different correction methods are available for various experimental conditions:
#' 
#' - Schildkraut2010: Updated salt correction method that accounts for monovalent and divalent cations
#' - Wetmur1991: Classic salt correction method for monovalent cations
#' - SantaLucia1996: DNA-specific salt correction
#' - SantaLucia1998-1: Improved DNA salt correction
#' - SantaLucia1998-2: Alternative DNA salt correction (requires sequence information)
#' - Owczarzy2004: Comprehensive salt correction including effects of divalent cations (requires sequence information)
#' - Owczarzy2008: Updated comprehensive salt correction (requires sequence information)
#' 
#' @references
#' 
#' Schildkraut C, Lifson S. Dependence of the melting temperature of DNA on salt concentration.
#' Biopolymers. 1965;3(2):195-208.
#' 
#' Wetmur JG. DNA Probes: Applications of the Principles of Nucleic Acid Hybridization.
#' Critical Reviews in Biochemistry and Molecular Biology. 1991;26(3-4):227-259.
#' 
#' SantaLucia J. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor
#' thermodynamics. Proceedings of the National Academy of Sciences. 1998;95(4):1460-1465.
#' 
#' Owczarzy R, Moreira BG, Manthey JA, et al. Predicting stability of DNA duplexes in solutions
#' containing magnesium and monovalent cations. Biochemistry. 2008;47(19):5336-5353.
#' 
#' @author Junhui Li
#' 
#' @examples
#' 
#' salt_correct(Na = 50, Mg = 1.5, method = "Owczarzy2008", 
#'                input_seq = "ATGCGATGCG")
#' 
#' @export salt_correct

salt_correct <- function(Na=0,
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
                            input_seq,
                            ambiguous = FALSE){
  method <- match.arg(method)
  if (method %in% c("SantaLucia1998-2","Owczarzy2004","Owczarzy2008")){
    if(is.null(input_seq)){
      stop("'input_seq' should not be NULL when method is one of 'SantaLucia1998-2','Owczarzy2004','Owczarzy2008'")
    }else{
      if(length(input_seq)==1){
        input_seq <- s2c(input_seq)
      }
      mySeq <- toupper(input_seq)
      nSeq <- length(mySeq)
      ptGC <- gc(mySeq,ambiguous=ambiguous)
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
    corr <-  0
  } else {
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
  }
  return(corr)
}