\name{BSgenome.Ecoli.NCBI.ASM584v2}
\docType{package}

\alias{BSgenome.Ecoli.NCBI.ASM584v2-package}
\alias{BSgenome.Ecoli.NCBI.ASM584v2}
\alias{Ecoli}

\title{Full genomic sequences for Escherichia coli (NCBI assembly ASM584v2)}

\description{
  Full genomic sequences for Escherichia coli as provided by NCBI (assembly ASM584v2, assembly accession GCF_000005845.2). The sequences are stored in DNAString objects.
}

\author{Junhui Li <ljh.biostat@gmail.com>}

\seealso{
  \itemize{
    \item \link[BSgenome]{BSgenome} objects in the \pkg{BSgenome}
          software package.

    \item The \code{\link[GenomeInfoDb]{seqinfo}} getter and
          \link[GenomeInfoDb]{Seqinfo} objects in the \pkg{GenomeInfoDb}
          package.

    \item The \code{\link[GenomeInfoDb]{seqlevelsStyle}} getter and
          setter in the \pkg{GenomeInfoDb} package.

    \item \link[Biostrings]{DNAString} objects in the \pkg{Biostrings}
          package.

    \item The \code{\link[BSgenome]{available.genomes}} function
          in the \pkg{BSgenome} software package.

    \item The \link[BSgenomeForge]{forgeBSgenomeDataPkgFromNCBI}
          function in the \pkg{BSgenomeForge} package that was used to
          create this package.
  }
}

\examples{
BSgenome.Ecoli.NCBI.ASM584v2
bsg <- BSgenome.Ecoli.NCBI.ASM584v2
head(seqlengths(bsg))
seqinfo(bsg)

## ---------------------------------------------------------------------
## Genome-wide motif searching
## ---------------------------------------------------------------------

## See the GenomeSearching vignette in the BSgenome software
## package for some examples of genome-wide motif searching using
## Biostrings and the BSgenome data packages:
if (interactive())
    vignette("GenomeSearching", package="BSgenome")
}

\keyword{package}
\keyword{data}
