###
###

.pkgname <- "BSgenome.Ecoli.NCBI.ASM584v2"

.circ_seqs <- c("U00096.3")

.onLoad <- function(libname, pkgname)
{
    if (pkgname != .pkgname)
        stop("package name (", pkgname, ") is not ",
             "the expected name (", .pkgname, ")")
    extdata_dirpath <- system.file("extdata", package=pkgname,
                                   lib.loc=libname, mustWork=TRUE)

    ## Make and export BSgenome object.
    bsgenome <- BSgenome(
        organism="Escherichia coli",
        common_name=NA,
        genome="ASM584v2",
        provider="NCBI",
        release_date=NA,
        source_url=NA,
        seqnames=NULL,
        circ_seqs=.circ_seqs,
        mseqnames=NULL,
        seqs_pkgname=pkgname,
        seqs_dirpath=extdata_dirpath
    )

    ns <- asNamespace(pkgname)

    objname <- pkgname
    assign(objname, bsgenome, envir=ns)
    namespaceExport(ns, objname)

    old_objname <- "Ecoli"
    assign(old_objname, bsgenome, envir=ns)
    namespaceExport(ns, old_objname)
}

