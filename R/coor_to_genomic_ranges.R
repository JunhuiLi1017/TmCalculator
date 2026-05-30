#' Convert genomic coordinate strings to a GRanges object
#'
#' Fast conversion of genomic coordinate strings to a \code{GRanges} object with
#' reference sequences fetched from installed \code{BSgenome.*} packages.
#' Designed for large sliding-window inputs: genome packages are loaded once,
#' coordinate strings are parsed in one pass, and sequences are extracted with
#' vectorized \code{\link[Biostrings]{getSeq}} (or optional chromosome
#' preloading).
#'
#' @param input Coordinate input. Either:
#'   \itemize{
#'     \item A list with \code{pkg_name} (BSgenome package name) and \code{seq}
#'       (character vector of coordinate strings). Preferred for many windows
#'       on the same genome.
#'     \item A plain character vector of coordinate strings (legacy format;
#'       genome package name is read from field 4 when present).
#'   }
#'   Supported colon-separated formats:
#'   \itemize{
#'     \item \code{chr:start-end:strand:region_id} - requires \code{pkg_name}
#'       in the input list.
#'     \item \code{chr:start-end:strand:pkg_name:region_id} - as produced by
#'       \code{\link{make_genomiccoord}}.
#'   }
#' @param complement_seq Optional complement coordinates in the same format as
#'   \code{input$seq}. When \code{NULL}, complements are generated automatically
#'   from \code{sequence}.
#' @param method Sequence extraction strategy:
#'   \describe{
#'     \item{\code{"vectorized"}}{One \code{getSeq()} call per genome package
#'       (default).}
#'     \item{\code{"preload_chr"}}{Load each chromosome once and extract windows
#'       with \code{subseq()}. Faster for dense whole-chromosome tiling; uses
#'       more memory.}
#'   }
#'
#' @return A \code{GRanges} object with metadata columns:
#'   \describe{
#'     \item{\code{sequence}}{Reference sequence for each interval.}
#'     \item{\code{complement}}{Complementary sequence.}
#'     \item{\code{GC}}{GC fraction (0-1) per interval.}
#'     \item{\code{region_id}}{Region identifier from the coordinate string.}
#'     \item{\code{genome_pkg}}{BSgenome package name used.}
#'   }
#'
#' @details
#' For genome-wide tiling with thousands of windows, pass coordinates as
#' \code{list(pkg_name = "BSgenome.Hsapiens.UCSC.hg38", seq = ...)} so the
#' genome package is loaded once instead of per interval.
#'
#' @examples
#' \dontrun{
#' coords <- c(
#'   "chr1:1000-1199:+:win1",
#'   "chr1:1200-1399:+:win2"
#' )
#' gr <- coor_to_genomic_ranges(
#'   list(pkg_name = "BSgenome.Hsapiens.UCSC.hg38", seq = coords)
#' )
#' gr
#' }
#'
#' @seealso \code{\link{to_genomic_ranges}}, \code{\link{to_genomic_ranges_fast}}
#'
#' @importFrom GenomicRanges GRanges seqnames start end width strand mcols
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlengths seqlevels
#' @importFrom Biostrings getSeq letterFrequency complement reverseComplement
#'   DNAStringSet DNAString subseq
#' @importFrom S4Vectors mcols
#' @encoding UTF-8
#' @author Junhui Li
#' @export

coor_to_genomic_ranges <- function(
    input,
    complement_seq = NULL,
    method         = c("vectorized", "preload_chr")
) {
  method <- match.arg(method)

  if (is.list(input)) {
    if (!all(c("pkg_name", "seq") %in% names(input)))
      stop("'input' list must have elements 'pkg_name' and 'seq'.")
    pkg_name  <- input$pkg_name
    coord_vec <- input$seq
  } else if (is.character(input)) {
    coord_vec <- input
    pkg_name  <- NULL
  } else {
    stop("'input' must be a list(pkg_name=, seq=) or a character vector.")
  }

  n <- length(coord_vec)
  if (n == 0L) stop("No coordinate strings provided.")

  parsed <- .parse_coord_strings(coord_vec, pkg_name_override = pkg_name)

  unique_pkgs <- unique(parsed$pkg_name)
  genome_objs <- .load_genome_packages(unique_pkgs)

  gr_query <- GenomicRanges::GRanges(
    seqnames = parsed$chr,
    ranges   = IRanges::IRanges(start = parsed$win_start, end = parsed$win_end),
    strand   = parsed$strand
  )
  names(gr_query) <- parsed$region_id

  primary_genome <- genome_objs[[unique_pkgs[1]]]
  sl <- GenomeInfoDb::seqlengths(primary_genome)
  shared_chrs <- intersect(GenomeInfoDb::seqlevels(gr_query), names(sl))
  GenomeInfoDb::seqlengths(gr_query)[shared_chrs] <- sl[shared_chrs]

  message(sprintf(
    "[coor_to_genomic_ranges] Extracting sequences for %s intervals (method='%s') ...",
    format(n, big.mark = ","), method
  ))
  t0 <- proc.time()

  if (method == "preload_chr") {
    seqs <- .getseq_preload_chr(gr_query, genome_objs, parsed)
  } else {
    seqs <- .getseq_vectorized(gr_query, genome_objs, parsed)
  }

  elapsed <- (proc.time() - t0)[["elapsed"]]
  message(sprintf(
    "[coor_to_genomic_ranges] Done in %.1f sec (%.2f ms/interval)",
    elapsed, elapsed / n * 1000
  ))

  gr_out <- gr_query
  S4Vectors::mcols(gr_out)$sequence   <- seqs
  S4Vectors::mcols(gr_out)$region_id  <- parsed$region_id
  S4Vectors::mcols(gr_out)$genome_pkg <- parsed$pkg_name

  gc_counts <- Biostrings::letterFrequency(seqs, letters = "GC", as.prob = FALSE)[, 1]
  win_widths <- parsed$win_end - parsed$win_start + 1L
  S4Vectors::mcols(gr_out)$GC <- round(gc_counts / win_widths, 6)

  if (!is.null(complement_seq)) {
    comp_input <- if (is.list(complement_seq)) {
      complement_seq
    } else {
      list(pkg_name = pkg_name %||% parsed$pkg_name[1], seq = complement_seq)
    }
    comp_gr <- coor_to_genomic_ranges(comp_input, method = method)
    S4Vectors::mcols(gr_out)$complement <- S4Vectors::mcols(comp_gr)$sequence
  } else {
    S4Vectors::mcols(gr_out)$complement <- Biostrings::complement(seqs)
  }

  gr_out
}


#' Parse coordinate strings into a data frame
#'
#' @param coord_vec Character vector of colon-separated coordinate strings.
#' @param pkg_name_override Optional BSgenome package name for 4-field strings.
#' @return Data frame with columns \code{chr}, \code{win_start}, \code{win_end},
#'   \code{strand}, \code{pkg_name}, and \code{region_id}.
#' @keywords internal
.parse_coord_strings <- function(coord_vec, pkg_name_override = NULL) {
  parts_list <- strsplit(coord_vec, ":", fixed = TRUE)
  n_fields   <- lengths(parts_list)

  bad <- which(!n_fields %in% c(4L, 5L))
  if (length(bad) > 0) {
    stop(sprintf(
      "%d coordinate string(s) have unexpected format (must have 4 or 5 colon-separated fields).\nFirst bad string: '%s'",
      length(bad), coord_vec[bad[1]]
    ))
  }

  chr <- vapply(parts_list, `[`, "", 1L)

  loc_raw   <- vapply(parts_list, `[`, "", 2L)
  loc_parts <- strsplit(loc_raw, "-", fixed = TRUE)
  win_start <- as.integer(vapply(loc_parts, `[`, "", 1L))
  win_end   <- as.integer(vapply(loc_parts, `[`, "", 2L))

  bad_coords <- which(is.na(win_start) | is.na(win_end) | win_start > win_end)
  if (length(bad_coords) > 0) {
    stop(sprintf(
      "Invalid start/end in %d coordinate string(s). First: '%s'",
      length(bad_coords), coord_vec[bad_coords[1]]
    ))
  }

  strand <- vapply(parts_list, `[`, "", 3L)
  bad_strand <- which(!strand %in% c("+", "-", "*"))
  if (length(bad_strand) > 0) {
    warning(sprintf(
      "%d strand values not in {+, -, *}. Defaulting to '*'.",
      length(bad_strand)
    ))
  }
  strand[!strand %in% c("+", "-", "*")] <- "*"

  pkg_name  <- character(length(coord_vec))
  region_id <- character(length(coord_vec))

  idx4 <- which(n_fields == 4L)
  idx5 <- which(n_fields == 5L)

  if (length(idx4) > 0) {
    if (is.null(pkg_name_override)) {
      stop("4-field coordinate strings require 'pkg_name' to be specified in the input list.")
    }
    pkg_name[idx4]  <- pkg_name_override
    region_id[idx4] <- vapply(parts_list[idx4], `[`, "", 4L)
  }

  if (length(idx5) > 0) {
    extracted_pkgs  <- vapply(parts_list[idx5], `[`, "", 4L)
    pkg_name[idx5]  <- if (!is.null(pkg_name_override)) pkg_name_override else extracted_pkgs
    region_id[idx5] <- vapply(parts_list[idx5], `[`, "", 5L)
  }

  data.frame(
    chr       = chr,
    win_start = win_start,
    win_end   = win_end,
    strand    = strand,
    pkg_name  = pkg_name,
    region_id = region_id,
    stringsAsFactors = FALSE
  )
}


#' Load installed BSgenome packages
#'
#' @param pkg_names Character vector of BSgenome package names.
#' @return Named list of genome objects.
#' @keywords internal
.load_genome_packages <- function(pkg_names) {
  genome_objs <- list()
  for (pkg in pkg_names) {
    genome_objs[[pkg]] <- .get_bsgenome_from_pkg(pkg)
    message(sprintf("[coor_to_genomic_ranges] Loaded genome: %s", pkg))
  }
  genome_objs
}


#' Load the BSgenome object from an installed BSgenome.* data package
#'
#' The genome object name is given by the \code{BSgenomeObjname} field in
#' \code{DESCRIPTION} (e.g. \code{Hsapiens}), not necessarily the package name.
#'
#' @param pkg Character scalar: installed BSgenome package name.
#' @return A \code{BSgenome} object.
#' @keywords internal
.get_bsgenome_from_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(
      "BSgenome package '%s' is not installed.\nInstall with: BiocManager::install('%s')",
      pkg, pkg
    ), call. = FALSE)
  }
  suppressPackageStartupMessages(
    suppressWarnings(library(pkg, character.only = TRUE))
  )
  objname <- utils::packageDescription(pkg, fields = "BSgenomeObjname")
  if (is.na(objname) || !nzchar(objname)) objname <- pkg
  ns <- asNamespace(pkg)
  if (exists(objname, envir = ns, inherits = FALSE)) {
    return(get(objname, envir = ns))
  }
  stop(sprintf(
    paste0(
      "BSgenome package '%s' is installed but genome object '%s' is missing.\n",
      "Re-forge and reinstall the data package (see BSgenomeForge::forgeBSgenomeDataPkgFromNCBI)."
    ),
    pkg, objname
  ), call. = FALSE)
}


#' Vectorized sequence extraction with getSeq
#'
#' @param gr_query Query \code{GRanges}.
#' @param genome_objs Named list of BSgenome objects.
#' @param parsed Parsed coordinate data frame.
#' @return \code{DNAStringSet} of extracted sequences.
#' @keywords internal
.getseq_vectorized <- function(gr_query, genome_objs, parsed) {
  n    <- length(gr_query)
  seqs <- Biostrings::DNAStringSet(rep("", n))

  for (pkg in names(genome_objs)) {
    genome  <- genome_objs[[pkg]]
    idx_pkg <- which(parsed$pkg_name == pkg)
    if (length(idx_pkg) == 0L) next

    gr_pkg <- gr_query[idx_pkg]
    seqs_pkg <- tryCatch(
      Biostrings::getSeq(genome, gr_pkg),
      error = function(e) {
        message(sprintf("  ERROR in getSeq for %s: %s", pkg, conditionMessage(e)))
        Biostrings::DNAStringSet(rep(NA_character_, length(gr_pkg)))
      }
    )
    seqs[idx_pkg] <- seqs_pkg
  }

  seqs
}


#' Sequence extraction by preloading whole chromosomes
#'
#' @param gr_query Query \code{GRanges}.
#' @param genome_objs Named list of BSgenome objects.
#' @param parsed Parsed coordinate data frame.
#' @return \code{DNAStringSet} of extracted sequences.
#' @keywords internal
.getseq_preload_chr <- function(gr_query, genome_objs, parsed) {
  n    <- length(gr_query)
  seqs <- vector("list", n)

  for (pkg in names(genome_objs)) {
    genome  <- genome_objs[[pkg]]
    idx_pkg <- which(parsed$pkg_name == pkg)
    if (length(idx_pkg) == 0L) next

    chrs_needed <- unique(parsed$chr[idx_pkg])

    for (chr in chrs_needed) {
      idx_chr <- idx_pkg[parsed$chr[idx_pkg] == chr]

      message(sprintf(
        "  Preloading %s from %s (%s intervals) ...",
        chr, pkg, format(length(idx_chr), big.mark = ",")
      ))

      chr_seq <- tryCatch(
        genome[[chr]],
        error = function(e) {
          message(sprintf("  ERROR loading %s: %s", chr, conditionMessage(e)))
          return(NULL)
        }
      )
      if (is.null(chr_seq)) {
        seqs[idx_chr] <- list(Biostrings::DNAString(""))
        next
      }

      chr_len  <- length(chr_seq)
      starts_i <- parsed$win_start[idx_chr]
      ends_i   <- parsed$win_end[idx_chr]
      starts_i <- pmax(1L, starts_i)
      ends_i   <- pmin(chr_len, ends_i)

      seqs_chr <- Biostrings::subseq(chr_seq, start = starts_i, end = ends_i)

      minus_idx <- which(parsed$strand[idx_chr] == "-")
      if (length(minus_idx) > 0) {
        seqs_chr[minus_idx] <- Biostrings::reverseComplement(seqs_chr[minus_idx])
      }

      seqs[idx_chr] <- as.list(seqs_chr)
      rm(chr_seq)
      base::gc(verbose = FALSE)
    }
  }

  Biostrings::DNAStringSet(unlist(seqs))
}


`%||%` <- function(a, b) if (!is.null(a)) a else b


#' Convert input sequences to a GRanges object (fast backend)
#'
#' Drop-in replacement for \code{\link{to_genomic_ranges}} that uses the fast
#' coordinate backend in \code{\link{coor_to_genomic_ranges}} for genomic
#' coordinate input. Accepts the same \code{input_seq} and \code{complement_seq}
#' arguments as \code{\link{to_genomic_ranges}}, plus a list input
#' \code{list(pkg_name = ..., seq = ...)} for large tiling jobs.
#'
#' @param method Sequence extraction method passed to
#'   \code{\link{coor_to_genomic_ranges}}. One of \code{"vectorized"} (default)
#'   or \code{"preload_chr"}.
#'
#' @return A \code{GRanges} object. See \code{\link{coor_to_genomic_ranges}} for
#'   metadata columns when coordinate input is used.
#'
#' @examples
#' \dontrun{
#' gr <- to_genomic_ranges_fast(
#'   list(
#'     pkg_name = "BSgenome.Hsapiens.UCSC.hg38",
#'     seq = c("chr1:1000-1199:+:win1", "chr1:1200-1399:+:win2")
#'   )
#' )
#' }
#'
#' @rdname to_genomic_ranges
#' @export
to_genomic_ranges_fast <- function(
    input_seq,
    complement_seq = NULL,
    method = c("vectorized", "preload_chr")
) {
  method <- match.arg(method)

  if (is.null(input_seq) || length(input_seq) == 0) {
    stop("Input sequence cannot be NULL or empty")
  }

  if (is.list(input_seq) && all(c("pkg_name", "seq") %in% names(input_seq))) {
    input_gr <- coor_to_genomic_ranges(input_seq, method = method)
  } else if (is.character(input_seq) && length(input_seq) == 1 && file.exists(input_seq)) {
    input_gr <- fa_to_genomic_ranges(input_seq)
  } else if (is.character(input_seq) && all(grepl(":", input_seq, fixed = TRUE))) {
    input_gr <- coor_to_genomic_ranges(input_seq, method = method)
  } else if (is.character(input_seq)) {
    input_gr <- vec_to_genomic_ranges(input_seq)
  } else {
    stop("Unrecognised input format.")
  }

  if (!is.null(complement_seq)) {
    comp_gr <- to_genomic_ranges_fast(complement_seq, method = method)
    S4Vectors::mcols(input_gr)$complement <- S4Vectors::mcols(comp_gr)$sequence
  } else if (!"complement" %in% names(S4Vectors::mcols(input_gr))) {
    S4Vectors::mcols(input_gr)$complement <- Biostrings::complement(
      S4Vectors::mcols(input_gr)$sequence
    )
  }

  input_gr
}
