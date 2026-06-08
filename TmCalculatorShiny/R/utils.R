# Internal helpers for TmCalculatorShiny -------------------------------------

#' NULL-coalescing operator
#' @noRd
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

#' Nearest-neighbour thermodynamic tables offered in the UI
#' @noRd
.tcs_nn_tables <- function() {
  c("DNA_NN_SantaLucia_2004",
    "DNA_NN_Breslauer_1986",
    "DNA_NN_Sugimoto_1996",
    "DNA_NN_Allawi_1998",
    "RNA_NN_Freier_1986",
    "RNA_NN_Xia_1998",
    "RNA_NN_Chen_2012",
    "RNA_DNA_NN_Sugimoto_1995")
}

#' GC empirical variants offered in the UI
#' @noRd
.tcs_gc_variants <- function() {
  c("Primer3Plus",
    "Chester1993",
    "QuikChange",
    "Schildkraut1965",
    "Wetmur1991_MELTING",
    "Wetmur1991_RNA",
    "Wetmur1991_RNA/DNA",
    "vonAhsen2001")
}

#' Salt-correction methods offered in the UI
#' @noRd
.tcs_salt_methods <- function() {
  c("Schildkraut2010",
    "Wetmur1991",
    "SantaLucia1996",
    "SantaLucia1998-1",
    "Owczarzy2004",
    "Owczarzy2008")
}

#' Discrete colour palettes for plots
#' @noRd
.tcs_palettes <- function() {
  list(
    Set1   = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
               "#ffff33", "#a65628", "#f781bf"),
    Dark2  = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
               "#e6ab02", "#a6761d", "#666666"),
    Viridis = c("#440154", "#3b528b", "#21908c", "#5dc963", "#fde725")
  )
}

#' Is a package installed (without attaching it)?
#' @noRd
.tcs_has_pkg <- function(pkg) {
  isTRUE(requireNamespace(pkg, quietly = TRUE))
}

#' Names of installed BSgenome data packages
#' @noRd
.tcs_installed_genomes <- function() {
  pkgs <- tryCatch(rownames(utils::installed.packages()),
                   error = function(e) character(0))
  genomes <- grep("^BSgenome\\.", pkgs, value = TRUE)
  sort(genomes)
}

#' Coerce a stored result to a data.frame for display
#' @noRd
.tcs_result_df <- function(x) {
  if (is.null(x)) return(NULL)
  gr <- .tcs_result_gr(x)
  if (is.null(gr)) return(NULL)
  as.data.frame(gr)
}

#' Coerce a stored result entry to a GRanges
#' @noRd
.tcs_result_gr <- function(x) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "GRanges")) return(x)
  if (is.list(x) && !is.null(x$gr)) return(x$gr)
  NULL
}

#' Parse a free-text sequence box into a character vector
#'
#' Accepts one sequence per line, optionally with a leading FASTA-style
#' \code{>name} header that is used as the region label.
#' @noRd
.tcs_parse_sequence_box <- function(txt) {
  if (is.null(txt) || !nzchar(trimws(txt))) return(NULL)
  lines <- strsplit(txt, "\r?\n")[[1]]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) return(NULL)

  if (any(startsWith(lines, ">"))) {
    seqs <- character(0)
    names_v <- character(0)
    cur_name <- NA_character_
    cur_seq <- ""
    flush <- function() {
      if (nzchar(cur_seq)) {
        seqs[[length(seqs) + 1L]] <<- cur_seq
        names_v[[length(names_v) + 1L]] <<- cur_name %||% paste0("seq", length(seqs))
      }
    }
    for (ln in lines) {
      if (startsWith(ln, ">")) {
        flush()
        cur_name <- trimws(sub("^>", "", ln))
        cur_seq <- ""
      } else {
        cur_seq <- paste0(cur_seq, gsub("\\s+", "", ln))
      }
    }
    flush()
    if (length(seqs) == 0) return(NULL)
    out <- toupper(seqs)
    names(out) <- names_v
    return(out)
  }

  # Plain sequences, one per line
  out <- toupper(gsub("\\s+", "", lines))
  names(out) <- paste0("seq", seq_along(out))
  out
}

#' Read a simple BED-like / CSV feature table into a GRanges
#'
#' Expects columns that can be mapped to seqnames/start/end; remaining
#' columns become metadata. Supports tab- or comma-separated files with a
#' header, or headerless 3+ column BED.
#' @noRd
.tcs_read_features <- function(path) {
  first <- readLines(path, n = 1L, warn = FALSE)
  sep <- if (grepl("\t", first)) "\t" else ","
  has_header <- grepl("[A-Za-z]", strsplit(first, sep)[[1]][1]) &&
    !grepl("^chr", tolower(first))
  df <- utils::read.csv(path, sep = sep, header = TRUE,
                        stringsAsFactors = FALSE, check.names = TRUE)
  nms <- tolower(names(df))
  chr_i <- which(nms %in% c("seqnames", "chrom", "chr", "chromosome"))[1]
  start_i <- which(nms %in% c("start", "chromstart", "begin"))[1]
  end_i <- which(nms %in% c("end", "chromend", "stop"))[1]

  if (is.na(chr_i) || is.na(start_i) || is.na(end_i)) {
    # Assume headerless BED: col1=chr, col2=start(0-based), col3=end
    df <- utils::read.csv(path, sep = sep, header = FALSE,
                          stringsAsFactors = FALSE)
    if (ncol(df) < 3)
      stop("Feature file needs at least chrom/start/end columns.")
    chr <- as.character(df[[1]])
    start <- as.integer(df[[2]]) + 1L  # BED is 0-based half-open
    end <- as.integer(df[[3]])
    meta <- if (ncol(df) > 3) df[, 4:ncol(df), drop = FALSE] else NULL
    if (!is.null(meta)) names(meta) <- paste0("V", 4:ncol(df))
  } else {
    chr <- as.character(df[[chr_i]])
    start <- as.integer(df[[start_i]])
    end <- as.integer(df[[end_i]])
    drop <- c(chr_i, start_i, end_i)
    strand_i <- which(nms == "strand")[1]
    meta_cols <- setdiff(seq_len(ncol(df)), c(drop, strand_i[!is.na(strand_i)]))
    meta <- if (length(meta_cols)) df[, meta_cols, drop = FALSE] else NULL
  }

  gr <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = start, end = end)
  )
  if (!is.null(meta) && ncol(meta) > 0) {
    for (cn in names(meta)) {
      v <- meta[[cn]]
      num <- suppressWarnings(as.numeric(v))
      GenomicRanges::mcols(gr)[[cn]] <- if (all(is.na(num) == is.na(v))) num else v
    }
  }
  gr
}
