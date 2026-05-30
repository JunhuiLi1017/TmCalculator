#' Filter invalid bases in nucleotide sequences
#' 
#' This function processes nucleotide sequences by converting characters to uppercase and replacing invalid bases with "". 
#' based on the specified method. The function preserves the sequence length and attributes (name and Tm) of each sequence.
#'
#' @param seq_list Input sequence in 5' to 3' direction. Must be provided as:
#'   - A list of sequences with attributes (name and Tm)
#'   
#' @param method Method to determine valid bases:
#' 
#' TM_Wallace: Valid bases are "A","B","C","D","G","H","I","K","M","N","R","S","T","V","W" and "Y"
#' 
#' TM_GC: Valid bases are "A","B","C","D","G","H","I","K","M","N","R","S","T","V","W", "X" and "Y"
#' 
#' TM_NN: Valid bases are "A","C","G","I" and "T" 
#' 
#' @returns Returns a list of sequences with the same structure as input, where invalid bases are replaced with ""
#' 
#' @author Junhui Li
#' 
#' @references 
#' 
#' \code{citation("TmCalculator")}
#' 

check_filter_seq <- function(seq_list, method) {
  if (is.null(seq_list) || length(seq_list) == 0) {
    stop("Input sequence list cannot be NULL or empty")
  }

  # Build a regex of characters to REMOVE (complement of valid set)
  if (method == "tm_wallace") {
    invalid_pattern <- "[^ABCDGHIKMNRSTVWY]"
  } else if (method == "tm_nn") {
    invalid_pattern <- "[^ACGTI]"
  } else if (method == "tm_gc") {
    invalid_pattern <- "[^ABCDGHIKMNRSTVWXY]"
  } else {
    stop("Invalid method specified")
  }

  # Vectorized single-sequence filter: uppercase + strip invalid chars in one pass
  filter_vec <- function(x) {
    gsub(invalid_pattern, "", toupper(as.character(x)), perl = TRUE)
  }

  # tm_nn special mode: pairwise filtering for sequence/complement with N handling
  if (method == "tm_nn" &&
      is.list(seq_list) &&
      all(c("sequence", "complement") %in% names(seq_list))) {

    seq_vec    <- seq_list$sequence
    comp_vec   <- seq_list$complement
    region_ids <- seq_list$region_ids

    if (is.null(region_ids)) {
      region_ids <- as.character(seq_along(seq_vec))
    }
    if (length(seq_vec) != length(comp_vec)) {
      stop("'sequence' and 'complement' must have the same length")
    }
    if (length(region_ids) != length(seq_vec)) {
      stop("'region_ids' must have the same length as 'sequence'")
    }

    # Fast N detection: works on DNAStringSet directly (no character coercion)
    if (inherits(seq_vec, "DNAStringSet")) {
      has_n <- (Biostrings::vcountPattern("N", seq_vec)  > 0) |
               (Biostrings::vcountPattern("N", comp_vec) > 0)
    } else {
      # character fallback
      has_n <- grepl("N", toupper(seq_vec),  fixed = FALSE) |
               grepl("N", toupper(comp_vec), fixed = FALSE)
    }

    skipped_regions <- character(0)
    if (any(has_n)) {
      skipped_regions <- region_ids[has_n]
      warning(
        sprintf(
          "Skipped %d region(s) because sequence or complement contains 'N'. See your_output$options$'Skipped regions containing N' for details",
          sum(has_n)
        ),
        call. = FALSE
      )
      seq_vec  <- seq_vec[!has_n]
      comp_vec <- comp_vec[!has_n]
    }
    if (length(seq_vec) == 0) {
      stop("No valid regions left for tm_nn calculation after filtering sequences with 'N'.")
    }

    # Vectorized filtering - one gsub call for all sequences
    filtered_seq  <- filter_vec(seq_vec)
    filtered_comp <- filter_vec(comp_vec)

    if (any(nchar(filtered_seq) < 2) || any(nchar(filtered_comp) < 2)) {
      stop("Invalid region or sequence in your input: too many Ns or region contains only Ns")
    }

    return(list(
      sequence        = filtered_seq,
      complement      = filtered_comp,
      kept            = !has_n,
      skipped_regions = skipped_regions
    ))
  }

  # Default mode: single sequence vector filtering - also vectorized now
  result <- filter_vec(seq_list)
  if (any(nchar(result) < 2)) {
    stop("Invalid region or sequence in your input: too many Ns or region contains only Ns")
  }
  result
}

#' convert a vector of characters into a string
#' 
#' Simply convert a vector of characters such as c("H","e","l","l","o","W","o","r","l","d") into a single string "HelloWorld".
#'
#' @param characters A vector of characters 
#' 
#' @returns Retrun a strings
#' 
#' @author Junhui Li
#' 
#' @references 
#' 
#' \code{citation("TmCalculator")}
#' 
#' @keywords internal

c2s <- function(characters){
  strings <- paste0(characters,collapse = "")
  return(strings)
}

#' convert a string into a vector of characters
#' 
#' Simply convert a single string such as "HelloWorld" into a vector of characters such as c("H","e","l","l","o","W","o","r","l","d")
#'
#' @param strings A single string such as "HelloWorld" 
#' 
#' @returns Retrun a vector of characters
#' 
#' @author Junhui Li
#' 
#' @references 
#' 
#' \code{citation("TmCalculator")}
#' 
#' @keywords internal

s2c <- function(strings){
  vec_char <- vector()
  for (i in nchar(strings):1){
    vec_char <- append(unlist(strsplit(strings,''))[i],vec_char)
  }
  return(vec_char)
}

#' Normalize Tm/GC metadata column names on a GRanges object
#'
#' Renames legacy lowercase \code{tm}/\code{gc} columns to \code{Tm}/\code{GC}.
#' @keywords internal
.normalize_tm_gc_metadata <- function(gr) {
  if (!inherits(gr, "GRanges")) {
    return(gr)
  }
  mc <- names(GenomicRanges::mcols(gr))
  if ("tm" %in% mc && !"Tm" %in% mc) {
    GenomicRanges::mcols(gr)$Tm <- GenomicRanges::mcols(gr)$tm
    GenomicRanges::mcols(gr)$tm <- NULL
  }
  if ("gc" %in% mc && !"GC" %in% mc) {
    GenomicRanges::mcols(gr)$GC <- GenomicRanges::mcols(gr)$gc
    GenomicRanges::mcols(gr)$gc <- NULL
  }
  gr
}

