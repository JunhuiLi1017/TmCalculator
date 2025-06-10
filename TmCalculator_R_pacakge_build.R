
shiny::runApp("~/dropbox/Project/UMMS/Github/JunhuiLi1017/TmCalculator/TmCalculatorApp.R")


load_all("~/dropbox/Project/UMMS/Github/JunhuiLi1017/TmCalculator")

pak::pkg_install("JunhuiLi1017/TmCalculator@dev")


library(TmCalculator)
?tm_calculate
?process_seq
?TmCalculator

primer = "GTGCCAGCAGCCGCGGTCAAAC"
complement_seq = as.character(generate_complement("GTGCCAGCCACCGCGGTTATAC"))
nchar(primer)
nchar(template)

setwd("~/Dropbox (UMass Medical School)/Project/UMMS/Github/JunhuiLi1017/TmCalculator/")
library(usethis)
library(devtools)
#install.packages("roxygen2")
##------
##git
##------
##
use_git()
library(shiny)
runApp(appDir = system.file('shiny', package = 'TmCalculator'))

##------
##add lecense to  to description
##------
use_mit_license()
#usethis::use_package_doc()

##------
##convert all R with roxygen2 format to Rd file in man folder
##------
document()
#roxygen2::roxygenise("~/dropbox/Project/UMMS/Github/JunhuiLi1017/StepReg/")
#?StepReg
##------
##check if there are some warnings
##------
warnings()

tools::showNonASCIIfile()

##------
##build and check R package
##------
rm(list = c("plot_tm_genome_tracks", "plot_tm_heatmap", "plot_tm_karyotype"))
devtools::load_all("~/dropbox/Project/UMMS/Github/JunhuiLi1017/TmCalculator/")
devtools::build()
#devtools::test()
devtools::check()
devtools::install()


input_seq <- c("ATGCGATGCG", "ATGCGATGCG")
result <- tm_calculate(input_seq,method=c("tm_gc"))
result$tm_gc$Options


input_seq <- "~/dropbox/Project/UMMS/Github/JunhuiLi1017/TmCalculator/inst/extdata/example1.fasta"
result <- tm_calculate(input_seq)
result$tm$Options
# Calculate Tm with specific method parameters
generate_complement("GCTAGCCGACAATGGGCAGATAGTAGAAA")





generate_random_dna_sequences <- function(num_sequences = 200, min_len = 60, max_len = 150) {
  # Define the possible DNA bases
  bases <- c("A", "T", "C", "G")
  
  # Initialize an empty vector to store the DNA sequences
  dna_sequences <- character(num_sequences)
  
  # Loop to generate each sequence
  for (i in 1:num_sequences) {
    # Randomly choose a length for the current sequence
    seq_length <- sample(min_len:max_len, 1)
    
    # Randomly select bases for the chosen length
    random_bases <- sample(bases, seq_length, replace = TRUE)
    
    # Collapse the bases into a single string
    dna_sequences[i] <- paste0(random_bases, collapse = "")
  }
  
  return(dna_sequences)
}




library(GenomicRanges) # For GRanges object
library(IRanges)       # For IRanges within GRanges

generate_dna_sequences_with_named_coords_by_chr <- function(
    chromosomes = c("chr1", "chr3", "chr13"),
    num_sequences_per_chr = c(100, 200, 100), # Must match length of 'chromosomes'
    min_len_per_chr = c(60, 70, 60),      # Must match length of 'chromosomes'
    max_len_per_chr = c(120, 130, 150),   # Must match length of 'chromosomes'
    # Optional: Named vector of chromosome lengths. Defaults to human approx.
    approx_chr_lengths = NULL
) {
  
  # --- Input Validation ---
  if (!is.character(chromosomes) || length(chromosomes) == 0) {
    stop("Input 'chromosomes' must be a non-empty character vector.")
  }
  if (!all(length(chromosomes) == length(num_sequences_per_chr),
           length(chromosomes) == length(min_len_per_chr),
           length(chromosomes) == length(max_len_per_chr))) {
    stop("Lengths of 'chromosomes', 'num_sequences_per_chr', 'min_len_per_chr', and 'max_len_per_chr' must all be the same.")
  }
  if (!is.numeric(num_sequences_per_chr) || any(num_sequences_per_chr < 0)) {
    stop("'num_sequences_per_chr' must be a numeric vector with non-negative values.")
  }
  if (!is.numeric(min_len_per_chr) || any(min_len_per_chr < 1)) {
    stop("'min_len_per_chr' must be a numeric vector with values >= 1.")
  }
  if (!is.numeric(max_len_per_chr) || any(max_len_per_chr < min_len_per_chr)) {
    stop("'max_len_per_chr' must be a numeric vector with values >= corresponding min_len_per_chr.")
  }
  
  
  # Define approximate chromosome lengths if not provided.
  if (is.null(approx_chr_lengths)) {
    # Using a comprehensive set of human chromosome lengths (hg38/GRCh38)
    approx_chr_lengths <- c(
      chr1 = 248956422, chr2 = 242193529, chr3 = 198295559, chr4 = 190214555,
      chr5 = 181538259, chr6 = 170805979, chr7 = 159345973, chr8 = 145138636,
      chr9 = 138394717, chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
      chr13 = 114364328, chr14 = 107043718, chr15 = 101991189, chr16 = 90338345,
      chr17 = 83257441, chr18 = 80373285, chr19 = 58617616, chr20 = 64441676,
      chr21 = 46709983, chr22 = 50818468, chrX = 156040895, chrY = 57227415,
      chrM = 16569 # Mitochondrial chromosome
    )
  }
  
  # Ensure all specified chromosomes have lengths
  missing_lengths <- setdiff(chromosomes, names(approx_chr_lengths))
  if (length(missing_lengths) > 0) {
    stop(paste("Approximate lengths missing for chromosomes:",
               paste(missing_lengths, collapse = ", "),
               ". Please provide them in 'approx_chr_lengths' or add them to the default set."))
  }
  
  # Define the possible DNA bases
  bases <- c("A", "T", "C", "G")
  
  # Initialize lists to store all generated sequences and their names
  all_dna_sequences_list <- list()
  all_names_list <- list()
  
  # --- Loop through each specified chromosome ---
  for (j in seq_along(chromosomes)) {
    current_chr <- chromosomes[j]
    num_to_generate <- num_sequences_per_chr[j]
    current_min_len <- min_len_per_chr[j]
    current_max_len <- max_len_per_chr[j]
    max_chr_coord <- approx_chr_lengths[current_chr]
    
    # Initialize temporary lists for sequences on the current chromosome
    chr_dna_sequences <- vector("list", num_to_generate)
    chr_names <- character(num_to_generate)
    
    message(paste0("Generating ", num_to_generate, " sequences for ", current_chr,
                   " (length range: ", current_min_len, "-", current_max_len, ")"))
    
    # --- Inner loop to generate sequences for the current chromosome ---
    for (i in 1:num_to_generate) {
      # 1. Generate DNA sequence length
      seq_length <- sample(current_min_len:current_max_len, 1)
      
      # 2. Generate random DNA bases
      random_bases <- sample(bases, seq_length, replace = TRUE)
      dna_seq <- paste0(random_bases, collapse = "")
      
      # 3. Randomly assign start and calculate end
      max_possible_start <- max_chr_coord - seq_length + 1
      
      if (max_possible_start < 1) {
        warning(paste0("Sequence length (", seq_length, ") is too long for ", current_chr,
                       " (length ", max_chr_coord, "). Skipping this sequence."))
        next # Skip to next iteration if sequence is too long
      }
      
      current_start <- sample(1:max_possible_start, 1)
      current_end <- current_start + seq_length - 1
      
      # 4. Store the sequence and its name
      chr_dna_sequences[[i]] <- dna_seq
      chr_names[i] <- paste0(current_chr, ":", current_start, "-", current_end)
    }
    
    # Add generated sequences and names for the current chromosome to the main lists
    all_dna_sequences_list <- c(all_dna_sequences_list, chr_dna_sequences)
    all_names_list <- c(all_names_list, chr_names)
  }
  
  # --- Combine all generated sequences into a single named character vector ---
  # Filter out any NULL entries if sequences were skipped due to length warnings
  final_sequences <- unlist(all_dna_sequences_list[!sapply(all_dna_sequences_list, is.null)])
  names(final_sequences) <- unlist(all_names_list[!sapply(all_dna_sequences_list, is.null)])
  
  return(final_sequences)
}

# --- Example Usage with your specified parameters ---
my_mixed_dna_vector <- generate_dna_sequences_with_named_coords_by_chr(
  chromosomes = c("chr1", "chr3", "chr13"),
  num_sequences_per_chr = c(100, 200, 100),
  min_len_per_chr = c(60, 70, 60),
  max_len_per_chr = c(120, 130, 150)
)
names(my_mixed_dna_vector) <- NULL
# You can access the sequences via mcols(my_gr_sequences)$sequence
# head(mcols(my_gr_sequences)$sequence)

# And check the structure of the GRanges object
# class(my_gr_sequences)
# 
# 
result <- tm_calculate(my_mixed_dna_vector)



