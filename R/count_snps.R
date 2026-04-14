#' Count SNPs in Consensus Sequences
#'
#' Counts polymorphic sites (masked as ambiguity codes) in consensus sequences
#' and calculates sequence statistics for each sample and locus.
#'
#' @param path_to_output_folder Path to HybPhaser output folder containing
#'   consensus sequences
#' @param fasta_file_with_targets Path to FASTA file with target sequences
#' @param targets_file_format Format of target file: "DNA" or "AA" (amino
#'   acid)
#' @param path_to_namelist Path to file containing sample names (one per
#'   line)
#' @param intronerated_contig Logical; whether to use intronerated contigs
#'   (TRUE) or normal contigs (FALSE)
#' @param subset_name Optional name for dataset optimization subset
#'
#' @return List containing two data frames:
#'   \item{tab_snps}{Proportion of SNPs per sample and locus}
#'   \item{tab_length}{Sequence length per sample and locus}
#' @export
#'
#' @examples
#' \dontrun{
#' results <- count_snps(
#'   path_to_output_folder = "hybphaser_output",
#'   fasta_file_with_targets = "targets.fasta",
#'   targets_file_format = "DNA",
#'   path_to_namelist = "samples.txt",
#'   intronerated_contig = FALSE
#' )
#' }
count_snps <- function(
  path_to_output_folder,
  fasta_file_with_targets,
  targets_file_format = c("DNA", "AA"),
  path_to_namelist,
  intronerated_contig = FALSE,
  subset_name = ""
) {
  # Validate inputs
  targets_file_format <- match.arg(targets_file_format)

  validate_paths(
    list(
      "Output folder" = path_to_output_folder,
      "Target FASTA" = fasta_file_with_targets,
      "Namelist" = path_to_namelist
    ),
    must_exist = TRUE,
    type = c("dir", "file", "file")
  )

  # Load required packages
  requireNamespace("seqinr", quietly = TRUE)
  requireNamespace("stringr", quietly = TRUE)

  # Create output directory for R objects
  output_robjects <- file.path(
    path_to_output_folder,
    "00_R_objects",
    subset_name
  )
  dir.create(output_robjects, showWarnings = FALSE, recursive = TRUE)

  # Read target sequences and extract target names
  targets <- seqinr::read.fasta(
    fasta_file_with_targets,
    seqtype = ifelse(targets_file_format == "AA", "AA", "DNA"),
    as.string = TRUE,
    set.attributes = FALSE
  )
  targets_name <- unique(gsub(".*-", "", labels(targets)))

  # Read sample names
  samples <- readLines(path_to_namelist)

  # Set intronerated file naming
  if (intronerated_contig) {
    intronerated_name <- "intronerated"
    intronerated_underscore <- "_"
  } else {
    intronerated_name <- ""
    intronerated_underscore <- ""
  }

  # Initialize result tables
  tab_snps <- data.frame(loci = targets_name)
  tab_length <- data.frame(loci = targets_name)

  # Process each sample
  message("Counting SNPs across samples...")
  start_time <- Sys.time()

  for (sample in samples) {
    elapsed <- round(difftime(Sys.time(), start_time, units = "secs"), 0)
    message(sprintf("%s (%ds)", sample, elapsed))

    tab_sample <- data.frame(
      targets = targets_name,
      seq_length = NA_real_,
      ambis = NA_real_,
      ambi_prop = NA_real_
    )

    consensus_dir <- file.path(
      path_to_output_folder,
      "01_data",
      sample,
      paste0(intronerated_name, intronerated_underscore, "consensus")
    )

    if (!dir.exists(consensus_dir)) {
      warning("Consensus directory not found for sample: ", sample)
      next
    }

    consensus_files <- list.files(consensus_dir, full.names = TRUE)

    for (consensus_file in consensus_files) {
      gene <- gsub(
        "(_intronerated|).fasta$",
        "",
        basename(consensus_file)
      )

      if (file.info(consensus_file)$size > 1) {
        stats <- .seq_stats(consensus_file)
      } else {
        stats <- c(NA_real_, NA_real_)
      }

      target_idx <- match(gene, tab_sample$targets)
      tab_sample$seq_length[target_idx] <- stats[1]
      tab_sample$ambis[target_idx] <- stats[2]
    }

    tab_sample$ambi_prop <- tab_sample$ambis / tab_sample$seq_length
    tab_snps[[sample]] <- tab_sample$ambi_prop
    tab_length[[sample]] <- tab_sample$seq_length
  }

  # Set row/column names
  rownames(tab_snps) <- targets_name
  rownames(tab_length) <- targets_name

  # Save results
  saveRDS(tab_snps, file = file.path(output_robjects, "Table_SNPs.Rds"))
  saveRDS(
    tab_length,
    file = file.path(output_robjects, "Table_consensus_length.Rds")
  )

  message("SNP counting completed")

  list(
    tab_snps = tab_snps,
    tab_length = tab_length
  )
}


#' Calculate Sequence Statistics from FASTA File
#'
#' Internal helper function to count ambiguity codes in a sequence
#'
#' @param file Path to FASTA file
#'
#' @return Numeric vector with sequence length and ambiguity count
#' @keywords internal
.seq_stats <- function(file) {
  # Try to read FASTA file, handle errors for empty/invalid files
  fasta <- tryCatch(
    seqinr::read.fasta(file, as.string = TRUE, set.attributes = FALSE),
    error = function(e) list()
  )

  if (length(fasta) == 0) {
    return(c(NA_real_, NA_real_))
  }

  # Remove N, ?, and - (case-insensitive)
  seq <- gsub("[Nn?-]", "", fasta[[1]])
  seq_length <- stringr::str_length(seq)

  # Count ambiguity codes (2-way ambiguity codes = 1, 3-way = 2)
  ambis_2way <- stringr::str_count(seq, "[YKRSMWykrsmw]")
  ambis_3way <- stringr::str_count(seq, "[DHBVdhbv]")
  total_ambis <- ambis_2way + (ambis_3way * 2)

  c(round(seq_length, 0), total_ambis)
}
