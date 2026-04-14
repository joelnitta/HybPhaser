#' Extract Mapped Reads
#'
#' Extracts reads mapped to target sequences from HybPhaser or HybPiper
#' output, concatenates them per sample, and optionally deduplicates by
#' sequence.
#'
#' This is an R implementation of the `2_extract_mapped_reads.sh` workflow
#' step.
#'
#' @param base_dir Path to base folder of either:
#'   \itemize{
#'   \item HybPhaser output (contains `01_data/`)
#'   \item HybPiper-RBGV output (contains `04_processed_gene_directories/`)
#'   \item HybPiper output (contains sample directories)
#'   }
#' @param output_dir Path to output directory for mapped read files.
#'   Default is `./mapped_reads/`.
#' @param namelist Optional path to a file containing sample names, one per
#'   line. If `NULL`, all sample directories are processed.
#' @param remove_duplicate_sequences Logical; if `TRUE`, deduplicate reads by
#'   sequence (not just duplicated names).
#'
#' @return Invisibly returns a list with:
#' \describe{
#'   \item{output_dir}{Output directory path}
#'   \item{samples}{Character vector of processed samples}
#'   \item{files}{Named character vector of output files by sample}
#'   \item{n_samples}{Number of processed samples}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' result <- extract_mapped_reads(
#'   base_dir = "hybphaser_output",
#'   output_dir = "mapped_reads"
#' )
#'
#' result_dedup <- extract_mapped_reads(
#'   base_dir = "hybphaser_output",
#'   output_dir = "mapped_reads",
#'   namelist = "samples.txt",
#'   remove_duplicate_sequences = TRUE
#' )
#' }
extract_mapped_reads <- function(
  base_dir = "./",
  output_dir = "./mapped_reads/",
  namelist = NULL,
  remove_duplicate_sequences = FALSE
) {
  if (!dir.exists(base_dir)) {
    stop("base_dir does not exist: ", base_dir)
  }

  if (!is.null(namelist) && !file.exists(namelist)) {
    stop("namelist does not exist: ", namelist)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  data_info <- get_samples_dir_for_extract(base_dir)
  samples_dir <- data_info$samples_dir
  input_type <- data_info$input_type

  if (!is.null(namelist)) {
    samples <- readLines(namelist, warn = FALSE)
    samples <- samples[nzchar(samples)]
  } else {
    samples <- list_sample_dirs(samples_dir)
  }

  output_files <- stats::setNames(character(length(samples)), samples)

  for (sample in samples) {
    message("Generating on-target reads files for ", sample)

    fasta_files <- get_mapped_read_files(
      samples_dir = samples_dir,
      sample = sample,
      input_type = input_type
    )

    output_file <- file.path(output_dir, paste0(sample, "_mapped_reads.fasta"))

    concatenate_fasta_files(fasta_files, output_file)

    if (remove_duplicate_sequences) {
      message("Deduplicating reads files for ", sample)
      dedup_file <- file.path(
        output_dir,
        paste0(sample, "_mapped_reads_dedup-seq.fasta")
      )
      deduplicate_fasta_by_sequence(output_file, dedup_file)
      file.remove(output_file)
      output_files[[sample]] <- dedup_file
    } else {
      message("Keeping all reads for ", sample)
      output_files[[sample]] <- output_file
    }
  }

  invisible(list(
    output_dir = output_dir,
    samples = samples,
    files = output_files,
    n_samples = as.numeric(length(samples))
  ))
}


# Identify expected sample directory layout for extraction.
get_samples_dir_for_extract <- function(base_dir) {
  if (dir.exists(file.path(base_dir, "01_data"))) {
    return(list(
      samples_dir = file.path(base_dir, "01_data"),
      input_type = "phaser"
    ))
  }

  if (dir.exists(file.path(base_dir, "04_processed_gene_directories"))) {
    return(list(
      samples_dir = file.path(base_dir, "04_processed_gene_directories"),
      input_type = "piper"
    ))
  }

  list(samples_dir = base_dir, input_type = "piper")
}


# Return sample folder names in a directory.
list_sample_dirs <- function(samples_dir) {
  entries <- list.files(samples_dir, full.names = TRUE, recursive = FALSE)
  entries <- entries[file.info(entries)$isdir]
  basename(entries)
}


# Collect mapped read FASTA files for a sample by input type.
get_mapped_read_files <- function(samples_dir, sample, input_type) {
  if (identical(input_type, "phaser")) {
    reads_dir <- file.path(samples_dir, sample, "reads")
    if (!dir.exists(reads_dir)) {
      return(character(0))
    }

    unpaired <- list.files(
      reads_dir,
      pattern = "_unpaired\\.fasta$",
      full.names = TRUE
    )
    combined <- list.files(
      reads_dir,
      pattern = "_combined\\.fasta$",
      full.names = TRUE
    )
    return(c(unpaired, combined))
  }

  sample_dir <- file.path(samples_dir, sample)
  if (!dir.exists(sample_dir)) {
    return(character(0))
  }

  unpaired <- list.files(
    sample_dir,
    pattern = "_unpaired\\.fasta$",
    recursive = TRUE,
    full.names = TRUE
  )
  interleaved <- list.files(
    sample_dir,
    pattern = "_interleaved\\.fasta$",
    recursive = TRUE,
    full.names = TRUE
  )

  c(unpaired, interleaved)
}


# Concatenate FASTA files to one output FASTA.
concatenate_fasta_files <- function(fasta_files, output_file) {
  if (length(fasta_files) == 0) {
    file.create(output_file)
    return(invisible(NULL))
  }

  lines <- unlist(lapply(fasta_files, readLines, warn = FALSE))
  writeLines(lines, output_file)
  invisible(NULL)
}


# Deduplicate FASTA entries by sequence while keeping first header seen.
deduplicate_fasta_by_sequence <- function(input_file, output_file) {
  if (!file.exists(input_file)) {
    file.create(output_file)
    return(invisible(NULL))
  }

  fasta <- seqinr::read.fasta(
    file = input_file,
    as.string = TRUE,
    set.attributes = FALSE,
    seqtype = "DNA"
  )

  if (length(fasta) == 0) {
    file.create(output_file)
    return(invisible(NULL))
  }

  seqs <- vapply(fasta, function(x) x[[1]], character(1))
  headers <- names(fasta)
  keep_idx <- !duplicated(seqs)

  seqinr::write.fasta(
    sequences = as.list(seqs[keep_idx]),
    names = headers[keep_idx],
    file.out = output_file,
    nbchar = 1e6
  )

  invisible(NULL)
}
