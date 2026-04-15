#' Write Clade Reference CSV
#'
#' Writes a clade-reference CSV file with `samples` and `abb` columns.
#'
#' @param reference_samples Character vector of reference sample names.
#' @param output_file Path to output CSV file.
#'
#' @return Character scalar. Path to created CSV file.
#' @export
#'
#' @examples
#' \dontrun{
#' write_clade_reference_csv(
#'   reference_samples = c("sampleA", "sampleB"),
#'   output_file = "04_clade_association/clade_references_example.csv"
#' )
#' }
write_clade_reference_csv <- function(reference_samples, output_file) {
  if (!is.character(reference_samples) || length(reference_samples) == 0) {
    stop("'reference_samples' must be a non-empty character vector")
  }

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    data.frame(
      samples = reference_samples,
      abb = paste0("C", seq_along(reference_samples))
    ),
    output_file,
    row.names = FALSE
  )

  output_file
}


#' Write Phasing Prep CSV
#'
#' Writes a minimal phasing-prep CSV file with columns `sample`, `ref1`,
#' and `abb1`.
#'
#' @param reference_samples Character vector of reference sample names.
#' @param output_file Path to output CSV file.
#' @param abbreviation Character scalar for `abb1` value. Default is "A".
#'
#' @return Character scalar. Path to created CSV file.
#' @export
#'
#' @examples
#' \dontrun{
#' write_phasing_prep_csv(
#'   reference_samples = c("sampleA", "sampleB"),
#'   output_file = "05_phasing/phasing_prep_example.csv"
#' )
#' }
write_phasing_prep_csv <- function(
  reference_samples,
  output_file,
  abbreviation = "A"
) {
  if (!is.character(reference_samples) || length(reference_samples) == 0) {
    stop("'reference_samples' must be a non-empty character vector")
  }

  if (!is.character(abbreviation) || length(abbreviation) != 1) {
    stop("'abbreviation' must be a single character value")
  }

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    data.frame(
      sample = reference_samples[[1]],
      ref1 = reference_samples[[1]],
      abb1 = abbreviation
    ),
    output_file,
    row.names = FALSE
  )

  output_file
}


#' Generate Sequence Lists After Assessment
#'
#' Wrapper around [generate_sequence_lists()] that takes an assessment object
#' to make dependency ordering explicit in `{targets}` pipelines.
#'
#' @param assessment Output from [assess_dataset()].
#' @param ... Arguments passed to [generate_sequence_lists()].
#'
#' @return List returned by [generate_sequence_lists()].
#' @export
#'
#' @examples
#' \dontrun{
#' seq_lists <- generate_sequence_lists_after_assessment(
#'   assessment = assessment_result,
#'   path_to_output_folder = "hybphaser_output",
#'   fasta_file_with_targets = "targets.fasta",
#'   targets_file_format = "DNA",
#'   path_to_namelist = "namelist.txt"
#' )
#' }
generate_sequence_lists_after_assessment <- function(assessment, ...) {
  if (is.null(assessment)) {
    stop("'assessment' must be provided.")
  }

  # Force evaluation of assessment to ensure targets recognizes the dependency
  force(assessment)

  generate_sequence_lists(...)
}


#' Get Reference Samples from Consensus Sequence Files
#'
#' Returns reference sample names from consensus FASTA filenames. Includes a
#' `sequence_lists` argument so `{targets}` pipelines can enforce ordering after
#' [generate_sequence_lists()].
#'
#' @param sequence_lists Output from [generate_sequence_lists()].
#' @param reference_dir Directory containing `*_consensus.fasta` files.
#' @param pattern Regex pattern used to match consensus files.
#'
#' @return Character vector of reference sample names.
#' @export
#'
#' @examples
#' \\dontrun{
#' refs <- get_reference_samples(
#'   sequence_lists = seq_lists,
#'   reference_dir = "hybphaser_output/03_sequence_lists/samples_consensus"
#' )
#' }
get_reference_samples <- function(
  sequence_lists,
  reference_dir,
  pattern = "_consensus\\\\.fasta$"
) {
  if (is.null(sequence_lists)) {
    stop("'sequence_lists' must be provided.")
  }

  # Force evaluation to ensure targets recognizes the dependency
  force(sequence_lists)

  if (!dir.exists(reference_dir)) {
    stop("Reference directory not found: ", reference_dir)
  }

  reference_files <- list.files(reference_dir, pattern = pattern)
  reference_samples <- sub(pattern, "", reference_files)

  if (length(reference_samples) == 0) {
    stop("No consensus reference files found in: ", reference_dir)
  }

  reference_samples
}


#' Run HybPiper Test Dataset with Clean Output Directory
#'
#' Removes an existing output directory before calling
#' [run_hybpiper_test_dataset()]. This prevents stale files from previous
#' runs from triggering HybPiper overwrite errors in `{targets}` workflows.
#'
#' @param output_dir Output directory passed to [run_hybpiper_test_dataset()].
#' @param ... Additional arguments passed to [run_hybpiper_test_dataset()].
#'
#' @return List returned by [run_hybpiper_test_dataset()].
#' @export
#'
#' @examples
#' \\dontrun{
#' hp_data <- run_hybpiper_test_dataset_clean(
#'   output_dir = "targets_hybpiper_test",
#'   samples = c("EG30"),
#'   conda_env = "hybpiper_env",
#'   cpu = 1,
#'   mapper = "bwa"
#' )
#' }
run_hybpiper_test_dataset_clean <- function(output_dir, ...) {
  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("'output_dir' must be a single character path")
  }

  # Remove only FASTQ files (not sample or output directories).
  # This stops HybPiper from seeing duplicate read files on rerun while
  # keeping any pre-existing hybphaser output directories in place.
  # Keeping those directories avoids a Docker Desktop / VirtioFS timing
  # issue on macOS: Docker cannot bind-mount a directory that was created
  # in the same moment the container starts, so we preserve the output dir
  # across runs rather than deleting and immediately recreating it.
  if (dir.exists(output_dir)) {
    fastq_files <- list.files(
      output_dir,
      pattern = "\\.fastq(\\.gz)?$",
      full.names = TRUE
    )
    if (length(fastq_files) > 0L) {
      file.remove(fastq_files)
    }
  }

  # Pass --force_overwrite so HybPiper can re-assemble into existing
  # sample directories without erroring on "directory already exists".
  run_hybpiper_test_dataset(
    output_dir = output_dir,
    ...,
    other_args = "--force_overwrite"
  )
}
