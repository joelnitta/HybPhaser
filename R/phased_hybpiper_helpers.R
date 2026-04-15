#' Run HybPiper on Phased Read Files
#'
#' Runs `hybpiper assemble` for each phased read file in a directory.
#' This is useful after [run_phasing()] in workflows where phased reads are
#' assembled separately.
#'
#' @param phased_reads_dir Path to directory containing phased read files.
#' @param phased_hp_dir Path to output directory for phased HybPiper results.
#' @param targets_file Path to DNA targets FASTA file.
#' @param conda_env Conda environment containing HybPiper.
#' @param cpu Number of CPU threads per HybPiper run.
#'
#' @return Character scalar. The `phased_hp_dir` path.
#' @export
#'
#' @examples
#' \dontrun{
#' out_dir <- run_phased_hybpiper(
#'   phased_reads_dir = "hybphaser_output/05_phasing/phased_reads",
#'   phased_hp_dir = "hybpiper_output_phased",
#'   targets_file = "test_targets.fasta"
#' )
#' }
run_phased_hybpiper <- function(
  phased_reads_dir,
  phased_hp_dir,
  targets_file,
  conda_env = "hybpiper_env",
  cpu = 1
) {
  if (!dir.exists(phased_reads_dir)) {
    stop("Phased reads directory not found: ", phased_reads_dir)
  }

  if (!file.exists(targets_file)) {
    stop("Targets file not found: ", targets_file)
  }

  phased_read_files <- list.files(
    phased_reads_dir,
    pattern = "\\.fastq(\\.gz)?$",
    full.names = TRUE
  )

  if (length(phased_read_files) == 0) {
    stop("No phased read files found in: ", phased_reads_dir)
  }

  phased_samples <- sub(
    "\\.fastq(\\.gz)?$",
    "",
    basename(phased_read_files)
  )

  dir.create(phased_hp_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_along(phased_samples)) {
    status <- run_conda(
      command = "hybpiper",
      env_name = conda_env,
      wd = phased_hp_dir,
      args = c(
        "assemble",
        "--readfiles",
        phased_read_files[i],
        "--targetfile_dna",
        targets_file,
        "--prefix",
        phased_samples[i],
        "--bwa",
        "--cpu",
        as.character(cpu)
      ),
      stdout = "",
      stderr = ""
    )

    if (!identical(as.integer(status), 0L)) {
      stop(
        "HybPiper assemble failed for phased sample '",
        phased_samples[i],
        "' with exit code: ",
        status
      )
    }
  }

  phased_hp_dir
}


#' Create Namelist from Phased Read Files
#'
#' Creates a namelist text file from phased read file names in a directory.
#'
#' @param phased_reads_dir Path to directory containing phased read files.
#' @param phased_hp_dir Directory where the namelist file should be written.
#' @param output_filename Name of output namelist file.
#'
#' @return Character scalar. Path to created namelist file.
#' @export
#'
#' @examples
#' \dontrun{
#' namelist_path <- make_phased_namelist(
#'   phased_reads_dir = "hybphaser_output/05_phasing/phased_reads",
#'   phased_hp_dir = "hybpiper_output_phased"
#' )
#' }
make_phased_namelist <- function(
  phased_reads_dir,
  phased_hp_dir,
  output_filename = "namelist_phased.txt"
) {
  if (!dir.exists(phased_reads_dir)) {
    stop("Phased reads directory not found: ", phased_reads_dir)
  }

  phased_read_files <- list.files(
    phased_reads_dir,
    pattern = "\\.fastq(\\.gz)?$",
    full.names = TRUE
  )

  if (length(phased_read_files) == 0) {
    stop("No phased read files found in: ", phased_reads_dir)
  }

  phased_samples <- sub(
    "\\.fastq(\\.gz)?$",
    "",
    basename(phased_read_files)
  )

  dir.create(phased_hp_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(phased_hp_dir, output_filename)
  writeLines(phased_samples, out_file)
  out_file
}
