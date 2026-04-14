#' Check if Conda is Available
#'
#' Verifies that conda is installed and can run commands.
#'
#' @param quiet Logical; if TRUE, suppresses messages
#'
#' @return Logical; TRUE if conda is available, FALSE otherwise
#' @export
check_conda <- function(quiet = FALSE) {
  conda_path <- Sys.which("conda")

  if (conda_path == "") {
    if (!quiet) {
      warning("Conda command not found. Please install Conda.")
    }
    return(FALSE)
  }

  result <- system2(
    conda_path,
    args = c("info", "--json"),
    stdout = FALSE,
    stderr = FALSE
  )

  if (result != 0) {
    if (!quiet) {
      warning("Conda command found but failed to run.")
    }
    return(FALSE)
  }

  if (!quiet) {
    message("Conda is available")
  }

  TRUE
}


#' Run a Command in a Conda Environment
#'
#' Runs a command via `conda run`, optionally within a named conda
#' environment.
#'
#' @param command Character; command to run (for example, "hybpiper")
#' @param args Character vector of arguments passed to command
#' @param env_name Optional conda environment name. If NULL, conda uses
#'   the default environment.
#' @param conda_executable Character; path or name of conda executable
#' @param wd Working directory for the command
#' @param ... Additional arguments passed to [system2()]
#'
#' @return Exit status code from [system2()]
#' @export
run_conda <- function(
  command,
  args = character(0),
  env_name = NULL,
  conda_executable = "conda",
  wd = getwd(),
  ...
) {
  if (!is.character(command) || length(command) != 1 || nchar(command) == 0) {
    stop("'command' must be a non-empty character string")
  }

  if (!is.character(args)) {
    stop("'args' must be a character vector")
  }

  if (
    !is.null(env_name) &&
      (!is.character(env_name) || length(env_name) != 1 || nchar(env_name) == 0)
  ) {
    stop("'env_name' must be NULL or a non-empty character string")
  }

  if (!check_conda(quiet = TRUE)) {
    stop("Conda is not available. Please install and configure Conda.")
  }

  conda_args <- c("run", "--no-capture-output")
  if (!is.null(env_name)) {
    conda_args <- c(conda_args, "-n", env_name)
  }
  conda_args <- c(conda_args, command, args)

  old_wd <- setwd(wd)
  on.exit(setwd(old_wd), add = TRUE)

  system2(conda_executable, args = conda_args, ...)
}


# Download and extract HybPiper official test dataset from upstream.
.download_hybpiper_test_dataset <- function(dest_dir) {
  dataset_url <- paste0(
    "https://raw.githubusercontent.com/mossmatters/",
    "HybPiper/master/test_dataset.tar.gz"
  )

  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  archive_path <- file.path(dest_dir, "test_dataset.tar.gz")
  dataset_dir <- file.path(dest_dir, "test_dataset")

  if (dir.exists(dataset_dir)) {
    unlink(dataset_dir, recursive = TRUE, force = TRUE)
  }

  status <- tryCatch(
    utils::download.file(
      url = dataset_url,
      destfile = archive_path,
      mode = "wb",
      quiet = TRUE
    ),
    error = function(e) e
  )

  if (inherits(status, "error") || !identical(status, 0L)) {
    stop(
      "Failed to download HybPiper test dataset from ",
      dataset_url,
      ". Check internet access and try again."
    )
  }

  utils::untar(archive_path, exdir = dest_dir)

  targets_file <- file.path(dataset_dir, "test_targets.fasta")
  namelist_file <- file.path(dataset_dir, "namelist.txt")
  reads_archive <- file.path(dataset_dir, "test_reads.fastq.tar.gz")

  if (!file.exists(targets_file) || !file.exists(reads_archive)) {
    stop("Downloaded HybPiper test dataset is missing required files")
  }

  list(
    dataset_dir = dataset_dir,
    targets_file = targets_file,
    namelist_file = namelist_file,
    reads_archive = reads_archive
  )
}


#' Check if HybPiper is Available in a Conda Environment
#'
#' Verifies that `hybpiper` can be executed from the specified conda
#' environment.
#'
#' @param env_name Character; conda environment name
#' @param quiet Logical; if TRUE, suppresses messages
#'
#' @return Logical; TRUE if HybPiper is available in the environment
#' @export
check_hybpiper_conda <- function(env_name = "hybpiper_env", quiet = FALSE) {
  if (!check_conda(quiet = TRUE)) {
    if (!quiet) {
      warning("Conda is not available")
    }
    return(FALSE)
  }

  status <- suppressWarnings(system2(
    "conda",
    args = c(
      "run",
      "--no-capture-output",
      "-n",
      env_name,
      "hybpiper",
      "--version"
    ),
    stdout = FALSE,
    stderr = FALSE
  ))

  if (status == 0) {
    if (!quiet) {
      message("HybPiper is available in conda env '", env_name, "'")
    }
    return(TRUE)
  }

  if (!quiet) {
    warning(
      "HybPiper is not available in conda env '",
      env_name,
      "'. Install with: conda install -n ",
      env_name,
      " -c bioconda hybpiper"
    )
  }
  FALSE
}


#' Run HybPiper Assemble via Conda
#'
#' Runs `hybpiper assemble` for one paired-end sample in a conda
#' environment.
#'
#' @param sample_name Character; sample prefix for HybPiper output
#' @param f_read Character; path to forward read file (R1)
#' @param r_read Character; path to reverse read file (R2)
#' @param targets_file Character; path to DNA targets FASTA file
#' @param wd Character; working directory where HybPiper output is written
#' @param mapper Character; mapping method, one of "bwa" or "blast"
#' @param cpu Integer; number of CPUs for HybPiper
#' @param conda_env Character; conda environment name containing HybPiper
#' @param other_args Character vector; additional CLI args passed to
#'   `hybpiper assemble`
#'
#' @return Character path to expected sample output directory
#' @export
hybpiper_assemble <- function(
  sample_name,
  f_read,
  r_read,
  targets_file,
  wd,
  mapper = c("bwa", "blast"),
  cpu = 1,
  conda_env = "hybpiper_env",
  other_args = character(0)
) {
  mapper <- match.arg(mapper)

  if (!file.exists(f_read)) {
    stop("Forward read file not found: ", f_read)
  }
  if (!file.exists(r_read)) {
    stop("Reverse read file not found: ", r_read)
  }
  if (!file.exists(targets_file)) {
    stop("Targets file not found: ", targets_file)
  }

  if (!dir.exists(wd)) {
    dir.create(wd, recursive = TRUE)
  }

  if (!is.numeric(cpu) || length(cpu) != 1 || cpu < 1) {
    stop("'cpu' must be a positive number")
  }

  if (!check_hybpiper_conda(conda_env, quiet = TRUE)) {
    stop(
      "HybPiper not found in conda environment '",
      conda_env,
      "'. Install it with: conda install -n ",
      conda_env,
      " -c bioconda hybpiper"
    )
  }

  mapper_flag <- if (mapper == "bwa") "--bwa" else "--blast"

  message("Running HybPiper assemble for ", sample_name)

  status <- run_conda(
    command = "hybpiper",
    env_name = conda_env,
    wd = wd,
    args = c(
      "assemble",
      "--readfiles",
      normalizePath(f_read),
      normalizePath(r_read),
      "--targetfile_dna",
      normalizePath(targets_file),
      "--prefix",
      sample_name,
      "--cpu",
      as.character(as.integer(cpu)),
      mapper_flag,
      other_args
    ),
    stdout = "",
    stderr = ""
  )

  if (status != 0) {
    stop("HybPiper assemble failed for sample ", sample_name)
  }

  file.path(wd, sample_name)
}


#' Run HybPiper Stats via Conda
#'
#' Runs `hybpiper stats` across one or more assembled samples.
#'
#' @param targets_file Character; path to target FASTA file
#' @param namelist Character; path to text file with one sample per line
#' @param wd Character; working directory for output files
#' @param mode Character; one of "gene" or "supercontig"
#' @param dna Logical; TRUE for DNA targets, FALSE for amino-acid targets
#' @param seq_lengths_filename Character; basename for sequence lengths TSV
#' @param stats_filename Character; basename for stats TSV
#' @param conda_env Character; conda environment name containing HybPiper
#' @param other_args Character vector; additional CLI args passed to
#'   `hybpiper stats`
#'
#' @return Named list with paths to output TSV files
#' @export
hybpiper_stats <- function(
  targets_file,
  namelist,
  wd = getwd(),
  mode = c("gene", "supercontig"),
  dna = TRUE,
  seq_lengths_filename = "seq_lengths",
  stats_filename = "hybpiper_stats",
  conda_env = "hybpiper_env",
  other_args = character(0)
) {
  mode <- match.arg(mode)

  if (!file.exists(targets_file)) {
    stop("Targets file not found: ", targets_file)
  }
  if (!file.exists(namelist)) {
    stop("Namelist file not found: ", namelist)
  }

  if (!dir.exists(wd)) {
    dir.create(wd, recursive = TRUE)
  }

  if (!check_hybpiper_conda(conda_env, quiet = TRUE)) {
    stop(
      "HybPiper not found in conda environment '",
      conda_env,
      "'. Install it with: conda install -n ",
      conda_env,
      " -c bioconda hybpiper"
    )
  }

  target_flag <- if (isTRUE(dna)) "--targetfile_dna" else "--targetfile_aa"

  message("Running HybPiper stats")

  status <- run_conda(
    command = "hybpiper",
    env_name = conda_env,
    wd = wd,
    args = c(
      "stats",
      target_flag,
      normalizePath(targets_file),
      "--seq_lengths_filename",
      seq_lengths_filename,
      "--stats_filename",
      stats_filename,
      mode,
      normalizePath(namelist),
      other_args
    ),
    stdout = "",
    stderr = ""
  )

  if (status != 0) {
    stop("HybPiper stats failed")
  }

  list(
    seq_lengths = file.path(wd, paste0(seq_lengths_filename, ".tsv")),
    stats = file.path(wd, paste0(stats_filename, ".tsv"))
  )
}


#' Run Official HybPiper Test Dataset via Conda
#'
#' Downloads and runs the official HybPiper test dataset via conda. The
#' dataset is fetched from the upstream HybPiper repository each run.
#'
#' @param output_dir Character; directory where HybPiper output is written
#' @param samples Character vector of sample names from the test dataset
#' @param conda_env Character; conda environment name containing HybPiper
#' @param cpu Integer; number of CPUs for HybPiper
#' @param keep_fastq Logical; if TRUE, keep extracted FASTQ files in
#'   `output_dir`; otherwise remove them after assembly
#' @param mapper Character; one of "bwa" or "blast"
#'
#' @return Named list with paths for downstream HybPhaser functions:
#'   `hybpiper_dir`, `targets_file`, and `namelist`
#' @export
run_hybpiper_test_dataset <- function(
  output_dir,
  samples = c("EG30", "EG98", "MWL2"),
  conda_env = "hybpiper_env",
  cpu = 1,
  keep_fastq = FALSE,
  mapper = c("bwa", "blast")
) {
  mapper <- match.arg(mapper)

  if (
    !is.character(samples) || length(samples) == 0 || any(nchar(samples) == 0)
  ) {
    stop("'samples' must be a non-empty character vector")
  }

  downloaded <- .download_hybpiper_test_dataset(output_dir)
  targets_file <- downloaded$targets_file
  tar_reads <- downloaded$reads_archive

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  targets_out <- file.path(output_dir, "test_targets.fasta")
  file.copy(targets_file, targets_out, overwrite = TRUE)

  namelist_out <- file.path(output_dir, "namelist.txt")
  writeLines(samples, namelist_out)

  tar_entries <- utils::untar(tar_reads, list = TRUE)

  for (sample_name in samples) {
    f_read <- paste0(sample_name, "_R1_test.fastq")
    r_read <- paste0(sample_name, "_R2_test.fastq")

    if (!all(c(f_read, r_read) %in% tar_entries)) {
      stop(
        "Sample not found in downloaded HybPiper test dataset: ",
        sample_name
      )
    }

    utils::untar(
      tarfile = tar_reads,
      files = c(f_read, r_read),
      exdir = output_dir
    )

    f_path <- file.path(output_dir, f_read)
    r_path <- file.path(output_dir, r_read)

    if (!file.exists(f_path) || !file.exists(r_path)) {
      stop("Could not extract FASTQ files for sample: ", sample_name)
    }

    hybpiper_assemble(
      sample_name = sample_name,
      f_read = f_path,
      r_read = r_path,
      targets_file = targets_out,
      wd = output_dir,
      mapper = mapper,
      cpu = cpu,
      conda_env = conda_env
    )

    if (!keep_fastq) {
      unlink(c(f_path, r_path))
    }
  }

  message("HybPiper test dataset run complete")

  list(
    hybpiper_dir = output_dir,
    targets_file = targets_out,
    namelist = namelist_out,
    samples = samples
  )
}
