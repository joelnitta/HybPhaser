#' Generate Consensus Sequences via Docker
#'
#' Generates consensus sequences from HybPiper output by remapping reads to
#' contigs and calling variants. Runs the HybPhaser bash script in a Docker
#' container.
#'
#' @param hybpiper_dir Path to HybPiper output directory (on host machine)
#' @param output_dir Path to output directory for HybPhaser results (on host
#'   machine). Will be created if it doesn't exist.
#' @param namelist Path to file containing sample names, one per line. If
#'   NULL, processes all samples in hybpiper_dir.
#' @param sample Single sample name to process. Ignored if namelist is
#'   provided.
#' @param intronerate Logical; if TRUE, use intronerated supercontigs
#' @param cleanup Logical; if TRUE, remove intermediate BAM and VCF files
#' @param threads Integer; number of CPU threads to use. Default is 1.
#' @param min_depth Integer; minimum read depth for variant calling.
#'   Default is 10.
#' @param min_allele_freq Numeric; minimum allele frequency (0-1) for
#'   calling heterozygous sites. Default is 0.15.
#' @param min_allele_count Integer; minimum number of reads supporting an
#'   allele. Default is 4.
#' @param docker_image Docker image name. Default is
#'   "joelnitta/hybphaser:latest"
#' @param pull_image Logical; if TRUE, pull Docker image if not found
#'
#' @return Exit code (0 for success, non-zero for failure)
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate consensus sequences for all samples
#' run_generate_consensus_sequences(
#'   hybpiper_dir = "path/to/hybpiper_output",
#'   output_dir = "path/to/hybphaser_output",
#'   namelist = "samples.txt",
#'   threads = 4
#' )
#'
#' # Process a single sample
#' run_generate_consensus_sequences(
#'   hybpiper_dir = "hybpiper_out",
#'   output_dir = "hybphaser_out",
#'   sample = "sample001"
#' )
#' }
run_generate_consensus_sequences <- function(
  hybpiper_dir,
  output_dir,
  namelist = NULL,
  sample = NULL,
  intronerate = FALSE,
  cleanup = FALSE,
  threads = 1,
  min_depth = 10,
  min_allele_freq = 0.15,
  min_allele_count = 4,
  docker_image = "joelnitta/hybphaser:latest",
  pull_image = FALSE
) {
  # Validate Docker
  if (!check_docker(quiet = TRUE)) {
    stop("Docker is not available. Please install and start Docker.")
  }

  if (!check_docker_image(docker_image, pull = pull_image)) {
    stop(
      "Docker image '",
      docker_image,
      "' not found. ",
      "Set pull_image = TRUE to download it."
    )
  }

  # Validate input directories exist
  if (!dir.exists(hybpiper_dir)) {
    stop("HybPiper directory not found: ", hybpiper_dir)
  }

  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }

  # Validate namelist if provided
  if (!is.null(namelist) && !file.exists(namelist)) {
    stop("Namelist file not found: ", namelist)
  }

  # Get absolute paths
  hybpiper_abs <- .normalize_docker_path(hybpiper_dir)
  output_abs <- .normalize_docker_path(output_dir)

  # Set up volume mounts
  volumes <- c(
    "/data/hybpiper" = hybpiper_abs,
    "/data/output" = output_abs
  )

  # Build command arguments
  cmd <- c("/opt/hybphaser/1_generate_consensus_sequences.sh")

  # Add options
  cmd <- c(cmd, "-p", "/data/hybpiper")
  cmd <- c(cmd, "-o", "/data/output")
  cmd <- c(cmd, "-t", as.character(threads))
  cmd <- c(cmd, "-d", as.character(min_depth))
  cmd <- c(cmd, "-f", as.character(min_allele_freq))
  cmd <- c(cmd, "-a", as.character(min_allele_count))

  if (intronerate) {
    cmd <- c(cmd, "-i")
  }

  if (cleanup) {
    cmd <- c(cmd, "-c")
  }

  # Handle sample/namelist
  if (!is.null(namelist)) {
    # Copy namelist to a temp location and mount it
    namelist_abs <- .normalize_docker_path(namelist)
    volumes["/data/namelist.txt"] <- namelist_abs
    cmd <- c(cmd, "-n", "/data/namelist.txt")
  } else if (!is.null(sample)) {
    cmd <- c(cmd, "-s", sample)
  }

  # Run Docker command
  message("Running consensus sequence generation in Docker...")
  message("Command: ", paste(cmd, collapse = " "))

  result <- .run_docker(
    cmd = cmd,
    volumes = volumes,
    image = docker_image,
    stdout = "",
    stderr = ""
  )

  if (result == 0) {
    message("Consensus sequences generated successfully")
    message("Output directory: ", output_dir)
  } else {
    warning("Consensus sequence generation failed with exit code: ", result)
  }

  invisible(result)
}


#' Extract Mapped Reads via Docker
#'
#' Extracts and concatenates reads that mapped to target sequences from
#' HybPiper or HybPhaser output. Runs the HybPhaser bash script in a Docker
#' container.
#'
#' @param base_dir Path to base directory containing sample data (HybPiper
#'   or HybPhaser output)
#' @param output_dir Path to output directory for extracted reads
#' @param namelist Path to file containing sample names, one per line. If
#'   NULL, processes all samples found.
#' @param remove_duplicates Logical; if TRUE, remove duplicate sequences
#'   (not just duplicate names)
#' @param docker_image Docker image name. Default is
#'   "joelnitta/hybphaser:latest"
#' @param pull_image Logical; if TRUE, pull Docker image if not found
#'
#' @return Exit code (0 for success, non-zero for failure)
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract mapped reads from all samples
#' run_extract_mapped_reads(
#'   base_dir = "path/to/hybphaser_output",
#'   output_dir = "path/to/mapped_reads"
#' )
#'
#' # Extract with duplicate removal
#' run_extract_mapped_reads(
#'   base_dir = "hybphaser_output",
#'   output_dir = "mapped_reads",
#'   remove_duplicates = TRUE,
#'   namelist = "samples_to_process.txt"
#' )
#' }
run_extract_mapped_reads <- function(
  base_dir,
  output_dir,
  namelist = NULL,
  remove_duplicates = FALSE,
  docker_image = "joelnitta/hybphaser:latest",
  pull_image = FALSE
) {
  # Validate Docker
  if (!check_docker(quiet = TRUE)) {
    stop("Docker is not available. Please install and start Docker.")
  }

  if (!check_docker_image(docker_image, pull = pull_image)) {
    stop(
      "Docker image '",
      docker_image,
      "' not found. ",
      "Set pull_image = TRUE to download it."
    )
  }

  # Validate input directory exists
  if (!dir.exists(base_dir)) {
    stop("Base directory not found: ", base_dir)
  }

  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }

  # Validate namelist if provided
  if (!is.null(namelist) && !file.exists(namelist)) {
    stop("Namelist file not found: ", namelist)
  }

  # Get absolute paths
  base_abs <- .normalize_docker_path(base_dir)
  output_abs <- .normalize_docker_path(output_dir)

  # Set up volume mounts
  volumes <- c(
    "/data/base" = base_abs,
    "/data/output" = output_abs
  )

  # Build command
  cmd <- c("/opt/hybphaser/2_extract_mapped_reads.sh")
  cmd <- c(cmd, "-b", "/data/base")
  cmd <- c(cmd, "-o", "/data/output")

  if (remove_duplicates) {
    cmd <- c(cmd, "-s")
  }

  # Handle namelist
  if (!is.null(namelist)) {
    namelist_abs <- .normalize_docker_path(namelist)
    volumes["/data/namelist.txt"] <- namelist_abs
    cmd <- c(cmd, "-n", "/data/namelist.txt")
  }

  # Run Docker command
  message("Extracting mapped reads in Docker...")
  message("Command: ", paste(cmd, collapse = " "))

  result <- .run_docker(
    cmd = cmd,
    volumes = volumes,
    image = docker_image,
    stdout = "",
    stderr = ""
  )

  if (result == 0) {
    message("Mapped reads extracted successfully")
    message("Output directory: ", output_dir)
  } else {
    warning("Extraction failed with exit code: ", result)
  }

  invisible(result)
}
