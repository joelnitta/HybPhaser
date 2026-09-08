#' Generate Consensus Sequences
#'
#' Generates consensus sequences from HybPiper output by remapping reads to
#' contigs and calling variants, by running a bundled bash script. By default
#' the script runs inside the `rhybphaser` Docker container; set
#' `engine = "local"` to run it directly against tools on `PATH`
#' (BWA, SAMtools, BCFtools, BBMap, GNU parallel).
#'
#' @param hybpiper_dir Path to HybPiper output directory (on host machine)
#' @param output_dir Path to output directory for HybPhaser results (on host
#'   machine). Will be created if it doesn't exist.
#' @param engine One of `"docker"` (default) or `"local"`.
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
#'   "joelnitta/rhybphaser:latest"
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
  engine = c("docker", "local"),
  docker_image = "joelnitta/rhybphaser:latest",
  pull_image = FALSE
) {
  engine <- match.arg(engine)

  # Validate inputs first so path mistakes are reported without requiring
  # Docker or local tools.
  if (!dir.exists(hybpiper_dir)) {
    stop("HybPiper directory not found: ", hybpiper_dir)
  }

  if (!is.null(namelist) && !file.exists(namelist)) {
    stop("Namelist file not found: ", namelist)
  }

  # Flags common to both engines.
  common_flags <- c(
    "-t", as.character(threads),
    "-d", as.character(min_depth),
    "-f", as.character(min_allele_freq),
    "-a", as.character(min_allele_count)
  )
  if (intronerate) {
    common_flags <- c(common_flags, "-i")
  }
  if (cleanup) {
    common_flags <- c(common_flags, "-c")
  }

  if (engine == "local") {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      message("Created output directory: ", output_dir)
    }
    script <- hybphaser_scripts("1_generate_consensus_sequences.sh")
    args <- c(
      "-p", normalizePath(hybpiper_dir),
      "-o", normalizePath(output_dir),
      common_flags
    )
    if (!is.null(namelist)) {
      args <- c(args, "-n", normalizePath(namelist))
    } else if (!is.null(sample)) {
      args <- c(args, "-s", sample)
    }

    message("Running consensus sequence generation locally...")
    message("Command: ", script, " ", paste(args, collapse = " "))
    result <- system2(script, args = args, stdout = "", stderr = "")

    if (result == 0) {
      message("Consensus sequences generated successfully")
      message("Output directory: ", output_dir)
    } else {
      warning(
        "Consensus sequence generation failed with exit code: ",
        result
      )
    }
    return(invisible(result))
  }

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

  # Get absolute paths
  hybpiper_abs <- .normalize_docker_path(hybpiper_dir)
  output_abs <- .normalize_docker_path(output_dir)
  hybpiper_norm <- normalizePath(hybpiper_abs, mustWork = FALSE)
  output_norm <- normalizePath(output_abs, mustWork = FALSE)

  # Determine Docker volume mounts and the container-side output path.
  #
  # Docker Desktop on macOS (VirtioFS) has a bug: a directory created on the
  # host immediately before a container starts appears in the parent directory
  # listing inside the container, but the directory itself cannot be opened
  # for writes by the container process (ENOENT / "No such file or directory").
  #
  # The reliable fix when output_dir is NESTED inside hybpiper_dir:
  #   - mount only hybpiper_dir as /data/hybpiper
  #   - reference the output as /data/hybpiper/<rel>
  #   - let Docker's own mkdir -p create the subdirectory (Docker-created
  #     directories are reliably accessible)
  #
  # When output_dir is OUTSIDE hybpiper_dir a separate bind mount is used and
  # the directory is pre-created on the host first (standard approach).
  rel <- .relative_path(output_norm, hybpiper_norm)

  if (!is.null(rel)) {
    # output inside hybpiper_dir: single bind mount, Docker creates subdir
    volumes <- c("/data/hybpiper" = hybpiper_abs)
    container_output <- paste0("/data/hybpiper/", rel)
  } else {
    # output outside hybpiper_dir: pre-create + separate bind mount
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      message("Created output directory: ", output_dir)
    }
    volumes <- c(
      "/data/hybpiper" = hybpiper_abs,
      "/data/output" = output_abs
    )
    container_output <- "/data/output"
  }

  # Build command
  cmd <- c(
    "/opt/rhybphaser/1_generate_consensus_sequences.sh",
    "-p", "/data/hybpiper",
    "-o", container_output,
    common_flags
  )

  if (!is.null(namelist)) {
    namelist_abs <- .normalize_docker_path(namelist)
    volumes["/data/namelist.txt"] <- namelist_abs
    cmd <- c(cmd, "-n", "/data/namelist.txt")
  } else if (!is.null(sample)) {
    cmd <- c(cmd, "-s", sample)
  }

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
    warning(
      "Consensus sequence generation failed with exit code: ",
      result
    )
  }

  invisible(result)
}

#' Extract Mapped Reads
#'
#' Extracts and concatenates reads that mapped to target sequences from
#' HybPiper or HybPhaser output, by running a bundled bash script. By default
#' the script runs inside the `rhybphaser` Docker container; set
#' `engine = "local"` to run it directly.
#'
#' @param base_dir Path to base directory containing sample data (HybPiper
#'   or HybPhaser output)
#' @param output_dir Path to output directory for extracted reads
#' @param namelist Path to file containing sample names, one per line. If
#'   NULL, processes all samples found.
#' @param remove_duplicates Logical; if TRUE, remove duplicate sequences
#'   (not just duplicate names)
#' @param engine One of `"docker"` (default) or `"local"`.
#' @param docker_image Docker image name. Default is
#'   "joelnitta/rhybphaser:latest"
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
  engine = c("docker", "local"),
  docker_image = "joelnitta/rhybphaser:latest",
  pull_image = FALSE
) {
  engine <- match.arg(engine)

  # Validate inputs first so path mistakes are reported without requiring
  # Docker or local tools.
  if (!dir.exists(base_dir)) {
    stop("Base directory not found: ", base_dir)
  }

  if (!is.null(namelist) && !file.exists(namelist)) {
    stop("Namelist file not found: ", namelist)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }

  if (engine == "local") {
    script <- hybphaser_scripts("2_extract_mapped_reads.sh")
    args <- c(
      "-b", normalizePath(base_dir),
      "-o", normalizePath(output_dir)
    )
    if (remove_duplicates) {
      args <- c(args, "-s")
    }
    if (!is.null(namelist)) {
      args <- c(args, "-n", normalizePath(namelist))
    }

    message("Extracting mapped reads locally...")
    message("Command: ", script, " ", paste(args, collapse = " "))
    result <- system2(script, args = args, stdout = "", stderr = "")

    if (result == 0) {
      message("Mapped reads extracted successfully")
      message("Output directory: ", output_dir)
    } else {
      warning("Extraction failed with exit code: ", result)
    }
    return(invisible(result))
  }

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

  # Get absolute paths
  base_abs <- .normalize_docker_path(base_dir)
  output_abs <- .normalize_docker_path(output_dir)

  # Set up volume mounts
  volumes <- c(
    "/data/base" = base_abs,
    "/data/output" = output_abs
  )

  # Build command
  cmd <- c("/opt/rhybphaser/2_extract_mapped_reads.sh")
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
