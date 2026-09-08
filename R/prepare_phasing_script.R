#' Run BBSplit Phasing
#'
#' Generates an executable bash script for the BBSplit phasing step and
#' executes it.
#'
#' @param path_to_phasing_folder Path to phasing output folder.
#' @param csv_file_with_phasing_prep_info CSV with sample and reference
#'   information. Expected columns are `sample` or `samples`, then
#'   `ref1, abb1, ref2, abb2, ...`.
#' @param path_to_read_files_phasing Folder containing read files for phasing.
#' @param read_type_4phasing Read type: `"paired-end"` or `"single-end"`.
#' @param ID_read_pair1 Identifier for read pair 1 (paired-end only).
#' @param ID_read_pair2 Identifier for read pair 2 (paired-end only).
#' @param reference_sequence_folder Folder containing reference sequences.
#' @param folder_for_phased_reads Optional output folder for phased reads.
#'   If empty, uses `path_to_phasing_folder/phased_reads`.
#' @param folder_for_phasing_stats Optional output folder for phasing stats.
#'   If empty, uses `path_to_phasing_folder/phasing_stats`.
#' @param path_to_bbmap_executables Optional path to BBMap executables.
#'   If empty, `bbsplit.sh` is expected on `PATH`.
#' @param no_of_threads_phasing Number of threads for BBSplit. Use `0` or
#'   `"auto"` to omit thread argument.
#' @param java_memory_usage_phasing Optional Java memory (e.g., `"2G"`).
#' @param engine One of `"docker"` (default) or `"local"`. With `"docker"`,
#'   BBSplit runs in the `rhybphaser` container. With `"local"`, the generated
#'   script runs `bbsplit.sh` from `PATH` (or `path_to_bbmap_executables`);
#'   this needs a POSIX shell and is not available on native Windows.
#' @param docker_image Docker image used when `engine = "docker"`.
#' @param pull_image Logical; if `TRUE`, pull the Docker image when missing.
#'
#' @return Invisibly returns a list with script path, commands, output folders,
#'   selected samples, and run status.
#' @export
run_phasing <- function(
  path_to_phasing_folder,
  csv_file_with_phasing_prep_info,
  path_to_read_files_phasing,
  read_type_4phasing = c("paired-end", "single-end"),
  ID_read_pair1 = "",
  ID_read_pair2 = "",
  reference_sequence_folder,
  folder_for_phased_reads = "",
  folder_for_phasing_stats = "",
  path_to_bbmap_executables = "",
  no_of_threads_phasing = 1,
  java_memory_usage_phasing = "",
  engine = c("docker", "local"),
  docker_image = "joelnitta/rhybphaser:latest",
  pull_image = FALSE
) {
  engine <- match.arg(engine)

  validate_paths(
    paths = list(
      "csv_file_with_phasing_prep_info" = csv_file_with_phasing_prep_info,
      "path_to_read_files_phasing" = path_to_read_files_phasing,
      "reference_sequence_folder" = reference_sequence_folder
    ),
    must_exist = TRUE,
    type = c("file", "dir", "dir")
  )

  read_type_4phasing <- match.arg(read_type_4phasing)

  if (read_type_4phasing == "paired-end") {
    if (ID_read_pair1 == "" || ID_read_pair2 == "") {
      stop("ID_read_pair1 and ID_read_pair2 are required for paired-end reads")
    }
  }

  prep_phasing <- utils::read.csv(
    csv_file_with_phasing_prep_info,
    header = TRUE,
    stringsAsFactors = FALSE
  )
  prep_phasing[is.na(prep_phasing)] <- ""

  sample_col <- detect_phasing_sample_column(prep_phasing)
  samples_to_phase <- prep_phasing[[sample_col]]

  refseqs_fullpath <- list.files(reference_sequence_folder, full.names = TRUE)
  refseqs <- basename(refseqs_fullpath)
  refseqs_samplenames <- gsub(
    "(_intronerated)*_consensus\\.fasta|(_intronerated)*_contig\\.fasta",
    "",
    refseqs
  )

  if (folder_for_phased_reads == "") {
    folder_for_phased_reads <- file.path(path_to_phasing_folder, "phased_reads")
  }
  dir.create(folder_for_phased_reads, recursive = TRUE, showWarnings = FALSE)

  if (folder_for_phasing_stats == "") {
    folder_for_phasing_stats <- file.path(
      path_to_phasing_folder,
      "phasing_stats"
    )
  }
  dir.create(folder_for_phasing_stats, recursive = TRUE, showWarnings = FALSE)

  ref_commands <- build_phasing_reference_commands(
    prep_phasing = prep_phasing,
    sample_col = sample_col,
    refseqs_fullpath = refseqs_fullpath,
    refseqs_samplenames = refseqs_samplenames
  )

  threadtext <- ""
  if (
    !(identical(no_of_threads_phasing, 0) ||
      identical(no_of_threads_phasing, "auto"))
  ) {
    threadtext <- paste0(" threads=", no_of_threads_phasing)
  }

  pXmx <- ""
  if (java_memory_usage_phasing != "") {
    pXmx <- paste0(" -Xmx", java_memory_usage_phasing)
  }

  # Resolve a local bbsplit.sh: used to run the script when engine = "local",
  # and baked into the generated script as a record of what was run.
  bbsplit_sh <- "bbsplit.sh"
  if (path_to_bbmap_executables != "") {
    bbsplit_sh <- file.path(path_to_bbmap_executables, "bbsplit.sh")
    if (!file.exists(bbsplit_sh)) {
      stop("Could not find bbsplit.sh at: ", bbsplit_sh)
    }
  } else {
    local_bbsplit <- Sys.which("bbsplit.sh")
    if (nzchar(local_bbsplit)) {
      bbsplit_sh <- unname(local_bbsplit)
    }
  }

  reads <- list.files(path_to_read_files_phasing, full.names = TRUE)

  phasing_commands <- build_phasing_commands(
    prep_phasing = prep_phasing,
    sample_col = sample_col,
    ref_commands = ref_commands,
    reads = reads,
    read_type_4phasing = read_type_4phasing,
    ID_read_pair1 = ID_read_pair1,
    ID_read_pair2 = ID_read_pair2,
    bbsplit_sh = bbsplit_sh,
    threadtext = threadtext,
    folder_for_phased_reads = folder_for_phased_reads,
    folder_for_phasing_stats = folder_for_phasing_stats,
    pXmx = pXmx
  )

  phasing_script_file <- file.path(
    path_to_phasing_folder,
    "run_bbsplit4phasing.sh"
  )
  writeLines(c("#!/bin/bash", phasing_commands), phasing_script_file)
  Sys.chmod(phasing_script_file, mode = "0755")

  if (engine == "local") {
    if (identical(bbsplit_sh, "bbsplit.sh") || !file.exists(bbsplit_sh)) {
      stop(
        "engine = \"local\" but bbsplit.sh was not found. Install BBMap and ",
        "put it on PATH, set path_to_bbmap_executables, or use ",
        "engine = \"docker\"."
      )
    }
    run_status <- system2(phasing_script_file, stdout = "", stderr = "")
  } else {
    if (!check_docker(quiet = TRUE)) {
      stop("engine = \"docker\" but Docker is not available")
    }
    if (!check_docker_image(docker_image, pull = pull_image)) {
      stop(
        "Docker image '",
        docker_image,
        "' not found. Set pull_image = TRUE to download it."
      )
    }

    run_status <- run_phasing_commands_in_docker(
      commands = phasing_commands,
      bbsplit_sh = bbsplit_sh,
      reference_sequence_folder = reference_sequence_folder,
      path_to_read_files_phasing = path_to_read_files_phasing,
      folder_for_phased_reads = folder_for_phased_reads,
      folder_for_phasing_stats = folder_for_phasing_stats,
      docker_image = docker_image
    )
  }

  invisible(list(
    script_path = phasing_script_file,
    commands = phasing_commands,
    folder_for_phased_reads = folder_for_phased_reads,
    folder_for_phasing_stats = folder_for_phasing_stats,
    samples = samples_to_phase,
    run_status = run_status
  ))
}


#' Run BBSplit Phasing Using a HybPhaser Config File
#'
#' Reads required phasing settings from `config.txt` and runs
#' `run_phasing()`.
#'
#' @param config_file Path to HybPhaser configuration file.
#' @param engine One of `"docker"` (default) or `"local"`; passed to
#'   [run_phasing()].
#' @param pull_image Logical; if `TRUE`, pull the Docker image when missing.
#'
#' @return Invisibly returns the same object as `run_phasing()`.
#' @export
run_phasing_from_config <- function(
  config_file = "./config.txt",
  engine = c("docker", "local"),
  pull_image = FALSE
) {
  engine <- match.arg(engine)
  required <- c(
    "path_to_phasing_folder",
    "csv_file_with_phasing_prep_info",
    "path_to_read_files_phasing",
    "read_type_4phasing",
    "ID_read_pair1",
    "ID_read_pair2",
    "reference_sequence_folder",
    "folder_for_phased_reads",
    "folder_for_phasing_stats",
    "path_to_bbmap_executables",
    "no_of_threads_phasing",
    "java_memory_usage_phasing"
  )

  cfg <- read_config(config_file, required_vars = required)

  run_phasing(
    path_to_phasing_folder = cfg$path_to_phasing_folder,
    csv_file_with_phasing_prep_info = cfg$csv_file_with_phasing_prep_info,
    path_to_read_files_phasing = cfg$path_to_read_files_phasing,
    read_type_4phasing = cfg$read_type_4phasing,
    ID_read_pair1 = cfg$ID_read_pair1,
    ID_read_pair2 = cfg$ID_read_pair2,
    reference_sequence_folder = cfg$reference_sequence_folder,
    folder_for_phased_reads = cfg$folder_for_phased_reads,
    folder_for_phasing_stats = cfg$folder_for_phasing_stats,
    path_to_bbmap_executables = cfg$path_to_bbmap_executables,
    no_of_threads_phasing = cfg$no_of_threads_phasing,
    java_memory_usage_phasing = cfg$java_memory_usage_phasing,
    engine = engine,
    pull_image = pull_image
  )
}


# Detect sample column name in phasing prep table.
detect_phasing_sample_column <- function(prep_phasing) {
  cols <- colnames(prep_phasing)
  if ("samples" %in% cols) {
    return("samples")
  }
  if ("sample" %in% cols) {
    return("sample")
  }
  cols[[1]]
}


# Build reference argument fragments per phasing sample.
build_phasing_reference_commands <- function(
  prep_phasing,
  sample_col,
  refseqs_fullpath,
  refseqs_samplenames
) {
  ref_cols <- setdiff(colnames(prep_phasing), sample_col)
  if (length(ref_cols) < 2) {
    stop("Phasing prep CSV must include ref/abb column pairs")
  }

  ref_commands <- character(nrow(prep_phasing))

  for (i in seq_len(nrow(prep_phasing))) {
    ref_co <- character(0)

    for (j in seq(1, length(ref_cols), by = 2)) {
      ref_col <- ref_cols[[j]]
      abb_col <- if (j + 1 <= length(ref_cols)) ref_cols[[j + 1]] else ""

      ref_name <- prep_phasing[[ref_col]][[i]]
      if (is.null(ref_name) || ref_name == "") {
        next
      }

      if (identical(abb_col, "")) {
        stop("Missing abbreviation column for reference column: ", ref_col)
      }

      ref_abb <- prep_phasing[[abb_col]][[i]]
      if (is.null(ref_abb) || ref_abb == "") {
        stop("Missing abbreviation for sample row ", i, " reference ", ref_name)
      }

      ref_file <- refseqs_fullpath[match(ref_name, refseqs_samplenames)]
      if (is.na(ref_file)) {
        stop("Reference sequence not found for: ", ref_name)
      }

      ref_co <- c(ref_co, paste0("ref_", ref_abb, "=", ref_file))
    }

    if (length(ref_co) == 0) {
      stop("No references found for sample row ", i)
    }

    ref_commands[[i]] <- paste(ref_co, collapse = " ")
  }

  ref_commands
}


# Build BBSplit phasing commands.
build_phasing_commands <- function(
  prep_phasing,
  sample_col,
  ref_commands,
  reads,
  read_type_4phasing,
  ID_read_pair1,
  ID_read_pair2,
  bbsplit_sh,
  threadtext,
  folder_for_phased_reads,
  folder_for_phasing_stats,
  pXmx
) {
  commands <- character(nrow(prep_phasing))
  read_names <- basename(reads)

  for (i in seq_len(nrow(prep_phasing))) {
    sample_name <- prep_phasing[[sample_col]][[i]]

    if (read_type_4phasing == "paired-end") {
      read_1_name <- paste0(sample_name, ID_read_pair1)
      read_2_name <- paste0(sample_name, ID_read_pair2)
      read_file_1 <- reads[match(read_1_name, read_names)]
      read_file_2 <- reads[match(read_2_name, read_names)]

      if (is.na(read_file_1) || is.na(read_file_2)) {
        stop("Could not find paired-end read files for sample: ", sample_name)
      }

      cmd <- paste0(
        bbsplit_sh,
        " ambiguous=all ambiguous2=all",
        threadtext,
        " ",
        ref_commands[[i]],
        " in=",
        read_file_1,
        " in2=",
        read_file_2,
        " basename=",
        file.path(folder_for_phased_reads, paste0(sample_name, "_to_%.fastq")),
        " refstats=",
        file.path(
          folder_for_phasing_stats,
          paste0(sample_name, "_phasing-stats.txt")
        ),
        pXmx
      )
    } else {
      read_file <- reads[grep(sample_name, reads, fixed = TRUE)]
      if (length(read_file) == 0) {
        stop("Could not find single-end read file for sample: ", sample_name)
      }
      if (length(read_file) > 1) {
        stop("Multiple single-end read files matched sample: ", sample_name)
      }

      cmd <- paste0(
        bbsplit_sh,
        " ambiguous=all ambiguous2=all",
        threadtext,
        " ",
        ref_commands[[i]],
        " in=",
        read_file,
        " basename=",
        file.path(folder_for_phased_reads, paste0(sample_name, "_to_%.fastq")),
        " refstats=",
        file.path(
          folder_for_phasing_stats,
          paste0(sample_name, "_phasing-stats.txt")
        ),
        pXmx
      )
    }

    commands[[i]] <- cmd
  }

  commands
}


# Execute prepared phasing commands in Docker by remapping host paths.
run_phasing_commands_in_docker <- function(
  commands,
  bbsplit_sh,
  reference_sequence_folder,
  path_to_read_files_phasing,
  folder_for_phased_reads,
  folder_for_phasing_stats,
  docker_image
) {
  volumes <- c(
    "/data/refs" = reference_sequence_folder,
    "/data/reads" = path_to_read_files_phasing,
    "/data/phased" = folder_for_phased_reads,
    "/data/stats" = folder_for_phasing_stats
  )

  run_status <- 0L
  for (cmd in commands) {
    docker_cmd <- cmd
    docker_cmd <- sub("^\\S+", "bbsplit.sh", docker_cmd)
    if (bbsplit_sh != "bbsplit.sh") {
      docker_cmd <- gsub(bbsplit_sh, "bbsplit.sh", docker_cmd, fixed = TRUE)
    }
    docker_cmd <- gsub(
      reference_sequence_folder,
      "/data/refs",
      docker_cmd,
      fixed = TRUE
    )
    docker_cmd <- gsub(
      path_to_read_files_phasing,
      "/data/reads",
      docker_cmd,
      fixed = TRUE
    )
    docker_cmd <- gsub(
      folder_for_phased_reads,
      "/data/phased",
      docker_cmd,
      fixed = TRUE
    )
    docker_cmd <- gsub(
      folder_for_phasing_stats,
      "/data/stats",
      docker_cmd,
      fixed = TRUE
    )

    cmd_args <- strsplit(trimws(docker_cmd), "[[:space:]]+")[[1]]
    status_i <- .run_docker(
      cmd = cmd_args,
      volumes = volumes,
      image = docker_image,
      stdout = "",
      stderr = ""
    )

    if (status_i != 0) {
      run_status <- as.integer(status_i)
      break
    }
  }

  run_status
}
