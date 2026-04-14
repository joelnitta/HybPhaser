#' Prepare BBSplit Script for Clade Association
#'
#' Generates an executable bash script to run BBSplit for clade association.
#' Optionally executes the generated script directly from R.
#'
#' @param path_to_clade_association_folder Path to clade association output
#'   folder.
#' @param csv_file_with_clade_reference_names CSV file with reference sample
#'   names and abbreviations (first two columns).
#' @param path_to_reference_sequences Folder with reference sequence files.
#' @param path_to_read_files_cladeassociation Folder with read files for clade
#'   association.
#' @param read_type_cladeassociation Read type: `"single-end"` or
#'   `"paired-end"`.
#' @param ID_read_pair1 Identifier for read 1 files (paired-end only).
#' @param ID_read_pair2 Identifier for read 2 files (paired-end only).
#' @param file_with_samples_included Optional text file with sample names to
#'   include. Use `""` or `"none"` to include all samples.
#' @param path_to_bbmap Optional path to BBMap binaries. If empty, `bbsplit.sh`
#'   is expected on `PATH`.
#' @param no_of_threads Number of threads for BBSplit. Use `0` or `"auto"`
#'   to omit thread argument.
#' @param run_clade_association_mapping_in_R Logical; if `TRUE`, run the
#'   generated script in R.
#' @param java_memory_usage_clade_association Optional Java memory text for
#'   BBSplit (e.g., `"2G"`, `"512m"`).
#' @param docker_fallback Logical; if `TRUE` and local `bbsplit.sh` is not
#'   available, execute BBSplit commands via Docker. Default `TRUE`.
#' @param docker_image Docker image used for fallback execution. Default is
#'   `"joelnitta/hybphaser:latest"`.
#' @param pull_image Logical; if `TRUE`, pull Docker image when not available.
#'
#' @return Invisibly returns a list with script path, generated commands,
#'   selected read files, stats folder, and run status.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- prepare_bbsplit_script(
#'   path_to_clade_association_folder = "04_clade_association",
#'   csv_file_with_clade_reference_names = "clade_references.csv",
#'   path_to_reference_sequences = "03_sequence_lists/samples_consensus",
#'   path_to_read_files_cladeassociation = "mapped_reads",
#'   read_type_cladeassociation = "single-end"
#' )
#' }
prepare_bbsplit_script <- function(
  path_to_clade_association_folder,
  csv_file_with_clade_reference_names,
  path_to_reference_sequences,
  path_to_read_files_cladeassociation,
  read_type_cladeassociation = c("single-end", "paired-end"),
  ID_read_pair1 = "",
  ID_read_pair2 = "",
  file_with_samples_included = "",
  path_to_bbmap = "",
  no_of_threads = 1,
  run_clade_association_mapping_in_R = FALSE,
  java_memory_usage_clade_association = "",
  docker_fallback = TRUE,
  docker_image = "joelnitta/hybphaser:latest",
  pull_image = FALSE
) {
  validate_paths(
    paths = list(
      "csv_file_with_clade_reference_names" = csv_file_with_clade_reference_names,
      "path_to_reference_sequences" = path_to_reference_sequences,
      "path_to_read_files_cladeassociation" = path_to_read_files_cladeassociation
    ),
    must_exist = TRUE,
    type = c("file", "dir", "dir")
  )

  if (
    !(file_with_samples_included %in% c("", "none")) &&
      !file.exists(file_with_samples_included)
  ) {
    stop(
      "file_with_samples_included not found: ",
      file_with_samples_included
    )
  }

  read_type_cladeassociation <- match.arg(read_type_cladeassociation)

  if (read_type_cladeassociation == "paired-end") {
    if (ID_read_pair1 == "" || ID_read_pair2 == "") {
      stop("ID_read_pair1 and ID_read_pair2 are required for paired-end reads")
    }
  }

  dir.create(
    path_to_clade_association_folder,
    recursive = TRUE,
    showWarnings = FALSE
  )
  folder_bbsplit_stats <- file.path(
    path_to_clade_association_folder,
    "bbsplit_stats"
  )
  dir.create(folder_bbsplit_stats, recursive = TRUE, showWarnings = FALSE)

  read_files <- list_read_files_for_clade(
    path_to_read_files_cladeassociation = path_to_read_files_cladeassociation,
    read_type_cladeassociation = read_type_cladeassociation,
    ID_read_pair1 = ID_read_pair1
  )

  read_files <- filter_read_files_by_sample_list(
    read_files = read_files,
    file_with_samples_included = file_with_samples_included,
    ID_read_pair1 = ID_read_pair1
  )

  ref_command <- build_bbsplit_reference_command(
    csv_file_with_clade_reference_names = csv_file_with_clade_reference_names,
    path_to_reference_sequences = path_to_reference_sequences
  )

  threadtext <- ""
  if (!(identical(no_of_threads, 0) || identical(no_of_threads, "auto"))) {
    threadtext <- paste0(" threads=", no_of_threads)
  }

  caXmx <- ""
  if (java_memory_usage_clade_association != "") {
    caXmx <- paste0(" -Xmx", java_memory_usage_clade_association)
  }

  bbsplit_sh <- "bbsplit.sh"
  if (path_to_bbmap != "") {
    bbsplit_sh <- file.path(path_to_bbmap, "bbsplit.sh")
    if (!file.exists(bbsplit_sh)) {
      stop("Could not find bbsplit.sh at: ", bbsplit_sh)
    }
  }

  local_bbsplit <- ""
  if (path_to_bbmap == "") {
    local_bbsplit <- Sys.which("bbsplit.sh")
    if (local_bbsplit != "") {
      bbsplit_sh <- local_bbsplit
    }
  }

  script_path <- file.path(
    path_to_clade_association_folder,
    "run_bbsplit4clade_association.sh"
  )

  commands <- build_bbsplit_commands(
    read_files = read_files,
    read_type_cladeassociation = read_type_cladeassociation,
    path_to_read_files_cladeassociation = path_to_read_files_cladeassociation,
    ID_read_pair1 = ID_read_pair1,
    ID_read_pair2 = ID_read_pair2,
    bbsplit_sh = bbsplit_sh,
    ref_command = ref_command,
    threadtext = threadtext,
    folder_bbsplit_stats = folder_bbsplit_stats,
    caXmx = caXmx
  )

  writeLines(c("#!/bin/bash", commands), script_path)
  Sys.chmod(script_path, mode = "0755")

  run_status <- NA_integer_
  if (isTRUE(run_clade_association_mapping_in_R)) {
    local_exec_available <- path_to_bbmap != "" || local_bbsplit != ""

    if (local_exec_available) {
      run_status <- system2(script_path, stdout = "", stderr = "")
    } else if (isTRUE(docker_fallback)) {
      if (!check_docker(quiet = TRUE)) {
        stop("Docker is not available for BBSplit fallback execution")
      }
      if (!check_docker_image(docker_image, pull = pull_image)) {
        stop(
          "Docker image '",
          docker_image,
          "' not found. Set pull_image = TRUE to download it."
        )
      }

      run_status <- run_bbsplit_commands_in_docker(
        commands = commands,
        bbsplit_sh = bbsplit_sh,
        path_to_reference_sequences = path_to_reference_sequences,
        path_to_read_files_cladeassociation = path_to_read_files_cladeassociation,
        folder_bbsplit_stats = folder_bbsplit_stats,
        docker_image = docker_image
      )
    } else {
      stop(
        "bbsplit.sh not found on PATH. Install BBMap, provide path_to_bbmap,",
        " or set docker_fallback = TRUE."
      )
    }
  }

  invisible(list(
    script_path = script_path,
    commands = commands,
    read_files = read_files,
    stats_folder = folder_bbsplit_stats,
    run_status = run_status
  ))
}


# Execute prepared BBSplit commands in Docker by remapping host paths.
run_bbsplit_commands_in_docker <- function(
  commands,
  bbsplit_sh,
  path_to_reference_sequences,
  path_to_read_files_cladeassociation,
  folder_bbsplit_stats,
  docker_image
) {
  volumes <- c(
    "/data/refs" = path_to_reference_sequences,
    "/data/reads" = path_to_read_files_cladeassociation,
    "/data/stats" = folder_bbsplit_stats
  )

  run_status <- 0L
  for (cmd in commands) {
    docker_cmd <- cmd
    docker_cmd <- sub("^\\S+", "bbsplit.sh", docker_cmd)
    if (bbsplit_sh != "bbsplit.sh") {
      docker_cmd <- gsub(bbsplit_sh, "bbsplit.sh", docker_cmd, fixed = TRUE)
    }
    docker_cmd <- gsub(
      path_to_reference_sequences,
      "/data/refs",
      docker_cmd,
      fixed = TRUE
    )
    docker_cmd <- gsub(
      path_to_read_files_cladeassociation,
      "/data/reads",
      docker_cmd,
      fixed = TRUE
    )
    docker_cmd <- gsub(
      folder_bbsplit_stats,
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


#' Prepare BBSplit Script Using a HybPhaser Config File
#'
#' Reads clade-association variables from a HybPhaser config file and prepares
#' the BBSplit script.
#'
#' @param config_file Path to HybPhaser configuration file.
#'
#' @return Invisibly returns a list with script path, commands, selected read
#'   files, stats folder, and run status.
#' @export
prepare_bbsplit_script_from_config <- function(config_file = "./config.txt") {
  required <- c(
    "path_to_clade_association_folder",
    "csv_file_with_clade_reference_names",
    "path_to_reference_sequences",
    "path_to_read_files_cladeassociation",
    "read_type_cladeassociation",
    "ID_read_pair1",
    "ID_read_pair2",
    "file_with_samples_included",
    "path_to_bbmap",
    "no_of_threads_clade_association",
    "run_clade_association_mapping_in_R",
    "java_memory_usage_clade_association"
  )

  cfg <- read_config(config_file, required_vars = required)

  prepare_bbsplit_script(
    path_to_clade_association_folder = cfg$path_to_clade_association_folder,
    csv_file_with_clade_reference_names = cfg$csv_file_with_clade_reference_names,
    path_to_reference_sequences = cfg$path_to_reference_sequences,
    path_to_read_files_cladeassociation = cfg$path_to_read_files_cladeassociation,
    read_type_cladeassociation = cfg$read_type_cladeassociation,
    ID_read_pair1 = cfg$ID_read_pair1,
    ID_read_pair2 = cfg$ID_read_pair2,
    file_with_samples_included = cfg$file_with_samples_included,
    path_to_bbmap = cfg$path_to_bbmap,
    no_of_threads = cfg$no_of_threads_clade_association,
    run_clade_association_mapping_in_R = tolower(
      cfg$run_clade_association_mapping_in_R
    ) ==
      "yes",
    java_memory_usage_clade_association = cfg$java_memory_usage_clade_association
  )
}


# List read files according to single-end or paired-end mode.
list_read_files_for_clade <- function(
  path_to_read_files_cladeassociation,
  read_type_cladeassociation,
  ID_read_pair1
) {
  read_files <- list.files(
    path_to_read_files_cladeassociation,
    full.names = FALSE,
    include.dirs = FALSE
  )

  if (read_type_cladeassociation == "paired-end") {
    read_files <- read_files[grepl(ID_read_pair1, read_files, fixed = TRUE)]
  }

  read_files
}


# Filter read files to a subset of samples when a list is provided.
filter_read_files_by_sample_list <- function(
  read_files,
  file_with_samples_included,
  ID_read_pair1
) {
  if (file_with_samples_included %in% c("", "none")) {
    return(read_files)
  }

  samples_in <- readLines(file_with_samples_included, warn = FALSE)
  samples_in <- trimws(samples_in)
  samples_in <- samples_in[nzchar(samples_in)]

  read_files_noID <- gsub(ID_read_pair1, "", read_files, fixed = TRUE)
  read_files_noend <- gsub("[.]fast.*", "", read_files_noID)
  read_files[match(samples_in, read_files_noend)]
}


# Build BBSplit reference arguments from clade reference CSV.
build_bbsplit_reference_command <- function(
  csv_file_with_clade_reference_names,
  path_to_reference_sequences
) {
  sample_sequences <- list.files(path_to_reference_sequences)

  ref_samples <- utils::read.csv(
    csv_file_with_clade_reference_names,
    header = TRUE,
    stringsAsFactors = FALSE
  )

  if (ncol(ref_samples) < 2) {
    stop("csv_file_with_clade_reference_names must have at least 2 columns")
  }

  colnames(ref_samples)[1:2] <- c("samples", "abb")

  ref_samples_files <- sample_sequences[match(
    ref_samples$samples,
    gsub("(_intronerated|_consensus|_contig).*", "", sample_sequences)
  )]

  if (any(is.na(ref_samples_files))) {
    missing_samples <- ref_samples$samples[is.na(ref_samples_files)]
    stop(
      "Could not find reference sequence files for: ",
      paste(missing_samples, collapse = ", ")
    )
  }

  abb_values <- trimws(ref_samples$abb)
  if (all(abb_values == "")) {
    return(paste0(
      "ref=",
      paste(
        file.path(path_to_reference_sequences, ref_samples_files),
        collapse = ","
      )
    ))
  }

  paste(
    paste0(
      "ref_",
      ref_samples$abb,
      "=",
      file.path(path_to_reference_sequences, ref_samples_files)
    ),
    collapse = " "
  )
}


# Build per-sample BBSplit command lines.
build_bbsplit_commands <- function(
  read_files,
  read_type_cladeassociation,
  path_to_read_files_cladeassociation,
  ID_read_pair1,
  ID_read_pair2,
  bbsplit_sh,
  ref_command,
  threadtext,
  folder_bbsplit_stats,
  caXmx
) {
  commands <- character(0)

  for (read_file in stats::na.omit(read_files)) {
    input1 <- file.path(path_to_read_files_cladeassociation, read_file)

    if (read_type_cladeassociation == "paired-end") {
      input2 <- sub(ID_read_pair1, ID_read_pair2, input1, fixed = TRUE)
      sample_tag <- sub(ID_read_pair1, "", read_file, fixed = TRUE)

      cmd <- paste0(
        bbsplit_sh,
        " ",
        ref_command,
        " in=",
        input1,
        " in2=",
        input2,
        threadtext,
        " ambiguous=random ambiguous2=all refstats=",
        file.path(
          folder_bbsplit_stats,
          paste0(sample_tag, "_bbsplit-stats.txt")
        ),
        caXmx
      )
    } else {
      cmd <- paste0(
        bbsplit_sh,
        " ",
        ref_command,
        " in=",
        input1,
        threadtext,
        " ambiguous=random ambiguous2=all refstats=",
        file.path(
          folder_bbsplit_stats,
          paste0(read_file, "_bbsplit-stats.txt")
        ),
        caXmx
      )
    }

    commands <- c(commands, cmd)
  }

  commands
}
