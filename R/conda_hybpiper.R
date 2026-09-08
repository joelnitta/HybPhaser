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


#' Command hint shown when HybPiper is missing from a conda environment
#'
#' conda-forge is listed before bioconda (the bioconda-recommended channel
#' order); getting it wrong pulls mis-built dependencies such as GNU
#' `parallel`.
#'
#' @param env_name Character; conda environment name
#' @return Character; a one-line `conda install` suggestion
#' @noRd
.hybpiper_install_hint <- function(env_name) {
  paste0(
    "Install it with: conda install -n ",
    env_name,
    " -c conda-forge -c bioconda hybpiper"
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
      "'. ",
      .hybpiper_install_hint(env_name)
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
      "'. ",
      .hybpiper_install_hint(conda_env)
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


.read_fasta_records <- function(path) {
  lines <- readLines(path, warn = FALSE)

  if (length(lines) == 0) {
    return(data.frame(id = character(0), sequence = character(0)))
  }

  record_ids <- character(0)
  record_sequences <- character(0)
  current_id <- NULL
  current_sequence <- character(0)

  for (line in lines) {
    if (startsWith(line, ">")) {
      if (!is.null(current_id)) {
        record_ids <- c(record_ids, current_id)
        record_sequences <- c(
          record_sequences,
          paste0(current_sequence, collapse = "")
        )
      }

      current_id <- strsplit(sub("^>", "", line), "\\s+")[[1]][1]
      current_sequence <- character(0)
    } else {
      current_sequence <- c(current_sequence, trimws(line))
    }
  }

  if (!is.null(current_id)) {
    record_ids <- c(record_ids, current_id)
    record_sequences <- c(
      record_sequences,
      paste0(current_sequence, collapse = "")
    )
  }

  data.frame(
    id = record_ids,
    sequence = record_sequences,
    stringsAsFactors = FALSE
  )
}


.sequence_length_without_ns <- function(sequence) {
  nchar(gsub("[Nn[:space:]]", "", sequence))
}


.count_nonempty_lines <- function(path) {
  if (!file.exists(path)) {
    return(0L)
  }

  lines <- trimws(readLines(path, warn = FALSE))
  sum(nzchar(lines))
}


.discover_supercontig_genes <- function(sample_names, wd) {
  genes <- character(0)

  for (sample_name in sample_names) {
    sample_dir <- file.path(wd, sample_name)
    if (!dir.exists(sample_dir)) {
      next
    }

    sample_files <- list.files(
      sample_dir,
      pattern = "_supercontig\\.fasta$",
      recursive = TRUE,
      full.names = TRUE
    )

    if (length(sample_files) == 0) {
      next
    }

    genes <- c(
      genes,
      sub("_supercontig\\.fasta$", "", basename(sample_files))
    )
  }

  sort(unique(genes))
}


.parse_target_gene_name <- function(record_id, known_genes = character(0)) {
  if (length(known_genes) > 0) {
    separators <- c("_", "-")
    matches <- known_genes[
      vapply(
        known_genes,
        function(gene_name) {
          if (identical(record_id, gene_name)) {
            return(TRUE)
          }

          any(vapply(
            separators,
            function(separator) {
              startsWith(record_id, paste0(gene_name, separator)) ||
                endsWith(record_id, paste0(separator, gene_name))
            },
            logical(1)
          ))
        },
        logical(1)
      )
    ]

    if (length(matches) > 0) {
      return(matches[which.max(nchar(matches))])
    }
  }

  if (grepl("-", record_id, fixed = TRUE)) {
    return(sub(".*-", "", record_id))
  }

  if (grepl("_", record_id, fixed = TRUE)) {
    return(sub("_[^_]+$", "", record_id))
  }

  record_id
}


.read_target_mean_lengths <- function(
  targets_file,
  dna = TRUE,
  known_genes = character(0)
) {
  records <- .read_fasta_records(targets_file)

  reference_lengths <- list()

  for (index in seq_len(nrow(records))) {
    gene_name <- .parse_target_gene_name(
      records$id[[index]],
      known_genes = known_genes
    )
    record_length <- .sequence_length_without_ns(records$sequence[[index]])

    if (!isTRUE(dna)) {
      record_length <- record_length * 3L
    }

    reference_lengths[[gene_name]] <- c(
      reference_lengths[[gene_name]],
      record_length
    )
  }

  if (length(reference_lengths) == 0) {
    return(setNames(integer(0), character(0)))
  }

  mean_lengths <- vapply(
    reference_lengths,
    function(lengths) as.integer(round(mean(lengths))),
    integer(1)
  )
  mean_lengths[order(names(mean_lengths))]
}


.collect_supercontig_sample_stats <- function(sample_name, wd, gene_names) {
  sample_dir <- file.path(wd, sample_name)
  seq_lengths <- stats::setNames(integer(length(gene_names)), gene_names)

  if (!dir.exists(sample_dir)) {
    warning("Sample directory not found: ", sample_dir, call. = FALSE)
    return(list(
      seq_lengths = seq_lengths,
      total_bases = 0L,
      genes_mapped = 0L,
      genes_with_contigs = 0L,
      genes_with_seqs = 0L,
      paralog_warnings_long = 0L,
      paralog_warnings_depth = 0L,
      genes_without_stitched_contigs = 0L,
      genes_with_stitched_contigs = 0L,
      genes_with_stitched_contigs_skipped = 0L,
      genes_with_chimera_warning = 0L
    ))
  }

  supercontig_files <- list.files(
    sample_dir,
    pattern = "_supercontig\\.fasta$",
    recursive = TRUE,
    full.names = TRUE
  )

  for (supercontig_file in supercontig_files) {
    gene_name <- sub("_supercontig\\.fasta$", "", basename(supercontig_file))
    if (!gene_name %in% names(seq_lengths)) {
      seq_lengths[gene_name] <- 0L
    }

    records <- .read_fasta_records(supercontig_file)
    if (nrow(records) == 0) {
      next
    }

    seq_lengths[gene_name] <- sum(vapply(
      records$sequence,
      .sequence_length_without_ns,
      integer(1)
    ))
  }

  genes_with_seqs <- sum(seq_lengths > 0)

  genes_mapped <- .count_nonempty_lines(file.path(
    sample_dir,
    paste0(sample_name, "_genes_with_mapped_reads.txt")
  ))
  if (genes_mapped == 0L) {
    genes_mapped <- genes_with_seqs
  }

  genes_with_contigs <- .count_nonempty_lines(file.path(
    sample_dir,
    paste0(sample_name, "_genes_with_contigs.txt")
  ))
  if (genes_with_contigs == 0L) {
    genes_with_contigs <- genes_with_seqs
  }

  paralog_warnings_long <- .count_nonempty_lines(file.path(
    sample_dir,
    paste0(sample_name, "_genes_with_long_paralog_warnings.txt")
  ))
  paralog_warnings_depth <- .count_nonempty_lines(file.path(
    sample_dir,
    paste0(sample_name, "_genes_with_paralog_warnings_by_contig_depth.csv")
  ))

  stitched_contig_file <- file.path(
    sample_dir,
    paste0(sample_name, "_genes_with_stitched_contig.csv")
  )
  genes_without_stitched_contigs <- 0L
  genes_with_stitched_contigs <- 0L
  genes_with_stitched_contigs_skipped <- 0L

  if (file.exists(stitched_contig_file)) {
    stitched_stats <- readLines(stitched_contig_file, warn = FALSE)
    for (line in stitched_stats) {
      parts <- strsplit(line, ",", fixed = TRUE)[[1]]
      if (length(parts) < 3) {
        next
      }

      stat <- trimws(parts[[3]])
      if (grepl("single Exonerate hit", stat, fixed = TRUE)) {
        genes_without_stitched_contigs <- genes_without_stitched_contigs + 1L
      } else if (grepl("Stitched contig produced", stat, fixed = TRUE)) {
        genes_with_stitched_contigs <- genes_with_stitched_contigs + 1L
      } else if (grepl("Stitched contig step skipped", stat, fixed = TRUE)) {
        genes_with_stitched_contigs_skipped <- genes_with_stitched_contigs_skipped +
          1L
      }
    }
  } else {
    genes_with_stitched_contigs <- genes_with_seqs
  }

  chimera_file <- file.path(
    sample_dir,
    paste0(
      sample_name,
      "_genes_derived_from_putative_chimeric_stitched_contig.csv"
    )
  )
  genes_with_chimera_warning <- 0L
  if (file.exists(chimera_file)) {
    chimera_stats <- readLines(chimera_file, warn = FALSE)
    for (line in chimera_stats) {
      parts <- strsplit(line, ",", fixed = TRUE)[[1]]
      if (length(parts) < 3) {
        next
      }

      if (
        grepl("Chimera WARNING for stitched_contig.", parts[[3]], fixed = TRUE)
      ) {
        genes_with_chimera_warning <- genes_with_chimera_warning + 1L
      }
    }
  }

  list(
    seq_lengths = seq_lengths,
    total_bases = sum(seq_lengths),
    genes_mapped = genes_mapped,
    genes_with_contigs = genes_with_contigs,
    genes_with_seqs = genes_with_seqs,
    paralog_warnings_long = paralog_warnings_long,
    paralog_warnings_depth = paralog_warnings_depth,
    genes_without_stitched_contigs = genes_without_stitched_contigs,
    genes_with_stitched_contigs = genes_with_stitched_contigs,
    genes_with_stitched_contigs_skipped = genes_with_stitched_contigs_skipped,
    genes_with_chimera_warning = genes_with_chimera_warning
  )
}


.compute_recovery_thresholds <- function(seq_lengths, mean_lengths) {
  valid_genes <- mean_lengths > 0

  c(
    GenesAt25pct = sum(valid_genes & seq_lengths > mean_lengths * 0.25),
    GenesAt50pct = sum(valid_genes & seq_lengths > mean_lengths * 0.50),
    GenesAt75pct = sum(valid_genes & seq_lengths > mean_lengths * 0.75),
    GenesAt150pct = sum(valid_genes & seq_lengths > mean_lengths * 1.50)
  )
}


.run_supercontig_stats_native <- function(
  targets_file,
  namelist,
  wd,
  dna,
  seq_lengths_filename,
  stats_filename
) {
  sample_names <- trimws(readLines(namelist, warn = FALSE))
  sample_names <- sample_names[nzchar(sample_names)]

  discovered_genes <- .discover_supercontig_genes(sample_names, wd)
  mean_lengths <- .read_target_mean_lengths(
    targets_file,
    dna = dna,
    known_genes = discovered_genes
  )

  gene_names <- sort(unique(c(names(mean_lengths), discovered_genes)))
  sample_stats <- stats::setNames(
    vector("list", length(sample_names)),
    sample_names
  )

  for (sample_name in sample_names) {
    sample_stats[[sample_name]] <- .collect_supercontig_sample_stats(
      sample_name = sample_name,
      wd = wd,
      gene_names = gene_names
    )
  }

  if (length(gene_names) > 0) {
    for (gene_name in setdiff(gene_names, names(mean_lengths))) {
      observed_lengths <- vapply(
        sample_stats,
        function(sample_stat) sample_stat$seq_lengths[[gene_name]],
        integer(1)
      )
      observed_lengths <- observed_lengths[observed_lengths > 0]

      mean_lengths[[gene_name]] <- if (length(observed_lengths) > 0) {
        round(mean(observed_lengths))
      } else {
        0L
      }
    }

    mean_lengths <- mean_lengths[gene_names]
  }

  seq_lengths_path <- file.path(wd, paste0(seq_lengths_filename, ".tsv"))
  stats_path <- file.path(wd, paste0(stats_filename, ".tsv"))

  seq_length_lines <- c(
    paste(c("Species", gene_names), collapse = "\t"),
    paste(c("MeanLength", as.character(mean_lengths)), collapse = "\t")
  )

  stats_rows <- vector("list", length(sample_names))
  names(stats_rows) <- sample_names

  for (sample_name in sample_names) {
    sample_stat <- sample_stats[[sample_name]]
    sample_seq_lengths <- sample_stat$seq_lengths[gene_names]
    thresholds <- .compute_recovery_thresholds(sample_seq_lengths, mean_lengths)

    seq_length_lines <- c(
      seq_length_lines,
      paste(c(sample_name, as.character(sample_seq_lengths)), collapse = "\t")
    )

    stats_rows[[sample_name]] <- data.frame(
      Name = sample_name,
      NumReads = 0L,
      ReadsMapped = 0L,
      PctOnTarget = sprintf("%.1f", 0),
      GenesMapped = sample_stat$genes_mapped,
      GenesWithContigs = sample_stat$genes_with_contigs,
      GenesWithSeqs = sample_stat$genes_with_seqs,
      GenesAt25pct = unname(thresholds[["GenesAt25pct"]]),
      GenesAt50pct = unname(thresholds[["GenesAt50pct"]]),
      GenesAt75pct = unname(thresholds[["GenesAt75pct"]]),
      GenesAt150pct = unname(thresholds[["GenesAt150pct"]]),
      ParalogWarningsLong = sample_stat$paralog_warnings_long,
      ParalogWarningsDepth = sample_stat$paralog_warnings_depth,
      GenesWithoutStitchedContigs = sample_stat$genes_without_stitched_contigs,
      GenesWithStitchedContigs = sample_stat$genes_with_stitched_contigs,
      GenesWithStitchedContigsSkipped = sample_stat$genes_with_stitched_contigs_skipped,
      GenesWithChimeraWarning = sample_stat$genes_with_chimera_warning,
      TotalBasesRecovered = sample_stat$total_bases,
      stringsAsFactors = FALSE
    )
  }

  writeLines(seq_length_lines, con = seq_lengths_path)
  utils::write.table(
    do.call(rbind, stats_rows),
    file = stats_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  warning(
    "Used native R supercontig stats fallback after 'hybpiper stats' failed. ",
    "GenesAt*pct columns are approximate for supercontigs.",
    call. = FALSE
  )

  list(
    seq_lengths = seq_lengths_path,
    stats = stats_path
  )
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
#' @details
#' When `mode = "supercontig"`, HybPhaser retries with a native R fallback if
#' `hybpiper stats` fails. The fallback writes HybPiper-compatible
#' `seq_lengths.tsv` and `hybpiper_stats.tsv` files directly from the assembled
#' supercontig FASTA files.
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
      "'. ",
      .hybpiper_install_hint(conda_env)
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
    if (identical(mode, "supercontig")) {
      message("Falling back to native R supercontig stats")

      return(.run_supercontig_stats_native(
        targets_file = targets_file,
        namelist = namelist,
        wd = wd,
        dna = dna,
        seq_lengths_filename = seq_lengths_filename,
        stats_filename = stats_filename
      ))
    }

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
#' @param other_args Character vector; additional CLI args passed to
#'   `hybpiper assemble`
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
  mapper = c("bwa", "blast"),
  other_args = character(0)
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
      conda_env = conda_env,
      other_args = other_args
    )

    if (!keep_fastq) {
      unlink(c(f_path, r_path))
    }
  }

  message("HybPiper test dataset run complete")

  list(
    hybpiper_dir = normalizePath(output_dir, mustWork = TRUE),
    targets_file = normalizePath(targets_out, mustWork = TRUE),
    namelist = normalizePath(namelist_out, mustWork = TRUE),
    samples = samples
  )
}
