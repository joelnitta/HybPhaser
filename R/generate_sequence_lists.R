#' Generate sequence lists for loci and samples
#'
#' Creates FASTA sequence lists organized by locus and by sample, applying
#' quality filters from dataset optimization. Generates four sets of output
#' files (locus-level and sample-level, for both consensus and contig
#' sequences) that can be used for downstream phylogenetic analyses and
#' phasing.
#'
#' @param path_to_output_folder Path to the base HybPiper output folder
#'   containing the "01_data" subfolder with sample results.
#' @param fasta_file_with_targets Path to FASTA file with target sequences
#'   used in HybPiper run.
#' @param targets_file_format Character string, either "DNA" or "AA"
#'   (amino acid).
#' @param path_to_namelist Path to text file with sample names (one per line)
#'   to include in sequence lists.
#' @param intronerated_contig Logical. If TRUE, processes intronerated contigs
#'   instead of regular contigs. Default FALSE.
#' @param subset_name Optional name for dataset optimization subset. If
#'   provided, creates output in subfolder "03_sequence_lists_<subset_name>"
#'   and loads R objects from corresponding assessment subfolder. If empty
#'   string (default), uses "03_sequence_lists" and loads from
#'   "00_R_objects".
#'
#' @details
#' This function performs the following steps:
#'
#' 1. **Load filtering results**: Reads R objects saved by [assess_dataset()]:
#'    - `outsamples_missing`: Samples removed for missing data
#'    - `outloci_missing`: Loci removed for missing data
#'    - `outloci_para_all`: Loci removed as paralogs across all samples
#'    - `outloci_para_each`: Loci removed as paralogs for individual samples
#'    - `Table_SNPs_cleaned.Rds`: Cleaned SNP table
#'
#' 2. **Generate locus-level sequence lists**: For each target locus,
#'    concatenates sequences from all samples into a single FASTA file
#'    (one file per locus). Creates both consensus and contig versions.
#'
#' 3. **Generate sample-level sequence lists**: For each sample, concatenates
#'    sequences from all loci into a single FASTA file (one file per sample).
#'    Creates both consensus and contig versions.
#'
#' 4. **Apply quality filters**: Removes sequences based on filters from
#'    `assess_dataset()`:
#'    - Samples removed for missing data
#'    - Loci removed for missing data or as paralogs (all samples)
#'    - Loci removed as paralogs for specific samples
#'
#' The function handles both Linux and non-Linux systems, using bash commands
#' (`cat`, `sed`) on Linux for faster processing, and R file operations
#' otherwise.
#'
#' Generated sequences are in non-interleaved FASTA format.
#'
#' @section Output directories:
#' Four directories are created under `03_sequence_lists` (or
#' `03_sequence_lists_<subset_name>`):
#'
#' - `loci_consensus/`: One file per locus with consensus sequences from all
#'   samples
#' - `loci_contigs/`: One file per locus with contig sequences from all
#'   samples
#' - `samples_consensus/`: One file per sample with consensus sequences from
#'   all loci
#' - `samples_contigs/`: One file per sample with contig sequences from all
#'   loci
#'
#' @section File naming:
#' - Locus files: `<locus>_consensus.fasta` or `<locus>_intronerated_consensus.fasta`
#' - Sample files: `<sample>_consensus.fasta` or `<sample>_intronerated_consensus.fasta`
#' - Similar naming for contig files (with `_contig` or `_contigs`)
#'
#' @return Invisibly returns a list with components:
#'   - `output_dir`: Path to main output directory
#'   - `loci_consensus_dir`: Path to locus consensus sequence directory
#'   - `loci_contigs_dir`: Path to locus contig sequence directory
#'   - `samples_consensus_dir`: Path to sample consensus sequence directory
#'   - `samples_contigs_dir`: Path to sample contig sequence directory
#'   - `n_loci_files`: Number of locus files created
#'   - `n_sample_files`: Number of sample files created
#'   - `loci_removed`: Character vector of loci removed from sequence lists
#'   - `samples_removed`: Character vector of samples removed from sequence lists
#'
#' @examples
#' \dontrun{
#' # Generate sequence lists after running assess_dataset
#' result <- generate_sequence_lists(
#'   path_to_output_folder = "output",
#'   fasta_file_with_targets = "targets.fasta",
#'   targets_file_format = "DNA",
#'   path_to_namelist = "namelist.txt"
#' )
#'
#' # With subset name
#' result <- generate_sequence_lists(
#'   path_to_output_folder = "output",
#'   fasta_file_with_targets = "targets.fasta",
#'   targets_file_format = "DNA",
#'   path_to_namelist = "namelist.txt",
#'   subset_name = "filtered"
#' )
#'
#' # View results
#' cat("Created", result$n_loci_files, "locus files\n")
#' cat("Created", result$n_sample_files, "sample files\n")
#' }
#'
#' @seealso [assess_dataset()], [count_snps()]
#'
#' @export
generate_sequence_lists <- function(
  path_to_output_folder,
  fasta_file_with_targets,
  targets_file_format = "DNA",
  path_to_namelist,
  intronerated_contig = FALSE,
  subset_name = ""
) {
  # Validate inputs
  if (!dir.exists(path_to_output_folder)) {
    stop("path_to_output_folder does not exist: ", path_to_output_folder)
  }

  if (!file.exists(fasta_file_with_targets)) {
    stop("fasta_file_with_targets does not exist: ", fasta_file_with_targets)
  }

  if (!file.exists(path_to_namelist)) {
    stop("path_to_namelist does not exist: ", path_to_namelist)
  }

  if (!targets_file_format %in% c("DNA", "AA")) {
    stop("targets_file_format must be 'DNA' or 'AA'")
  }

  # Set up paths based on subset name
  if (subset_name != "") {
    folder_subset_add <- paste0("_", subset_name)
  } else {
    folder_subset_add <- ""
  }

  output_Robjects <- file.path(
    path_to_output_folder,
    "00_R_objects",
    subset_name
  )

  output_sequences <- file.path(
    path_to_output_folder,
    paste0("03_sequence_lists", folder_subset_add)
  )

  # Check that R objects directory exists
  if (!dir.exists(output_Robjects)) {
    stop(
      "R objects directory does not exist: ",
      output_Robjects,
      "\nPlease run assess_dataset() first."
    )
  }

  # Set naming conventions for intronerated contigs
  if (intronerated_contig) {
    intronerated_name <- "intronerated"
    intronerated_underscore <- "_"
  } else {
    intronerated_name <- ""
    intronerated_underscore <- ""
  }

  # Read targets and sample names
  targets <- seqinr::read.fasta(
    fasta_file_with_targets,
    as.string = TRUE,
    set.attributes = FALSE
  )
  targets_name <- unique(gsub(".*-", "", labels(targets)))
  samples <- readLines(path_to_namelist)

  # Load R objects from assess_dataset
  message("Loading filtering results from assess_dataset...")

  outsamples_missing <- readRDS(
    file.path(output_Robjects, "outsamples_missing.Rds")
  )
  outloci_missing <- readRDS(
    file.path(output_Robjects, "outloci_missing.Rds")
  )
  outloci_para_all <- readRDS(
    file.path(output_Robjects, "outloci_para_all.Rds")
  )
  outloci_para_each <- readRDS(
    file.path(output_Robjects, "outloci_para_each.Rds")
  )
  tab_snps_cl2b <- readRDS(
    file.path(output_Robjects, "Table_SNPs_cleaned.Rds")
  )

  # Prepare data matrices
  tab_snps <- as.matrix(tab_snps_cl2b)
  loci <- t(tab_snps)
  failed_loci <- which(colSums(is.na(loci)) == nrow(loci))
  failed_samples <- which(colSums(is.na(tab_snps)) == nrow(tab_snps))

  # Create output directories
  message("Creating output directories...")

  folder4seq_consensus_loci <- file.path(output_sequences, "loci_consensus")
  folder4seq_consensus_samples <- file.path(
    output_sequences,
    "samples_consensus"
  )
  folder4seq_contig_loci <- file.path(output_sequences, "loci_contigs")
  folder4seq_contig_samples <- file.path(output_sequences, "samples_contigs")

  # Delete existing directories to prevent errors
  unlink(
    c(
      folder4seq_consensus_loci,
      folder4seq_consensus_samples,
      folder4seq_contig_loci,
      folder4seq_contig_samples
    ),
    recursive = TRUE
  )

  # Create fresh directories
  dir.create(output_sequences, showWarnings = FALSE)
  dir.create(folder4seq_consensus_loci, showWarnings = FALSE)
  dir.create(folder4seq_consensus_samples, showWarnings = FALSE)
  dir.create(folder4seq_contig_loci, showWarnings = FALSE)
  dir.create(folder4seq_contig_samples, showWarnings = FALSE)

  # Generate locus-level sequence lists
  message("Generating locus-level sequence lists...")

  generate_locus_lists(
    path_to_output_folder = path_to_output_folder,
    targets_name = targets_name,
    folder4seq_consensus_loci = folder4seq_consensus_loci,
    folder4seq_contig_loci = folder4seq_contig_loci,
    intronerated_contig = intronerated_contig,
    intronerated_name = intronerated_name,
    intronerated_underscore = intronerated_underscore
  )

  # Convert interleaved to non-interleaved format for contigs
  message("Converting contig files to non-interleaved format...")
  convert_to_noninterleaved(folder4seq_contig_loci)

  # Remove samples not in namelist
  message("Filtering samples not in namelist...")
  remove_unlisted_samples(
    path_to_output_folder = path_to_output_folder,
    samples = samples,
    folder4seq_consensus_loci = folder4seq_consensus_loci,
    folder4seq_contig_loci = folder4seq_contig_loci
  )

  # Remove loci from dataset optimization (missing data and paralogs)
  message("Removing loci based on dataset optimization filters...")

  loci_to_remove <- c(names(failed_loci), outloci_missing, outloci_para_all)

  remove_loci_from_locus_files(
    loci_to_remove = loci_to_remove,
    folder4seq_consensus_loci = folder4seq_consensus_loci,
    folder4seq_contig_loci = folder4seq_contig_loci,
    intronerated_contig = intronerated_contig
  )

  # Remove samples and per-sample paralogs from locus files
  message("Removing samples and per-sample paralogs from locus files...")

  samples_to_remove_4all <- c()
  if (length(failed_samples) > 0) {
    samples_to_remove_4all <- names(failed_samples)
  }
  if (length(outsamples_missing) > 0) {
    samples_to_remove_4all <- unique(c(
      samples_to_remove_4all,
      outsamples_missing
    ))
  }

  remove_samples_from_locus_files(
    tab_snps_cl2b = tab_snps_cl2b,
    failed_loci = failed_loci,
    loci_to_remove = loci_to_remove,
    samples_to_remove_4all = samples_to_remove_4all,
    outloci_para_each = outloci_para_each,
    folder4seq_consensus_loci = folder4seq_consensus_loci,
    folder4seq_contig_loci = folder4seq_contig_loci,
    intronerated_name = intronerated_name,
    intronerated_underscore = intronerated_underscore
  )

  # Generate sample-level sequence lists
  message("Generating sample-level sequence lists...")

  samples_in <- samples
  if (length(failed_samples) != 0) {
    samples_in <- samples_in[!(samples %in% names(failed_samples))]
  }

  generate_sample_lists(
    path_to_output_folder = path_to_output_folder,
    samples = samples_in,
    folder4seq_consensus_samples = folder4seq_consensus_samples,
    folder4seq_contig_samples = folder4seq_contig_samples,
    intronerated_contig = intronerated_contig,
    intronerated_name = intronerated_name,
    intronerated_underscore = intronerated_underscore
  )

  # Convert interleaved to non-interleaved format for sample contigs
  message("Converting sample contig files to non-interleaved format...")
  convert_to_noninterleaved(folder4seq_contig_samples)

  # Remove samples based on missing data
  message("Removing samples based on dataset optimization filters...")

  remove_samples_from_sample_files(
    outsamples_missing = outsamples_missing,
    folder4seq_consensus_samples = folder4seq_consensus_samples,
    folder4seq_contig_samples = folder4seq_contig_samples
  )

  # Remove loci from sample files
  message("Removing loci from sample files...")

  samples_in <- samples
  if (length(outsamples_missing) != 0) {
    samples_in <- samples_in[!(samples %in% outsamples_missing)]
  }

  loci_to_remove_4all <- c()
  if (length(failed_loci) > 0) {
    loci_to_remove_4all <- names(failed_loci)
  }
  if (length(outloci_missing) > 0) {
    loci_to_remove_4all <- c(loci_to_remove_4all, outloci_missing)
  }
  if (length(outloci_para_all) > 0) {
    loci_to_remove_4all <- c(loci_to_remove_4all, outloci_para_all)
  }

  remove_loci_from_sample_files(
    samples_in = samples_in,
    loci_to_remove_4all = loci_to_remove_4all,
    outloci_para_each = outloci_para_each,
    folder4seq_consensus_samples = folder4seq_consensus_samples,
    folder4seq_contig_samples = folder4seq_contig_samples,
    intronerated_name = intronerated_name,
    intronerated_underscore = intronerated_underscore
  )

  # Count output files
  n_loci_consensus <- as.numeric(length(list.files(folder4seq_consensus_loci)))
  n_loci_contigs <- as.numeric(length(list.files(folder4seq_contig_loci)))
  n_samples_consensus <- as.numeric(
    length(list.files(folder4seq_consensus_samples))
  )
  n_samples_contigs <- as.numeric(length(list.files(folder4seq_contig_samples)))

  message("Sequence list generation complete!")
  message(
    "  Locus files created: ",
    n_loci_consensus,
    " consensus, ",
    n_loci_contigs,
    " contigs"
  )
  message(
    "  Sample files created: ",
    n_samples_consensus,
    " consensus, ",
    n_samples_contigs,
    " contigs"
  )

  # Prepare return value
  result <- list(
    output_dir = output_sequences,
    loci_consensus_dir = folder4seq_consensus_loci,
    loci_contigs_dir = folder4seq_contig_loci,
    samples_consensus_dir = folder4seq_consensus_samples,
    samples_contigs_dir = folder4seq_contig_samples,
    n_loci_files = n_loci_consensus,
    n_sample_files = n_samples_consensus,
    loci_removed = loci_to_remove,
    samples_removed = unique(c(names(failed_samples), outsamples_missing))
  )

  invisible(result)
}


# Helper function: Generate locus-level sequence lists
generate_locus_lists <- function(
  path_to_output_folder,
  targets_name,
  folder4seq_consensus_loci,
  folder4seq_contig_loci,
  intronerated_contig,
  intronerated_name,
  intronerated_underscore
) {
  if (Sys.info()["sysname"] == "Linux") {
    # Use bash commands on Linux for speed
    if (intronerated_contig) {
      for (locus in targets_name) {
        command_cat_consensus <- paste(
          "cat",
          file.path(
            path_to_output_folder,
            "01_data/*/intronerated_consensus/",
            paste0(locus, "_intronerated.fasta")
          ),
          ">",
          file.path(
            folder4seq_consensus_loci,
            paste0(locus, "_intronerated_consensus.fasta")
          )
        )
        system(command_cat_consensus, ignore.stderr = TRUE)

        command_cat_contig <- paste(
          "cat",
          file.path(
            path_to_output_folder,
            "01_data/*/intronerated_contigs/",
            paste0(locus, "_intronerated.fasta")
          ),
          ">",
          file.path(
            folder4seq_contig_loci,
            paste0(locus, "_intronerated_contig.fasta")
          )
        )
        system(command_cat_contig, ignore.stderr = TRUE)

        command_remove_locus_in_seqnames_consensus <- paste0(
          "sed -i 's/-",
          locus,
          "//g' ",
          file.path(
            folder4seq_consensus_loci,
            paste0(locus, "_intronerated_consensus.fasta")
          )
        )
        system(command_remove_locus_in_seqnames_consensus)

        command_remove_locus_in_seqnames_contig <- paste0(
          "sed -i 's/-",
          locus,
          "//g' ",
          file.path(
            folder4seq_contig_loci,
            paste0(locus, "_intronerated_contig.fasta")
          )
        )
        system(command_remove_locus_in_seqnames_contig)
      }
    } else {
      for (locus in targets_name) {
        command_cat_consensus <- paste(
          "cat",
          file.path(
            path_to_output_folder,
            "01_data/*/consensus/",
            paste0(locus, ".fasta")
          ),
          ">",
          file.path(
            folder4seq_consensus_loci,
            paste0(locus, "_consensus.fasta")
          )
        )
        command_cat_contig <- paste(
          "cat",
          file.path(
            path_to_output_folder,
            "01_data/*/contigs/",
            paste0(locus, ".fasta")
          ),
          ">",
          file.path(
            folder4seq_contig_loci,
            paste0(locus, "_contig.fasta")
          )
        )
        system(command_cat_consensus, ignore.stderr = TRUE)
        system(command_cat_contig, ignore.stderr = TRUE)

        command_remove_locus_in_seqnames_consensus <- paste0(
          "sed -i 's/-",
          locus,
          "//g' ",
          file.path(
            folder4seq_consensus_loci,
            paste0(locus, "_consensus.fasta")
          )
        )
        command_remove_locus_in_seqnames_contig <- paste0(
          "sed -i 's/-",
          locus,
          "//g' ",
          file.path(
            folder4seq_contig_loci,
            paste0(locus, "_contig.fasta")
          )
        )
        system(command_remove_locus_in_seqnames_consensus)
        system(command_remove_locus_in_seqnames_contig)
      }
    }
  } else {
    # Use R file operations on non-Linux systems
    if (intronerated_contig) {
      for (locus in targets_name) {
        # List all fasta files from that locus for all samples
        fasta_files <- list.files(
          path = file.path(path_to_output_folder, "01_data/"),
          pattern = paste0(locus, "_intronerated.fasta"),
          recursive = TRUE,
          full.names = TRUE
        )

        # Select consensus/contig files
        consensus_files <- grep("consensus", fasta_files, value = TRUE)
        contigs_files <- grep("contigs", fasta_files, value = TRUE)

        # Define output files
        output_file_consensus <- file.path(
          folder4seq_consensus_loci,
          paste0(locus, "_intronerated_consensus.fasta")
        )
        output_file_contigs <- file.path(
          folder4seq_contig_loci,
          paste0(locus, "_intronerated_contig.fasta")
        )

        # Generate output files
        file.create(output_file_consensus)
        file.create(output_file_contigs)

        # Append sample fastas to the empty output file
        file.append(output_file_consensus, consensus_files)
        file.append(output_file_contigs, contigs_files)

        # Read lines of each file and remove the "-locus" from sequence names
        lines_consensus <- gsub(
          paste0("-", locus),
          "",
          readLines(output_file_consensus)
        )
        lines_contigs <- gsub(
          paste0("-", locus),
          "",
          readLines(output_file_contigs)
        )

        # Write lines into files
        write(lines_consensus, file = output_file_consensus)
        write(lines_contigs, file = output_file_contigs)
      }
    } else {
      for (locus in targets_name) {
        # List all fasta files from that locus for all samples
        fasta_files <- list.files(
          path = file.path(path_to_output_folder, "01_data/"),
          pattern = paste0(locus, ".fasta"),
          recursive = TRUE,
          full.names = TRUE
        )

        # Select consensus/contig files
        consensus_files <- grep("consensus", fasta_files, value = TRUE)
        contigs_files <- grep("contigs", fasta_files, value = TRUE)

        # Define output files
        output_file_consensus <- file.path(
          folder4seq_consensus_loci,
          paste0(locus, "_consensus.fasta")
        )
        output_file_contigs <- file.path(
          folder4seq_contig_loci,
          paste0(locus, "_contig.fasta")
        )

        # Generate output files
        file.create(output_file_consensus)
        file.create(output_file_contigs)

        # Append sample fastas to the empty output file
        file.append(output_file_consensus, consensus_files)
        file.append(output_file_contigs, contigs_files)

        # Read lines of each file and remove the "-locus" from sequence names
        lines_consensus <- gsub(
          paste0("-", locus),
          "",
          readLines(output_file_consensus)
        )
        lines_contigs <- gsub(
          paste0("-", locus),
          "",
          readLines(output_file_contigs)
        )

        # Write lines into files
        write(lines_consensus, file = output_file_consensus)
        write(lines_contigs, file = output_file_contigs)
      }
    }
  }
}


# Helper function: Convert interleaved FASTA to non-interleaved
convert_to_noninterleaved <- function(folder) {
  for (file_path in list.files(folder, full.names = TRUE)) {
    lines <- readLines(file_path)
    conx <- file(file_path)
    writeLines(
      stringr::str_split(
        paste(gsub("(>.*)", ":\\1:", lines), collapse = ""),
        pattern = ":"
      )[[1]][-1],
      conx
    )
    close(conx)
  }
}


# Helper function: Remove samples not in namelist
remove_unlisted_samples <- function(
  path_to_output_folder,
  samples,
  folder4seq_consensus_loci,
  folder4seq_contig_loci
) {
  hybpiper_result_dirs <- list.dirs(
    file.path(path_to_output_folder, "01_data"),
    full.names = FALSE,
    recursive = FALSE
  )
  dirs_not_in_sample_list <- hybpiper_result_dirs[
    !(hybpiper_result_dirs %in% samples)
  ]

  if (length(dirs_not_in_sample_list) != 0) {
    # Consensus
    for (raw_consensus_file in list.files(
      folder4seq_consensus_loci,
      full.names = TRUE
    )) {
      locus_consensus <- readLines(raw_consensus_file)
      lines_with_samplename <- which(
        gsub(">", "", locus_consensus) %in% dirs_not_in_sample_list
      )
      if (length(lines_with_samplename) != 0) {
        lines_to_remove <- c(lines_with_samplename, lines_with_samplename + 1)
        locus_file_red <- locus_consensus[-lines_to_remove]
        conn <- file(raw_consensus_file)
        writeLines(locus_file_red, conn)
        close(conn)
      }
    }

    # Contig
    for (raw_contig_file in list.files(
      folder4seq_contig_loci,
      full.names = TRUE
    )) {
      locus_contig <- readLines(raw_contig_file)
      lines_with_samplename <- which(
        gsub(">", "", locus_contig) %in% dirs_not_in_sample_list
      )
      if (length(lines_with_samplename) != 0) {
        lines_to_remove <- c(lines_with_samplename, lines_with_samplename + 1)
        locus_file_hp_red <- locus_contig[-lines_to_remove]
        conn <- file(raw_contig_file)
        writeLines(locus_file_hp_red, conn)
        close(conn)
      }
    }
  }
}


# Helper function: Remove loci from locus files
remove_loci_from_locus_files <- function(
  loci_to_remove,
  folder4seq_consensus_loci,
  folder4seq_contig_loci,
  intronerated_contig
) {
  if (length(loci_to_remove) == 0) {
    return(invisible(NULL))
  }

  loci_files_consensus <- list.files(
    path = folder4seq_consensus_loci,
    full.names = TRUE
  )
  loci_files_contig <- list.files(
    path = folder4seq_contig_loci,
    full.names = TRUE
  )

  if (!intronerated_contig) {
    loci_files_to_remove_consensus <- loci_files_consensus[
      which(
        gsub(".*/(.*)_consensus.fasta", "\\1", loci_files_consensus) %in%
          loci_to_remove
      )
    ]
    loci_files_to_remove_contig <- loci_files_contig[
      which(
        gsub(".*/(.*)_contig.fasta", "\\1", loci_files_contig) %in%
          loci_to_remove
      )
    ]
  } else {
    loci_files_to_remove_consensus <- loci_files_consensus[
      which(
        gsub(
          ".*/(.*)_intronerated_consensus.fasta",
          "\\1",
          loci_files_consensus
        ) %in%
          loci_to_remove
      )
    ]
    loci_files_to_remove_contig <- loci_files_contig[
      which(
        gsub(".*/(.*)_intronerated_contig.fasta", "\\1", loci_files_contig) %in%
          loci_to_remove
      )
    ]
  }

  file.remove(loci_files_to_remove_consensus)
  file.remove(loci_files_to_remove_contig)
}


# Helper function: Remove samples from locus files
remove_samples_from_locus_files <- function(
  tab_snps_cl2b,
  failed_loci,
  loci_to_remove,
  samples_to_remove_4all,
  outloci_para_each,
  folder4seq_consensus_loci,
  folder4seq_contig_loci,
  intronerated_name,
  intronerated_underscore
) {
  for (locus in rownames(tab_snps_cl2b)) {
    # Skip loci that were already removed (failed, missing data, or paralogs)
    if (locus %in% c(names(failed_loci), loci_to_remove)) {
      next
    }

    # ADDITIONAL CHECK: Skip if the locus files don't exist
    # (defensive check in case they were removed but still in rownames)
    locus_consensus_file <- file.path(
      folder4seq_consensus_loci,
      paste0(
        locus,
        "_",
        intronerated_name,
        intronerated_underscore,
        "consensus.fasta"
      )
    )
    locus_contig_file <- file.path(
      folder4seq_contig_loci,
      paste0(
        locus,
        "_",
        intronerated_name,
        intronerated_underscore,
        "contig.fasta"
      )
    )

    if (!file.exists(locus_consensus_file) && !file.exists(locus_contig_file)) {
      next
    }

    if (length(grep(paste0("\\b", locus, "\\b"), outloci_para_each)) > 0) {
      samples_to_remove <- c(
        samples_to_remove_4all,
        names(outloci_para_each[grep(
          paste0("\\b", locus, "\\b"),
          outloci_para_each
        )])
      )
    } else {
      samples_to_remove <- samples_to_remove_4all
    }

    if (length(samples_to_remove) != 0) {
      # Consensus
      if (file.exists(locus_consensus_file)) {
        locus_consensus <- readLines(locus_consensus_file)
        lines_with_samplename <- which(
          gsub(">", "", locus_consensus) %in% samples_to_remove
        )
        if (length(lines_with_samplename) != 0) {
          lines_to_remove <- c(lines_with_samplename, lines_with_samplename + 1)
          locus_file_red <- locus_consensus[-lines_to_remove]
          conn <- file(locus_consensus_file)
          writeLines(locus_file_red, conn)
          close(conn)
        }
      }

      # Contig
      if (file.exists(locus_contig_file)) {
        locus_contig <- readLines(locus_contig_file)
        lines_with_samplename <- which(
          gsub(">", "", locus_contig) %in% samples_to_remove
        )
        if (length(lines_with_samplename) != 0) {
          lines_to_remove <- c(lines_with_samplename, lines_with_samplename + 1)
          locus_file_hp_red <- locus_contig[-lines_to_remove]
          conn <- file(locus_contig_file)
          writeLines(locus_file_hp_red, conn)
          close(conn)
        }
      }
    }
  }
}


# Helper function: Generate sample-level sequence lists
generate_sample_lists <- function(
  path_to_output_folder,
  samples,
  folder4seq_consensus_samples,
  folder4seq_contig_samples,
  intronerated_contig,
  intronerated_name,
  intronerated_underscore
) {
  if (Sys.info()["sysname"] == "Linux") {
    # Use bash commands on Linux for speed
    if (intronerated_contig) {
      for (sample in samples) {
        command_cat_loci_consensus <- paste(
          "cat",
          file.path(
            path_to_output_folder,
            "01_data/",
            sample,
            "/intronerated_consensus/*.fasta"
          ),
          ">",
          file.path(
            folder4seq_consensus_samples,
            paste0(sample, "_intronerated_consensus.fasta")
          )
        )
        system(command_cat_loci_consensus)

        command_cat_loci_contig <- paste(
          "cat",
          file.path(
            path_to_output_folder,
            "01_data/",
            sample,
            "/intronerated_contigs/*.fasta"
          ),
          ">",
          file.path(
            folder4seq_contig_samples,
            paste0(sample, "_intronerated_contig.fasta")
          )
        )
        system(command_cat_loci_contig)
      }
    } else {
      for (sample in samples) {
        command_cat_loci_consensus <- paste(
          "cat",
          file.path(
            path_to_output_folder,
            "01_data/",
            sample,
            "/consensus/*.fasta"
          ),
          ">",
          file.path(
            folder4seq_consensus_samples,
            paste0(sample, "_consensus.fasta")
          )
        )
        system(command_cat_loci_consensus)

        command_cat_loci_contig <- paste(
          "cat",
          file.path(
            path_to_output_folder,
            "01_data/",
            sample,
            "/contigs/*.fasta"
          ),
          ">",
          file.path(
            folder4seq_contig_samples,
            paste0(sample, "_contig.fasta")
          )
        )
        system(command_cat_loci_contig)
      }
    }
  } else {
    # Use R file operations on non-Linux systems
    if (intronerated_contig) {
      for (sample in samples) {
        # List all fasta files from the sample for all loci
        consensus_files_samples <- list.files(
          path = file.path(
            path_to_output_folder,
            "01_data/",
            sample,
            "/intronerated_consensus/"
          ),
          pattern = "*.fasta",
          full.names = TRUE
        )
        contigs_files_samples <- list.files(
          path = file.path(
            path_to_output_folder,
            "01_data/",
            sample,
            "/intronerated_contigs/"
          ),
          pattern = "*.fasta",
          full.names = TRUE
        )

        # Define output files
        output_file_consensus_samples <- file.path(
          folder4seq_consensus_samples,
          paste0(sample, "_intronerated_consensus.fasta")
        )
        output_file_contigs_samples <- file.path(
          folder4seq_contig_samples,
          paste0(sample, "_intronerated_contig.fasta")
        )

        # Create output files
        file.create(output_file_consensus_samples)
        file.create(output_file_contigs_samples)

        # Append fasta files to empty output file
        file.append(output_file_consensus_samples, consensus_files_samples)
        file.append(output_file_contigs_samples, contigs_files_samples)
      }
    } else {
      for (sample in samples) {
        # List all fasta files from the sample for all loci
        consensus_files_samples <- list.files(
          path = file.path(
            path_to_output_folder,
            "01_data/",
            sample,
            "/consensus/"
          ),
          pattern = "*.fasta",
          full.names = TRUE
        )
        contigs_files_samples <- list.files(
          path = file.path(
            path_to_output_folder,
            "01_data/",
            sample,
            "/contigs/"
          ),
          pattern = "*.fasta",
          full.names = TRUE
        )

        # Define output files
        output_file_consensus_samples <- file.path(
          folder4seq_consensus_samples,
          paste0(sample, "_consensus.fasta")
        )
        output_file_contigs_samples <- file.path(
          folder4seq_contig_samples,
          paste0(sample, "_contig.fasta")
        )

        # Create output files
        file.create(output_file_consensus_samples)
        file.create(output_file_contigs_samples)

        # Append fasta files to empty output file
        file.append(output_file_consensus_samples, consensus_files_samples)
        file.append(output_file_contigs_samples, contigs_files_samples)
      }
    }
  }
}


# Helper function: Remove samples from sample files
remove_samples_from_sample_files <- function(
  outsamples_missing,
  folder4seq_consensus_samples,
  folder4seq_contig_samples
) {
  if (length(outsamples_missing) == 0) {
    return(invisible(NULL))
  }

  sample_files_consensus <- list.files(
    path = folder4seq_consensus_samples,
    full.names = TRUE
  )
  sample_files_contig <- list.files(
    path = folder4seq_contig_samples,
    full.names = TRUE
  )

  samples_files_to_remove_consensus <- sample_files_consensus[
    which(
      gsub(".*/(.*)_consensus.fasta", "\\1", sample_files_consensus) %in%
        outsamples_missing
    )
  ]
  samples_files_to_remove_contig <- sample_files_contig[
    which(
      gsub(".*/(.*)_contig.fasta", "\\1", sample_files_contig) %in%
        outsamples_missing
    )
  ]

  file.remove(samples_files_to_remove_consensus)
  file.remove(samples_files_to_remove_contig)
}


# Helper function: Remove loci from sample files
remove_loci_from_sample_files <- function(
  samples_in,
  loci_to_remove_4all,
  outloci_para_each,
  folder4seq_consensus_samples,
  folder4seq_contig_samples,
  intronerated_name,
  intronerated_underscore
) {
  for (sample in samples_in) {
    if (length(which(names(outloci_para_each) %in% sample)) > 0) {
      loci_to_remove <- c(
        loci_to_remove_4all,
        names(outloci_para_each[[which(names(outloci_para_each) %in% sample)]])
      )
    } else {
      loci_to_remove <- loci_to_remove_4all
    }

    if (length(loci_to_remove) != 0) {
      # Consensus
      consensus_file2clean <- file.path(
        folder4seq_consensus_samples,
        paste0(
          sample,
          intronerated_underscore,
          intronerated_name,
          "_consensus.fasta"
        )
      )

      # Contig
      contig_file2clean <- file.path(
        folder4seq_contig_samples,
        paste0(
          sample,
          intronerated_underscore,
          intronerated_name,
          "_contig.fasta"
        )
      )

      # Skip only when neither file exists
      if (
        !file.exists(consensus_file2clean) && !file.exists(contig_file2clean)
      ) {
        next
      }

      # Process consensus file
      if (file.exists(consensus_file2clean)) {
        samples_consensus <- readLines(consensus_file2clean)
        lines_with_lociname <- which(
          gsub(">.*-", "", samples_consensus) %in% loci_to_remove
        )
        if (length(lines_with_lociname) != 0) {
          lines_to_remove <- c(lines_with_lociname, lines_with_lociname + 1)
          sample_file_consensus_red <- samples_consensus[-lines_to_remove]
          conn <- file(consensus_file2clean)
          writeLines(sample_file_consensus_red, conn)
          close(conn)
        }
      }

      # Process contig file
      if (file.exists(contig_file2clean)) {
        samples_contig <- readLines(contig_file2clean)
        lines_with_lociname <- which(
          gsub(">.*-", "", samples_contig) %in% loci_to_remove
        )
        if (length(lines_with_lociname) != 0) {
          lines_to_remove <- c(lines_with_lociname, lines_with_lociname + 1)
          sample_file_contig_red <- samples_contig[-lines_to_remove]
          conn <- file(contig_file2clean)
          writeLines(sample_file_contig_red, conn)
          close(conn)
        }
      }
    }
  }
}
