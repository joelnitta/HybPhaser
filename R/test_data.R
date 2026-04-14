#' Create Test HybPiper Output Using Real Test Data
#'
#' Creates a minimal but realistic HybPiper output directory structure using
#' test data downloaded from the upstream HybPiper project. This creates a
#' structure suitable for testing HybPhaser functions with realistic gene names
#' and sequences.
#'
#' @param output_dir Path where test HybPiper output should be created
#' @param samples Character vector of sample names. If NULL (default), uses
#'   samples from HybPiper's test dataset namelist.
#' @param n_genes Integer; number of genes to include (default: 5). Will use
#'   the first n genes from the HybPiper test targets file.
#'
#' @return Path to created test directory (invisibly)
#' @export
#'
#' @examples
#' \dontrun{
#' # Create test data from real HybPiper test files
#' test_dir <- create_hybpiper_test_output(tempdir())
#'
#' # Use with HybPhaser functions
#' run_generate_consensus_sequences(
#'   hybpiper_dir = test_dir,
#'   output_dir = file.path(tempdir(), "hybphaser_out")
#' )
#' }
create_hybpiper_test_output <- function(
  output_dir,
  samples = NULL,
  n_genes = 5
) {
  download_dir <- tempfile("hybpiper_test_dataset_")
  on.exit(unlink(download_dir, recursive = TRUE, force = TRUE), add = TRUE)

  downloaded <- .download_hybpiper_test_dataset(download_dir)
  targets_file <- downloaded$targets_file
  namelist_file <- downloaded$namelist_file

  if (!file.exists(targets_file)) {
    stop("Test targets file not found: ", targets_file)
  }

  # Read samples from namelist if not provided
  if (is.null(samples)) {
    if (file.exists(namelist_file)) {
      samples <- readLines(namelist_file)
    } else {
      samples <- c("EG30", "EG98", "MWL2") # Default samples
    }
  }

  # Read target sequences
  targets <- seqinr::read.fasta(
    targets_file,
    seqtype = "DNA",
    as.string = TRUE,
    set.attributes = FALSE
  )

  # Extract gene names (everything after the dash in header like "Artocarpus-gene001")
  gene_names <- gsub(".*-", "", names(targets))

  # Limit to n_genes
  if (length(gene_names) > n_genes) {
    gene_names <- gene_names[seq_len(n_genes)]
    targets <- targets[seq_len(n_genes)]
  }

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Create structure for each sample
  for (sample in samples) {
    sample_dir <- file.path(output_dir, sample)

    for (i in seq_along(gene_names)) {
      gene <- gene_names[i]
      gene_dir <- file.path(sample_dir, gene)
      dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

      # Create contig file using real target sequence
      contig_file <- file.path(gene_dir, paste0(gene, "_contigs.fasta"))
      target_seq <- as.character(targets[[i]])

      # Add some variation to simulate real contigs (occasional N)
      seq_chars <- strsplit(target_seq, "")[[1]]
      # Randomly replace ~1% of bases with N
      n_snps <- max(1, round(length(seq_chars) * 0.01))
      snp_positions <- sample(seq_along(seq_chars), n_snps)
      seq_chars[snp_positions] <- "N"
      varied_seq <- paste(seq_chars, collapse = "")

      contig_content <- paste0(">", gene, "_contig\n", varied_seq)
      writeLines(contig_content, contig_file)

      # Create interleaved reads file
      reads_file <- file.path(gene_dir, paste0(gene, "_interleaved.fasta"))

      # Generate reads from the target sequence
      n_reads <- 50
      read_length <- 100
      target_length <- nchar(target_seq)
      read_seqs <- character(n_reads * 2) # Header + sequence

      for (j in seq_len(n_reads)) {
        read_seqs[(j - 1) * 2 + 1] <- paste0(">read_", j, "/1")

        # Extract random substring from target
        if (target_length > read_length) {
          start_pos <- sample(1:(target_length - read_length + 1), 1)
          read_seq <- substr(target_seq, start_pos, start_pos + read_length - 1)
        } else {
          read_seq <- target_seq
        }

        # Add occasional SNPs (convert some bases to N)
        if (j %% 10 == 0) {
          read_chars <- strsplit(read_seq, "")[[1]]
          if (length(read_chars) > 0) {
            snp_pos <- sample(seq_along(read_chars), 1)
            read_chars[snp_pos] <- "N"
            read_seq <- paste(read_chars, collapse = "")
          }
        }

        read_seqs[(j - 1) * 2 + 2] <- read_seq
      }

      writeLines(read_seqs, reads_file)
    }
  }

  message("Created HybPiper test output at: ", output_dir)
  message("  Samples: ", paste(samples, collapse = ", "))
  message("  Genes: ", paste(gene_names, collapse = ", "))
  message("  Using real HybPiper test target sequences")

  invisible(output_dir)
}


#' Create Test HybPiper Output
#'
#' Creates a minimal but realistic HybPiper output directory structure for
#' testing and examples. Includes mock sequence data.
#'
#' @param base_dir Path where test data should be created
#' @param samples Character vector of sample names. Default is
#'   c("sample1", "sample2")
#' @param genes Character vector of gene names. Default is
#'   c("gene001", "gene002")
#'
#' @return Path to created test directory (invisibly)
#' @export
#'
#' @examples
#' \dontrun{
#' # Create test data
#' test_dir <- create_test_hybpiper_output(tempdir())
#'
#' # Use with HybPhaser functions
#' run_generate_consensus_sequences(
#'   hybpiper_dir = test_dir,
#'   output_dir = file.path(tempdir(), "hybphaser_out"),
#'   sample = "sample1"
#' )
#' }
create_test_hybpiper_output <- function(
  base_dir,
  samples = c("sample1", "sample2"),
  genes = c("gene001", "gene002")
) {
  # Create base directory
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE)
  }

  # Create structure for each sample
  for (sample in samples) {
    sample_dir <- file.path(base_dir, sample)

    for (gene in genes) {
      gene_dir <- file.path(sample_dir, gene)
      dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

      # Create a contig file
      contig_file <- file.path(gene_dir, paste0(gene, "_contigs.fasta"))
      contig_seq <- paste0(
        ">",
        gene,
        "_contig\n",
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCG"
      )
      writeLines(contig_seq, contig_file)

      # Create interleaved reads file (simulating mapped reads)
      reads_file <- file.path(gene_dir, paste0(gene, "_interleaved.fasta"))

      # Generate reads with some variation
      n_reads <- 20
      read_seqs <- character(n_reads * 2) # Header + sequence for each

      for (i in seq_len(n_reads)) {
        read_seqs[(i - 1) * 2 + 1] <- paste0(">read_", i, "/1")
        # Create reads with occasional SNPs (N for simplicity)
        base_seq <- "ATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        if (i %% 5 == 0) {
          # Add a variant every 5th read
          base_seq <- gsub("A", "N", base_seq, fixed = TRUE)[1]
        }
        read_seqs[(i - 1) * 2 + 2] <- base_seq
      }

      writeLines(read_seqs, reads_file)
    }
  }

  message("Created test HybPiper output at: ", base_dir)
  message("  Samples: ", paste(samples, collapse = ", "))
  message("  Genes: ", paste(genes, collapse = ", "))

  invisible(base_dir)
}


#' Create Test Target File
#'
#' Creates a FASTA file with target sequences (baits) for testing
#'
#' @param file_path Path where target file should be created
#' @param genes Character vector of gene names
#' @param seq_length Integer; length of sequences to generate
#'
#' @return Path to created file (invisibly)
#' @export
#'
#' @examples
#' \dontrun{
#' target_file <- create_test_targets(
#'   file.path(tempdir(), "targets.fasta")
#' )
#' }
create_test_targets <- function(
  file_path,
  genes = c("gene001", "gene002"),
  seq_length = 150
) {
  # Ensure directory exists
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)

  # Generate target sequences
  target_seqs <- character(length(genes) * 2)

  for (i in seq_along(genes)) {
    target_seqs[(i - 1) * 2 + 1] <- paste0(">species1-", genes[i])

    # Generate a sequence
    bases <- c("A", "T", "C", "G")
    seq <- paste(sample(bases, seq_length, replace = TRUE), collapse = "")
    target_seqs[(i - 1) * 2 + 2] <- seq
  }

  writeLines(target_seqs, file_path)
  message("Created target file: ", file_path)

  invisible(file_path)
}


#' Create Test Namelist File
#'
#' Creates a namelist file for testing
#'
#' @param file_path Path where namelist should be created
#' @param samples Character vector of sample names
#'
#' @return Path to created file (invisibly)
#' @export
#'
#' @examples
#' \dontrun{
#' namelist <- create_test_namelist(
#'   file.path(tempdir(), "namelist.txt"),
#'   samples = c("sample1", "sample2")
#' )
#' }
create_test_namelist <- function(file_path, samples = c("sample1", "sample2")) {
  # Ensure directory exists
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)

  writeLines(samples, file_path)
  message("Created namelist: ", file_path)

  invisible(file_path)
}


#' Create Complete Test Dataset
#'
#' Creates a complete test dataset including HybPiper output, targets,
#' and namelist
#'
#' @param base_dir Base directory for test data
#' @param samples Character vector of sample names
#' @param genes Character vector of gene names
#'
#' @return Named list with paths to created files
#' @export
#'
#' @examples
#' \dontrun{
#' test_data <- create_test_dataset(tempdir())
#'
#' # Use with HybPhaser
#' run_generate_consensus_sequences(
#'   hybpiper_dir = test_data$hybpiper_dir,
#'   output_dir = file.path(tempdir(), "output"),
#'   namelist = test_data$namelist
#' )
#' }
create_test_dataset <- function(
  base_dir,
  samples = c("sample1", "sample2"),
  genes = c("gene001", "gene002")
) {
  # Create directories
  hybpiper_dir <- file.path(base_dir, "test_hybpiper")

  # Create HybPiper output
  create_test_hybpiper_output(hybpiper_dir, samples, genes)

  # Create targets file
  targets_file <- file.path(base_dir, "test_targets.fasta")
  create_test_targets(targets_file, genes)

  # Create namelist
  namelist_file <- file.path(base_dir, "test_namelist.txt")
  create_test_namelist(namelist_file, samples)

  message("\nComplete test dataset created!")

  list(
    base_dir = base_dir,
    hybpiper_dir = hybpiper_dir,
    targets_file = targets_file,
    namelist = namelist_file,
    samples = samples,
    genes = genes
  )
}


#' Create Complete Test Dataset Using Real HybPiper Data
#'
#' Creates a complete test dataset using real target sequences and gene names
#' from the HybPiper project's test dataset. The dataset is downloaded from
#' upstream when this function runs. This provides more realistic
#' testing data than synthetic sequences.
#'
#' @param base_dir Base directory for test data
#' @param samples Character vector of sample names. If NULL (default), uses
#'   samples from HybPiper test dataset.
#' @param n_genes Integer; number of genes to include (default: 5)
#'
#' @return Named list with paths to created files and metadata
#' @export
#'
#' @examples
#' \dontrun{
#' # Create test data using real HybPiper sequences
#' test_data <- create_real_test_dataset(tempdir())
#'
#' # Use with HybPhaser
#' run_generate_consensus_sequences(
#'   hybpiper_dir = test_data$hybpiper_dir,
#'   output_dir = file.path(tempdir(), "output"),
#'   namelist = test_data$namelist
#' )
#' }
create_real_test_dataset <- function(
  base_dir,
  samples = NULL,
  n_genes = 5
) {
  download_dir <- tempfile("hybpiper_test_dataset_")
  on.exit(unlink(download_dir, recursive = TRUE, force = TRUE), add = TRUE)

  downloaded <- .download_hybpiper_test_dataset(download_dir)
  real_targets <- downloaded$targets_file
  real_namelist <- downloaded$namelist_file

  # Read samples if not provided
  if (is.null(samples)) {
    if (file.exists(real_namelist)) {
      all_samples <- readLines(real_namelist)
      # Use first 3 samples for faster testing
      samples <- all_samples[seq_len(min(3, length(all_samples)))]
    } else {
      samples <- c("EG30", "EG98", "MWL2")
    }
  }

  # Get gene names from real targets
  targets <- seqinr::read.fasta(
    real_targets,
    seqtype = "DNA",
    as.string = FALSE,
    set.attributes = FALSE
  )
  gene_names <- gsub(".*-", "", names(targets))
  gene_names <- gene_names[seq_len(min(n_genes, length(gene_names)))]

  # Create directories
  hybpiper_dir <- file.path(base_dir, "test_hybpiper")

  # Create HybPiper output using real test data
  create_hybpiper_test_output(hybpiper_dir, samples, n_genes)

  # Copy real targets file to base_dir
  targets_file <- file.path(base_dir, "test_targets.fasta")
  file.copy(real_targets, targets_file, overwrite = TRUE)

  # Create namelist with selected samples
  namelist_file <- file.path(base_dir, "test_namelist.txt")
  writeLines(samples, namelist_file)

  message("\nComplete test dataset created using real HybPiper data!")
  message("  Targets from: downloaded HybPiper test dataset")
  message("  Gene names: ", paste(gene_names, collapse = ", "))

  list(
    base_dir = base_dir,
    hybpiper_dir = hybpiper_dir,
    targets_file = targets_file,
    namelist = namelist_file,
    samples = samples,
    genes = gene_names,
    real_data = TRUE
  )
}
