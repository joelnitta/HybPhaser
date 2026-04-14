test_that("generate_sequence_lists validates input paths", {
  # Create minimal test data
  base_dir <- tempfile()
  dir.create(base_dir, recursive = TRUE)

  target_file <- tempfile(fileext = ".fasta")
  writeLines(c(">species1-gene1", "ATCGATCGATCGATCG"), target_file)

  namelist_file <- tempfile()
  writeLines("sample1", namelist_file)

  # Test missing output folder
  expect_error(
    generate_sequence_lists(
      path_to_output_folder = file.path(base_dir, "nonexistent"),
      fasta_file_with_targets = target_file,
      path_to_namelist = namelist_file
    ),
    "path_to_output_folder does not exist"
  )

  # Test missing targets file
  expect_error(
    generate_sequence_lists(
      path_to_output_folder = base_dir,
      fasta_file_with_targets = file.path(base_dir, "nonexistent.fasta"),
      path_to_namelist = namelist_file
    ),
    "fasta_file_with_targets does not exist"
  )

  # Test missing namelist file
  expect_error(
    generate_sequence_lists(
      path_to_output_folder = base_dir,
      fasta_file_with_targets = target_file,
      path_to_namelist = file.path(base_dir, "nonexistent.txt")
    ),
    "path_to_namelist does not exist"
  )

  # Clean up
  unlink(base_dir, recursive = TRUE)
  unlink(target_file)
  unlink(namelist_file)
})


test_that("generate_sequence_lists validates targets_file_format", {
  base_dir <- tempfile()
  dir.create(base_dir, recursive = TRUE)

  target_file <- tempfile(fileext = ".fasta")
  writeLines(c(">species1-gene1", "ATCGATCGATCGATCG"), target_file)

  namelist_file <- tempfile()
  writeLines("sample1", namelist_file)

  # Test invalid format
  expect_error(
    generate_sequence_lists(
      path_to_output_folder = base_dir,
      fasta_file_with_targets = target_file,
      targets_file_format = "INVALID",
      path_to_namelist = namelist_file
    ),
    "targets_file_format must be 'DNA' or 'AA'"
  )

  # Clean up
  unlink(base_dir, recursive = TRUE)
  unlink(target_file)
  unlink(namelist_file)
})


test_that("generate_sequence_lists requires assess_dataset results", {
  # Create test directory structure without assess_dataset results
  base_dir <- tempfile()
  dir.create(file.path(base_dir, "01_data"), recursive = TRUE)

  target_file <- tempfile(fileext = ".fasta")
  writeLines(c(">species1-gene1", "ATCGATCGATCGATCG"), target_file)

  namelist_file <- tempfile()
  writeLines("sample1", namelist_file)

  # Should fail because assess_dataset hasn't been run
  expect_error(
    generate_sequence_lists(
      path_to_output_folder = base_dir,
      fasta_file_with_targets = target_file,
      path_to_namelist = namelist_file
    ),
    "R objects directory does not exist.*Please run assess_dataset"
  )

  # Clean up
  unlink(base_dir, recursive = TRUE)
  unlink(target_file)
  unlink(namelist_file)
})


test_that("generate_sequence_lists creates output directories", {
  # Create comprehensive test structure
  base_dir <- tempfile()

  # Create HybPiper output structure
  dir.create(
    file.path(base_dir, "01_data", "sample1", "consensus"),
    recursive = TRUE
  )
  dir.create(
    file.path(base_dir, "01_data", "sample1", "contigs"),
    recursive = TRUE
  )
  dir.create(
    file.path(base_dir, "01_data", "sample2", "consensus"),
    recursive = TRUE
  )
  dir.create(
    file.path(base_dir, "01_data", "sample2", "contigs"),
    recursive = TRUE
  )

  # Create R objects directory with required files
  robjects_dir <- file.path(base_dir, "00_R_objects", "")
  dir.create(robjects_dir, recursive = TRUE)

  # Create minimal R objects
  saveRDS(character(0), file.path(robjects_dir, "outsamples_missing.Rds"))
  saveRDS(character(0), file.path(robjects_dir, "outloci_missing.Rds"))
  saveRDS(character(0), file.path(robjects_dir, "outloci_para_all.Rds"))
  saveRDS(list(), file.path(robjects_dir, "outloci_para_each.Rds"))

  # Create test SNP table
  snp_table <- matrix(0.01, nrow = 2, ncol = 2)
  rownames(snp_table) <- c("gene1", "gene2")
  colnames(snp_table) <- c("sample1", "sample2")
  saveRDS(snp_table, file.path(robjects_dir, "Table_SNPs_cleaned.Rds"))

  # Create test FASTA files
  writeLines(
    c(">sample1-gene1", "ATCGATCGRATCGYATCG"),
    file.path(base_dir, "01_data", "sample1", "consensus", "gene1.fasta")
  )
  writeLines(
    c(">sample1-gene1", "ATCGATCGRATCGYATCG"),
    file.path(base_dir, "01_data", "sample1", "contigs", "gene1.fasta")
  )
  writeLines(
    c(">sample1-gene2", "ATCGATCGWATCGNATCG"),
    file.path(base_dir, "01_data", "sample1", "consensus", "gene2.fasta")
  )
  writeLines(
    c(">sample1-gene2", "ATCGATCGWATCGNATCG"),
    file.path(base_dir, "01_data", "sample1", "contigs", "gene2.fasta")
  )

  writeLines(
    c(">sample2-gene1", "ATCGATCGMATCGKATCG"),
    file.path(base_dir, "01_data", "sample2", "consensus", "gene1.fasta")
  )
  writeLines(
    c(">sample2-gene1", "ATCGATCGMATCGKATCG"),
    file.path(base_dir, "01_data", "sample2", "contigs", "gene1.fasta")
  )
  writeLines(
    c(">sample2-gene2", "ATCGATCGRATCGYATCG"),
    file.path(base_dir, "01_data", "sample2", "consensus", "gene2.fasta")
  )
  writeLines(
    c(">sample2-gene2", "ATCGATCGRATCGYATCG"),
    file.path(base_dir, "01_data", "sample2", "contigs", "gene2.fasta")
  )

  # Create targets file
  target_file <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">species1-gene1",
      "ATCGATCGATCGATCG",
      ">species1-gene2",
      "ATCGATCGATCGATCG"
    ),
    target_file
  )

  # Create namelist
  namelist_file <- tempfile()
  writeLines(c("sample1", "sample2"), namelist_file)

  # Run function
  result <- generate_sequence_lists(
    path_to_output_folder = base_dir,
    fasta_file_with_targets = target_file,
    targets_file_format = "DNA",
    path_to_namelist = namelist_file
  )

  # Check that output directories were created
  expect_true(dir.exists(result$output_dir))
  expect_true(dir.exists(result$loci_consensus_dir))
  expect_true(dir.exists(result$loci_contigs_dir))
  expect_true(dir.exists(result$samples_consensus_dir))
  expect_true(dir.exists(result$samples_contigs_dir))

  # Check directory names
  expect_match(result$output_dir, "03_sequence_lists$")
  expect_match(result$loci_consensus_dir, "loci_consensus$")
  expect_match(result$loci_contigs_dir, "loci_contigs$")
  expect_match(result$samples_consensus_dir, "samples_consensus$")
  expect_match(result$samples_contigs_dir, "samples_contigs$")

  # Clean up
  unlink(base_dir, recursive = TRUE)
  unlink(target_file)
  unlink(namelist_file)
})


test_that("generate_sequence_lists creates locus files", {
  # Create comprehensive test structure
  base_dir <- tempfile()

  # Create HybPiper output structure
  dir.create(
    file.path(base_dir, "01_data", "sample1", "consensus"),
    recursive = TRUE
  )
  dir.create(
    file.path(base_dir, "01_data", "sample1", "contigs"),
    recursive = TRUE
  )
  dir.create(
    file.path(base_dir, "01_data", "sample2", "consensus"),
    recursive = TRUE
  )
  dir.create(
    file.path(base_dir, "01_data", "sample2", "contigs"),
    recursive = TRUE
  )

  # Create R objects directory with required files
  robjects_dir <- file.path(base_dir, "00_R_objects", "")
  dir.create(robjects_dir, recursive = TRUE)

  saveRDS(character(0), file.path(robjects_dir, "outsamples_missing.Rds"))
  saveRDS(character(0), file.path(robjects_dir, "outloci_missing.Rds"))
  saveRDS(character(0), file.path(robjects_dir, "outloci_para_all.Rds"))
  saveRDS(list(), file.path(robjects_dir, "outloci_para_each.Rds"))

  snp_table <- matrix(0.01, nrow = 2, ncol = 2)
  rownames(snp_table) <- c("gene1", "gene2")
  colnames(snp_table) <- c("sample1", "sample2")
  saveRDS(snp_table, file.path(robjects_dir, "Table_SNPs_cleaned.Rds"))

  # Create test FASTA files
  writeLines(
    c(">sample1-gene1", "ATCGATCGRATCGYATCG"),
    file.path(base_dir, "01_data", "sample1", "consensus", "gene1.fasta")
  )
  writeLines(
    c(">sample1-gene1", "ATCGATCGRATCGYATCG"),
    file.path(base_dir, "01_data", "sample1", "contigs", "gene1.fasta")
  )
  writeLines(
    c(">sample1-gene2", "ATCGATCGWATCGNATCG"),
    file.path(base_dir, "01_data", "sample1", "consensus", "gene2.fasta")
  )
  writeLines(
    c(">sample1-gene2", "ATCGATCGWATCGNATCG"),
    file.path(base_dir, "01_data", "sample1", "contigs", "gene2.fasta")
  )

  writeLines(
    c(">sample2-gene1", "ATCGATCGMATCGKATCG"),
    file.path(base_dir, "01_data", "sample2", "consensus", "gene1.fasta")
  )
  writeLines(
    c(">sample2-gene1", "ATCGATCGMATCGKATCG"),
    file.path(base_dir, "01_data", "sample2", "contigs", "gene1.fasta")
  )
  writeLines(
    c(">sample2-gene2", "ATCGATCGRATCGYATCG"),
    file.path(base_dir, "01_data", "sample2", "consensus", "gene2.fasta")
  )
  writeLines(
    c(">sample2-gene2", "ATCGATCGRATCGYATCG"),
    file.path(base_dir, "01_data", "sample2", "contigs", "gene2.fasta")
  )

  target_file <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">species1-gene1",
      "ATCGATCGATCGATCG",
      ">species1-gene2",
      "ATCGATCGATCGATCG"
    ),
    target_file
  )

  namelist_file <- tempfile()
  writeLines(c("sample1", "sample2"), namelist_file)

  # Run function
  result <- generate_sequence_lists(
    path_to_output_folder = base_dir,
    fasta_file_with_targets = target_file,
    targets_file_format = "DNA",
    path_to_namelist = namelist_file
  )

  # Check that locus files were created
  locus_consensus_files <- list.files(result$loci_consensus_dir)
  locus_contig_files <- list.files(result$loci_contigs_dir)

  expect_true(length(locus_consensus_files) > 0)
  expect_true(length(locus_contig_files) > 0)

  # Check file naming convention
  expect_true(any(grepl("gene1_consensus.fasta", locus_consensus_files)))
  expect_true(any(grepl("gene2_consensus.fasta", locus_consensus_files)))
  expect_true(any(grepl("gene1_contig.fasta", locus_contig_files)))
  expect_true(any(grepl("gene2_contig.fasta", locus_contig_files)))

  # Check that locus files contain sequences from multiple samples
  gene1_consensus <- readLines(
    file.path(result$loci_consensus_dir, "gene1_consensus.fasta")
  )
  expect_true(any(grepl("sample1", gene1_consensus)))
  expect_true(any(grepl("sample2", gene1_consensus)))

  # Check that locus name is removed from sequence headers
  expect_false(any(grepl("-gene1", gene1_consensus)))

  # Clean up
  unlink(base_dir, recursive = TRUE)
  unlink(target_file)
  unlink(namelist_file)
})


test_that("generate_sequence_lists creates sample files", {
  # Create comprehensive test structure
  base_dir <- tempfile()

  # Create HybPiper output structure
  dir.create(
    file.path(base_dir, "01_data", "sample1", "consensus"),
    recursive = TRUE
  )
  dir.create(
    file.path(base_dir, "01_data", "sample1", "contigs"),
    recursive = TRUE
  )
  dir.create(
    file.path(base_dir, "01_data", "sample2", "consensus"),
    recursive = TRUE
  )
  dir.create(
    file.path(base_dir, "01_data", "sample2", "contigs"),
    recursive = TRUE
  )

  # Create R objects directory with required files
  robjects_dir <- file.path(base_dir, "00_R_objects", "")
  dir.create(robjects_dir, recursive = TRUE)

  saveRDS(character(0), file.path(robjects_dir, "outsamples_missing.Rds"))
  saveRDS(character(0), file.path(robjects_dir, "outloci_missing.Rds"))
  saveRDS(character(0), file.path(robjects_dir, "outloci_para_all.Rds"))
  saveRDS(list(), file.path(robjects_dir, "outloci_para_each.Rds"))

  snp_table <- matrix(0.01, nrow = 2, ncol = 2)
  rownames(snp_table) <- c("gene1", "gene2")
  colnames(snp_table) <- c("sample1", "sample2")
  saveRDS(snp_table, file.path(robjects_dir, "Table_SNPs_cleaned.Rds"))

  # Create test FASTA files
  writeLines(
    c(">sample1-gene1", "ATCGATCGRATCGYATCG"),
    file.path(base_dir, "01_data", "sample1", "consensus", "gene1.fasta")
  )
  writeLines(
    c(">sample1-gene1", "ATCGATCGRATCGYATCG"),
    file.path(base_dir, "01_data", "sample1", "contigs", "gene1.fasta")
  )
  writeLines(
    c(">sample1-gene2", "ATCGATCGWATCGNATCG"),
    file.path(base_dir, "01_data", "sample1", "consensus", "gene2.fasta")
  )
  writeLines(
    c(">sample1-gene2", "ATCGATCGWATCGNATCG"),
    file.path(base_dir, "01_data", "sample1", "contigs", "gene2.fasta")
  )

  writeLines(
    c(">sample2-gene1", "ATCGATCGMATCGKATCG"),
    file.path(base_dir, "01_data", "sample2", "consensus", "gene1.fasta")
  )
  writeLines(
    c(">sample2-gene1", "ATCGATCGMATCGKATCG"),
    file.path(base_dir, "01_data", "sample2", "contigs", "gene1.fasta")
  )
  writeLines(
    c(">sample2-gene2", "ATCGATCGRATCGYATCG"),
    file.path(base_dir, "01_data", "sample2", "consensus", "gene2.fasta")
  )
  writeLines(
    c(">sample2-gene2", "ATCGATCGRATCGYATCG"),
    file.path(base_dir, "01_data", "sample2", "contigs", "gene2.fasta")
  )

  target_file <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">species1-gene1",
      "ATCGATCGATCGATCG",
      ">species1-gene2",
      "ATCGATCGATCGATCG"
    ),
    target_file
  )

  namelist_file <- tempfile()
  writeLines(c("sample1", "sample2"), namelist_file)

  # Run function
  result <- generate_sequence_lists(
    path_to_output_folder = base_dir,
    fasta_file_with_targets = target_file,
    targets_file_format = "DNA",
    path_to_namelist = namelist_file
  )

  # Check that sample files were created
  sample_consensus_files <- list.files(result$samples_consensus_dir)
  sample_contig_files <- list.files(result$samples_contigs_dir)

  expect_true(length(sample_consensus_files) > 0)
  expect_true(length(sample_contig_files) > 0)

  # Check file naming convention
  expect_true(any(grepl("sample1_consensus.fasta", sample_consensus_files)))
  expect_true(any(grepl("sample2_consensus.fasta", sample_consensus_files)))
  expect_true(any(grepl("sample1_contig", sample_contig_files)))
  expect_true(any(grepl("sample2_contig", sample_contig_files)))

  # Check that sample files contain sequences from multiple loci
  sample1_consensus <- readLines(
    file.path(result$samples_consensus_dir, "sample1_consensus.fasta")
  )
  expect_true(any(grepl("gene1", sample1_consensus)))
  expect_true(any(grepl("gene2", sample1_consensus)))

  # Clean up
  unlink(base_dir, recursive = TRUE)
  unlink(target_file)
  unlink(namelist_file)
})


test_that("generate_sequence_lists handles subset_name parameter", {
  # Create comprehensive test structure
  base_dir <- tempfile()

  # Create HybPiper output structure
  dir.create(
    file.path(base_dir, "01_data", "sample1", "consensus"),
    recursive = TRUE
  )
  dir.create(
    file.path(base_dir, "01_data", "sample1", "contigs"),
    recursive = TRUE
  )

  # Create R objects directory with subset name
  robjects_dir <- file.path(base_dir, "00_R_objects", "filtered")
  dir.create(robjects_dir, recursive = TRUE)

  saveRDS(character(0), file.path(robjects_dir, "outsamples_missing.Rds"))
  saveRDS(character(0), file.path(robjects_dir, "outloci_missing.Rds"))
  saveRDS(character(0), file.path(robjects_dir, "outloci_para_all.Rds"))
  saveRDS(list(), file.path(robjects_dir, "outloci_para_each.Rds"))

  snp_table <- matrix(0.01, nrow = 1, ncol = 1)
  rownames(snp_table) <- "gene1"
  colnames(snp_table) <- "sample1"
  saveRDS(snp_table, file.path(robjects_dir, "Table_SNPs_cleaned.Rds"))

  # Create test FASTA file
  writeLines(
    c(">sample1-gene1", "ATCGATCGRATCGYATCG"),
    file.path(base_dir, "01_data", "sample1", "consensus", "gene1.fasta")
  )
  writeLines(
    c(">sample1-gene1", "ATCGATCGRATCGYATCG"),
    file.path(base_dir, "01_data", "sample1", "contigs", "gene1.fasta")
  )

  target_file <- tempfile(fileext = ".fasta")
  writeLines(c(">species1-gene1", "ATCGATCGATCGATCG"), target_file)

  namelist_file <- tempfile()
  writeLines("sample1", namelist_file)

  # Run function with subset name
  result <- generate_sequence_lists(
    path_to_output_folder = base_dir,
    fasta_file_with_targets = target_file,
    targets_file_format = "DNA",
    path_to_namelist = namelist_file,
    subset_name = "filtered"
  )

  # Check that output directory includes subset name
  expect_match(result$output_dir, "03_sequence_lists_filtered$")

  # Check that directories were created
  expect_true(dir.exists(result$output_dir))
  expect_true(dir.exists(result$loci_consensus_dir))

  # Clean up
  unlink(base_dir, recursive = TRUE)
  unlink(target_file)
  unlink(namelist_file)
})


test_that("generate_sequence_lists applies missing data filters", {
  # Create comprehensive test structure
  base_dir <- tempfile()

  # Create HybPiper output structure for 3 samples
  for (sample in c("sample1", "sample2", "sample3")) {
    dir.create(
      file.path(base_dir, "01_data", sample, "consensus"),
      recursive = TRUE
    )
    dir.create(
      file.path(base_dir, "01_data", sample, "contigs"),
      recursive = TRUE
    )
  }

  # Create R objects directory with filters
  robjects_dir <- file.path(base_dir, "00_R_objects", "")
  dir.create(robjects_dir, recursive = TRUE)

  # Mark sample3 as removed for missing data
  saveRDS("sample3", file.path(robjects_dir, "outsamples_missing.Rds"))
  # Mark gene3 as removed for missing data
  saveRDS("gene3", file.path(robjects_dir, "outloci_missing.Rds"))
  saveRDS(character(0), file.path(robjects_dir, "outloci_para_all.Rds"))
  saveRDS(list(), file.path(robjects_dir, "outloci_para_each.Rds"))

  # SNP table with only retained samples and loci
  snp_table <- matrix(0.01, nrow = 2, ncol = 2)
  rownames(snp_table) <- c("gene1", "gene2")
  colnames(snp_table) <- c("sample1", "sample2")
  saveRDS(snp_table, file.path(robjects_dir, "Table_SNPs_cleaned.Rds"))

  # Create test FASTA files for all samples and genes
  for (sample in c("sample1", "sample2", "sample3")) {
    for (gene in c("gene1", "gene2", "gene3")) {
      writeLines(
        c(paste0(">", sample, "-", gene), "ATCGATCGRATCGYATCG"),
        file.path(
          base_dir,
          "01_data",
          sample,
          "consensus",
          paste0(gene, ".fasta")
        )
      )
      writeLines(
        c(paste0(">", sample, "-", gene), "ATCGATCGRATCGYATCG"),
        file.path(
          base_dir,
          "01_data",
          sample,
          "contigs",
          paste0(gene, ".fasta")
        )
      )
    }
  }

  target_file <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">species1-gene1",
      "ATCGATCGATCGATCG",
      ">species1-gene2",
      "ATCGATCGATCGATCG",
      ">species1-gene3",
      "ATCGATCGATCGATCG"
    ),
    target_file
  )

  namelist_file <- tempfile()
  writeLines(c("sample1", "sample2", "sample3"), namelist_file)

  # Run function
  result <- generate_sequence_lists(
    path_to_output_folder = base_dir,
    fasta_file_with_targets = target_file,
    targets_file_format = "DNA",
    path_to_namelist = namelist_file
  )

  # Check that gene3 locus file was NOT created
  locus_files <- list.files(result$loci_consensus_dir)
  expect_false(any(grepl("gene3", locus_files)))
  expect_true(any(grepl("gene1", locus_files)))
  expect_true(any(grepl("gene2", locus_files)))

  # Check that sample3 file was NOT created
  sample_files <- list.files(result$samples_consensus_dir)
  expect_false(any(grepl("sample3", sample_files)))
  expect_true(any(grepl("sample1", sample_files)))
  expect_true(any(grepl("sample2", sample_files)))

  # Check that gene1 locus file doesn't contain sample3
  gene1_consensus <- readLines(
    file.path(result$loci_consensus_dir, "gene1_consensus.fasta")
  )
  expect_false(any(grepl("sample3", gene1_consensus)))
  expect_true(any(grepl("sample1", gene1_consensus)))

  # Check that sample1 file doesn't contain gene3
  sample1_consensus <- readLines(
    file.path(result$samples_consensus_dir, "sample1_consensus.fasta")
  )
  expect_false(any(grepl("gene3", sample1_consensus)))
  expect_true(any(grepl("gene1", sample1_consensus)))

  # Clean up
  unlink(base_dir, recursive = TRUE)
  unlink(target_file)
  unlink(namelist_file)
})


test_that("generate_sequence_lists applies paralog filters", {
  # Create comprehensive test structure
  base_dir <- tempfile()

  # Create HybPiper output structure for 2 samples
  for (sample in c("sample1", "sample2")) {
    dir.create(
      file.path(base_dir, "01_data", sample, "consensus"),
      recursive = TRUE
    )
    dir.create(
      file.path(base_dir, "01_data", sample, "contigs"),
      recursive = TRUE
    )
  }

  # Create R objects directory with paralog filters
  robjects_dir <- file.path(base_dir, "00_R_objects", "")
  dir.create(robjects_dir, recursive = TRUE)

  saveRDS(character(0), file.path(robjects_dir, "outsamples_missing.Rds"))
  saveRDS(character(0), file.path(robjects_dir, "outloci_missing.Rds"))
  # Mark gene3 as paralog for all samples
  saveRDS("gene3", file.path(robjects_dir, "outloci_para_all.Rds"))
  # Mark gene2 as paralog for sample2 only
  outloci_para_each <- list(sample2 = c(gene2 = "gene2"))
  saveRDS(outloci_para_each, file.path(robjects_dir, "outloci_para_each.Rds"))

  # SNP table
  snp_table <- matrix(0.01, nrow = 3, ncol = 2)
  rownames(snp_table) <- c("gene1", "gene2", "gene3")
  colnames(snp_table) <- c("sample1", "sample2")
  saveRDS(snp_table, file.path(robjects_dir, "Table_SNPs_cleaned.Rds"))

  # Create test FASTA files for all samples and genes
  for (sample in c("sample1", "sample2")) {
    for (gene in c("gene1", "gene2", "gene3")) {
      writeLines(
        c(paste0(">", sample, "-", gene), "ATCGATCGRATCGYATCG"),
        file.path(
          base_dir,
          "01_data",
          sample,
          "consensus",
          paste0(gene, ".fasta")
        )
      )
      writeLines(
        c(paste0(">", sample, "-", gene), "ATCGATCGRATCGYATCG"),
        file.path(
          base_dir,
          "01_data",
          sample,
          "contigs",
          paste0(gene, ".fasta")
        )
      )
    }
  }

  target_file <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">species1-gene1",
      "ATCGATCGATCGATCG",
      ">species1-gene2",
      "ATCGATCGATCGATCG",
      ">species1-gene3",
      "ATCGATCGATCGATCG"
    ),
    target_file
  )

  namelist_file <- tempfile()
  writeLines(c("sample1", "sample2"), namelist_file)

  # Run function
  result <- generate_sequence_lists(
    path_to_output_folder = base_dir,
    fasta_file_with_targets = target_file,
    targets_file_format = "DNA",
    path_to_namelist = namelist_file
  )

  # Check that gene3 locus file was NOT created (paralog for all)
  locus_files <- list.files(result$loci_consensus_dir)
  expect_false(any(grepl("gene3", locus_files)))

  # Check that gene2 locus file was created but doesn't contain sample2
  gene2_consensus <- readLines(
    file.path(result$loci_consensus_dir, "gene2_consensus.fasta")
  )
  expect_false(any(grepl("sample2", gene2_consensus)))
  expect_true(any(grepl("sample1", gene2_consensus)))

  # Check that sample2 file doesn't contain gene2 (paralog for this sample)
  sample2_consensus <- readLines(
    file.path(result$samples_consensus_dir, "sample2_consensus.fasta")
  )
  expect_false(any(grepl("gene2", sample2_consensus)))
  expect_true(any(grepl("gene1", sample2_consensus)))

  # Check that sample1 file contains gene2 (not a paralog for sample1)
  sample1_consensus <- readLines(
    file.path(result$samples_consensus_dir, "sample1_consensus.fasta")
  )
  expect_true(any(grepl("gene2", sample1_consensus)))

  # Check that neither sample file contains gene3
  expect_false(any(grepl("gene3", sample1_consensus)))
  expect_false(any(grepl("gene3", sample2_consensus)))

  # Clean up
  unlink(base_dir, recursive = TRUE)
  unlink(target_file)
  unlink(namelist_file)
})


test_that("generate_sequence_lists returns correct structure", {
  # Create minimal test structure
  base_dir <- tempfile()

  dir.create(
    file.path(base_dir, "01_data", "sample1", "consensus"),
    recursive = TRUE
  )
  dir.create(
    file.path(base_dir, "01_data", "sample1", "contigs"),
    recursive = TRUE
  )

  robjects_dir <- file.path(base_dir, "00_R_objects", "")
  dir.create(robjects_dir, recursive = TRUE)

  saveRDS(character(0), file.path(robjects_dir, "outsamples_missing.Rds"))
  saveRDS(character(0), file.path(robjects_dir, "outloci_missing.Rds"))
  saveRDS(character(0), file.path(robjects_dir, "outloci_para_all.Rds"))
  saveRDS(list(), file.path(robjects_dir, "outloci_para_each.Rds"))

  snp_table <- matrix(0.01, nrow = 1, ncol = 1)
  rownames(snp_table) <- "gene1"
  colnames(snp_table) <- "sample1"
  saveRDS(snp_table, file.path(robjects_dir, "Table_SNPs_cleaned.Rds"))

  writeLines(
    c(">sample1-gene1", "ATCGATCGRATCGYATCG"),
    file.path(base_dir, "01_data", "sample1", "consensus", "gene1.fasta")
  )
  writeLines(
    c(">sample1-gene1", "ATCGATCGRATCGYATCG"),
    file.path(base_dir, "01_data", "sample1", "contigs", "gene1.fasta")
  )

  target_file <- tempfile(fileext = ".fasta")
  writeLines(c(">species1-gene1", "ATCGATCGATCGATCG"), target_file)

  namelist_file <- tempfile()
  writeLines("sample1", namelist_file)

  # Run function
  result <- generate_sequence_lists(
    path_to_output_folder = base_dir,
    fasta_file_with_targets = target_file,
    targets_file_format = "DNA",
    path_to_namelist = namelist_file
  )

  # Check return structure
  expect_type(result, "list")
  expect_named(
    result,
    c(
      "output_dir",
      "loci_consensus_dir",
      "loci_contigs_dir",
      "samples_consensus_dir",
      "samples_contigs_dir",
      "n_loci_files",
      "n_sample_files",
      "loci_removed",
      "samples_removed"
    )
  )

  # Check types
  expect_type(result$output_dir, "character")
  expect_type(result$n_loci_files, "double")
  expect_type(result$n_sample_files, "double")
  expect_type(result$loci_removed, "character")
  expect_type(result$samples_removed, "character")

  # Check counts are positive
  expect_true(result$n_loci_files > 0)
  expect_true(result$n_sample_files > 0)

  # Clean up
  unlink(base_dir, recursive = TRUE)
  unlink(target_file)
  unlink(namelist_file)
})


test_that("generate_sequence_lists handles missing files gracefully", {
  # This test simulates the real-world scenario where assess_dataset()
  # deletes files from disk but they're still referenced in the SNP table.
  # The function should skip missing files without error.

  base_dir <- tempfile()

  # Create HybPiper output structure for 3 samples
  for (sample in c("sample1", "sample2", "sample3")) {
    dir.create(
      file.path(base_dir, "01_data", sample, "consensus"),
      recursive = TRUE
    )
    dir.create(
      file.path(base_dir, "01_data", sample, "contigs"),
      recursive = TRUE
    )
  }

  # Create R objects directory
  robjects_dir <- file.path(base_dir, "00_R_objects", "")
  dir.create(robjects_dir, recursive = TRUE)

  # Mark sample3 as removed (simulating missing data filter)
  saveRDS("sample3", file.path(robjects_dir, "outsamples_missing.Rds"))
  saveRDS(character(0), file.path(robjects_dir, "outloci_missing.Rds"))
  # Mark gene3 as paralog for all (simulating paralog filter)
  saveRDS("gene3", file.path(robjects_dir, "outloci_para_all.Rds"))
  saveRDS(list(), file.path(robjects_dir, "outloci_para_each.Rds"))

  # SNP table includes ALL samples and loci (even those marked for
  # removal)
  # This simulates the real state after assess_dataset() where the table
  # contains items that have been filtered but files were already deleted
  snp_table <- matrix(0.01, nrow = 3, ncol = 3)
  rownames(snp_table) <- c("gene1", "gene2", "gene3")
  colnames(snp_table) <- c("sample1", "sample2", "sample3")
  saveRDS(snp_table, file.path(robjects_dir, "Table_SNPs_cleaned.Rds"))

  # Create FASTA files ONLY for non-removed items
  # Simulate that assess_dataset() already deleted gene3 and sample3
  # files
  for (sample in c("sample1", "sample2")) {
    # NOT sample3
    for (gene in c("gene1", "gene2")) {
      # NOT gene3
      writeLines(
        c(paste0(">", sample, "-", gene), "ATCGATCGRATCGYATCG"),
        file.path(
          base_dir,
          "01_data",
          sample,
          "consensus",
          paste0(gene, ".fasta")
        )
      )
      writeLines(
        c(paste0(">", sample, "-", gene), "ATCGATCGRATCGYATCG"),
        file.path(
          base_dir,
          "01_data",
          sample,
          "contigs",
          paste0(gene, ".fasta")
        )
      )
    }
  }

  target_file <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">species1-gene1",
      "ATCGATCGATCGATCG",
      ">species1-gene2",
      "ATCGATCGATCGATCG",
      ">species1-gene3",
      "ATCGATCGATCGATCG"
    ),
    target_file
  )

  namelist_file <- tempfile()
  writeLines(c("sample1", "sample2", "sample3"), namelist_file)

  # Run function - should NOT fail even though:
  # 1. gene3 files don't exist (paralog, deleted)
  # 2. sample3 files don't exist (removed sample, deleted)
  # 3. Both are still in the SNP table
  expect_no_error({
    result <- generate_sequence_lists(
      path_to_output_folder = base_dir,
      fasta_file_with_targets = target_file,
      targets_file_format = "DNA",
      path_to_namelist = namelist_file
    )
  })

  # Verify output was created successfully
  expect_true(dir.exists(result$output_dir))
  expect_true(result$n_loci_files > 0)
  expect_true(result$n_sample_files > 0)

  # Verify gene3 was not included (paralog for all)
  locus_files <- list.files(result$loci_consensus_dir)
  expect_false(any(grepl("gene3", locus_files)))
  expect_true(any(grepl("gene1", locus_files)))
  expect_true(any(grepl("gene2", locus_files)))

  # Verify sample3 was not included (removed sample)
  sample_files <- list.files(result$samples_consensus_dir)
  expect_false(any(grepl("sample3", sample_files)))
  expect_true(any(grepl("sample1", sample_files)))
  expect_true(any(grepl("sample2", sample_files)))

  # Verify the files that were created don't reference removed items
  gene1_consensus <- readLines(
    file.path(result$loci_consensus_dir, "gene1_consensus.fasta")
  )
  expect_false(any(grepl("sample3", gene1_consensus)))

  sample1_consensus <- readLines(
    file.path(result$samples_consensus_dir, "sample1_consensus.fasta")
  )
  expect_false(any(grepl("gene3", sample1_consensus)))

  # Clean up
  unlink(base_dir, recursive = TRUE)
  unlink(target_file)
  unlink(namelist_file)
})
