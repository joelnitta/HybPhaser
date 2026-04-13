test_that(".seq_stats counts ambiguity codes correctly", {
  # Create a temporary FASTA file with ambiguity codes
  temp_fasta <- tempfile(fileext = ".fasta")

  # Sequence with: 2 2-way ambigs (R, Y) and 1 3-way ambig (V)
  # Should count as 2 + 2 = 4
  writeLines(
    c(
      ">test_seq",
      "ATCGRATCGYATCGVATCG" # R and Y are 2-way, V is 3-way
    ),
    temp_fasta
  )

  result <- HybPhaser:::.seq_stats(temp_fasta)

  # Length should be 19 (all characters)
  expect_equal(result[1], 19)
  # Ambiguities: R(1) + Y(1) + V(1)*2 = 4
  expect_equal(result[2], 4)

  unlink(temp_fasta)
})

test_that(".seq_stats handles sequences with gaps and Ns", {
  temp_fasta <- tempfile(fileext = ".fasta")

  # Sequence with gaps, Ns - these should be removed from length
  writeLines(
    c(
      ">test_seq",
      "ATCGN-NATCG"
    ),
    temp_fasta
  )

  result <- HybPhaser:::.seq_stats(temp_fasta)

  # Length should be 8 (11 - 2 N's - 1 gap)
  expect_equal(result[1], 8)
  # No ambiguity codes
  expect_equal(result[2], 0)

  unlink(temp_fasta)
})

test_that(".seq_stats handles empty files", {
  temp_fasta <- tempfile(fileext = ".fasta")
  writeLines("", temp_fasta)

  result <- HybPhaser:::.seq_stats(temp_fasta)

  expect_true(is.na(result[1]))
  expect_true(is.na(result[2]))

  unlink(temp_fasta)
})

test_that("count_snps validates input paths", {
  expect_error(
    count_snps(
      path_to_output_folder = "/nonexistent/path",
      fasta_file_with_targets = "targets.fasta",
      targets_file_format = "DNA",
      path_to_namelist = "samples.txt"
    ),
    "Output folder not found"
  )
})

test_that("count_snps validates targets_file_format argument", {
  # Create minimal temp files
  temp_dir <- tempfile()
  dir.create(temp_dir)
  temp_targets <- tempfile(fileext = ".fasta")
  writeLines(">gene1\nATCG", temp_targets)
  temp_namelist <- tempfile()
  writeLines("sample1", temp_namelist)

  expect_error(
    count_snps(
      path_to_output_folder = temp_dir,
      fasta_file_with_targets = temp_targets,
      targets_file_format = "INVALID",
      path_to_namelist = temp_namelist
    ),
    "'arg' should be one of"
  )

  unlink(temp_dir, recursive = TRUE)
  unlink(temp_targets)
  unlink(temp_namelist)
})

test_that("count_snps processes a minimal dataset", {
  # Set up minimal test directory structure
  base_dir <- tempfile()
  data_dir <- file.path(base_dir, "01_data", "sample1", "consensus")
  dir.create(data_dir, recursive = TRUE)

  # Create a target file
  target_file <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">species1-gene1",
      "ATCGATCGATCG"
    ),
    target_file
  )

  # Create a consensus file for gene1
  consensus_file <- file.path(data_dir, "gene1.fasta")
  writeLines(
    c(
      ">gene1",
      "ATCGATRYATCG" # R and Y are ambiguity codes
    ),
    consensus_file
  )

  # Create namelist
  namelist_file <- tempfile()
  writeLines("sample1", namelist_file)

  # Run count_snps
  result <- count_snps(
    path_to_output_folder = base_dir,
    fasta_file_with_targets = target_file,
    targets_file_format = "DNA",
    path_to_namelist = namelist_file,
    intronerated_contig = FALSE
  )

  # Check results
  expect_type(result, "list")
  expect_named(result, c("tab_snps", "tab_length"))
  expect_s3_class(result$tab_snps, "data.frame")
  expect_s3_class(result$tab_length, "data.frame")
  expect_true("sample1" %in% colnames(result$tab_snps))
  expect_equal(nrow(result$tab_snps), 1) # One locus

  # Check that saved files exist
  expect_true(file.exists(file.path(
    base_dir,
    "00_R_objects",
    "Table_SNPs.Rds"
  )))
  expect_true(file.exists(file.path(
    base_dir,
    "00_R_objects",
    "Table_consensus_length.Rds"
  )))

  # Clean up
  unlink(base_dir, recursive = TRUE)
  unlink(target_file)
  unlink(namelist_file)
})
