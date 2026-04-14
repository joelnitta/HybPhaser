test_that("assess_dataset validates input files", {
  # Non-existent targets file
  expect_error(
    assess_dataset(
      snp_table = matrix(0.1, nrow = 5, ncol = 3),
      length_table = matrix(100, nrow = 5, ncol = 3),
      targets_file = "/nonexistent/targets.fasta",
      targets_type = "DNA",
      output_dir = tempdir()
    ),
    "targets_file does not exist"
  )
})

test_that("assess_dataset validates targets_type argument", {
  # Create minimal targets file
  temp_targets <- tempfile(fileext = ".fasta")
  writeLines(c(">gene1", "ATCG"), temp_targets)

  expect_error(
    assess_dataset(
      snp_table = matrix(0.1, nrow = 5, ncol = 3),
      length_table = matrix(100, nrow = 5, ncol = 3),
      targets_file = temp_targets,
      targets_type = "INVALID",
      output_dir = tempdir()
    ),
    "should be one of"
  )

  unlink(temp_targets)
})

test_that("assess_dataset loads SNP table from file path", {
  temp_dir <- tempfile()
  dir.create(temp_dir, recursive = TRUE)

  # Create test data
  snp_mat <- matrix(
    c(0.01, 0.02, 0.03, 0.015, 0.025, 0.035, 0.02, 0.03, 0.04),
    nrow = 3,
    ncol = 3,
    dimnames = list(
      c("gene1", "gene2", "gene3"),
      c("sample1", "sample2", "sample3")
    )
  )
  length_mat <- matrix(
    c(500, 600, 550, 480, 590, 540, 510, 610, 560),
    nrow = 3,
    ncol = 3,
    dimnames = list(
      c("gene1", "gene2", "gene3"),
      c("sample1", "sample2", "sample3")
    )
  )

  snp_file <- file.path(temp_dir, "snps.rds")
  length_file <- file.path(temp_dir, "lengths.rds")
  saveRDS(snp_mat, snp_file)
  saveRDS(length_mat, length_file)

  # Create targets file
  temp_targets <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">species1-gene1",
      "ATCGATCGATCGATCGATCGATCGATCGATCGATCG",
      ">species1-gene2",
      "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA",
      ">species1-gene3",
      "TGACTGACTGACTGACTGACTGACTGACTGACTGAC"
    ),
    temp_targets
  )

  # Run assessment
  result <- assess_dataset(
    snp_table = snp_file,
    length_table = length_file,
    targets_file = temp_targets,
    targets_type = "DNA",
    output_dir = temp_dir,
    min_loci_per_sample_prop = 0.5,
    min_samples_per_locus_prop = 0.5
  )

  expect_type(result, "list")
  expect_named(
    result,
    c(
      "snp_table_cleaned",
      "length_table_cleaned",
      "samples_removed",
      "loci_removed_missing",
      "loci_removed_paralogs_all",
      "loci_removed_paralogs_each",
      "summary_table"
    )
  )

  unlink(temp_dir, recursive = TRUE)
  unlink(temp_targets)
})

test_that("assess_dataset works with matrix input directly", {
  temp_dir <- tempfile()
  dir.create(temp_dir, recursive = TRUE)

  # Create test data
  snp_mat <- matrix(
    c(0.01, 0.02, 0.03, 0.015, 0.025, 0.035, 0.02, 0.03),
    nrow = 4,
    ncol = 2,
    dimnames = list(
      c("gene1", "gene2", "gene3", "gene4"),
      c("sample1", "sample2")
    )
  )
  length_mat <- matrix(
    c(500, 600, 550, 520, 480, 590, 540, 510),
    nrow = 4,
    ncol = 2,
    dimnames = list(
      c("gene1", "gene2", "gene3", "gene4"),
      c("sample1", "sample2")
    )
  )

  # Create targets file
  temp_targets <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">species1-gene1",
      paste(rep("A", 600), collapse = ""),
      ">species1-gene2",
      paste(rep("T", 650), collapse = ""),
      ">species1-gene3",
      paste(rep("G", 580), collapse = ""),
      ">species1-gene4",
      paste(rep("C", 550), collapse = "")
    ),
    temp_targets
  )

  # Run assessment with matrices directly
  result <- assess_dataset(
    snp_table = snp_mat,
    length_table = length_mat,
    targets_file = temp_targets,
    targets_type = "DNA",
    output_dir = temp_dir,
    min_loci_per_sample_prop = 0.3,
    min_samples_per_locus_prop = 0.3
  )

  expect_type(result, "list")
  expect_s3_class(result$summary_table, "data.frame")
  expect_equal(nrow(result$summary_table), 2) # 2 samples

  unlink(temp_dir, recursive = TRUE)
  unlink(temp_targets)
})

test_that("assess_dataset handles missing data correctly", {
  temp_dir <- tempfile()
  dir.create(temp_dir, recursive = TRUE)

  # Create data with NAs (missing loci)
  # gene3 is completely missing (all NAs)
  snp_mat <- matrix(
    c(
      0.01,
      0.025,
      NA, # sample1: gene1, gene2, gene3
      NA,
      0.035,
      NA, # sample2
      0.03,
      NA,
      NA, # sample3
      0.015,
      0.03,
      NA # sample4
    ),
    nrow = 3,
    ncol = 4,
    dimnames = list(
      c("gene1", "gene2", "gene3"),
      c("sample1", "sample2", "sample3", "sample4")
    )
  )
  length_mat <- matrix(
    c(
      500,
      480,
      NA, # sample1
      NA,
      590,
      NA, # sample2
      550,
      NA,
      NA, # sample3
      520,
      510,
      NA # sample4
    ),
    nrow = 3,
    ncol = 4,
    dimnames = list(
      c("gene1", "gene2", "gene3"),
      c("sample1", "sample2", "sample3", "sample4")
    )
  )

  temp_targets <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">species1-gene1",
      paste(rep("A", 600), collapse = ""),
      ">species1-gene2",
      paste(rep("T", 600), collapse = ""),
      ">species1-gene3",
      paste(rep("G", 600), collapse = "")
    ),
    temp_targets
  )

  result <- assess_dataset(
    snp_table = snp_mat,
    length_table = length_mat,
    targets_file = temp_targets,
    targets_type = "DNA",
    output_dir = temp_dir,
    min_loci_per_sample_prop = 0.3,
    min_samples_per_locus_prop = 0.3
  )

  # Check that completely failed locus was removed
  expect_true("gene3" %in% result$loci_removed_missing)

  # Check output files exist
  expect_true(file.exists(file.path(
    temp_dir,
    "02_assessment",
    "1_Summary_missing_data.txt"
  )))
  expect_true(file.exists(file.path(
    temp_dir,
    "02_assessment",
    "0_Table_SNPs.csv"
  )))

  unlink(temp_dir, recursive = TRUE)
  unlink(temp_targets)
})

test_that("assess_dataset applies paralog threshold correctly", {
  temp_dir <- tempfile()
  dir.create(temp_dir, recursive = TRUE)

  # Create data with one high-SNP locus (putative paralog)
  snp_mat <- matrix(
    c(
      0.01,
      0.02,
      0.015,
      0.025, # normal loci
      0.15,
      0.18,
      0.16,
      0.20
    ), # high SNP locus (paralog)
    nrow = 2,
    ncol = 4,
    dimnames = list(
      c("gene1", "gene2"),
      c("sample1", "sample2", "sample3", "sample4")
    )
  )
  length_mat <- matrix(
    c(500, 600, 550, 520, 480, 590, 540, 510),
    nrow = 2,
    ncol = 4,
    dimnames = list(
      c("gene1", "gene2"),
      c("sample1", "sample2", "sample3", "sample4")
    )
  )

  temp_targets <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">species1-gene1",
      paste(rep("A", 600), collapse = ""),
      ">species1-gene2",
      paste(rep("T", 600), collapse = "")
    ),
    temp_targets
  )

  # Use numeric threshold
  result <- assess_dataset(
    snp_table = snp_mat,
    length_table = length_mat,
    targets_file = temp_targets,
    targets_type = "DNA",
    output_dir = temp_dir,
    min_loci_per_sample_prop = 0.3,
    min_samples_per_locus_prop = 0.3,
    paralog_threshold = 0.10 # Remove loci with >10% mean SNPs
  )

  # gene2 should be flagged as paralog
  expect_true("gene2" %in% result$loci_removed_paralogs_all)
  expect_false("gene1" %in% result$loci_removed_paralogs_all)

  unlink(temp_dir, recursive = TRUE)
  unlink(temp_targets)
})

test_that("assess_dataset handles subset_name parameter", {
  temp_dir <- tempfile()
  dir.create(temp_dir, recursive = TRUE)

  snp_mat <- matrix(
    0.02,
    nrow = 2,
    ncol = 2,
    dimnames = list(c("gene1", "gene2"), c("sample1", "sample2"))
  )
  length_mat <- matrix(
    500,
    nrow = 2,
    ncol = 2,
    dimnames = list(c("gene1", "gene2"), c("sample1", "sample2"))
  )

  temp_targets <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">sp1-gene1",
      paste(rep("A", 600), collapse = ""),
      ">sp1-gene2",
      paste(rep("T", 600), collapse = "")
    ),
    temp_targets
  )

  result <- assess_dataset(
    snp_table = snp_mat,
    length_table = length_mat,
    targets_file = temp_targets,
    targets_type = "DNA",
    output_dir = temp_dir,
    subset_name = "optimized"
  )

  # Check that subset folder was created
  expect_true(dir.exists(file.path(temp_dir, "02_assessment_optimized")))
  expect_true(dir.exists(file.path(temp_dir, "00_R_objects", "optimized")))

  unlink(temp_dir, recursive = TRUE)
  unlink(temp_targets)
})

test_that("assess_dataset removes sample outliers when enabled", {
  temp_dir <- tempfile()
  dir.create(temp_dir, recursive = TRUE)

  # Create data where sample2 has one extreme outlier locus
  snp_mat <- matrix(
    c(0.01, 0.02, 0.015, 0.018, 0.02, 0.30, 0.022, 0.025), # gene2-sample2 is outlier
    nrow = 4,
    ncol = 2,
    dimnames = list(
      c("gene1", "gene2", "gene3", "gene4"),
      c("sample1", "sample2")
    )
  )
  length_mat <- matrix(
    500,
    nrow = 4,
    ncol = 2,
    dimnames = list(
      c("gene1", "gene2", "gene3", "gene4"),
      c("sample1", "sample2")
    )
  )

  temp_targets <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">sp1-gene1",
      paste(rep("A", 600), collapse = ""),
      ">sp1-gene2",
      paste(rep("T", 600), collapse = ""),
      ">sp1-gene3",
      paste(rep("G", 600), collapse = ""),
      ">sp1-gene4",
      paste(rep("C", 600), collapse = "")
    ),
    temp_targets
  )

  result <- assess_dataset(
    snp_table = snp_mat,
    length_table = length_mat,
    targets_file = temp_targets,
    targets_type = "DNA",
    output_dir = temp_dir,
    remove_sample_outliers = TRUE
  )

  # Check that per-sample outliers were detected
  expect_true(length(result$loci_removed_paralogs_each) > 0)

  unlink(temp_dir, recursive = TRUE)
  unlink(temp_targets)
})

test_that("assess_dataset generates expected output files", {
  temp_dir <- tempfile()
  dir.create(temp_dir, recursive = TRUE)

  snp_mat <- matrix(
    runif(12, 0.01, 0.05),
    nrow = 4,
    ncol = 3,
    dimnames = list(paste0("gene", 1:4), paste0("sample", 1:3))
  )
  length_mat <- matrix(
    runif(12, 400, 600),
    nrow = 4,
    ncol = 3,
    dimnames = list(paste0("gene", 1:4), paste0("sample", 1:3))
  )

  temp_targets <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">sp1-gene1",
      paste(rep("A", 600), collapse = ""),
      ">sp1-gene2",
      paste(rep("T", 600), collapse = ""),
      ">sp1-gene3",
      paste(rep("G", 600), collapse = ""),
      ">sp1-gene4",
      paste(rep("C", 600), collapse = "")
    ),
    temp_targets
  )

  result <- assess_dataset(
    snp_table = snp_mat,
    length_table = length_mat,
    targets_file = temp_targets,
    targets_type = "DNA",
    output_dir = temp_dir
  )

  assess_dir <- file.path(temp_dir, "02_assessment")

  # Check that all expected files exist
  expect_true(file.exists(file.path(
    assess_dir,
    "1_Data_recovered_overview.pdf"
  )))
  expect_true(file.exists(file.path(
    assess_dir,
    "1_Data_recovered_overview.png"
  )))
  expect_true(file.exists(file.path(
    assess_dir,
    "1_Data_recovered_per_sample.csv"
  )))
  expect_true(file.exists(file.path(
    assess_dir,
    "1_Data_recovered_per_locus.csv"
  )))
  expect_true(file.exists(file.path(assess_dir, "1_Summary_missing_data.txt")))
  expect_true(file.exists(file.path(
    assess_dir,
    "2a_Paralogs_for_all_samples.pdf"
  )))
  expect_true(file.exists(file.path(
    assess_dir,
    "2b_Paralogs_for_each_sample.pdf"
  )))
  expect_true(file.exists(file.path(assess_dir, "2_Summary_Paralogs.txt")))
  expect_true(file.exists(file.path(assess_dir, "3_LH_vs_AD.pdf")))
  expect_true(file.exists(file.path(assess_dir, "3_varLH_vs_AD.pdf")))
  expect_true(file.exists(file.path(assess_dir, "4_Summary_table.csv")))
  expect_true(file.exists(file.path(assess_dir, "0_Table_SNPs.csv")))
  expect_true(file.exists(file.path(
    assess_dir,
    "0_Table_consensus_length.csv"
  )))
  expect_true(file.exists(file.path(
    assess_dir,
    "0_namelist_included_samples.txt"
  )))

  # Check R objects directory
  robjects_dir <- file.path(temp_dir, "00_R_objects")
  expect_true(file.exists(file.path(robjects_dir, "Table_SNPs_cleaned.Rds")))
  expect_true(file.exists(file.path(
    robjects_dir,
    "Table_consensus_length_cleaned.Rds"
  )))
  expect_true(file.exists(file.path(robjects_dir, "Summary_table.Rds")))

  unlink(temp_dir, recursive = TRUE)
  unlink(temp_targets)
})

test_that("assess_dataset summary table has correct structure", {
  temp_dir <- tempfile()
  dir.create(temp_dir, recursive = TRUE)

  snp_mat <- matrix(
    c(0.01, 0.02, 0.00, 0.015, 0.025, 0.00, 0.022, 0.03),
    nrow = 4,
    ncol = 2,
    dimnames = list(paste0("gene", 1:4), c("sample1", "sample2"))
  )
  length_mat <- matrix(
    c(500, 600, 550, 520, 480, 590, 540, 510),
    nrow = 4,
    ncol = 2,
    dimnames = list(paste0("gene", 1:4), c("sample1", "sample2"))
  )

  temp_targets <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">sp1-gene1",
      paste(rep("A", 600), collapse = ""),
      ">sp1-gene2",
      paste(rep("T", 600), collapse = ""),
      ">sp1-gene3",
      paste(rep("G", 600), collapse = ""),
      ">sp1-gene4",
      paste(rep("C", 600), collapse = "")
    ),
    temp_targets
  )

  result <- assess_dataset(
    snp_table = snp_mat,
    length_table = length_mat,
    targets_file = temp_targets,
    targets_type = "DNA",
    output_dir = temp_dir
  )

  summary <- result$summary_table

  # Check columns
  expect_true(all(
    c(
      "sample",
      "bp",
      "bpoftarget",
      "paralogs_all",
      "paralogs_each",
      "nloci",
      "allele_divergence",
      "locus_heterozygosity"
    ) %in%
      colnames(summary)
  ))

  # Check number of rows
  expect_equal(nrow(summary), 2)

  # Check sample names
  expect_true(all(c("sample1", "sample2") %in% summary$sample))

  # Check that metrics are numeric and reasonable
  expect_true(all(summary$bp > 0))
  expect_true(all(summary$bpoftarget >= 0 & summary$bpoftarget <= 100))
  expect_true(all(summary$allele_divergence >= 0))
  expect_true(all(
    summary$locus_heterozygosity >= 0 &
      summary$locus_heterozygosity <= 100
  ))

  unlink(temp_dir, recursive = TRUE)
  unlink(temp_targets)
})

test_that("assess_dataset handles data frames from count_snps", {
  # Regression test: ensure assess_dataset can handle the data frame output
  # from count_snps without colSums errors
  temp_dir <- tempfile()
  dir.create(temp_dir, recursive = TRUE)

  # Simulate count_snps output: data frames with only numeric columns
  snp_df <- data.frame(
    sample1 = c(0.01, 0.02, 0.015),
    sample2 = c(0.025, 0.03, 0.02),
    row.names = c("gene1", "gene2", "gene3")
  )

  length_df <- data.frame(
    sample1 = c(500, 600, 550),
    sample2 = c(520, 610, 540),
    row.names = c("gene1", "gene2", "gene3")
  )

  # Verify input is data frame (not matrix)
  expect_s3_class(snp_df, "data.frame")
  expect_s3_class(length_df, "data.frame")

  # Verify all columns are numeric
  expect_true(all(sapply(snp_df, is.numeric)))
  expect_true(all(sapply(length_df, is.numeric)))

  temp_targets <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">sp1-gene1",
      paste(rep("A", 600), collapse = ""),
      ">sp1-gene2",
      paste(rep("T", 600), collapse = ""),
      ">sp1-gene3",
      paste(rep("G", 600), collapse = "")
    ),
    temp_targets
  )

  # This should not produce colSums error
  expect_no_error(
    result <- assess_dataset(
      snp_table = snp_df,
      length_table = length_df,
      targets_file = temp_targets,
      targets_type = "DNA",
      output_dir = temp_dir,
      min_loci_per_sample_prop = 0.3,
      min_samples_per_locus_prop = 0.3
    )
  )

  # Verify output structure
  expect_type(result, "list")
  expect_s3_class(result$summary_table, "data.frame")

  # Verify summary table has expected columns and metrics
  expect_true("allele_divergence" %in% colnames(result$summary_table))
  expect_true("locus_heterozygosity" %in% colnames(result$summary_table))
  expect_equal(nrow(result$summary_table), 2) # 2 samples

  unlink(temp_dir, recursive = TRUE)
  unlink(temp_targets)
})
