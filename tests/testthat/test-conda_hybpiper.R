test_that("check_conda returns logical", {
  result <- check_conda(quiet = TRUE)
  expect_type(result, "logical")
  expect_length(result, 1)
})


test_that("check_hybpiper_conda returns logical", {
  result <- check_hybpiper_conda(env_name = "base", quiet = TRUE)
  expect_type(result, "logical")
  expect_length(result, 1)
})


test_that("run_conda validates inputs", {
  expect_error(run_conda(command = ""), "non-empty")
  expect_error(run_conda(command = 1), "non-empty")
  expect_error(run_conda(command = "echo", args = 1), "character vector")
  expect_error(
    run_conda(command = "echo", env_name = ""),
    "non-empty character"
  )
})


test_that("hybpiper_assemble validates paths", {
  temp_dir <- tempfile()
  dir.create(temp_dir)

  fake_file <- file.path(temp_dir, "missing.fastq")

  expect_error(
    hybpiper_assemble(
      sample_name = "sample1",
      f_read = fake_file,
      r_read = fake_file,
      targets_file = fake_file,
      wd = temp_dir
    ),
    "Forward read file not found"
  )

  unlink(temp_dir, recursive = TRUE)
})


test_that("hybpiper_stats validates paths", {
  fake_file <- tempfile(fileext = ".txt")

  expect_error(
    hybpiper_stats(
      targets_file = fake_file,
      namelist = fake_file
    ),
    "Targets file not found"
  )
})


test_that("hybpiper_stats falls back to native supercontig stats", {
  temp_dir <- tempfile()
  dir.create(temp_dir)

  targets_file <- file.path(temp_dir, "targets.fasta")
  writeLines(
    c(
      ">geneA_OG0001_species1",
      "ATGCATGC",
      ">geneA_OG0001_species2",
      "ATGCATGCATGC",
      ">geneB_OG0002_species1",
      "ATGCATGCAT"
    ),
    targets_file
  )

  namelist <- file.path(temp_dir, "namelist.txt")
  writeLines(c("sample1", "sample2"), namelist)

  sample1_dir <- file.path(
    temp_dir,
    "sample1",
    "geneA_OG0001",
    "sample1",
    "sequences",
    "intron"
  )
  dir.create(sample1_dir, recursive = TRUE)
  writeLines(
    c(
      ">geneA_OG0001",
      "ATGCNNNNAT"
    ),
    file.path(sample1_dir, "geneA_OG0001_supercontig.fasta")
  )
  writeLines(
    "geneA_OG0001,x,Stitched contig produced",
    file.path(temp_dir, "sample1", "sample1_genes_with_stitched_contig.csv")
  )

  sample2_dir <- file.path(
    temp_dir,
    "sample2",
    "geneB_OG0002",
    "sample2",
    "sequences",
    "intron"
  )
  dir.create(sample2_dir, recursive = TRUE)
  writeLines(
    c(
      ">geneB_OG0002",
      "ATGCATGCATGC"
    ),
    file.path(sample2_dir, "geneB_OG0002_supercontig.fasta")
  )

  local_mocked_bindings(
    check_hybpiper_conda = function(env_name = "hybpiper_env", quiet = FALSE) {
      TRUE
    },
    run_conda = function(...) 1,
    .package = "rhybphaser"
  )

  expect_warning(
    result <- hybpiper_stats(
      targets_file = targets_file,
      namelist = namelist,
      wd = temp_dir,
      mode = "supercontig"
    ),
    "Used native R supercontig stats fallback"
  )

  expect_true(file.exists(result$seq_lengths))
  expect_true(file.exists(result$stats))

  seq_lengths <- read.delim(result$seq_lengths, sep = "\t", check.names = FALSE)
  stats <- read.delim(result$stats, sep = "\t", check.names = FALSE)

  expect_identical(
    names(seq_lengths),
    c("Species", "geneA_OG0001", "geneB_OG0002")
  )
  expect_identical(seq_lengths$Species, c("MeanLength", "sample1", "sample2"))
  expect_identical(as.integer(seq_lengths$geneA_OG0001), c(10L, 6L, 0L))
  expect_identical(as.integer(seq_lengths$geneB_OG0002), c(10L, 0L, 12L))

  expect_identical(
    names(stats),
    c(
      "Name",
      "NumReads",
      "ReadsMapped",
      "PctOnTarget",
      "GenesMapped",
      "GenesWithContigs",
      "GenesWithSeqs",
      "GenesAt25pct",
      "GenesAt50pct",
      "GenesAt75pct",
      "GenesAt150pct",
      "ParalogWarningsLong",
      "ParalogWarningsDepth",
      "GenesWithoutStitchedContigs",
      "GenesWithStitchedContigs",
      "GenesWithStitchedContigsSkipped",
      "GenesWithChimeraWarning",
      "TotalBasesRecovered"
    )
  )

  expect_identical(stats$Name, c("sample1", "sample2"))
  expect_identical(stats$GenesWithSeqs, c(1L, 1L))
  expect_identical(stats$GenesMapped, c(1L, 1L))
  expect_identical(stats$GenesWithContigs, c(1L, 1L))
  expect_identical(stats$GenesWithStitchedContigs, c(1L, 1L))
  expect_identical(stats$TotalBasesRecovered, c(6L, 12L))

  unlink(temp_dir, recursive = TRUE)
})


test_that(".read_target_mean_lengths keeps gene names aligned with values", {
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  targets_file <- file.path(temp_dir, "targets.fasta")
  # geneZ is listed before geneA, and the two genes have different mean
  # lengths, so a name/value misalignment would be visible in the result.
  writeLines(
    c(
      ">geneZ_OG0002_species1",
      "ATGCATGCAT", # 10
      ">geneZ_OG0002_species2",
      "ATGCATGCATGCATGCAT", # 18 -> mean 14
      ">geneA_OG0001_species1",
      "ATGC" # 4 -> mean 4
    ),
    targets_file
  )

  result <- rhybphaser:::.read_target_mean_lengths(
    targets_file,
    dna = TRUE,
    known_genes = c("geneA_OG0001", "geneZ_OG0002")
  )

  expect_identical(names(result), c("geneA_OG0001", "geneZ_OG0002"))
  expect_identical(unname(result), c(4L, 14L))
})


test_that("run_hybpiper_test_dataset validates sample input", {
  temp_dir <- tempfile()
  dir.create(temp_dir)

  expect_error(
    run_hybpiper_test_dataset(
      output_dir = temp_dir,
      samples = "",
      keep_fastq = FALSE
    ),
    "must be a non-empty character vector"
  )

  unlink(temp_dir, recursive = TRUE)
})
