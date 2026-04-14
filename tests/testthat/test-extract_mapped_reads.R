test_that("extract_mapped_reads validates inputs", {
  expect_error(
    extract_mapped_reads(base_dir = "/nonexistent/path"),
    "base_dir does not exist"
  )

  temp_base <- tempfile()
  dir.create(temp_base)

  expect_error(
    extract_mapped_reads(base_dir = temp_base, namelist = "/missing/list.txt"),
    "namelist does not exist"
  )

  unlink(temp_base, recursive = TRUE)
})


test_that("extract_mapped_reads concatenates reads for HybPhaser layout", {
  base_dir <- tempfile()
  output_dir <- tempfile()

  reads_dir <- file.path(base_dir, "01_data", "sample1", "reads")
  dir.create(reads_dir, recursive = TRUE)

  writeLines(
    c(">read1", "AAAA", ">read2", "CCCC"),
    file.path(reads_dir, "gene1_unpaired.fasta")
  )
  writeLines(
    c(">read3", "GGGG", ">read4", "TTTT"),
    file.path(reads_dir, "gene1_combined.fasta")
  )

  result <- extract_mapped_reads(
    base_dir = base_dir,
    output_dir = output_dir
  )

  out_file <- file.path(output_dir, "sample1_mapped_reads.fasta")
  expect_true(file.exists(out_file))
  expect_equal(result$samples, "sample1")
  expect_equal(result$n_samples, 1)

  lines <- readLines(out_file)
  expect_true(any(lines == ">read1"))
  expect_true(any(lines == ">read4"))

  unlink(base_dir, recursive = TRUE)
  unlink(output_dir, recursive = TRUE)
})


test_that("extract_mapped_reads deduplicates by sequence", {
  base_dir <- tempfile()
  output_dir <- tempfile()

  reads_dir <- file.path(base_dir, "01_data", "sample1", "reads")
  dir.create(reads_dir, recursive = TRUE)

  # Two reads share sequence AAAA but have different names.
  writeLines(
    c(">read1", "AAAA", ">read2", "CCCC"),
    file.path(reads_dir, "gene1_unpaired.fasta")
  )
  writeLines(
    c(">read3", "AAAA", ">read4", "GGGG"),
    file.path(reads_dir, "gene1_combined.fasta")
  )

  extract_mapped_reads(
    base_dir = base_dir,
    output_dir = output_dir,
    remove_duplicate_sequences = TRUE
  )

  dedup_file <- file.path(output_dir, "sample1_mapped_reads_dedup-seq.fasta")
  raw_file <- file.path(output_dir, "sample1_mapped_reads.fasta")

  expect_true(file.exists(dedup_file))
  expect_false(file.exists(raw_file))

  dedup <- seqinr::read.fasta(
    file = dedup_file,
    as.string = TRUE,
    set.attributes = FALSE,
    seqtype = "DNA"
  )
  seqs <- toupper(vapply(dedup, function(x) x[[1]], character(1)))
  expect_equal(length(unique(seqs)), length(seqs))
  expect_equal(sort(unique(seqs)), sort(c("AAAA", "CCCC", "GGGG")))

  unlink(base_dir, recursive = TRUE)
  unlink(output_dir, recursive = TRUE)
})


test_that("extract_mapped_reads respects namelist", {
  base_dir <- tempfile()
  output_dir <- tempfile()

  reads_dir1 <- file.path(base_dir, "01_data", "sample1", "reads")
  reads_dir2 <- file.path(base_dir, "01_data", "sample2", "reads")
  dir.create(reads_dir1, recursive = TRUE)
  dir.create(reads_dir2, recursive = TRUE)

  writeLines(c(">r1", "AAAA"), file.path(reads_dir1, "gene_unpaired.fasta"))
  writeLines(c(">r2", "CCCC"), file.path(reads_dir2, "gene_unpaired.fasta"))

  namelist <- tempfile(fileext = ".txt")
  writeLines("sample2", namelist)

  result <- extract_mapped_reads(
    base_dir = base_dir,
    output_dir = output_dir,
    namelist = namelist
  )

  expect_equal(result$samples, "sample2")
  expect_true(file.exists(file.path(output_dir, "sample2_mapped_reads.fasta")))
  expect_false(file.exists(file.path(output_dir, "sample1_mapped_reads.fasta")))

  unlink(base_dir, recursive = TRUE)
  unlink(output_dir, recursive = TRUE)
  unlink(namelist)
})
