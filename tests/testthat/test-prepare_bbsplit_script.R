test_that("prepare_bbsplit_script validates inputs", {
  expect_error(
    prepare_bbsplit_script(
      path_to_clade_association_folder = tempdir(),
      csv_file_with_clade_reference_names = "/missing/ref.csv",
      path_to_reference_sequences = tempdir(),
      path_to_read_files_cladeassociation = tempdir()
    ),
    "csv_file_with_clade_reference_names not found"
  )

  tmp <- tempfile()
  dir.create(tmp)
  refs <- file.path(tmp, "refs")
  reads <- file.path(tmp, "reads")
  out <- file.path(tmp, "out")
  dir.create(refs)
  dir.create(reads)

  csv <- file.path(tmp, "clade_refs.csv")
  write.csv(data.frame(samples = "s1", abb = "A"), csv, row.names = FALSE)
  writeLines(c(">s1", "ATCG"), file.path(refs, "s1_consensus.fasta"))

  expect_error(
    prepare_bbsplit_script(
      path_to_clade_association_folder = out,
      csv_file_with_clade_reference_names = csv,
      path_to_reference_sequences = refs,
      path_to_read_files_cladeassociation = reads,
      read_type_cladeassociation = "paired-end"
    ),
    "ID_read_pair1 and ID_read_pair2 are required"
  )

  unlink(tmp, recursive = TRUE)
})


test_that("prepare_bbsplit_script creates single-end commands", {
  tmp <- tempfile()
  dir.create(tmp)

  refs <- file.path(tmp, "refs")
  reads <- file.path(tmp, "reads")
  out <- file.path(tmp, "out")
  dir.create(refs)
  dir.create(reads)

  writeLines(c(">ref1", "ATCG"), file.path(refs, "ref1_consensus.fasta"))
  writeLines(c(">ref2", "ATCG"), file.path(refs, "ref2_consensus.fasta"))

  clade_refs <- file.path(tmp, "clade_refs.csv")
  utils::write.csv(
    data.frame(samples = c("ref1", "ref2"), abb = c("R1", "R2")),
    clade_refs,
    row.names = FALSE
  )

  writeLines(c("@r1", "AAAA"), file.path(reads, "sample1.fastq"))
  writeLines(c("@r2", "CCCC"), file.path(reads, "sample2.fastq"))

  result <- prepare_bbsplit_script(
    path_to_clade_association_folder = out,
    csv_file_with_clade_reference_names = clade_refs,
    path_to_reference_sequences = refs,
    path_to_read_files_cladeassociation = reads,
    read_type_cladeassociation = "single-end",
    no_of_threads = 2,
    java_memory_usage_clade_association = "1G"
  )

  expect_true(file.exists(result$script_path))
  expect_equal(length(result$commands), 2)
  expect_true(any(grepl("ref_R1=", result$commands, fixed = TRUE)))
  expect_true(all(grepl("threads=2", result$commands, fixed = TRUE)))
  expect_true(all(grepl("-Xmx1G", result$commands, fixed = TRUE)))

  script_lines <- readLines(result$script_path)
  expect_equal(script_lines[[1]], "#!/bin/bash")
  expect_true(any(grepl("sample1.fastq_bbsplit-stats.txt", script_lines)))

  unlink(tmp, recursive = TRUE)
})


test_that("prepare_bbsplit_script filters paired-end reads by namelist", {
  tmp <- tempfile()
  dir.create(tmp)

  refs <- file.path(tmp, "refs")
  reads <- file.path(tmp, "reads")
  out <- file.path(tmp, "out")
  dir.create(refs)
  dir.create(reads)

  writeLines(c(">ref1", "ATCG"), file.path(refs, "ref1_consensus.fasta"))
  clade_refs <- file.path(tmp, "clade_refs.csv")
  utils::write.csv(
    data.frame(samples = "ref1", abb = "R1"),
    clade_refs,
    row.names = FALSE
  )

  writeLines(c("@r1", "AAAA"), file.path(reads, "s1_R1.fastq"))
  writeLines(c("@r1", "TTTT"), file.path(reads, "s1_R2.fastq"))
  writeLines(c("@r2", "CCCC"), file.path(reads, "s2_R1.fastq"))
  writeLines(c("@r2", "GGGG"), file.path(reads, "s2_R2.fastq"))

  include <- file.path(tmp, "samples.txt")
  writeLines("s2", include)

  result <- prepare_bbsplit_script(
    path_to_clade_association_folder = out,
    csv_file_with_clade_reference_names = clade_refs,
    path_to_reference_sequences = refs,
    path_to_read_files_cladeassociation = reads,
    read_type_cladeassociation = "paired-end",
    ID_read_pair1 = "_R1.fastq",
    ID_read_pair2 = "_R2.fastq",
    file_with_samples_included = include
  )

  expect_equal(length(result$commands), 1)
  expect_true(grepl("s2_R1.fastq", result$commands[[1]], fixed = TRUE))
  expect_true(grepl("s2_R2.fastq", result$commands[[1]], fixed = TRUE))
  expect_false(grepl("s1_R1.fastq", result$commands[[1]], fixed = TRUE))

  unlink(tmp, recursive = TRUE)
})


test_that("prepare_bbsplit_script can execute generated script", {
  tmp <- tempfile()
  dir.create(tmp)

  refs <- file.path(tmp, "refs")
  reads <- file.path(tmp, "reads")
  bbmap <- file.path(tmp, "bbmap")
  out <- file.path(tmp, "out")
  dir.create(refs)
  dir.create(reads)
  dir.create(bbmap)

  writeLines(c(">ref1", "ATCG"), file.path(refs, "ref1_consensus.fasta"))
  clade_refs <- file.path(tmp, "clade_refs.csv")
  utils::write.csv(
    data.frame(samples = "ref1", abb = "R1"),
    clade_refs,
    row.names = FALSE
  )

  writeLines(c("@r1", "AAAA"), file.path(reads, "sample1.fastq"))

  fake_bbsplit <- file.path(bbmap, "bbsplit.sh")
  writeLines(c("#!/bin/bash", "exit 0"), fake_bbsplit)
  Sys.chmod(fake_bbsplit, mode = "0755")

  result <- prepare_bbsplit_script(
    path_to_clade_association_folder = out,
    csv_file_with_clade_reference_names = clade_refs,
    path_to_reference_sequences = refs,
    path_to_read_files_cladeassociation = reads,
    read_type_cladeassociation = "single-end",
    path_to_bbmap = bbmap,
    run_clade_association_mapping_in_R = TRUE
  )

  expect_equal(result$run_status, 0)
  expect_true(file.exists(result$script_path))

  unlink(tmp, recursive = TRUE)
})


test_that("prepare_bbsplit_script_from_config reads config values", {
  tmp <- tempfile()
  dir.create(tmp)

  refs <- file.path(tmp, "refs")
  reads <- file.path(tmp, "reads")
  out <- file.path(tmp, "out")
  dir.create(refs)
  dir.create(reads)

  writeLines(c(">ref1", "ATCG"), file.path(refs, "ref1_consensus.fasta"))
  clade_refs <- file.path(tmp, "clade_refs.csv")
  utils::write.csv(
    data.frame(samples = "ref1", abb = "R1"),
    clade_refs,
    row.names = FALSE
  )
  writeLines(c("@r1", "AAAA"), file.path(reads, "sample1.fastq"))

  cfg <- file.path(tmp, "config.txt")
  writeLines(
    c(
      paste0("path_to_clade_association_folder = ", dQuote(out)),
      paste0("csv_file_with_clade_reference_names = ", dQuote(clade_refs)),
      paste0("path_to_reference_sequences = ", dQuote(refs)),
      paste0("path_to_read_files_cladeassociation = ", dQuote(reads)),
      'read_type_cladeassociation = "single-end"',
      'ID_read_pair1 = ""',
      'ID_read_pair2 = ""',
      'file_with_samples_included = ""',
      'path_to_bbmap = ""',
      'no_of_threads_clade_association = "auto"',
      'run_clade_association_mapping_in_R = "no"',
      'java_memory_usage_clade_association = ""'
    ),
    cfg
  )

  result <- prepare_bbsplit_script_from_config(cfg)
  expect_true(file.exists(result$script_path))
  expect_equal(length(result$commands), 1)

  unlink(tmp, recursive = TRUE)
})
