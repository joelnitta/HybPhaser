test_that("prepare_phasing_script validates paired-end IDs", {
  tmp <- tempfile()
  dir.create(tmp)

  refs <- file.path(tmp, "refs")
  reads <- file.path(tmp, "reads")
  out <- file.path(tmp, "out")
  dir.create(refs)
  dir.create(reads)

  writeLines(c(">ref1", "ATCG"), file.path(refs, "ref1_consensus.fasta"))
  prep <- file.path(tmp, "phasing_prep.csv")
  utils::write.csv(
    data.frame(sample = "S1", ref1 = "ref1", abb1 = "R1"),
    prep,
    row.names = FALSE
  )

  expect_error(
    prepare_phasing_script(
      path_to_phasing_folder = out,
      csv_file_with_phasing_prep_info = prep,
      path_to_read_files_phasing = reads,
      read_type_4phasing = "paired-end",
      reference_sequence_folder = refs
    ),
    "ID_read_pair1 and ID_read_pair2 are required"
  )

  unlink(tmp, recursive = TRUE)
})


test_that("prepare_phasing_script builds single-end commands", {
  tmp <- tempfile()
  dir.create(tmp)

  refs <- file.path(tmp, "refs")
  reads <- file.path(tmp, "reads")
  out <- file.path(tmp, "out")
  dir.create(refs)
  dir.create(reads)

  writeLines(c(">ref1", "ATCG"), file.path(refs, "ref1_consensus.fasta"))
  writeLines(c(">ref2", "ATCG"), file.path(refs, "ref2_consensus.fasta"))

  prep <- file.path(tmp, "phasing_prep.csv")
  utils::write.csv(
    data.frame(
      sample = c("S1", "S2"),
      ref1 = c("ref1", "ref2"),
      abb1 = c("R1", "R2")
    ),
    prep,
    row.names = FALSE
  )

  writeLines(c("@r1", "AAAA"), file.path(reads, "S1.fastq"))
  writeLines(c("@r2", "CCCC"), file.path(reads, "S2.fastq"))

  res <- prepare_phasing_script(
    path_to_phasing_folder = out,
    csv_file_with_phasing_prep_info = prep,
    path_to_read_files_phasing = reads,
    read_type_4phasing = "single-end",
    reference_sequence_folder = refs,
    no_of_threads_phasing = 2,
    java_memory_usage_phasing = "1G"
  )

  expect_true(file.exists(res$script_path))
  expect_equal(length(res$commands), 2)
  expect_true(all(grepl("ambiguous=all", res$commands, fixed = TRUE)))
  expect_true(all(grepl("threads=2", res$commands, fixed = TRUE)))
  expect_true(all(grepl("-Xmx1G", res$commands, fixed = TRUE)))
  expect_true(any(grepl("basename=.*S1_to_%.fastq", res$commands)))

  unlink(tmp, recursive = TRUE)
})


test_that("prepare_phasing_script builds paired-end commands", {
  tmp <- tempfile()
  dir.create(tmp)

  refs <- file.path(tmp, "refs")
  reads <- file.path(tmp, "reads")
  out <- file.path(tmp, "out")
  dir.create(refs)
  dir.create(reads)

  writeLines(c(">ref1", "ATCG"), file.path(refs, "ref1_consensus.fasta"))

  prep <- file.path(tmp, "phasing_prep.csv")
  utils::write.csv(
    data.frame(sample = "S1", ref1 = "ref1", abb1 = "R1"),
    prep,
    row.names = FALSE
  )

  writeLines(c("@r1", "AAAA"), file.path(reads, "S1_R1.fastq"))
  writeLines(c("@r1", "TTTT"), file.path(reads, "S1_R2.fastq"))

  res <- prepare_phasing_script(
    path_to_phasing_folder = out,
    csv_file_with_phasing_prep_info = prep,
    path_to_read_files_phasing = reads,
    read_type_4phasing = "paired-end",
    ID_read_pair1 = "_R1.fastq",
    ID_read_pair2 = "_R2.fastq",
    reference_sequence_folder = refs
  )

  expect_equal(length(res$commands), 1)
  expect_true(grepl("in=.*S1_R1.fastq", res$commands[[1]]))
  expect_true(grepl("in2=.*S1_R2.fastq", res$commands[[1]]))

  unlink(tmp, recursive = TRUE)
})


test_that("prepare_phasing_script can execute generated script", {
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

  prep <- file.path(tmp, "phasing_prep.csv")
  utils::write.csv(
    data.frame(sample = "S1", ref1 = "ref1", abb1 = "R1"),
    prep,
    row.names = FALSE
  )

  writeLines(c("@r1", "AAAA"), file.path(reads, "S1.fastq"))

  fake_bbsplit <- file.path(bbmap, "bbsplit.sh")
  writeLines(c("#!/bin/bash", "exit 0"), fake_bbsplit)
  Sys.chmod(fake_bbsplit, mode = "0755")

  res <- prepare_phasing_script(
    path_to_phasing_folder = out,
    csv_file_with_phasing_prep_info = prep,
    path_to_read_files_phasing = reads,
    read_type_4phasing = "single-end",
    reference_sequence_folder = refs,
    path_to_bbmap_executables = bbmap,
    run_bash_script_in_R = TRUE
  )

  expect_equal(res$run_status, 0)
  expect_true(file.exists(res$script_path))

  unlink(tmp, recursive = TRUE)
})


test_that("prepare_phasing_script_from_config reads config values", {
  tmp <- tempfile()
  dir.create(tmp)

  refs <- file.path(tmp, "refs")
  reads <- file.path(tmp, "reads")
  out <- file.path(tmp, "out")
  dir.create(refs)
  dir.create(reads)

  writeLines(c(">ref1", "ATCG"), file.path(refs, "ref1_consensus.fasta"))
  prep <- file.path(tmp, "phasing_prep.csv")
  utils::write.csv(
    data.frame(sample = "S1", ref1 = "ref1", abb1 = "R1"),
    prep,
    row.names = FALSE
  )
  writeLines(c("@r1", "AAAA"), file.path(reads, "S1.fastq"))

  cfg <- file.path(tmp, "config.txt")
  writeLines(
    c(
      paste0("path_to_phasing_folder = ", dQuote(out)),
      paste0("csv_file_with_phasing_prep_info = ", dQuote(prep)),
      paste0("path_to_read_files_phasing = ", dQuote(reads)),
      'read_type_4phasing = "single-end"',
      'ID_read_pair1 = ""',
      'ID_read_pair2 = ""',
      paste0("reference_sequence_folder = ", dQuote(refs)),
      'folder_for_phased_reads = ""',
      'folder_for_phasing_stats = ""',
      'path_to_bbmap_executables = ""',
      'no_of_threads_phasing = "auto"',
      'run_bash_script_in_R = "no"',
      'java_memory_usage_phasing = ""'
    ),
    cfg
  )

  res <- prepare_phasing_script_from_config(cfg)
  expect_true(file.exists(res$script_path))
  expect_equal(length(res$commands), 1)

  unlink(tmp, recursive = TRUE)
})
