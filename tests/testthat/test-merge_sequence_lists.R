test_that("merge_sequence_lists validates required paths", {
  expect_error(
    merge_sequence_lists(
      path_to_sequence_lists_normal = "/missing/normal",
      path_to_sequence_lists_phased = "/missing/phased",
      path_of_sequence_lists_output = tempfile("out_"),
      path_to_namelist_normal = "/missing/normal.txt",
      path_to_namelist_phased = "/missing/phased.txt"
    ),
    "path_to_sequence_lists_normal not found"
  )
})


test_that("merge_sequence_lists merges phased and non-phased sequence lists", {
  tmp <- tempfile()
  dir.create(tmp)

  normal_dir <- file.path(tmp, "normal")
  phased_dir <- file.path(tmp, "phased")
  out_dir <- file.path(tmp, "merged")
  dir.create(normal_dir)
  dir.create(phased_dir)

  normal_namelist <- file.path(tmp, "namelist_normal.txt")
  phased_namelist <- file.path(tmp, "namelist_phased.txt")

  writeLines(c("S1", "S2"), normal_namelist)
  writeLines(c("S2_to_A", "S2_to_B"), phased_namelist)

  writeLines(
    c(
      ">S1",
      "AAAA",
      ">S2",
      "CCCC"
    ),
    file.path(normal_dir, "locus1_consensus.fasta")
  )
  writeLines(
    c(
      ">S1",
      "TTTT"
    ),
    file.path(normal_dir, "locus2_consensus.fasta")
  )

  writeLines(
    c(
      ">S2_to_A",
      "GGGG",
      ">S2_to_B",
      "GGGA"
    ),
    file.path(phased_dir, "locus1_consensus.fasta")
  )
  writeLines(
    c(
      ">S2_to_A",
      "TTTA"
    ),
    file.path(phased_dir, "locus3_consensus.fasta")
  )

  res <- merge_sequence_lists(
    path_to_sequence_lists_normal = normal_dir,
    path_to_sequence_lists_phased = phased_dir,
    path_of_sequence_lists_output = out_dir,
    path_to_namelist_normal = normal_namelist,
    path_to_namelist_phased = phased_namelist,
    exchange_phased_with_not_phased_samples = "yes",
    include_phased_seqlists_when_non_phased_locus_absent = "yes"
  )

  expect_true(file.exists(file.path(out_dir, "locus1_consensus.fasta")))
  expect_true(file.exists(file.path(out_dir, "locus2_consensus.fasta")))
  expect_true(file.exists(file.path(out_dir, "locus3_consensus.fasta")))
  expect_true("S2" %in% res$samples_excluded)

  locus1 <- readLines(file.path(out_dir, "locus1_consensus.fasta"))
  expect_false(any(locus1 == ">S2"))
  expect_true(any(locus1 == ">S2_to_A"))
  expect_true(any(locus1 == ">S2_to_B"))

  locus3 <- readLines(file.path(out_dir, "locus3_consensus.fasta"))
  expect_true(any(locus3 == ">S2_to_A"))

  unlink(tmp, recursive = TRUE)
})


test_that("merge_sequence_lists applies locus include and sample exclude files", {
  tmp <- tempfile()
  dir.create(tmp)

  normal_dir <- file.path(tmp, "normal")
  phased_dir <- file.path(tmp, "phased")
  out_dir <- file.path(tmp, "merged")
  dir.create(normal_dir)
  dir.create(phased_dir)

  normal_namelist <- file.path(tmp, "namelist_normal.txt")
  phased_namelist <- file.path(tmp, "namelist_phased.txt")
  loci_include <- file.path(tmp, "loci_include.txt")
  samples_exclude <- file.path(tmp, "samples_exclude.txt")

  writeLines(c("S1", "S2"), normal_namelist)
  writeLines(c("S2_to_A"), phased_namelist)
  writeLines("locus2", loci_include)
  writeLines("S1", samples_exclude)

  writeLines(
    c(
      ">S1",
      "AAAA",
      ">S2",
      "CCCC"
    ),
    file.path(normal_dir, "locus1_consensus.fasta")
  )
  writeLines(
    c(
      ">S1",
      "TTTT",
      ">S2",
      "GGGG"
    ),
    file.path(normal_dir, "locus2_consensus.fasta")
  )
  writeLines(
    c(
      ">S2_to_A",
      "GGGA"
    ),
    file.path(phased_dir, "locus2_consensus.fasta")
  )

  res <- merge_sequence_lists(
    path_to_sequence_lists_normal = normal_dir,
    path_to_sequence_lists_phased = phased_dir,
    path_of_sequence_lists_output = out_dir,
    path_to_namelist_normal = normal_namelist,
    path_to_namelist_phased = phased_namelist,
    file_with_loci_included = loci_include,
    file_with_samples_excluded = samples_exclude,
    exchange_phased_with_not_phased_samples = "no",
    include_phased_seqlists_when_non_phased_locus_absent = "no"
  )

  expect_false(file.exists(file.path(out_dir, "locus1_consensus.fasta")))
  expect_true(file.exists(file.path(out_dir, "locus2_consensus.fasta")))
  expect_false("S1" %in% res$samples_included)

  locus2 <- readLines(file.path(out_dir, "locus2_consensus.fasta"))
  expect_false(any(locus2 == ">S1"))

  unlink(tmp, recursive = TRUE)
})


test_that("merge_sequence_lists_from_config reads config and runs", {
  tmp <- tempfile()
  dir.create(tmp)

  normal_dir <- file.path(tmp, "normal")
  phased_dir <- file.path(tmp, "phased")
  out_dir <- file.path(tmp, "merged")
  dir.create(normal_dir)
  dir.create(phased_dir)

  normal_namelist <- file.path(tmp, "namelist_normal.txt")
  phased_namelist <- file.path(tmp, "namelist_phased.txt")

  writeLines("S1", normal_namelist)
  writeLines("S1_to_A", phased_namelist)

  writeLines(c(">S1", "AAAA"), file.path(normal_dir, "locus1_consensus.fasta"))
  writeLines(
    c(">S1_to_A", "AAAT"),
    file.path(phased_dir, "locus1_consensus.fasta")
  )

  cfg <- file.path(tmp, "config.txt")
  writeLines(
    c(
      paste0("path_to_sequence_lists_normal = ", dQuote(normal_dir)),
      paste0("path_to_sequence_lists_phased = ", dQuote(phased_dir)),
      paste0("path_of_sequence_lists_output = ", dQuote(out_dir)),
      paste0("path_to_namelist_normal = ", dQuote(normal_namelist)),
      paste0("path_to_namelist_phased = ", dQuote(phased_namelist)),
      'file_with_samples_included = ""',
      'file_with_samples_excluded = ""',
      'file_with_loci_excluded = ""',
      'file_with_loci_included = ""',
      'exchange_phased_with_not_phased_samples = "yes"',
      'include_phased_seqlists_when_non_phased_locus_absent = "no"'
    ),
    cfg
  )

  res <- merge_sequence_lists_from_config(cfg)
  expect_true(file.exists(file.path(res$output_dir, "locus1_consensus.fasta")))

  unlink(tmp, recursive = TRUE)
})
