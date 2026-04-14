test_that("collate_bbsplit_results validates required paths", {
  expect_error(
    collate_bbsplit_results(
      path_to_clade_association_folder = "/missing/folder",
      csv_file_with_clade_reference_names = "/missing/clade_refs.csv"
    ),
    "path_to_clade_association_folder not found"
  )
})


test_that("collate_bbsplit_results collates and normalizes stats", {
  tmp <- tempfile()
  dir.create(tmp)

  ca_dir <- file.path(tmp, "04_clade_association")
  stats_dir <- file.path(ca_dir, "bbsplit_stats")
  dir.create(stats_dir, recursive = TRUE)

  clade_refs <- file.path(tmp, "clade_refs.csv")
  utils::write.csv(
    data.frame(sample = c("EG30", "EG98"), abbreviation = c("C1", "C2")),
    clade_refs,
    row.names = FALSE
  )

  # BBSplit-style stats files (first two columns used)
  writeLines(
    c(
      "name\tunambiguousReads",
      "C1\t100",
      "C2\t10"
    ),
    file.path(stats_dir, "EG30_bbsplit-stats.txt")
  )
  writeLines(
    c(
      "name\tunambiguousReads",
      "C1\t20",
      "C2\t80"
    ),
    file.path(stats_dir, "EG98_bbsplit-stats.txt")
  )

  res <- collate_bbsplit_results(
    path_to_clade_association_folder = ca_dir,
    csv_file_with_clade_reference_names = clade_refs
  )

  expect_true(file.exists(file.path(ca_dir, "Table_clade_association.csv")))
  expect_true(file.exists(file.path(
    ca_dir,
    "Table_clade_association_normalised.csv"
  )))

  expect_true(all(
    c("sample", "C1", "C2") %in%
      colnames(res$table_clade_association)
  ))

  # Self-normalized entries should be 1.
  idx_eg30 <- which(res$table_clade_association_normalized$sample == "EG30")
  idx_eg98 <- which(res$table_clade_association_normalized$sample == "EG98")
  expect_equal(
    as.numeric(res$table_clade_association_normalized$C1[idx_eg30]),
    1
  )
  expect_equal(
    as.numeric(res$table_clade_association_normalized$C2[idx_eg98]),
    1
  )

  unlink(tmp, recursive = TRUE)
})


test_that("collate_bbsplit_results merges summary table when available", {
  tmp <- tempfile()
  dir.create(tmp)

  ca_dir <- file.path(tmp, "04_clade_association")
  stats_dir <- file.path(ca_dir, "bbsplit_stats")
  dir.create(stats_dir, recursive = TRUE)

  clade_refs <- file.path(tmp, "clade_refs.csv")
  utils::write.csv(
    data.frame(sample = c("EG30", "EG98"), abbreviation = c("C1", "C2")),
    clade_refs,
    row.names = FALSE
  )

  writeLines(
    c("name\tunambiguousReads", "C1\t100", "C2\t10"),
    file.path(stats_dir, "EG30_bbsplit-stats.txt")
  )
  writeLines(
    c("name\tunambiguousReads", "C1\t20", "C2\t80"),
    file.path(stats_dir, "EG98_bbsplit-stats.txt")
  )

  robj_dir <- file.path(tmp, "00_R_objects", "")
  dir.create(robj_dir, recursive = TRUE)
  summary_table <- data.frame(sample = c("EG30", "EG98"), SNPs = c(0.1, 0.2))
  saveRDS(summary_table, file.path(robj_dir, "Summary_table.Rds"))

  res <- collate_bbsplit_results(
    path_to_clade_association_folder = ca_dir,
    csv_file_with_clade_reference_names = clade_refs,
    path_to_output_folder = tmp,
    subset_name = ""
  )

  expect_false(is.null(res$table_with_summary))
  expect_true("SNPs" %in% colnames(res$table_with_summary))
  expect_true(file.exists(file.path(
    ca_dir,
    "Table_clade_association_and_summary_table.csv"
  )))

  unlink(tmp, recursive = TRUE)
})


test_that("collate_bbsplit_results_from_config reads config and runs", {
  tmp <- tempfile()
  dir.create(tmp)

  ca_dir <- file.path(tmp, "04_clade_association")
  stats_dir <- file.path(ca_dir, "bbsplit_stats")
  dir.create(stats_dir, recursive = TRUE)

  clade_refs <- file.path(tmp, "clade_refs.csv")
  utils::write.csv(
    data.frame(sample = c("EG30", "EG98"), abbreviation = c("C1", "C2")),
    clade_refs,
    row.names = FALSE
  )

  writeLines(
    c("name\tunambiguousReads", "C1\t100", "C2\t10"),
    file.path(stats_dir, "EG30_bbsplit-stats.txt")
  )
  writeLines(
    c("name\tunambiguousReads", "C1\t20", "C2\t80"),
    file.path(stats_dir, "EG98_bbsplit-stats.txt")
  )

  robj_dir <- file.path(tmp, "00_R_objects", "")
  dir.create(robj_dir, recursive = TRUE)
  saveRDS(
    data.frame(sample = c("EG30", "EG98"), SNPs = c(0.1, 0.2)),
    file.path(robj_dir, "Summary_table.Rds")
  )

  cfg <- file.path(tmp, "config.txt")
  writeLines(
    c(
      paste0("path_to_clade_association_folder = ", dQuote(ca_dir)),
      paste0("csv_file_with_clade_reference_names = ", dQuote(clade_refs)),
      paste0("path_to_output_folder = ", dQuote(tmp)),
      'name_for_dataset_optimization_subset = ""'
    ),
    cfg
  )

  res <- collate_bbsplit_results_from_config(cfg)
  expect_true(file.exists(res$output_files[["table"]]))

  unlink(tmp, recursive = TRUE)
})
