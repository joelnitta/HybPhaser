test_that("collate_phasing_stats validates required paths", {
  expect_error(
    collate_phasing_stats(
      path_to_phasing_folder = "/missing/phasing",
      csv_file_with_phasing_prep_info = "/missing/phasing_prep.csv"
    ),
    "path_to_phasing_folder not found"
  )
})


test_that("collate_phasing_stats creates table from stats files", {
  tmp <- tempfile()
  dir.create(tmp)

  phasing_dir <- file.path(tmp, "05_phasing")
  stats_dir <- file.path(phasing_dir, "phasing_stats")
  dir.create(stats_dir, recursive = TRUE)

  phasing_prep <- file.path(phasing_dir, "phasing_prep.csv")
  utils::write.csv(
    data.frame(
      sample = c("EG30", "EG98"),
      ref1 = c("R1", "R1"),
      abb1 = c("A", "A"),
      ref2 = c("R2", "R2"),
      abb2 = c("B", "B")
    ),
    phasing_prep,
    row.names = FALSE
  )

  writeLines(
    c(
      "name\tpercent",
      "A\t50",
      "B\t50"
    ),
    file.path(stats_dir, "EG30_phasing-stats.txt")
  )
  writeLines(
    c(
      "name\tpercent",
      "A\t75",
      "B\t25"
    ),
    file.path(stats_dir, "EG98_phasing-stats.txt")
  )

  res <- collate_phasing_stats(
    path_to_phasing_folder = phasing_dir,
    csv_file_with_phasing_prep_info = phasing_prep
  )

  expect_true(file.exists(file.path(phasing_dir, "Table_phasing_stats.csv")))
  expect_true(all(c("sample", "A", "B") %in% colnames(res$table_phasing_stats)))

  idx_eg30 <- which(res$table_phasing_stats$sample == "EG30")
  idx_eg98 <- which(res$table_phasing_stats$sample == "EG98")
  expect_equal(as.numeric(res$table_phasing_stats$A[idx_eg30]), 0.5)
  expect_equal(as.numeric(res$table_phasing_stats$B[idx_eg30]), 0.5)
  expect_equal(as.numeric(res$table_phasing_stats$A[idx_eg98]), 0.75)
  expect_equal(as.numeric(res$table_phasing_stats$B[idx_eg98]), 0.25)

  unlink(tmp, recursive = TRUE)
})


test_that("collate_phasing_stats_from_config reads config and runs", {
  tmp <- tempfile()
  dir.create(tmp)

  phasing_dir <- file.path(tmp, "05_phasing")
  stats_dir <- file.path(phasing_dir, "phasing_stats")
  dir.create(stats_dir, recursive = TRUE)

  phasing_prep <- file.path(phasing_dir, "phasing_prep.csv")
  utils::write.csv(
    data.frame(
      samples = c("EG30"),
      ref1 = c("R1"),
      abb1 = c("A")
    ),
    phasing_prep,
    row.names = FALSE
  )

  writeLines(
    c(
      "name\tpercent",
      "A\t100"
    ),
    file.path(stats_dir, "EG30_phasing-stats.txt")
  )

  cfg <- file.path(tmp, "config.txt")
  writeLines(
    c(
      paste0("path_to_phasing_folder = ", dQuote(phasing_dir)),
      paste0("csv_file_with_phasing_prep_info = ", dQuote(phasing_prep)),
      paste0("folder_for_phasing_stats = ", dQuote(stats_dir))
    ),
    cfg
  )

  res <- collate_phasing_stats_from_config(cfg)
  expect_true(file.exists(res$output_file))

  unlink(tmp, recursive = TRUE)
})
