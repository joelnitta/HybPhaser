test_that("read_config reads valid configuration file", {
  # Create a temporary config file
  config_file <- tempfile(fileext = ".txt")
  writeLines(
    c(
      'path_to_output_folder = "/tmp/output"',
      'fasta_file_with_targets = "/tmp/targets.fasta"',
      'targets_file_format = "DNA"'
    ),
    config_file
  )

  config <- read_config(config_file)

  expect_type(config, "list")
  expect_equal(config$path_to_output_folder, "/tmp/output")
  expect_equal(config$targets_file_format, "DNA")

  # Clean up  unlink(config_file)
})

test_that("read_config fails when file doesn't exist", {
  expect_error(
    read_config("nonexistent_config.txt"),
    "Configuration file not found"
  )
})

test_that("read_config validates required variables", {
  config_file <- tempfile(fileext = ".txt")
  writeLines('path_to_output_folder = "/tmp/output"', config_file)

  expect_error(
    read_config(
      config_file,
      required_vars = c("path_to_output_folder", "missing_var")
    ),
    "Missing required config variables: missing_var"
  )

  unlink(config_file)
})

test_that("validate_paths checks file existence", {
  # Create temp file
  temp_file <- tempfile()
  writeLines("test", temp_file)

  expect_true(
    validate_paths(
      list("test file" = temp_file),
      must_exist = TRUE,
      type = "file"
    )
  )

  expect_error(
    validate_paths(
      list("missing file" = "/nonexistent/path"),
      must_exist = TRUE
    ),
    "missing file not found"
  )

  unlink(temp_file)
})

test_that("validate_paths distinguishes files and directories", {
  temp_dir <- tempfile()
  dir.create(temp_dir)
  temp_file <- tempfile()
  writeLines("test", temp_file)

  expect_error(
    validate_paths(
      list("should be dir" = temp_file),
      must_exist = TRUE,
      type = "dir"
    ),
    "must be a directory"
  )

  expect_error(
    validate_paths(
      list("should be file" = temp_dir),
      must_exist = TRUE,
      type = "file"
    ),
    "must be a file"
  )

  unlink(temp_dir, recursive = TRUE)
  unlink(temp_file)
})

test_that("create_output_dirs creates directory structure", {
  base_dir <- tempfile()

  dirs <- create_output_dirs(
    base_dir,
    subdirs = c("dir1", "dir2"),
    subset_name = "test"
  )

  expect_true(dir.exists(file.path(base_dir, "dir1_test")))
  expect_true(dir.exists(file.path(base_dir, "dir2_test")))

  unlink(base_dir, recursive = TRUE)
})

test_that("create_output_dirs works without subset name", {
  base_dir <- tempfile()

  dirs <- create_output_dirs(
    base_dir,
    subdirs = c("dir1", "dir2"),
    subset_name = ""
  )

  expect_true(dir.exists(file.path(base_dir, "dir1")))
  expect_true(dir.exists(file.path(base_dir, "dir2")))
  expect_false(dir.exists(file.path(base_dir, "dir1_")))

  unlink(base_dir, recursive = TRUE)
})
