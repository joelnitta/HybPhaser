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
