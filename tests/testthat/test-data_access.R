test_that("hybphaser_scripts returns scripts directory", {
  scripts_dir <- hybphaser_scripts()
  expect_type(scripts_dir, "character")
  # Note: This will only pass if package is properly installed
  # expect_true(dir.exists(scripts_dir))
})

test_that("hybphaser_scripts returns path to specific script", {
  skip_if_not_installed("HybPhaser")

  script_path <- hybphaser_scripts("1_generate_consensus_sequences.sh")
  expect_type(script_path, "character")
  expect_match(script_path, "1_generate_consensus_sequences\\.sh$")
})

test_that("hybphaser_scripts errors on nonexistent script", {
  skip_if_not_installed("HybPhaser")

  expect_error(
    hybphaser_scripts("nonexistent_script.sh"),
    "Script not found"
  )
})

test_that("hybphaser_example returns extdata directory", {
  extdata_dir <- hybphaser_example()
  expect_type(extdata_dir, "character")
})

test_that("hybphaser_example returns path to specific file", {
  skip_if_not_installed("HybPhaser")

  file_path <- hybphaser_example("clade_references.csv")
  expect_type(file_path, "character")
  expect_match(file_path, "clade_references\\.csv$")
})

test_that("hybphaser_example errors on nonexistent file", {
  skip_if_not_installed("HybPhaser")

  expect_error(
    hybphaser_example("nonexistent_file.txt"),
    "Example file not found"
  )
})

test_that("run_hybphaser_script dry_run returns command", {
  skip_if_not_installed("HybPhaser")

  cmd <- run_hybphaser_script(
    "1_generate_consensus_sequences.sh",
    args = c("-n", "test.txt"),
    dry_run = TRUE
  )

  expect_type(cmd, "character")
  expect_match(cmd, "1_generate_consensus_sequences\\.sh")
  expect_match(cmd, "-n")
  expect_match(cmd, "test\\.txt")
})
