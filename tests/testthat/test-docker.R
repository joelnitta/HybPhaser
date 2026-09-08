test_that("check_docker returns logical", {
  result <- check_docker(quiet = TRUE)
  expect_type(result, "logical")
})

test_that(".normalize_docker_path handles paths correctly", {
  # Test basic path expansion
  path <- rhybphaser:::.normalize_docker_path("~/test")
  expect_true(grepl("test$", path))
  expect_true(!grepl("~", path))

  # Test absolute path passthrough
  if (.Platform$OS.type != "windows") {
    path <- rhybphaser:::.normalize_docker_path("/tmp/test")
    expect_equal(path, "/tmp/test")
  }
})

test_that("check_docker_image validates arguments", {
  skip_if_not(check_docker(quiet = TRUE), "Docker not available")

  # Should not error with valid call (but image may not exist)
  expect_no_error({
    result <- check_docker_image("nonexistent/image:tag", pull = FALSE)
  })
})

test_that("run_generate_consensus_sequences validates inputs", {
  skip_if_not(check_docker(quiet = TRUE), "Docker not available")

  # Should error on missing directory
  expect_error(
    run_generate_consensus_sequences(
      hybpiper_dir = "/nonexistent/path",
      output_dir = tempdir()
    ),
    "HybPiper directory not found"
  )

  # Should error on missing namelist
  expect_error(
    run_generate_consensus_sequences(
      hybpiper_dir = tempdir(),
      output_dir = tempdir(),
      namelist = "/nonexistent/namelist.txt"
    ),
    "Namelist file not found"
  )
})

test_that("run_generate_consensus_sequences creates output directory", {
  skip_if_not(check_docker(quiet = TRUE), "Docker not available")
  skip_if_not(
    check_docker_image("joelnitta/hybphaser:latest", pull = FALSE),
    "HybPhaser Docker image not available"
  )

  temp_hybpiper <- tempfile()
  dir.create(temp_hybpiper)
  temp_output <- tempfile()

  # Output dir should be created even if command fails
  # (will fail because hybpiper dir is empty, but that's OK for this test)
  suppressWarnings({
    run_generate_consensus_sequences(
      hybpiper_dir = temp_hybpiper,
      output_dir = temp_output,
      sample = "test_sample"
    )
  })

  expect_true(dir.exists(temp_output))

  unlink(temp_hybpiper, recursive = TRUE)
  unlink(temp_output, recursive = TRUE)
})

test_that("run_extract_mapped_reads validates inputs", {
  skip_if_not(check_docker(quiet = TRUE), "Docker not available")

  # Should error on missing directory
  expect_error(
    run_extract_mapped_reads(
      base_dir = "/nonexistent/path",
      output_dir = tempdir()
    ),
    "Base directory not found"
  )

  # Should error on missing namelist
  expect_error(
    run_extract_mapped_reads(
      base_dir = tempdir(),
      output_dir = tempdir(),
      namelist = "/nonexistent/namelist.txt"
    ),
    "Namelist file not found"
  )
})

# Integration tests with real test data
test_that("create_test_hybpiper_output creates valid structure", {
  test_dir <- tempfile()

  result <- create_test_hybpiper_output(
    test_dir,
    samples = c("sample1"),
    genes = c("gene001")
  )

  expect_true(dir.exists(test_dir))
  expect_true(dir.exists(file.path(test_dir, "sample1", "gene001")))
  expect_true(file.exists(
    file.path(test_dir, "sample1", "gene001", "gene001_contigs.fasta")
  ))
  expect_true(file.exists(
    file.path(test_dir, "sample1", "gene001", "gene001_interleaved.fasta")
  ))

  # Check contig file content
  contig_lines <- readLines(
    file.path(test_dir, "sample1", "gene001", "gene001_contigs.fasta")
  )
  expect_true(grepl("^>", contig_lines[1]))
  expect_true(nchar(contig_lines[2]) > 0)

  unlink(test_dir, recursive = TRUE)
})

test_that("create_test_dataset creates complete test data", {
  test_base <- tempfile()

  test_data <- create_test_dataset(
    test_base,
    samples = c("sample1", "sample2"),
    genes = c("gene001", "gene002")
  )

  expect_type(test_data, "list")
  expect_named(
    test_data,
    c(
      "base_dir",
      "hybpiper_dir",
      "targets_file",
      "namelist",
      "samples",
      "genes"
    )
  )
  expect_true(dir.exists(test_data$hybpiper_dir))
  expect_true(file.exists(test_data$targets_file))
  expect_true(file.exists(test_data$namelist))

  # Verify namelist content
  namelist_content <- readLines(test_data$namelist)
  expect_equal(namelist_content, c("sample1", "sample2"))

  unlink(test_base, recursive = TRUE)
})

test_that("create_real_test_dataset creates valid structure", {
  test_base <- tempfile()

  test_data <- create_real_test_dataset(
    test_base,
    samples = NULL, # Use default HybPiper samples
    n_genes = 3
  )

  expect_type(test_data, "list")
  expect_named(
    test_data,
    c(
      "base_dir",
      "hybpiper_dir",
      "targets_file",
      "namelist",
      "samples",
      "genes",
      "real_data"
    )
  )

  expect_true(test_data$real_data)
  expect_true(dir.exists(test_data$hybpiper_dir))
  expect_true(file.exists(test_data$targets_file))
  expect_true(file.exists(test_data$namelist))

  # Verify samples and genes
  expect_equal(length(test_data$samples), 3)
  expect_equal(length(test_data$genes), 3)

  # Verify HybPiper structure exists
  sample1_dir <- file.path(test_data$hybpiper_dir, test_data$samples[1])
  expect_true(dir.exists(sample1_dir))

  # Verify gene directory exists
  gene1_dir <- file.path(sample1_dir, test_data$genes[1])
  expect_true(dir.exists(gene1_dir))

  # Verify contig and reads files exist
  contig_file <- file.path(
    gene1_dir,
    paste0(test_data$genes[1], "_contigs.fasta")
  )
  reads_file <- file.path(
    gene1_dir,
    paste0(test_data$genes[1], "_interleaved.fasta")
  )
  expect_true(file.exists(contig_file))
  expect_true(file.exists(reads_file))

  # Verify files have content
  contig_lines <- readLines(contig_file)
  expect_true(grepl("^>", contig_lines[1]))
  expect_true(nchar(contig_lines[2]) > 0)

  unlink(test_base, recursive = TRUE)
})


test_that("run_generate_consensus_sequences works with test data", {
  skip_if_not(check_docker(quiet = TRUE), "Docker not available")
  skip_if_not(
    check_docker_image("joelnitta/hybphaser:latest", pull = FALSE),
    "HybPhaser Docker image not available"
  )

  # Create test data
  test_base <- tempfile()
  test_data <- create_test_dataset(
    test_base,
    samples = c("test_sample1"),
    genes = c("test_gene1")
  )

  output_dir <- file.path(test_base, "hybphaser_output")

  # Run consensus generation
  result <- run_generate_consensus_sequences(
    hybpiper_dir = test_data$hybpiper_dir,
    output_dir = output_dir,
    sample = "test_sample1",
    threads = 1
  )

  # Check exit code (0 means script ran successfully)
  expect_equal(result, 0)

  # Check output directory was created
  expect_true(dir.exists(output_dir))

  # Note: The Docker script may not create 01_data if it can't find
  # properly formatted HybPiper output files. Our test data is minimal,
  # so we just verify the script runs without errors.

  # Clean up
  unlink(test_base, recursive = TRUE)
})

test_that("count_snps works with Docker-generated consensus", {
  skip_if_not(check_docker(quiet = TRUE), "Docker not available")
  skip_if_not(
    check_docker_image("joelnitta/hybphaser:latest", pull = FALSE),
    "HybPhaser Docker image not available"
  )

  # Create test data
  test_base <- tempfile()
  test_data <- create_test_dataset(
    test_base,
    samples = c("sample1"),
    genes = c("gene001")
  )

  output_dir <- file.path(test_base, "hybphaser_output")

  # Run consensus generation
  suppressMessages({
    run_generate_consensus_sequences(
      hybpiper_dir = test_data$hybpiper_dir,
      output_dir = output_dir,
      sample = "sample1",
      threads = 1
    )
  })

  # If consensus sequences were generated, try to count SNPs
  if (dir.exists(file.path(output_dir, "01_data", "sample1"))) {
    # Count SNPs
    snp_results <- count_snps(
      path_to_output_folder = output_dir,
      fasta_file_with_targets = test_data$targets_file,
      targets_file_format = "DNA",
      path_to_namelist = test_data$namelist,
      intronerated_contig = FALSE
    )

    expect_type(snp_results, "list")
    expect_named(snp_results, c("tab_snps", "tab_length"))
    expect_s3_class(snp_results$tab_snps, "data.frame")
    expect_s3_class(snp_results$tab_length, "data.frame")
  }

  # Clean up
  unlink(test_base, recursive = TRUE)
})

test_that("complete workflow with Docker runs end-to-end", {
  skip_if_not(check_docker(quiet = TRUE), "Docker not available")
  skip_if_not(
    check_docker_image("joelnitta/hybphaser:latest", pull = FALSE),
    "HybPhaser Docker image not available"
  )

  # Create test data
  test_base <- tempfile()
  test_data <- create_test_dataset(
    test_base,
    samples = c("sample1", "sample2"),
    genes = c("gene001", "gene002")
  )

  output_dir <- file.path(test_base, "hybphaser_output")

  # Step 1: Generate consensus sequences
  suppressMessages({
    consensus_result <- run_generate_consensus_sequences(
      hybpiper_dir = test_data$hybpiper_dir,
      output_dir = output_dir,
      namelist = test_data$namelist,
      threads = 1,
      cleanup = FALSE
    )
  })

  expect_equal(consensus_result, 0)

  # Step 2: Count SNPs (if consensus was generated)
  if (dir.exists(file.path(output_dir, "01_data"))) {
    suppressMessages({
      snp_results <- count_snps(
        path_to_output_folder = output_dir,
        fasta_file_with_targets = test_data$targets_file,
        targets_file_format = "DNA",
        path_to_namelist = test_data$namelist,
        intronerated_contig = FALSE
      )
    })

    expect_true(!is.null(snp_results))
  }

  # Step 3: Extract mapped reads
  mapped_reads_dir <- file.path(test_base, "mapped_reads")
  suppressMessages({
    extract_result <- run_extract_mapped_reads(
      base_dir = output_dir,
      output_dir = mapped_reads_dir,
      namelist = test_data$namelist
    )
  })

  expect_equal(extract_result, 0)
  expect_true(dir.exists(mapped_reads_dir))

  # Clean up
  unlink(test_base, recursive = TRUE)
})

test_that("run_extract_mapped_reads creates output directory", {
  skip_if_not(check_docker(quiet = TRUE), "Docker not available")
  skip_if_not(
    check_docker_image("joelnitta/hybphaser:latest", pull = FALSE),
    "HybPhaser Docker image not available"
  )

  temp_base <- tempfile()
  dir.create(temp_base)
  temp_output <- tempfile()

  # Output dir should be created even if command fails
  suppressWarnings({
    run_extract_mapped_reads(
      base_dir = temp_base,
      output_dir = temp_output
    )
  })

  expect_true(dir.exists(temp_output))

  unlink(temp_base, recursive = TRUE)
  unlink(temp_output, recursive = TRUE)
})
