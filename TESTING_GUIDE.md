# Testing the Package Conversion So Far

This document shows how to test the package conversion work completed so far.

## What's Been Converted

The first major function from the HybPhaser pipeline has been converted:
- **`count_snps()`** - Counts SNPs in consensus sequences (was `1a_count_snps.R`)

Plus supporting infrastructure:
- Utility functions for config reading, path validation, and directory creation
- Data access functions for bash scripts and example files
- Complete test suite using testthat
- Package structure (DESCRIPTION, NAMESPACE, etc.)

## Quick Test

Run these commands in R from the package directory:

```r
# Install development tools if needed
install.packages(c("devtools", "roxygen2", "testthat"))

# Load the package in development mode
devtools::load_all()

# Generate documentation
devtools::document()

# Run all tests
devtools::test()

# Check the package
devtools::check()
```

## Manual Testing

### Test Utility Functions

```r
# Test configuration reading
temp_config <- tempfile(fileext = ".txt")
writeLines('
path_to_output_folder = "/tmp/output"
fasta_file_with_targets = "targets.fasta"
targets_file_format = "DNA"
', temp_config)

config <- read_config(temp_config)
print(config)

# Test path validation
validate_paths(
  list("temp config" = temp_config),
  must_exist = TRUE,
  type = "file"
)
```

### Test Data Access Functions

```r
# Get path to bash scripts directory
scripts_dir <- hybphaser_scripts()
print(scripts_dir)

# Get path to a specific script
consensus_script <- hybphaser_scripts("1_generate_consensus_sequences.sh")
print(consensus_script)

# Get path to example data
clade_refs <- hybphaser_example("clade_references.csv")
print(read.csv(clade_refs))

# Preview a bash command (doesn't actually run it)
run_hybphaser_script(
  "1_generate_consensus_sequences.sh",
  args = c("-n", "namelist.txt", "-h"),
  dry_run = TRUE
)
```

### Test SNP Counting

To test `count_snps()`, you need a minimal HybPhaser output structure. Here's a test setup:

```r
# Create minimal test structure
base_dir <- tempfile()
data_dir <- file.path(base_dir, "01_data", "sample1", "consensus")
dir.create(data_dir, recursive = TRUE)

# Create target file
target_file <- tempfile(fileext = ".fasta")
writeLines(c(
  ">species1-gene1",
  "ATCGATCGATCG"
), target_file)

# Create consensus file
consensus_file <- file.path(data_dir, "gene1.fasta")
writeLines(c(
  ">gene1",
  "ATCGATRYATCG"  # R and Y are ambiguity codes (SNPs)
), consensus_file)

# Create namelist
namelist <- tempfile()
writeLines("sample1", namelist)

# Run count_snps
result <- count_snps(
  path_to_output_folder = base_dir,
  fasta_file_with_targets = target_file,
  targets_file_format = "DNA",
  path_to_namelist = namelist,
  intronerated_contig = FALSE
)

# Check results
print(result$tab_snps)
print(result$tab_length)

# Check that output files were created
list.files(file.path(base_dir, "00_R_objects"), recursive = TRUE)

# Clean up
unlink(base_dir, recursive = TRUE)
unlink(target_file)
unlink(namelist)
```

## What to Look For

1. **Documentation**: Are the roxygen2 comments clear and complete?
2. **Function interfaces**: Are the function arguments well-named and documented?
3. **Error handling**: Do functions validate inputs appropriately?
4. **Tests**: Do the tests cover the main use cases and edge cases?
5. **Code style**: Is the code readable and following R conventions?

## Known Issues

- The package needs to be properly installed for `system.file()` calls in `hybphaser_scripts()` and `hybphaser_example()` to work correctly
- Some tests use `skip_if_not_installed()` to handle this during development

## Next Steps

After you review this first function, we'll continue converting:
1. `1b_assess_dataset.R` → `assess_dataset()`
2. `1c_generate_sequence_lists.R` → `generate_sequence_lists()`
3. And so on for the remaining scripts...
