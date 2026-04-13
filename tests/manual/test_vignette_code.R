# Testing Vignette Code

This
document
tests
the
code
from
the
HybPhaser
vignette
to
ensure
it
runs
correctly
when
the
package
is
installed.

## Test 1: Access Scripts

library(HybPhaser)

# Get path to scripts
script <- hybphaser_scripts("1_generate_consensus_sequences.sh")
cat("Script path:", script, "\n")
stopifnot(file.exists(script))

# List all scripts
scripts <- list.files(hybphaser_scripts())
cat("Available scripts:", paste(scripts, collapse = ", "), "\n")
stopifnot(length(scripts) > 0)

## Test 2: Access Example Data

clade_refs <- read.csv(hybphaser_example("clade_references.csv"))
cat("Clade references loaded:", nrow(clade_refs), "rows\n")
stopifnot(nrow(clade_refs) > 0)

phasing_prep <- read.csv(hybphaser_example("phasing_prep.csv"))
cat("Phasing prep loaded:", nrow(phasing_prep), "rows\n")
stopifnot(nrow(phasing_prep) > 0)

## Test 3: Count SNPs

# Create temporary directory structure
temp_dir <- tempfile()
dir.create(
  file.path(temp_dir, "01_data", "sample1", "consensus"),
  recursive = TRUE
)

# Create a simple target file
target_file <- file.path(temp_dir, "targets.fasta")
writeLines(
  c(
    ">species1-gene1",
    "ATCGATCGATCGATCGATCG",
    ">species1-gene2",
    "GCTAGCTAGCTAGCTAGCTA"
  ),
  target_file
)

# Create consensus sequences
writeLines(
  c(">gene1", "ATCGATRYATCGATCGATCG"),
  file.path(temp_dir, "01_data", "sample1", "consensus", "gene1.fasta")
)

writeLines(
  c(">gene2", "GCTAGCTAGCTAGCTAGCTA"),
  file.path(temp_dir, "01_data", "sample1", "consensus", "gene2.fasta")
)

# Create namelist
namelist <- file.path(temp_dir, "namelist.txt")
writeLines("sample1", namelist)

# Count SNPs
results <- count_snps(
  path_to_output_folder = temp_dir,
  fasta_file_with_targets = target_file,
  targets_file_format = "DNA",
  path_to_namelist = namelist,
  intronerated_contig = FALSE
)

cat("\nResults structure:\n")
str(results)

cat("\nSNP proportions:\n")
print(results$tab_snps)

cat("\nSequence lengths:\n")
print(results$tab_length)

# Verify output files
output_files <- list.files(
  file.path(temp_dir, "00_R_objects"),
  recursive = TRUE
)
cat("\nSaved output files:", paste(output_files, collapse = ", "), "\n")
stopifnot(length(output_files) == 2)

# Clean up
unlink(temp_dir, recursive = TRUE)

cat("\n✓ All vignette code tests passed!\n")
