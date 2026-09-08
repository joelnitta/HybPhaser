# rhybphaser Example Workflow

``` r

library(rhybphaser)

# Check if Docker is available for this vignette
has_docker <- check_docker(quiet = TRUE)
has_image <- FALSE
if (has_docker) {
  has_image <- check_docker_image("joelnitta/rhybphaser:latest",
    pull = FALSE
  )
}
#> Docker image 'joelnitta/rhybphaser:latest' not found. Set pull = TRUE to download it.

# The end-to-end workflow below needs conda, HybPiper, Docker and network
# access, so it does not run during package checks or when building the
# site. Set RHYBPHASER_RUN_VIGNETTE=true to execute it (the vignettes.yaml
# workflow does this).
run_full_vignette <- identical(Sys.getenv("RHYBPHASER_RUN_VIGNETTE"), "true")
```

## Introduction

HybPhaser is designed to detect and phase hybrid accessions in target
capture datasets from polyploid species. It works as an extension to the
[HybPiper](https://github.com/mossmatters/HybPiper) assembly pipeline.

The core workflow involves:

1.  **SNP Assessment** - Detect hybrids by measuring heterozygosity in
    consensus sequences
2.  **Dataset Quality Assessment** - Filter samples and loci based on
    missing data and heterozygosity
3.  **Sequence List Generation** - Create organized FASTA files for loci
    and samples
4.  **Clade Association** - Map reads to divergent reference clades
5.  **Phasing** - Separate reads according to similarity with parental
    clades

This vignette demonstrates the package workflow using the functions
currently available. Additional functionality is being added as the
package conversion continues.

## Example with Test Data

HybPhaser can download the official HybPiper test dataset on demand. You
can run HybPiper on that data via conda, then use the resulting real
HybPiper output with HybPhaser.

The example below is the full end-to-end workflow:

*Part 1 — SNP assessment (normal run)*

1.  Run HybPiper on all samples (conda)
2.  Generate consensus sequences (Docker)
3.  Count SNPs (R)
4.  Assess dataset quality (R)
5.  Generate sequence lists (R)

*Part 2 — Clade association*

6.  Extract mapped reads for clade association (Docker)
7.  Run BBSplit clade association (R)
8.  Collate BBSplit results (R)

*Part 3 — Phasing*

9.  Run BBSplit phasing (R)
10. Collate phasing BBSplit stats (R)

*Part 1 repeated — SNP assessment on phased reads*

11. Re-run HybPiper on phased reads (conda)
12. Re-run HybPhaser Part 1 on phased HybPiper output (Docker + R)

*Part 4 — Merge sequence lists*

13. Merge phased and non-phased sequence lists (R)

If HybPiper is installed in a conda environment (default name:
`hybpiper_env`), you can download and run the official test dataset and
continue directly with HybPhaser:

``` r

# 0. Optional checks
check_conda()
check_hybpiper_conda("hybpiper_env")
check_docker()
check_docker_image("joelnitta/rhybphaser:latest", pull = TRUE)

# 1. Download and run official HybPiper test data using conda
hp_data <- run_hybpiper_test_dataset(
  output_dir = tempfile("hybpiper_test_"),
  samples = c("EG30"),
  conda_env = "hybpiper_env",
  cpu = 1,
  mapper = "bwa"
)

# 2. Generate consensus sequences with HybPhaser Docker workflow
run_generate_consensus_sequences(
  hybpiper_dir = hp_data$hybpiper_dir,
  output_dir = file.path(hp_data$hybpiper_dir, "hybphaser_output"),
  namelist = hp_data$namelist,
  threads = 1
)

# 3. Count SNPs
snp_results <- count_snps(
  path_to_output_folder = file.path(
    hp_data$hybpiper_dir,
    "hybphaser_output"
  ),
  fasta_file_with_targets = hp_data$targets_file,
  targets_file_format = "DNA",
  path_to_namelist = hp_data$namelist
)

# 4. Assess dataset quality
# Generate tables and plots to assess sequence recovery and identify
# putative paralogs
assessment <- assess_dataset(
  snp_table = snp_results$tab_snps,
  length_table = snp_results$tab_length,
  targets_file = hp_data$targets_file,
  targets_type = "DNA",
  output_dir = file.path(hp_data$hybpiper_dir, "hybphaser_output"),
  min_loci_per_sample_prop = 0.5,
  min_samples_per_locus_prop = 0.5,
  paralog_threshold = "outliers",
  remove_sample_outliers = TRUE
)

# 5. Generate sequence lists
# Create FASTA files organized by locus and by sample, applying quality
# filters from dataset assessment
seq_lists <- generate_sequence_lists(
  path_to_output_folder = file.path(
    hp_data$hybpiper_dir,
    "hybphaser_output"
  ),
  fasta_file_with_targets = hp_data$targets_file,
  targets_file_format = "DNA",
  path_to_namelist = hp_data$namelist
)

# 6. Extract mapped reads for clade association
mapped_reads <- extract_mapped_reads(
  base_dir = file.path(hp_data$hybpiper_dir, "hybphaser_output"),
  output_dir = file.path(
    hp_data$hybpiper_dir,
    "hybphaser_output",
    "mapped_reads"
  ),
  namelist = hp_data$namelist,
  remove_duplicate_sequences = FALSE
)

# 7. Inspect results
head(snp_results$tab_snps)
head(snp_results$tab_length)
head(assessment$summary_table)
cat("\nSequence lists created:\n")
cat("  Loci:", seq_lists$n_loci_files, "\n")
cat("  Samples:", seq_lists$n_sample_files, "\n")
cat("\nMapped reads files created:", mapped_reads$n_samples, "\n")

# 8. Run BBSplit clade association
# Requires a clade reference CSV with two columns: sample and abbreviation.
# Here we derive a minimal example directly from generated sample sequences.
reference_dir <- file.path(
  hp_data$hybpiper_dir,
  "hybphaser_output",
  "03_sequence_lists",
  "samples_consensus"
)

reference_files <- list.files(reference_dir, pattern = "_consensus\\.fasta$")
reference_samples <- sub("_consensus\\.fasta$", "", reference_files)

clade_reference_csv <- file.path(
  hp_data$hybpiper_dir,
  "hybphaser_output",
  "04_clade_association",
  "clade_references_example.csv"
)
dir.create(dirname(clade_reference_csv), recursive = TRUE, showWarnings = FALSE)

utils::write.csv(
  data.frame(
    samples = reference_samples,
    abb = paste0("C", seq_along(reference_samples))
  ),
  clade_reference_csv,
  row.names = FALSE
)

bbsplit_job <- run_clade_association(
  path_to_clade_association_folder = file.path(
    hp_data$hybpiper_dir,
    "hybphaser_output",
    "04_clade_association"
  ),
  csv_file_with_clade_reference_names = clade_reference_csv,
  path_to_reference_sequences = reference_dir,
  path_to_read_files_cladeassociation = file.path(
    hp_data$hybpiper_dir,
    "hybphaser_output",
    "mapped_reads"
  ),
  read_type_cladeassociation = "single-end"
)

cat("\nBBSplit clade association run completed:\n")
cat("  Script:", bbsplit_job$script_path, "\n")
cat("  Status:", bbsplit_job$run_status, "\n")

# 9. Collate BBSplit results
clade_assoc <- collate_bbsplit_results(
  path_to_clade_association_folder = file.path(
    hp_data$hybpiper_dir,
    "hybphaser_output",
    "04_clade_association"
  ),
  csv_file_with_clade_reference_names = clade_reference_csv,
  path_to_output_folder = file.path(
    hp_data$hybpiper_dir,
    "hybphaser_output"
  )
)

head(clade_assoc$table_clade_association)

# 10. Run phasing
# Create a minimal phasing prep CSV. In real analyses this table usually
# includes multiple references per sample (ref1/abb1, ref2/abb2, ...).
phasing_prep_csv <- file.path(
  hp_data$hybpiper_dir,
  "hybphaser_output",
  "05_phasing",
  "phasing_prep_example.csv"
)
dir.create(dirname(phasing_prep_csv), recursive = TRUE, showWarnings = FALSE)

utils::write.csv(
  data.frame(
    sample = reference_samples[[1]],
    ref1 = reference_samples[[1]],
    abb1 = "A"
  ),
  phasing_prep_csv,
  row.names = FALSE
)

phasing_job <- run_phasing(
  path_to_phasing_folder = file.path(
    hp_data$hybpiper_dir,
    "hybphaser_output",
    "05_phasing"
  ),
  csv_file_with_phasing_prep_info = phasing_prep_csv,
  path_to_read_files_phasing = file.path(
    hp_data$hybpiper_dir,
    "hybphaser_output",
    "mapped_reads"
  ),
  read_type_4phasing = "single-end",
  reference_sequence_folder = reference_dir
)

cat("\nPhasing run completed:\n")
cat("  Script:", phasing_job$script_path, "\n")
cat("  Status:", phasing_job$run_status, "\n")

# 11. Collate phasing stats
phasing_stats <- collate_phasing_stats(
  path_to_phasing_folder = file.path(
    hp_data$hybpiper_dir,
    "hybphaser_output",
    "05_phasing"
  ),
  csv_file_with_phasing_prep_info = phasing_prep_csv
)

head(phasing_stats$table_phasing_stats)

# 11. Re-run HybPiper on phased reads
# After phasing (steps 9-10), BBSplit produces per-sample per-reference
# read files inside the phasing folder:
#   single-end : {sample}_to_{ref}.fastq
#   paired-end : {sample}_to_{ref}_1.fastq / {sample}_to_{ref}_2.fastq
# Each file represents one phased haplotype; reassemble them individually
# with HybPiper, writing output to a separate directory.
phased_reads_dir <- file.path(
  hp_data$hybpiper_dir,
  "hybphaser_output",
  "05_phasing",
  "phased_reads"
)
phased_hp_dir <- file.path(hp_data$hybpiper_dir, "hybpiper_output_phased")
dir.create(phased_hp_dir, recursive = TRUE, showWarnings = FALSE)

# Collect single-end phased read files (adjust pattern for paired-end)
phased_read_files <- list.files(
  phased_reads_dir,
  pattern = "\\.fastq(\\.gz)?$",
  full.names = TRUE
)
phased_samples <- sub("\\.fastq(\\.gz)?$", "", basename(phased_read_files))

phased_namelist <- file.path(phased_hp_dir, "namelist_phased.txt")
writeLines(phased_samples, phased_namelist)

# Run HybPiper assemble for each phased sample (single-end example)
for (i in seq_along(phased_samples)) {
  run_conda(
    command = "hybpiper",
    env_name = "hybpiper_env",
    wd = phased_hp_dir,
    args = c(
      "assemble",
      "--readfiles", phased_read_files[i],
      "--targetfile_dna", hp_data$targets_file,
      "--prefix", phased_samples[i],
      "--bwa", "--cpu", "1"
    ),
    stdout = "", stderr = ""
  )
}

# 12. Re-run HybPhaser Part 1 on phased HybPiper output
# Steps 2-5 are repeated on the phased dataset to generate phased sequence
# lists that step 13 can merge with the normal sequence lists.
phased_hp_output <- file.path(phased_hp_dir, "hybphaser_output_phased")

run_generate_consensus_sequences(
  hybpiper_dir = phased_hp_dir,
  output_dir = phased_hp_output,
  namelist = phased_namelist,
  threads = 1
)

phased_snps <- count_snps(
  path_to_output_folder = phased_hp_output,
  fasta_file_with_targets = hp_data$targets_file,
  targets_file_format = "DNA",
  path_to_namelist = phased_namelist
)

assess_dataset(
  snp_table = phased_snps$tab_snps,
  length_table = phased_snps$tab_length,
  targets_file = hp_data$targets_file,
  targets_type = "DNA",
  output_dir = phased_hp_output,
  min_loci_per_sample_prop = 0.5,
  min_samples_per_locus_prop = 0.5
)

generate_sequence_lists(
  path_to_output_folder = phased_hp_output,
  fasta_file_with_targets = hp_data$targets_file,
  targets_file_format = "DNA",
  path_to_namelist = phased_namelist
)

# 13. Merge phased and non-phased sequence lists
phased_loci_dir <- file.path(
  phased_hp_output,
  "03_sequence_lists",
  "loci_consensus"
)

merged_seq_lists <- merge_sequence_lists(
  path_to_sequence_lists_normal = file.path(
    hp_data$hybpiper_dir,
    "hybphaser_output",
    "03_sequence_lists",
    "loci_consensus"
  ),
  path_to_sequence_lists_phased = phased_loci_dir,
  path_of_sequence_lists_output = file.path(
    hp_data$hybpiper_dir,
    "hybphaser_output",
    "03_sequence_lists_merged",
    "loci_consensus"
  ),
  path_to_namelist_normal = hp_data$namelist,
  path_to_namelist_phased = phased_namelist,
  exchange_phased_with_not_phased_samples = "yes",
  include_phased_seqlists_when_non_phased_locus_absent = "no"
)

cat("\nMerged sequence-list output:", merged_seq_lists$output_dir, "\n")
```

## Next Steps

After generating merged sequence lists, the combined dataset is ready
for alignment and phylogenetic analyses.

Additional utility functions will continue to be added in future
updates.

## Getting Help

- Function documentation:
  [`?count_snps`](https://joelnitta.github.io/rhybphaser/reference/count_snps.md),
  [`?assess_dataset`](https://joelnitta.github.io/rhybphaser/reference/assess_dataset.md),
  [`?generate_sequence_lists`](https://joelnitta.github.io/rhybphaser/reference/generate_sequence_lists.md),
  [`?extract_mapped_reads`](https://joelnitta.github.io/rhybphaser/reference/extract_mapped_reads.md),
  [`?run_clade_association`](https://joelnitta.github.io/rhybphaser/reference/run_clade_association.md),
  [`?run_clade_association_from_config`](https://joelnitta.github.io/rhybphaser/reference/run_clade_association_from_config.md),
  [`?collate_bbsplit_results`](https://joelnitta.github.io/rhybphaser/reference/collate_bbsplit_results.md),
  [`?run_phasing`](https://joelnitta.github.io/rhybphaser/reference/run_phasing.md),
  [`?run_phasing_from_config`](https://joelnitta.github.io/rhybphaser/reference/run_phasing_from_config.md),
  [`?collate_phasing_stats`](https://joelnitta.github.io/rhybphaser/reference/collate_phasing_stats.md),
  [`?collate_phasing_stats_from_config`](https://joelnitta.github.io/rhybphaser/reference/collate_phasing_stats_from_config.md),
  [`?merge_sequence_lists`](https://joelnitta.github.io/rhybphaser/reference/merge_sequence_lists.md),
  [`?merge_sequence_lists_from_config`](https://joelnitta.github.io/rhybphaser/reference/merge_sequence_lists_from_config.md)
- Report issues: <https://github.com/joelnitta/rhybphaser/issues>
- Original HybPhaser method and code:
  <https://github.com/LarsNauheimer/HybPhaser>
- Original HybPhaser paper: Nauheimer et al. (2021), *Applications in
  Plant Sciences* 9(7): e11441, <https://doi.org/10.1002/aps3.11441>
