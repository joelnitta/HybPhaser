# Extract Mapped Reads

Extracts reads mapped to target sequences from HybPhaser or HybPiper
output, concatenates them per sample, and optionally deduplicates by
sequence.

## Usage

``` r
extract_mapped_reads(
  base_dir = "./",
  output_dir = "./mapped_reads/",
  namelist = NULL,
  remove_duplicate_sequences = FALSE
)
```

## Arguments

- base_dir:

  Path to base folder of either:

  - HybPhaser output (contains `01_data/`)

  - HybPiper-RBGV output (contains `04_processed_gene_directories/`)

  - HybPiper output (contains sample directories)

- output_dir:

  Path to output directory for mapped read files. Default is
  `./mapped_reads/`.

- namelist:

  Optional path to a file containing sample names, one per line. If
  `NULL`, all sample directories are processed.

- remove_duplicate_sequences:

  Logical; if `TRUE`, deduplicate reads by sequence (not just duplicated
  names).

## Value

Invisibly returns a list with:

- output_dir:

  Output directory path

- samples:

  Character vector of processed samples

- files:

  Named character vector of output files by sample

- n_samples:

  Number of processed samples

## Details

This is an R implementation of the `2_extract_mapped_reads.sh` workflow
step.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- extract_mapped_reads(
  base_dir = "hybphaser_output",
  output_dir = "mapped_reads"
)

result_dedup <- extract_mapped_reads(
  base_dir = "hybphaser_output",
  output_dir = "mapped_reads",
  namelist = "samples.txt",
  remove_duplicate_sequences = TRUE
)
} # }
```
