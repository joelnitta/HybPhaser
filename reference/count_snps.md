# Count SNPs in Consensus Sequences

Counts polymorphic sites (masked as ambiguity codes) in consensus
sequences and calculates sequence statistics for each sample and locus.

## Usage

``` r
count_snps(
  path_to_output_folder,
  fasta_file_with_targets,
  targets_file_format = c("DNA", "AA"),
  path_to_namelist,
  intronerated_contig = FALSE,
  subset_name = ""
)
```

## Arguments

- path_to_output_folder:

  Path to HybPhaser output folder containing consensus sequences

- fasta_file_with_targets:

  Path to FASTA file with target sequences

- targets_file_format:

  Format of target file: "DNA" or "AA" (amino acid)

- path_to_namelist:

  Path to file containing sample names (one per line)

- intronerated_contig:

  Logical; whether to use intronerated contigs (TRUE) or normal contigs
  (FALSE)

- subset_name:

  Optional name for dataset optimization subset

## Value

List containing two data frames:

- tab_snps:

  Proportion of SNPs per sample and locus

- tab_length:

  Sequence length per sample and locus

## Examples

``` r
if (FALSE) { # \dontrun{
results <- count_snps(
  path_to_output_folder = "hybphaser_output",
  fasta_file_with_targets = "targets.fasta",
  targets_file_format = "DNA",
  path_to_namelist = "samples.txt",
  intronerated_contig = FALSE
)
} # }
```
