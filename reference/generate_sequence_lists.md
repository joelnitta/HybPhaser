# Generate sequence lists for loci and samples

Creates FASTA sequence lists organized by locus and by sample, applying
quality filters from dataset optimization. Generates four sets of output
files (locus-level and sample-level, for both consensus and contig
sequences) that can be used for downstream phylogenetic analyses and
phasing.

## Usage

``` r
generate_sequence_lists(
  path_to_output_folder,
  fasta_file_with_targets,
  targets_file_format = "DNA",
  path_to_namelist,
  intronerated_contig = FALSE,
  subset_name = ""
)
```

## Arguments

- path_to_output_folder:

  Path to the base HybPiper output folder containing the "01_data"
  subfolder with sample results.

- fasta_file_with_targets:

  Path to FASTA file with target sequences used in HybPiper run.

- targets_file_format:

  Character string, either "DNA" or "AA" (amino acid).

- path_to_namelist:

  Path to text file with sample names (one per line) to include in
  sequence lists.

- intronerated_contig:

  Logical. If TRUE, processes intronerated contigs instead of regular
  contigs. Default FALSE.

- subset_name:

  Optional name for dataset optimization subset. If provided, creates
  output in subfolder "03_sequence_lists\_\<subset_name\>" and loads R
  objects from corresponding assessment subfolder. If empty string
  (default), uses "03_sequence_lists" and loads from "00_R_objects".

## Value

Invisibly returns a list with components:

- `output_dir`: Path to main output directory

- `loci_consensus_dir`: Path to locus consensus sequence directory

- `loci_contigs_dir`: Path to locus contig sequence directory

- `samples_consensus_dir`: Path to sample consensus sequence directory

- `samples_contigs_dir`: Path to sample contig sequence directory

- `n_loci_files`: Number of locus files created

- `n_sample_files`: Number of sample files created

- `loci_removed`: Character vector of loci removed from sequence lists

- `samples_removed`: Character vector of samples removed from sequence
  lists

## Details

This function performs the following steps:

1.  **Load filtering results**: Reads R objects saved by
    [`assess_dataset()`](https://joelnitta.github.io/rhybphaser/reference/assess_dataset.md):

    - `outsamples_missing`: Samples removed for missing data

    - `outloci_missing`: Loci removed for missing data

    - `outloci_para_all`: Loci removed as paralogs across all samples

    - `outloci_para_each`: Loci removed as paralogs for individual
      samples

    - `Table_SNPs_cleaned.Rds`: Cleaned SNP table

2.  **Generate locus-level sequence lists**: For each target locus,
    concatenates sequences from all samples into a single FASTA file
    (one file per locus). Creates both consensus and contig versions.

3.  **Generate sample-level sequence lists**: For each sample,
    concatenates sequences from all loci into a single FASTA file (one
    file per sample). Creates both consensus and contig versions.

4.  **Apply quality filters**: Removes sequences based on filters from
    [`assess_dataset()`](https://joelnitta.github.io/rhybphaser/reference/assess_dataset.md):

    - Samples removed for missing data

    - Loci removed for missing data or as paralogs (all samples)

    - Loci removed as paralogs for specific samples

The function handles both Linux and non-Linux systems, using bash
commands (`cat`, `sed`) on Linux for faster processing, and R file
operations otherwise.

Generated sequences are in non-interleaved FASTA format.

## Output directories

Four directories are created under `03_sequence_lists` (or
`03_sequence_lists_<subset_name>`):

- `loci_consensus/`: One file per locus with consensus sequences from
  all samples

- `loci_contigs/`: One file per locus with contig sequences from all
  samples

- `samples_consensus/`: One file per sample with consensus sequences
  from all loci

- `samples_contigs/`: One file per sample with contig sequences from all
  loci

## File naming

- Locus files: `<locus>_consensus.fasta` or
  `<locus>_intronerated_consensus.fasta`

- Sample files: `<sample>_consensus.fasta` or
  `<sample>_intronerated_consensus.fasta`

- Similar naming for contig files (with `_contig` or `_contigs`)

## See also

[`assess_dataset()`](https://joelnitta.github.io/rhybphaser/reference/assess_dataset.md),
[`count_snps()`](https://joelnitta.github.io/rhybphaser/reference/count_snps.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate sequence lists after running assess_dataset
result <- generate_sequence_lists(
  path_to_output_folder = "output",
  fasta_file_with_targets = "targets.fasta",
  targets_file_format = "DNA",
  path_to_namelist = "namelist.txt"
)

# With subset name
result <- generate_sequence_lists(
  path_to_output_folder = "output",
  fasta_file_with_targets = "targets.fasta",
  targets_file_format = "DNA",
  path_to_namelist = "namelist.txt",
  subset_name = "filtered"
)

# View results
cat("Created", result$n_loci_files, "locus files\n")
cat("Created", result$n_sample_files, "sample files\n")
} # }
```
