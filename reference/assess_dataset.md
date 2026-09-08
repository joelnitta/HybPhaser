# Assess HybPhaser dataset quality

Generate tables and graphs to assess sequence recovery, missing data,
and identify putative paralogs based on SNP proportions. This function
applies quality filters and produces comprehensive assessment reports.

## Usage

``` r
assess_dataset(
  snp_table,
  length_table,
  targets_file,
  targets_type = c("DNA", "AA"),
  output_dir,
  subset_name = "",
  min_loci_per_sample_prop = 0.5,
  min_target_length_per_sample_prop = 0.5,
  min_samples_per_locus_prop = 0.5,
  min_target_length_per_locus_prop = 0.5,
  paralog_threshold = "none",
  paralog_file = NULL,
  remove_sample_outliers = FALSE
)
```

## Arguments

- snp_table:

  Path to RDS file containing SNP proportions table (output from
  [`count_snps()`](https://joelnitta.github.io/rhybphaser/reference/count_snps.md)),
  or a matrix/data.frame with loci as rows and samples as columns.

- length_table:

  Path to RDS file containing sequence length table (output from
  [`count_snps()`](https://joelnitta.github.io/rhybphaser/reference/count_snps.md)),
  or a matrix/data.frame with loci as rows and samples as columns.

- targets_file:

  Path to FASTA file with target sequences.

- targets_type:

  Character string, either "DNA" or "AA" (amino acid). If "AA", lengths
  will be multiplied by 3 to get nucleotide length.

- output_dir:

  Path to output directory where assessment results will be saved.

- subset_name:

  Optional name for dataset optimization subset. If provided, creates
  subfolder "02_assessment\_\<subset_name\>". If empty string (default),
  creates "02_assessment" folder.

- min_loci_per_sample_prop:

  Minimum proportion of loci that must be recovered per sample (default:
  0.5). Samples below this threshold are removed.

- min_target_length_per_sample_prop:

  Minimum proportion of target sequence length that must be recovered
  per sample (default: 0.5). Samples below this threshold are removed.

- min_samples_per_locus_prop:

  Minimum proportion of samples that must recover a locus (default:
  0.5). Loci below this threshold are removed.

- min_target_length_per_locus_prop:

  Minimum proportion of target sequence length that must be recovered
  per locus (default: 0.5). Loci below this threshold are removed.

- paralog_threshold:

  Threshold for removing putative paralog loci across all samples. Can
  be:

  - Numeric: remove loci with mean SNP proportion above this value

  - "outliers": use 1.5\*IQR method to identify outliers

  - "file": read list of loci to remove from `paralog_file`

  - "none" or NA: skip paralog removal (default)

- paralog_file:

  Path to text file with list of paralog loci to remove (one per line).
  Only used if `paralog_threshold = "file"`.

- remove_sample_outliers:

  Logical. If TRUE, removes outlier loci for each sample individually
  using 1.5\*IQR method (default: FALSE).

## Value

Invisibly returns a list with components:

- `snp_table_cleaned`: Cleaned SNP table after filtering

- `length_table_cleaned`: Cleaned length table after filtering

- `samples_removed`: Character vector of removed sample names

- `loci_removed_missing`: Character vector of loci removed for missing
  data

- `loci_removed_paralogs_all`: Character vector of loci removed as
  paralogs across all samples

- `loci_removed_paralogs_each`: Named list of loci removed per sample

- `summary_table`: Data frame with heterozygosity and allele divergence
  metrics per sample

## Details

The function performs three main assessment steps:

**Step 1: Missing Data Assessment**

- Calculates proportion of loci recovered per sample

- Calculates proportion of samples recovering each locus

- Calculates proportion of target sequence length recovered

- Removes samples and loci failing to meet minimum thresholds

- Generates plots and tables showing data recovery

**Step 2a: Paralog Detection (All Samples)**

- Identifies loci with unusually high SNP proportions across all samples

- Applies specified threshold method

- Generates barplots showing mean SNP proportions per locus

**Step 2b: Paralog Detection (Per Sample)**

- Optionally identifies outlier loci for each sample individually

- Uses 1.5\*IQR method on loci with SNPs

- Generates boxplots showing SNP distributions per sample

**Step 3: Heterozygosity and Allele Divergence**

- Calculates locus heterozygosity (proportion of loci with SNPs)

- Calculates allele divergence (weighted mean SNP proportion)

- Generates summary table and scatter plots

Output files are saved to `<output_dir>/02_assessment[_<subset_name>]/`:

- `1_Data_recovered_overview.pdf/png`: Missing data plots

- `1_Data_recovered_per_sample.csv`: Recovery metrics per sample

- `1_Data_recovered_per_locus.csv`: Recovery metrics per locus

- `1_Summary_missing_data.txt`: Text summary of removed samples/loci

- `2a_Paralogs_for_all_samples.pdf/png`: Paralog detection plots

- `2a_List_of_paralogs_removed_for_all_samples.txt`: List of removed
  loci

- `2b_Paralogs_for_each_sample.pdf/png`: Per-sample outlier plots

- `2_Summary_Paralogs.txt`: Text summary of paralog removal

- `3_LH_vs_AD.pdf/png`: Locus heterozygosity vs allele divergence plot

- `3_varLH_vs_AD.pdf/png`: Multiple heterozygosity threshold plots

- `4_Summary_table.csv`: Final metrics per sample

- `0_Table_SNPs.csv`: Cleaned SNP table

- `0_Table_consensus_length.csv`: Cleaned length table

- `0_namelist_included_samples.txt`: List of retained samples

R objects are also saved to `<output_dir>/00_R_objects/<subset_name>/`:

- `Table_SNPs_cleaned.Rds`

- `Table_consensus_length_cleaned.Rds`

- `Summary_table.Rds`

- Various vectors of removed samples/loci

## Examples

``` r
if (FALSE) { # \dontrun{
# Run count_snps first
snp_results <- count_snps(
  consensus_dir = "01_seqs",
  output_dir = "output"
)

# Then assess dataset
assess_dataset(
  snp_table = file.path("output", "00_R_objects", "Table_SNPs.Rds"),
  length_table = file.path("output", "00_R_objects",
    "Table_consensus_length.Rds"),
  targets_file = "targets.fasta",
  targets_type = "DNA",
  output_dir = "output",
  min_loci_per_sample_prop = 0.5,
  min_samples_per_locus_prop = 0.5,
  paralog_threshold = "outliers",
  remove_sample_outliers = TRUE
)
} # }
```
