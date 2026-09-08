# Run HybPiper on Phased Read Files

Runs `hybpiper assemble` for each phased read file in a directory. This
is useful after
[`run_phasing()`](https://joelnitta.github.io/rhybphaser/reference/run_phasing.md)
in workflows where phased reads are assembled separately.

## Usage

``` r
run_phased_hybpiper(
  phased_reads_dir,
  phased_hp_dir,
  targets_file,
  conda_env = "hybpiper_env",
  cpu = 1
)
```

## Arguments

- phased_reads_dir:

  Path to directory containing phased read files.

- phased_hp_dir:

  Path to output directory for phased HybPiper results.

- targets_file:

  Path to DNA targets FASTA file.

- conda_env:

  Conda environment containing HybPiper.

- cpu:

  Number of CPU threads per HybPiper run.

## Value

Character scalar. The `phased_hp_dir` path.

## Examples

``` r
if (FALSE) { # \dontrun{
out_dir <- run_phased_hybpiper(
  phased_reads_dir = "hybphaser_output/05_phasing/phased_reads",
  phased_hp_dir = "hybpiper_output_phased",
  targets_file = "test_targets.fasta"
)
} # }
```
