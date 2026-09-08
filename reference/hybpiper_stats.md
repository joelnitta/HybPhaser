# Run HybPiper Stats via Conda

Runs `hybpiper stats` across one or more assembled samples.

## Usage

``` r
hybpiper_stats(
  targets_file,
  namelist,
  wd = getwd(),
  mode = c("gene", "supercontig"),
  dna = TRUE,
  seq_lengths_filename = "seq_lengths",
  stats_filename = "hybpiper_stats",
  conda_env = "hybpiper_env",
  other_args = character(0)
)
```

## Arguments

- targets_file:

  Character; path to target FASTA file

- namelist:

  Character; path to text file with one sample per line

- wd:

  Character; working directory for output files

- mode:

  Character; one of "gene" or "supercontig"

- dna:

  Logical; TRUE for DNA targets, FALSE for amino-acid targets

- seq_lengths_filename:

  Character; basename for sequence lengths TSV

- stats_filename:

  Character; basename for stats TSV

- conda_env:

  Character; conda environment name containing HybPiper

- other_args:

  Character vector; additional CLI args passed to `hybpiper stats`

## Value

Named list with paths to output TSV files

## Details

When `mode = "supercontig"`, HybPhaser retries with a native R fallback
if `hybpiper stats` fails. The fallback writes HybPiper-compatible
`seq_lengths.tsv` and `hybpiper_stats.tsv` files directly from the
assembled supercontig FASTA files.
