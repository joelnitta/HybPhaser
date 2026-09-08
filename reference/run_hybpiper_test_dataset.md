# Run Official HybPiper Test Dataset via Conda

Downloads and runs the official HybPiper test dataset via conda. The
dataset is fetched from the upstream HybPiper repository each run.

## Usage

``` r
run_hybpiper_test_dataset(
  output_dir,
  samples = c("EG30", "EG98", "MWL2"),
  conda_env = "hybpiper_env",
  cpu = 1,
  keep_fastq = FALSE,
  mapper = c("bwa", "blast"),
  other_args = character(0)
)
```

## Arguments

- output_dir:

  Character; directory where HybPiper output is written

- samples:

  Character vector of sample names from the test dataset

- conda_env:

  Character; conda environment name containing HybPiper

- cpu:

  Integer; number of CPUs for HybPiper

- keep_fastq:

  Logical; if TRUE, keep extracted FASTQ files in `output_dir`;
  otherwise remove them after assembly

- mapper:

  Character; one of "bwa" or "blast"

- other_args:

  Character vector; additional CLI args passed to `hybpiper assemble`

## Value

Named list with paths for downstream HybPhaser functions:
`hybpiper_dir`, `targets_file`, and `namelist`
