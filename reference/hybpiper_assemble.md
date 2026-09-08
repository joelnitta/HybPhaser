# Run HybPiper Assemble via Conda

Runs `hybpiper assemble` for one paired-end sample in a conda
environment.

## Usage

``` r
hybpiper_assemble(
  sample_name,
  f_read,
  r_read,
  targets_file,
  wd,
  mapper = c("bwa", "blast"),
  cpu = 1,
  conda_env = "hybpiper_env",
  other_args = character(0)
)
```

## Arguments

- sample_name:

  Character; sample prefix for HybPiper output

- f_read:

  Character; path to forward read file (R1)

- r_read:

  Character; path to reverse read file (R2)

- targets_file:

  Character; path to DNA targets FASTA file

- wd:

  Character; working directory where HybPiper output is written

- mapper:

  Character; mapping method, one of "bwa" or "blast"

- cpu:

  Integer; number of CPUs for HybPiper

- conda_env:

  Character; conda environment name containing HybPiper

- other_args:

  Character vector; additional CLI args passed to `hybpiper assemble`

## Value

Character path to expected sample output directory
