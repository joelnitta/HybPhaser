# Create Complete Test Dataset

Creates a complete test dataset including HybPiper output, targets, and
namelist

## Usage

``` r
create_test_dataset(
  base_dir,
  samples = c("sample1", "sample2"),
  genes = c("gene001", "gene002")
)
```

## Arguments

- base_dir:

  Base directory for test data

- samples:

  Character vector of sample names

- genes:

  Character vector of gene names

## Value

Named list with paths to created files

## Examples

``` r
if (FALSE) { # \dontrun{
test_data <- create_test_dataset(tempdir())

# Use with HybPhaser
run_generate_consensus_sequences(
  hybpiper_dir = test_data$hybpiper_dir,
  output_dir = file.path(tempdir(), "output"),
  namelist = test_data$namelist
)
} # }
```
