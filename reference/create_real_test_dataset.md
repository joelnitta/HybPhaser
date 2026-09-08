# Create Complete Test Dataset Using Real HybPiper Data

Creates a complete test dataset using real target sequences and gene
names from the HybPiper project's test dataset. The dataset is
downloaded from upstream when this function runs. This provides more
realistic testing data than synthetic sequences.

## Usage

``` r
create_real_test_dataset(base_dir, samples = NULL, n_genes = 5)
```

## Arguments

- base_dir:

  Base directory for test data

- samples:

  Character vector of sample names. If NULL (default), uses samples from
  HybPiper test dataset.

- n_genes:

  Integer; number of genes to include (default: 5)

## Value

Named list with paths to created files and metadata

## Examples

``` r
if (FALSE) { # \dontrun{
# Create test data using real HybPiper sequences
test_data <- create_real_test_dataset(tempdir())

# Use with HybPhaser
run_generate_consensus_sequences(
  hybpiper_dir = test_data$hybpiper_dir,
  output_dir = file.path(tempdir(), "output"),
  namelist = test_data$namelist
)
} # }
```
