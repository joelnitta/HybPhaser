# Create Test HybPiper Output Using Real Test Data

Creates a minimal but realistic HybPiper output directory structure
using test data downloaded from the upstream HybPiper project. This
creates a structure suitable for testing HybPhaser functions with
realistic gene names and sequences.

## Usage

``` r
create_hybpiper_test_output(output_dir, samples = NULL, n_genes = 5)
```

## Arguments

- output_dir:

  Path where test HybPiper output should be created

- samples:

  Character vector of sample names. If NULL (default), uses samples from
  HybPiper's test dataset namelist.

- n_genes:

  Integer; number of genes to include (default: 5). Will use the first n
  genes from the HybPiper test targets file.

## Value

Path to created test directory (invisibly)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create test data from real HybPiper test files
test_dir <- create_hybpiper_test_output(tempdir())

# Use with HybPhaser functions
run_generate_consensus_sequences(
  hybpiper_dir = test_dir,
  output_dir = file.path(tempdir(), "hybphaser_out")
)
} # }
```
