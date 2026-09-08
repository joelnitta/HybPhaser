# Create Test HybPiper Output

Creates a minimal but realistic HybPiper output directory structure for
testing and examples. Includes mock sequence data.

## Usage

``` r
create_test_hybpiper_output(
  base_dir,
  samples = c("sample1", "sample2"),
  genes = c("gene001", "gene002")
)
```

## Arguments

- base_dir:

  Path where test data should be created

- samples:

  Character vector of sample names. Default is c("sample1", "sample2")

- genes:

  Character vector of gene names. Default is c("gene001", "gene002")

## Value

Path to created test directory (invisibly)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create test data
test_dir <- create_test_hybpiper_output(tempdir())

# Use with HybPhaser functions
run_generate_consensus_sequences(
  hybpiper_dir = test_dir,
  output_dir = file.path(tempdir(), "hybphaser_out"),
  sample = "sample1"
)
} # }
```
