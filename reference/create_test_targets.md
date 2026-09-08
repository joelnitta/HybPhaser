# Create Test Target File

Creates a FASTA file with target sequences (baits) for testing

## Usage

``` r
create_test_targets(
  file_path,
  genes = c("gene001", "gene002"),
  seq_length = 150
)
```

## Arguments

- file_path:

  Path where target file should be created

- genes:

  Character vector of gene names

- seq_length:

  Integer; length of sequences to generate

## Value

Path to created file (invisibly)

## Examples

``` r
if (FALSE) { # \dontrun{
target_file <- create_test_targets(
  file.path(tempdir(), "targets.fasta")
)
} # }
```
