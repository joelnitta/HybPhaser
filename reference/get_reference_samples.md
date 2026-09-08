# Get Reference Samples from Consensus Sequence Files

Returns reference sample names from consensus FASTA filenames. Includes
a `sequence_lists` argument so `{targets}` pipelines can enforce
ordering after
[`generate_sequence_lists()`](https://joelnitta.github.io/rhybphaser/reference/generate_sequence_lists.md).

## Usage

``` r
get_reference_samples(
  sequence_lists,
  reference_dir,
  pattern = "_consensus\\.fasta$"
)
```

## Arguments

- sequence_lists:

  Output from
  [`generate_sequence_lists()`](https://joelnitta.github.io/rhybphaser/reference/generate_sequence_lists.md).

- reference_dir:

  Directory containing `*_consensus.fasta` files.

- pattern:

  Regex pattern used to match consensus files.

## Value

Character vector of reference sample names.

## Examples

``` r
if (FALSE) { # \dontrun{
refs <- get_reference_samples(
  sequence_lists = seq_lists,
  reference_dir = "hybphaser_output/03_sequence_lists/samples_consensus"
)
} # }
```
