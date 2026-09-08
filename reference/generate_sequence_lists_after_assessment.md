# Generate Sequence Lists After Assessment

Wrapper around
[`generate_sequence_lists()`](https://joelnitta.github.io/rhybphaser/reference/generate_sequence_lists.md)
that takes an assessment object to make dependency ordering explicit in
`{targets}` pipelines.

## Usage

``` r
generate_sequence_lists_after_assessment(assessment, ...)
```

## Arguments

- assessment:

  Output from
  [`assess_dataset()`](https://joelnitta.github.io/rhybphaser/reference/assess_dataset.md).

- ...:

  Arguments passed to
  [`generate_sequence_lists()`](https://joelnitta.github.io/rhybphaser/reference/generate_sequence_lists.md).

## Value

List returned by
[`generate_sequence_lists()`](https://joelnitta.github.io/rhybphaser/reference/generate_sequence_lists.md).

## Examples

``` r
if (FALSE) { # \dontrun{
seq_lists <- generate_sequence_lists_after_assessment(
  assessment = assessment_result,
  path_to_output_folder = "hybphaser_output",
  fasta_file_with_targets = "targets.fasta",
  targets_file_format = "DNA",
  path_to_namelist = "namelist.txt"
)
} # }
```
