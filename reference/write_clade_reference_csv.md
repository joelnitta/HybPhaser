# Write Clade Reference CSV

Writes a clade-reference CSV file with `samples` and `abb` columns.

## Usage

``` r
write_clade_reference_csv(reference_samples, output_file)
```

## Arguments

- reference_samples:

  Character vector of reference sample names.

- output_file:

  Path to output CSV file.

## Value

Character scalar. Path to created CSV file.

## Examples

``` r
if (FALSE) { # \dontrun{
write_clade_reference_csv(
  reference_samples = c("sampleA", "sampleB"),
  output_file = "04_clade_association/clade_references_example.csv"
)
} # }
```
