# Write Phasing Prep CSV

Writes a minimal phasing-prep CSV file with columns `sample`, `ref1`,
and `abb1`.

## Usage

``` r
write_phasing_prep_csv(reference_samples, output_file, abbreviation = "A")
```

## Arguments

- reference_samples:

  Character vector of reference sample names.

- output_file:

  Path to output CSV file.

- abbreviation:

  Character scalar for `abb1` value. Default is "A".

## Value

Character scalar. Path to created CSV file.

## Examples

``` r
if (FALSE) { # \dontrun{
write_phasing_prep_csv(
  reference_samples = c("sampleA", "sampleB"),
  output_file = "05_phasing/phasing_prep_example.csv"
)
} # }
```
