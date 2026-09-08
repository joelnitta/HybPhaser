# Collate Phasing BBSplit Stats

Collates phasing `refstats` files into a sample-by-reference table with
values expressed as proportions of unambiguous reads.

## Usage

``` r
collate_phasing_stats(
  path_to_phasing_folder,
  csv_file_with_phasing_prep_info,
  folder_for_phasing_stats = ""
)
```

## Arguments

- path_to_phasing_folder:

  Path to phasing folder where output table is written.

- csv_file_with_phasing_prep_info:

  CSV used for phasing preparation. Expected columns are `sample` or
  `samples`, then `ref1, abb1, ref2, abb2, ...`.

- folder_for_phasing_stats:

  Optional folder containing phasing stats files (default:
  `path_to_phasing_folder/phasing_stats`).

## Value

Invisibly returns a list with the collated table and output path.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- collate_phasing_stats(
  path_to_phasing_folder = "05_phasing",
  csv_file_with_phasing_prep_info = "phasing_prep.csv"
)
} # }
```
