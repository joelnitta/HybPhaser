# Collate BBSplit Results

Collates BBSplit `refstats` files into clade-association tables,
including a normalized table and optional merges with the HybPhaser
summary table.

## Usage

``` r
collate_bbsplit_results(
  path_to_clade_association_folder,
  csv_file_with_clade_reference_names,
  path_to_output_folder = NULL,
  subset_name = ""
)
```

## Arguments

- path_to_clade_association_folder:

  Path to clade association folder containing `bbsplit_stats/`.

- csv_file_with_clade_reference_names:

  CSV with clade reference sample names and abbreviations (first two
  columns).

- path_to_output_folder:

  Optional HybPhaser output folder used to find
  `00_R_objects/.../Summary_table.Rds` for merged output tables.

- subset_name:

  Optional subset folder name under `00_R_objects`.

## Value

Invisibly returns a list with generated tables and output paths.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- collate_bbsplit_results(
  path_to_clade_association_folder = "04_clade_association",
  csv_file_with_clade_reference_names = "clade_references.csv"
)
} # }
```
