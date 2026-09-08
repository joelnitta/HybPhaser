# Create Namelist from Phased Read Files

Creates a namelist text file from phased read file names in a directory.

## Usage

``` r
make_phased_namelist(
  phased_reads_dir,
  phased_hp_dir,
  output_filename = "namelist_phased.txt"
)
```

## Arguments

- phased_reads_dir:

  Path to directory containing phased read files.

- phased_hp_dir:

  Directory where the namelist file should be written.

- output_filename:

  Name of output namelist file.

## Value

Character scalar. Path to created namelist file.

## Examples

``` r
if (FALSE) { # \dontrun{
namelist_path <- make_phased_namelist(
  phased_reads_dir = "hybphaser_output/05_phasing/phased_reads",
  phased_hp_dir = "hybpiper_output_phased"
)
} # }
```
