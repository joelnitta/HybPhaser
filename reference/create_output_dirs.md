# Create Output Directory Structure

Creates output directories with optional subset names

## Usage

``` r
create_output_dirs(base_path, subdirs, subset_name = "", recursive = TRUE)
```

## Arguments

- base_path:

  Base output directory path

- subdirs:

  Character vector of subdirectory names to create

- subset_name:

  Optional subset name to append to directories

- recursive:

  Create parent directories if needed

## Value

Character vector of created directory paths (invisibly)
