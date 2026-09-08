# Get Path to HybPhaser Example Data

Returns the full path to example data files included with the package

## Usage

``` r
hybphaser_example(file_name = NULL)
```

## Arguments

- file_name:

  Name of the example file. If NULL, returns the extdata directory.
  Available files:

  - "clade_references.csv"

  - "phasing_prep.csv"

## Value

Character path to the example file or extdata directory

## Examples

``` r
# Get path to extdata directory
data_dir <- hybphaser_example()

# Get path to specific example file
clade_ref <- hybphaser_example("clade_references.csv")
```
