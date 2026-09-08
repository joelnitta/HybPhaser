# Validate File and Directory Paths

Checks that specified paths exist and have correct type (file or
directory)

## Usage

``` r
validate_paths(paths, must_exist = TRUE, type = "any")
```

## Arguments

- paths:

  Named list where names are path descriptions and values are paths

- must_exist:

  Logical; if TRUE, paths must exist

- type:

  Character vector matching paths; "file", "dir", or "any"

## Value

Logical; TRUE if all validations pass, error otherwise
