# Compute the relative path from a base to a target directory

Returns the relative portion of `path` with respect to `base`, using
forward slashes, if `path` is a strict descendant of `base`. Returns
`NULL` when `path` is the same as `base` or is not nested inside it.

## Usage

``` r
.relative_path(path, base)
```

## Arguments

- path:

  Character; absolute path of the target.

- base:

  Character; absolute path of the base directory.

## Value

Character relative path (no leading `/`), or `NULL`.
