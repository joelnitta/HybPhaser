# Create Test Namelist File

Creates a namelist file for testing

## Usage

``` r
create_test_namelist(file_path, samples = c("sample1", "sample2"))
```

## Arguments

- file_path:

  Path where namelist should be created

- samples:

  Character vector of sample names

## Value

Path to created file (invisibly)

## Examples

``` r
if (FALSE) { # \dontrun{
namelist <- create_test_namelist(
  file.path(tempdir(), "namelist.txt"),
  samples = c("sample1", "sample2")
)
} # }
```
