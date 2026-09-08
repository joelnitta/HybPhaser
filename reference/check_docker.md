# Check if Docker is Available

Verifies that Docker is installed and running on the system

## Usage

``` r
check_docker(quiet = FALSE)
```

## Arguments

- quiet:

  Logical; if TRUE, suppresses messages

## Value

Logical; TRUE if Docker is available, FALSE otherwise

## Examples

``` r
if (FALSE) { # \dontrun{
if (check_docker()) {
  message("Docker is available!")
}
} # }
```
