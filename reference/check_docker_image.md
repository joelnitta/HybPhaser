# Check if the rhybphaser Docker Image Exists

Verifies that the rhybphaser Docker image is available locally

## Usage

``` r
check_docker_image(image = "joelnitta/rhybphaser:latest", pull = FALSE)
```

## Arguments

- image:

  Name of the Docker image. Default is "joelnitta/rhybphaser:latest"

- pull:

  Logical; if TRUE and image is not found, attempt to pull it

## Value

Logical; TRUE if image exists, FALSE otherwise

## Examples

``` r
if (FALSE) { # \dontrun{
check_docker_image("joelnitta/rhybphaser:latest")
} # }
```
