# Run BBSplit Phasing Using a HybPhaser Config File

Reads required phasing settings from `config.txt` and runs
[`run_phasing()`](https://joelnitta.github.io/rhybphaser/reference/run_phasing.md).

## Usage

``` r
run_phasing_from_config(
  config_file = "./config.txt",
  engine = c("docker", "local"),
  pull_image = FALSE
)
```

## Arguments

- config_file:

  Path to HybPhaser configuration file.

- engine:

  One of `"docker"` (default) or `"local"`; passed to
  [`run_phasing()`](https://joelnitta.github.io/rhybphaser/reference/run_phasing.md).

- pull_image:

  Logical; if `TRUE`, pull the Docker image when missing.

## Value

Invisibly returns the same object as
[`run_phasing()`](https://joelnitta.github.io/rhybphaser/reference/run_phasing.md).
