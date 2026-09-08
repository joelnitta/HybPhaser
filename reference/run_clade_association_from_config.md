# Run BBSplit Clade Association Using a HybPhaser Config File

Reads clade-association variables from a HybPhaser config file and runs
BBSplit clade association.

## Usage

``` r
run_clade_association_from_config(
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
  [`run_clade_association()`](https://joelnitta.github.io/rhybphaser/reference/run_clade_association.md).

- pull_image:

  Logical; if `TRUE`, pull the Docker image when not available.

## Value

Invisibly returns a list with script path, commands, selected read
files, stats folder, and run status.
