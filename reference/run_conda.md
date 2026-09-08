# Run a Command in a Conda Environment

Runs a command via `conda run`, optionally within a named conda
environment.

## Usage

``` r
run_conda(
  command,
  args = character(0),
  env_name = NULL,
  conda_executable = "conda",
  wd = getwd(),
  ...
)
```

## Arguments

- command:

  Character; command to run (for example, "hybpiper")

- args:

  Character vector of arguments passed to command

- env_name:

  Optional conda environment name. If NULL, conda uses the default
  environment.

- conda_executable:

  Character; path or name of conda executable

- wd:

  Working directory for the command

- ...:

  Additional arguments passed to
  [`system2()`](https://rdrr.io/r/base/system2.html)

## Value

Exit status code from [`system2()`](https://rdrr.io/r/base/system2.html)
