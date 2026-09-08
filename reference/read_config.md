# Read and Validate HybPhaser Configuration

Reads a HybPhaser `config.txt` file and validates required parameters.
The file is a list of `key = value` lines, where `value` is a quoted
string or a bare number (trailing `# comments` are ignored). Values are
taken literally, so Windows paths with backslashes are read correctly.

## Usage

``` r
read_config(config_file, required_vars = NULL)
```

## Arguments

- config_file:

  Path to configuration file

- required_vars:

  Character vector of required variable names

## Value

Named list of configuration parameters

## Examples

``` r
if (FALSE) { # \dontrun{
config <- read_config("config.txt")
} # }
```
