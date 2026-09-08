# Run HybPiper Test Dataset with Clean Output Directory

Removes an existing output directory before calling
[`run_hybpiper_test_dataset()`](https://joelnitta.github.io/rhybphaser/reference/run_hybpiper_test_dataset.md).
This prevents stale files from previous runs from triggering HybPiper
overwrite errors in `{targets}` workflows.

## Usage

``` r
run_hybpiper_test_dataset_clean(output_dir, ...)
```

## Arguments

- output_dir:

  Output directory passed to
  [`run_hybpiper_test_dataset()`](https://joelnitta.github.io/rhybphaser/reference/run_hybpiper_test_dataset.md).

- ...:

  Additional arguments passed to
  [`run_hybpiper_test_dataset()`](https://joelnitta.github.io/rhybphaser/reference/run_hybpiper_test_dataset.md).

## Value

List returned by
[`run_hybpiper_test_dataset()`](https://joelnitta.github.io/rhybphaser/reference/run_hybpiper_test_dataset.md).

## Examples

``` r
if (FALSE) { # \dontrun{
hp_data <- run_hybpiper_test_dataset_clean(
  output_dir = "targets_hybpiper_test",
  samples = c("EG30"),
  conda_env = "hybpiper_env",
  cpu = 1,
  mapper = "bwa"
)
} # }
```
