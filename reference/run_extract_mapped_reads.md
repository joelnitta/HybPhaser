# Extract Mapped Reads

Extracts and concatenates reads that mapped to target sequences from
HybPiper or HybPhaser output, by running a bundled bash script. By
default the script runs inside the `rhybphaser` Docker container; set
`engine = "local"` to run it directly.

## Usage

``` r
run_extract_mapped_reads(
  base_dir,
  output_dir,
  namelist = NULL,
  remove_duplicates = FALSE,
  engine = c("docker", "local"),
  docker_image = "joelnitta/rhybphaser:latest",
  pull_image = FALSE
)
```

## Arguments

- base_dir:

  Path to base directory containing sample data (HybPiper or HybPhaser
  output)

- output_dir:

  Path to output directory for extracted reads

- namelist:

  Path to file containing sample names, one per line. If NULL, processes
  all samples found.

- remove_duplicates:

  Logical; if TRUE, remove duplicate sequences (not just duplicate
  names)

- engine:

  One of `"docker"` (default) or `"local"`. `"local"` runs the bundled
  bash script directly and needs a POSIX shell (not available on native
  Windows).

- docker_image:

  Docker image name. Default is "joelnitta/rhybphaser:latest"

- pull_image:

  Logical; if TRUE, pull Docker image if not found

## Value

Exit code (0 for success, non-zero for failure)

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract mapped reads from all samples
run_extract_mapped_reads(
  base_dir = "path/to/hybphaser_output",
  output_dir = "path/to/mapped_reads"
)

# Extract with duplicate removal
run_extract_mapped_reads(
  base_dir = "hybphaser_output",
  output_dir = "mapped_reads",
  remove_duplicates = TRUE,
  namelist = "samples_to_process.txt"
)
} # }
```
