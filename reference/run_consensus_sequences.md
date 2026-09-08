# Generate Consensus Sequences and Return Output Directory

Wrapper around
[`run_generate_consensus_sequences()`](https://joelnitta.github.io/rhybphaser/reference/run_generate_consensus_sequences.md)
that returns `output_dir` on success. This is useful in `{targets}`
pipelines where returning a path can make dependencies easier to track.

## Usage

``` r
run_consensus_sequences(
  hybpiper_dir,
  output_dir,
  namelist = NULL,
  sample = NULL,
  intronerate = FALSE,
  cleanup = FALSE,
  threads = 1,
  min_depth = 10,
  min_allele_freq = 0.15,
  min_allele_count = 4,
  engine = c("docker", "local"),
  docker_image = "joelnitta/rhybphaser:latest",
  pull_image = FALSE
)
```

## Arguments

- hybpiper_dir:

  Path to HybPiper output directory (on host machine)

- output_dir:

  Path to output directory for HybPhaser results (on host machine). Will
  be created if it doesn't exist.

- namelist:

  Path to file containing sample names, one per line. If NULL, processes
  all samples in hybpiper_dir.

- sample:

  Single sample name to process. Ignored if namelist is provided.

- intronerate:

  Logical; if TRUE, use intronerated supercontigs

- cleanup:

  Logical; if TRUE, remove intermediate BAM and VCF files

- threads:

  Integer; number of CPU threads to use. Default is 1.

- min_depth:

  Integer; minimum read depth for variant calling. Default is 10.

- min_allele_freq:

  Numeric; minimum allele frequency (0-1) for calling heterozygous
  sites. Default is 0.15.

- min_allele_count:

  Integer; minimum number of reads supporting an allele. Default is 4.

- engine:

  One of `"docker"` (default) or `"local"`. `"local"` runs the bundled
  bash script directly and needs a POSIX shell (not available on native
  Windows).

- docker_image:

  Docker image name. Default is "joelnitta/rhybphaser:latest"

- pull_image:

  Logical; if TRUE, pull Docker image if not found

## Value

Character scalar. The `output_dir` path.

## Examples

``` r
if (FALSE) { # \dontrun{
out_dir <- run_consensus_sequences(
  hybpiper_dir = "path/to/hybpiper_output",
  output_dir = "path/to/hybphaser_output",
  namelist = "samples.txt",
  threads = 4
)
} # }
```
