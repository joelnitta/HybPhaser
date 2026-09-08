# Generate Consensus Sequences

Generates consensus sequences from HybPiper output by remapping reads to
contigs and calling variants, by running a bundled bash script. By
default the script runs inside the `rhybphaser` Docker container; set
`engine = "local"` to run it directly against tools on `PATH` (BWA,
SAMtools, BCFtools, BBMap, GNU parallel).

## Usage

``` r
run_generate_consensus_sequences(
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

Exit code (0 for success, non-zero for failure)

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate consensus sequences for all samples
run_generate_consensus_sequences(
  hybpiper_dir = "path/to/hybpiper_output",
  output_dir = "path/to/hybphaser_output",
  namelist = "samples.txt",
  threads = 4
)

# Process a single sample
run_generate_consensus_sequences(
  hybpiper_dir = "hybpiper_out",
  output_dir = "hybphaser_out",
  sample = "sample001"
)
} # }
```
