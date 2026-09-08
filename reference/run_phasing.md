# Run BBSplit Phasing

Generates an executable bash script for the BBSplit phasing step and
executes it.

## Usage

``` r
run_phasing(
  path_to_phasing_folder,
  csv_file_with_phasing_prep_info,
  path_to_read_files_phasing,
  read_type_4phasing = c("paired-end", "single-end"),
  ID_read_pair1 = "",
  ID_read_pair2 = "",
  reference_sequence_folder,
  folder_for_phased_reads = "",
  folder_for_phasing_stats = "",
  path_to_bbmap_executables = "",
  no_of_threads_phasing = 1,
  java_memory_usage_phasing = "",
  engine = c("docker", "local"),
  docker_image = "joelnitta/rhybphaser:latest",
  pull_image = FALSE
)
```

## Arguments

- path_to_phasing_folder:

  Path to phasing output folder.

- csv_file_with_phasing_prep_info:

  CSV with sample and reference information. Expected columns are
  `sample` or `samples`, then `ref1, abb1, ref2, abb2, ...`.

- path_to_read_files_phasing:

  Folder containing read files for phasing.

- read_type_4phasing:

  Read type: `"paired-end"` or `"single-end"`.

- ID_read_pair1:

  Identifier for read pair 1 (paired-end only).

- ID_read_pair2:

  Identifier for read pair 2 (paired-end only).

- reference_sequence_folder:

  Folder containing reference sequences.

- folder_for_phased_reads:

  Optional output folder for phased reads. If empty, uses
  `path_to_phasing_folder/phased_reads`.

- folder_for_phasing_stats:

  Optional output folder for phasing stats. If empty, uses
  `path_to_phasing_folder/phasing_stats`.

- path_to_bbmap_executables:

  Optional path to BBMap executables. If empty, `bbsplit.sh` is expected
  on `PATH`.

- no_of_threads_phasing:

  Number of threads for BBSplit. Use `0` or `"auto"` to omit thread
  argument.

- java_memory_usage_phasing:

  Optional Java memory (e.g., `"2G"`).

- engine:

  One of `"docker"` (default) or `"local"`. With `"docker"`, BBSplit
  runs in the `rhybphaser` container. With `"local"`, the generated
  script runs `bbsplit.sh` from `PATH` (or `path_to_bbmap_executables`);
  this needs a POSIX shell and is not available on native Windows.

- docker_image:

  Docker image used when `engine = "docker"`.

- pull_image:

  Logical; if `TRUE`, pull the Docker image when missing.

## Value

Invisibly returns a list with script path, commands, output folders,
selected samples, and run status.
