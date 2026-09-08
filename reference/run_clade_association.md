# Run BBSplit Clade Association

Generates an executable bash script to run BBSplit for clade association
and executes the generated commands from R.

## Usage

``` r
run_clade_association(
  path_to_clade_association_folder,
  csv_file_with_clade_reference_names,
  path_to_reference_sequences,
  path_to_read_files_cladeassociation,
  read_type_cladeassociation = c("single-end", "paired-end"),
  ID_read_pair1 = "",
  ID_read_pair2 = "",
  file_with_samples_included = "",
  path_to_bbmap = "",
  no_of_threads = 1,
  java_memory_usage_clade_association = "",
  engine = c("docker", "local"),
  docker_image = "joelnitta/rhybphaser:latest",
  pull_image = FALSE
)
```

## Arguments

- path_to_clade_association_folder:

  Path to clade association output folder.

- csv_file_with_clade_reference_names:

  CSV file with reference sample names and abbreviations (first two
  columns).

- path_to_reference_sequences:

  Folder with reference sequence files.

- path_to_read_files_cladeassociation:

  Folder with read files for clade association.

- read_type_cladeassociation:

  Read type: `"single-end"` or `"paired-end"`.

- ID_read_pair1:

  Identifier for read 1 files (paired-end only).

- ID_read_pair2:

  Identifier for read 2 files (paired-end only).

- file_with_samples_included:

  Optional text file with sample names to include. Use `""` or `"none"`
  to include all samples.

- path_to_bbmap:

  Optional path to BBMap binaries. If empty, `bbsplit.sh` is expected on
  `PATH`.

- no_of_threads:

  Number of threads for BBSplit. Use `0` or `"auto"` to omit thread
  argument.

- java_memory_usage_clade_association:

  Optional Java memory text for BBSplit (e.g., `"2G"`, `"512m"`).

- engine:

  One of `"docker"` (default) or `"local"`. With `"docker"`, BBSplit
  runs in the `rhybphaser` container. With `"local"`, the generated
  script runs `bbsplit.sh` from `PATH` (or `path_to_bbmap`); this needs
  a POSIX shell and is not available on native Windows.

- docker_image:

  Docker image used when `engine = "docker"`. Default is
  `"joelnitta/rhybphaser:latest"`.

- pull_image:

  Logical; if `TRUE`, pull the Docker image when not available.

## Value

Invisibly returns a list with script path, generated commands, selected
read files, stats folder, and run status.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_clade_association(
  path_to_clade_association_folder = "04_clade_association",
  csv_file_with_clade_reference_names = "clade_references.csv",
  path_to_reference_sequences = "03_sequence_lists/samples_consensus",
  path_to_read_files_cladeassociation = "mapped_reads",
  read_type_cladeassociation = "single-end"
)
} # }
```
