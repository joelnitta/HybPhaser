# Merge Phased and Non-phased Sequence Lists

Combines phased and non-phased sequence lists, optionally applies locus
and sample subsets, and can replace non-phased samples with phased
derivatives.

## Usage

``` r
merge_sequence_lists(
  path_to_sequence_lists_normal,
  path_to_sequence_lists_phased,
  path_of_sequence_lists_output,
  path_to_namelist_normal,
  path_to_namelist_phased,
  file_with_samples_included = "",
  file_with_samples_excluded = "",
  file_with_loci_excluded = "",
  file_with_loci_included = "",
  exchange_phased_with_not_phased_samples = c("yes", "no"),
  include_phased_seqlists_when_non_phased_locus_absent = c("no", "yes")
)
```

## Arguments

- path_to_sequence_lists_normal:

  Folder with non-phased sequence lists.

- path_to_sequence_lists_phased:

  Folder with phased sequence lists.

- path_of_sequence_lists_output:

  Output folder for merged sequence lists.

- path_to_namelist_normal:

  Path to namelist for non-phased samples.

- path_to_namelist_phased:

  Path to namelist for phased samples.

- file_with_samples_included:

  Optional file listing samples to include.

- file_with_samples_excluded:

  Optional file listing samples to exclude.

- file_with_loci_excluded:

  Optional file listing loci to exclude.

- file_with_loci_included:

  Optional file listing loci to include.

- exchange_phased_with_not_phased_samples:

  Whether to remove non-phased accessions when phased derivatives exist
  (`"yes"` or `"no"`).

- include_phased_seqlists_when_non_phased_locus_absent:

  Whether to copy loci present only in phased lists (`"yes"` or `"no"`).

## Value

Invisibly returns a list with output files and retained samples/loci.

## Examples

``` r
if (FALSE) { # \dontrun{
merge_sequence_lists(
  path_to_sequence_lists_normal = "03_sequence_lists/loci_consensus",
  path_to_sequence_lists_phased = "03_sequence_lists_phased/loci_consensus",
  path_of_sequence_lists_output = "03_sequence_lists_merged/loci_consensus",
  path_to_namelist_normal = "namelist_normal.txt",
  path_to_namelist_phased = "namelist_phased.txt"
)
} # }
```
