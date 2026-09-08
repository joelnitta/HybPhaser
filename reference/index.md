# Package index

## Analysis workflow

The main steps, in the order they are run. See
[`vignette("rhybphaser")`](https://joelnitta.github.io/rhybphaser/articles/rhybphaser.md)
for a worked example.

- [`run_generate_consensus_sequences()`](https://joelnitta.github.io/rhybphaser/reference/run_generate_consensus_sequences.md)
  : Generate Consensus Sequences
- [`count_snps()`](https://joelnitta.github.io/rhybphaser/reference/count_snps.md)
  : Count SNPs in Consensus Sequences
- [`assess_dataset()`](https://joelnitta.github.io/rhybphaser/reference/assess_dataset.md)
  : Assess HybPhaser dataset quality
- [`generate_sequence_lists()`](https://joelnitta.github.io/rhybphaser/reference/generate_sequence_lists.md)
  : Generate sequence lists for loci and samples
- [`generate_sequence_lists_after_assessment()`](https://joelnitta.github.io/rhybphaser/reference/generate_sequence_lists_after_assessment.md)
  : Generate Sequence Lists After Assessment
- [`get_reference_samples()`](https://joelnitta.github.io/rhybphaser/reference/get_reference_samples.md)
  : Get Reference Samples from Consensus Sequence Files
- [`extract_mapped_reads()`](https://joelnitta.github.io/rhybphaser/reference/extract_mapped_reads.md)
  : Extract Mapped Reads
- [`run_clade_association()`](https://joelnitta.github.io/rhybphaser/reference/run_clade_association.md)
  : Run BBSplit Clade Association
- [`collate_bbsplit_results()`](https://joelnitta.github.io/rhybphaser/reference/collate_bbsplit_results.md)
  : Collate BBSplit Results
- [`run_phasing()`](https://joelnitta.github.io/rhybphaser/reference/run_phasing.md)
  : Run BBSplit Phasing
- [`collate_phasing_stats()`](https://joelnitta.github.io/rhybphaser/reference/collate_phasing_stats.md)
  : Collate Phasing BBSplit Stats
- [`run_phased_hybpiper()`](https://joelnitta.github.io/rhybphaser/reference/run_phased_hybpiper.md)
  : Run HybPiper on Phased Read Files
- [`make_phased_namelist()`](https://joelnitta.github.io/rhybphaser/reference/make_phased_namelist.md)
  : Create Namelist from Phased Read Files
- [`merge_sequence_lists()`](https://joelnitta.github.io/rhybphaser/reference/merge_sequence_lists.md)
  : Merge Phased and Non-phased Sequence Lists

## Config-file entry points

`config.txt`-driven wrappers for compatibility with the original
HybPhaser script workflow.

- [`read_config()`](https://joelnitta.github.io/rhybphaser/reference/read_config.md)
  : Read and Validate HybPhaser Configuration
- [`run_clade_association_from_config()`](https://joelnitta.github.io/rhybphaser/reference/run_clade_association_from_config.md)
  : Run BBSplit Clade Association Using a HybPhaser Config File
- [`collate_bbsplit_results_from_config()`](https://joelnitta.github.io/rhybphaser/reference/collate_bbsplit_results_from_config.md)
  : Collate BBSplit Results Using a HybPhaser Config File
- [`run_phasing_from_config()`](https://joelnitta.github.io/rhybphaser/reference/run_phasing_from_config.md)
  : Run BBSplit Phasing Using a HybPhaser Config File
- [`collate_phasing_stats_from_config()`](https://joelnitta.github.io/rhybphaser/reference/collate_phasing_stats_from_config.md)
  : Collate Phasing BBSplit Stats Using a HybPhaser Config File
- [`merge_sequence_lists_from_config()`](https://joelnitta.github.io/rhybphaser/reference/merge_sequence_lists_from_config.md)
  : Merge Sequence Lists Using a HybPhaser Config File

## Running HybPiper

Run HybPiper assembly and stats through a conda environment.

- [`hybpiper_assemble()`](https://joelnitta.github.io/rhybphaser/reference/hybpiper_assemble.md)
  : Run HybPiper Assemble via Conda
- [`hybpiper_stats()`](https://joelnitta.github.io/rhybphaser/reference/hybpiper_stats.md)
  : Run HybPiper Stats via Conda
- [`run_hybpiper_test_dataset()`](https://joelnitta.github.io/rhybphaser/reference/run_hybpiper_test_dataset.md)
  : Run Official HybPiper Test Dataset via Conda
- [`run_hybpiper_test_dataset_clean()`](https://joelnitta.github.io/rhybphaser/reference/run_hybpiper_test_dataset_clean.md)
  : Run HybPiper Test Dataset with Clean Output Directory
- [`run_conda()`](https://joelnitta.github.io/rhybphaser/reference/run_conda.md)
  : Run a Command in a Conda Environment

## Environment checks

- [`check_conda()`](https://joelnitta.github.io/rhybphaser/reference/check_conda.md)
  : Check if Conda is Available
- [`check_hybpiper_conda()`](https://joelnitta.github.io/rhybphaser/reference/check_hybpiper_conda.md)
  : Check if HybPiper is Available in a Conda Environment
- [`check_docker()`](https://joelnitta.github.io/rhybphaser/reference/check_docker.md)
  : Check if Docker is Available
- [`check_docker_image()`](https://joelnitta.github.io/rhybphaser/reference/check_docker_image.md)
  : Check if the rhybphaser Docker Image Exists

## Input-file helpers

- [`write_clade_reference_csv()`](https://joelnitta.github.io/rhybphaser/reference/write_clade_reference_csv.md)
  : Write Clade Reference CSV
- [`write_phasing_prep_csv()`](https://joelnitta.github.io/rhybphaser/reference/write_phasing_prep_csv.md)
  : Write Phasing Prep CSV

## Bundled scripts and example data

- [`hybphaser_scripts()`](https://joelnitta.github.io/rhybphaser/reference/hybphaser_scripts.md)
  : Get Path to HybPhaser Bash Scripts
- [`hybphaser_example()`](https://joelnitta.github.io/rhybphaser/reference/hybphaser_example.md)
  : Get Path to HybPhaser Example Data
- [`run_hybphaser_script()`](https://joelnitta.github.io/rhybphaser/reference/run_hybphaser_script.md)
  : Run HybPhaser Bash Script
- [`run_consensus_sequences()`](https://joelnitta.github.io/rhybphaser/reference/run_consensus_sequences.md)
  : Generate Consensus Sequences and Return Output Directory
- [`run_extract_mapped_reads()`](https://joelnitta.github.io/rhybphaser/reference/run_extract_mapped_reads.md)
  : Extract Mapped Reads

## Test-data generators

- [`create_test_dataset()`](https://joelnitta.github.io/rhybphaser/reference/create_test_dataset.md)
  : Create Complete Test Dataset
- [`create_real_test_dataset()`](https://joelnitta.github.io/rhybphaser/reference/create_real_test_dataset.md)
  : Create Complete Test Dataset Using Real HybPiper Data
- [`create_test_hybpiper_output()`](https://joelnitta.github.io/rhybphaser/reference/create_test_hybpiper_output.md)
  : Create Test HybPiper Output
- [`create_hybpiper_test_output()`](https://joelnitta.github.io/rhybphaser/reference/create_hybpiper_test_output.md)
  : Create Test HybPiper Output Using Real Test Data
- [`create_test_namelist()`](https://joelnitta.github.io/rhybphaser/reference/create_test_namelist.md)
  : Create Test Namelist File
- [`create_test_targets()`](https://joelnitta.github.io/rhybphaser/reference/create_test_targets.md)
  : Create Test Target File
