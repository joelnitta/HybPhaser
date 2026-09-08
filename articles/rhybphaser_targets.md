# rhybphaser Workflow with targets

## Introduction

This vignette shows how to run the HybPhaser end-to-end workflow with
{targets}. The analysis steps are the same as the standard workflow
vignette, but dependencies are declared explicitly in a pipeline.

Using {targets} gives you:

- Reproducibility: each target records what code and data created it
- Caching: up-to-date steps are skipped automatically
- Visibility: the full pipeline graph can be inspected and debugged

## Inputs

This workflow is self-contained. It downloads and runs the official
HybPiper test dataset, then generates clade-association and phasing CSV
files inside the pipeline.

## Create Pipeline Script from the Vignette

[`targets::tar_script()`](https://docs.ropensci.org/targets/reference/tar_script.html)
writes a `_targets.R` file from R code. This lets you define and
maintain your pipeline directly in the vignette.

``` r

targets::tar_dir({
  targets::tar_script(
    {
      library(targets)
      library(tarchetypes)
      library(rhybphaser)

      tar_option_set(
        format = "rds"
      )

      tar_plan(
        # setup from official HybPiper test dataset ----
        tar_target(
          hp_data,
          run_hybpiper_test_dataset_clean(
            output_dir = "targets_hybpiper_test",
            samples = c("EG30"),
            conda_env = "hybpiper_env",
            cpu = 1,
            mapper = "bwa"
          )
        ),
        tar_target(hybpiper_dir, hp_data$hybpiper_dir),
        tar_target(namelist, hp_data$namelist),
        tar_target(targets_file, hp_data$targets_file),

        # part 1 normal run ----
        tar_target(
          hybphaser_output_dir,
          run_consensus_sequences(
            hybpiper_dir = hybpiper_dir,
            output_dir = file.path(hybpiper_dir, "hybphaser_output"),
            namelist = namelist,
            threads = 1
          )
        ),
        tar_target(
          snp_results,
          count_snps(
            path_to_output_folder = hybphaser_output_dir,
            fasta_file_with_targets = targets_file,
            targets_file_format = "DNA",
            path_to_namelist = namelist
          )
        ),
        tar_target(
          dataset_assessment,
          assess_dataset(
            snp_table = snp_results$tab_snps,
            length_table = snp_results$tab_length,
            targets_file = targets_file,
            targets_type = "DNA",
            output_dir = hybphaser_output_dir,
            min_loci_per_sample_prop = 0.5,
            min_samples_per_locus_prop = 0.5,
            paralog_threshold = "outliers",
            remove_sample_outliers = TRUE
          )
        ),
        tar_target(
          sequence_lists,
          generate_sequence_lists_after_assessment(
            assessment = dataset_assessment,
            path_to_output_folder = hybphaser_output_dir,
            fasta_file_with_targets = targets_file,
            targets_file_format = "DNA",
            path_to_namelist = namelist
          )
        ),

        # generate clade/phasing input tables ----
        tar_target(
          reference_dir,
          file.path(hybphaser_output_dir, "03_sequence_lists", "samples_consensus")
        ),
        tar_target(
          reference_samples,
          get_reference_samples(
            sequence_lists = sequence_lists,
            reference_dir = reference_dir
          )
        ),
        tar_target(
          clade_ref_csv,
          write_clade_reference_csv(
            reference_samples = reference_samples,
            output_file = file.path(
              hybphaser_output_dir,
              "04_clade_association",
              "clade_references_example.csv"
            )
          ),
          format = "file"
        ),
        tar_target(
          phasing_prep_csv,
          write_phasing_prep_csv(
            reference_samples = reference_samples,
            output_file = file.path(
              hybphaser_output_dir,
              "05_phasing",
              "phasing_prep_example.csv"
            )
          ),
          format = "file"
        ),

        # part 2 clade association ----
        tar_target(
          mapped_reads,
          extract_mapped_reads(
            base_dir = hybphaser_output_dir,
            output_dir = file.path(hybphaser_output_dir, "mapped_reads"),
            namelist = namelist,
            remove_duplicate_sequences = FALSE
          )
        ),
        tar_target(
          clade_job,
          run_clade_association(
            path_to_clade_association_folder = file.path(
              hybphaser_output_dir,
              "04_clade_association"
            ),
            csv_file_with_clade_reference_names = clade_ref_csv,
            path_to_reference_sequences = reference_dir,
            path_to_read_files_cladeassociation = mapped_reads$output_dir,
            read_type_cladeassociation = "single-end"
          )
        ),
        tar_target(
          clade_table,
          collate_bbsplit_results(
            path_to_clade_association_folder = file.path(
              hybphaser_output_dir,
              "04_clade_association"
            ),
            csv_file_with_clade_reference_names = clade_ref_csv,
            path_to_output_folder = hybphaser_output_dir
          )
        ),

        # part 3 phasing ----
        tar_target(
          phasing_job,
          run_phasing(
            path_to_phasing_folder = file.path(hybphaser_output_dir, "05_phasing"),
            csv_file_with_phasing_prep_info = phasing_prep_csv,
            path_to_read_files_phasing = mapped_reads$output_dir,
            read_type_4phasing = "single-end",
            reference_sequence_folder = file.path(
              hybphaser_output_dir,
              "03_sequence_lists",
              "samples_consensus"
            )
          )
        ),
        tar_target(
          phasing_stats,
          collate_phasing_stats(
            path_to_phasing_folder = file.path(hybphaser_output_dir, "05_phasing"),
            csv_file_with_phasing_prep_info = phasing_prep_csv
          )
        ),

        # part 1 repeated on phased data ----
        tar_target(
          phased_hp_dir,
          run_phased_hybpiper(
            phased_reads_dir = phasing_job$folder_for_phased_reads,
            phased_hp_dir = file.path(hybpiper_dir, "hybpiper_output_phased"),
            targets_file = targets_file,
            conda_env = "hybpiper_env",
            cpu = 1
          )
        ),
        tar_target(
          phased_namelist,
          make_phased_namelist(
            phased_reads_dir = phasing_job$folder_for_phased_reads,
            phased_hp_dir = phased_hp_dir
          )
        ),
        tar_target(
          phased_hybphaser_output,
          run_consensus_sequences(
            hybpiper_dir = phased_hp_dir,
            output_dir = file.path(phased_hp_dir, "hybphaser_output_phased"),
            namelist = phased_namelist,
            threads = 1
          )
        ),
        tar_target(
          phased_snps,
          count_snps(
            path_to_output_folder = phased_hybphaser_output,
            fasta_file_with_targets = targets_file,
            targets_file_format = "DNA",
            path_to_namelist = phased_namelist
          )
        ),
        tar_target(
          phased_assessment,
          assess_dataset(
            snp_table = phased_snps$tab_snps,
            length_table = phased_snps$tab_length,
            targets_file = targets_file,
            targets_type = "DNA",
            output_dir = phased_hybphaser_output,
            min_loci_per_sample_prop = 0.5,
            min_samples_per_locus_prop = 0.5
          )
        ),
        tar_target(
          phased_sequence_lists,
          generate_sequence_lists_after_assessment(
            assessment = phased_assessment,
            path_to_output_folder = phased_hybphaser_output,
            fasta_file_with_targets = targets_file,
            targets_file_format = "DNA",
            path_to_namelist = phased_namelist
          )
        ),

        # part 4 merge sequence lists ----
        tar_target(
          merged_sequence_lists,
          {
            # enforce that both normal and phased sequence lists are generated
            sequence_lists
            phased_sequence_lists
            merge_sequence_lists(
              path_to_sequence_lists_normal = file.path(
                hybphaser_output_dir,
                "03_sequence_lists",
                "loci_consensus"
              ),
              path_to_sequence_lists_phased = file.path(
                phased_hybphaser_output,
                "03_sequence_lists",
                "loci_consensus"
              ),
              path_of_sequence_lists_output = file.path(
                hybphaser_output_dir,
                "03_sequence_lists_merged",
                "loci_consensus"
              ),
              path_to_namelist_normal = namelist,
              path_to_namelist_phased = phased_namelist,
              exchange_phased_with_not_phased_samples = "yes",
              include_phased_seqlists_when_non_phased_locus_absent = "no"
            )
          }
        )
      )
    },
    script = "_targets.R",
    ask = FALSE
  )
  targets::tar_make()
})
```

## Notes

- By default, this vignette does not execute during package checks.
- To run it end-to-end while rendering, set environment variable
  `RHYBPHASER_RUN_TARGETS_VIGNETTE=true` before knitting.
- The pipeline auto-generates all required input files and does not
  depend on pre-existing files in `data/`.
- [`run_consensus_sequences()`](https://joelnitta.github.io/rhybphaser/reference/run_consensus_sequences.md)
  is used instead of
  [`run_generate_consensus_sequences()`](https://joelnitta.github.io/rhybphaser/reference/run_generate_consensus_sequences.md)
  because returning a path is more natural for
  [targets](https://docs.ropensci.org/targets/) dependency tracking.
- [`run_phased_hybpiper()`](https://joelnitta.github.io/rhybphaser/reference/run_phased_hybpiper.md)
  and
  [`make_phased_namelist()`](https://joelnitta.github.io/rhybphaser/reference/make_phased_namelist.md)
  are exported package functions, so no custom helper file is required.
- [`write_clade_reference_csv()`](https://joelnitta.github.io/rhybphaser/reference/write_clade_reference_csv.md)
  and
  [`write_phasing_prep_csv()`](https://joelnitta.github.io/rhybphaser/reference/write_phasing_prep_csv.md)
  are exported helper functions used to generate CSV inputs inside the
  pipeline.
- Docker and conda availability are still required when you actually
  execute the pipeline.
