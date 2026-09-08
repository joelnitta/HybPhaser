# rhybphaser Workflow with targets

## Introduction

This vignette runs the same end-to-end analysis as
[`vignette("rhybphaser")`](https://joelnitta.github.io/rhybphaser/articles/rhybphaser.md),
but as a [`targets`](https://docs.ropensci.org/targets/) pipeline.
Compared with a plain script, `targets` gives you:

- **Caching** — re-running only re-computes steps whose inputs changed.
- **Reproducibility** — every result records the code and data that made
  it.
- **Inspection** — the dependency graph can be visualised
  ([`tar_visnetwork()`](https://docs.ropensci.org/targets/reference/tar_visnetwork.html))
  and any intermediate result read back with
  [`tar_read()`](https://docs.ropensci.org/targets/reference/tar_read.html).

The pipeline is self-contained: it downloads the official HybPiper test
dataset, assembles it, and generates the clade-association and phasing
input CSVs inside the pipeline, so there are no files to prepare.

## The pipeline

Put the following in `_targets.R` at the root of your project. Each
[`tar_target()`](https://docs.ropensci.org/targets/reference/tar_target.html)
wraps one `rhybphaser` step; the comments mark the parts of the
HybPhaser workflow.

``` r

library(targets)
library(tarchetypes)
library(rhybphaser)

tar_option_set(format = "rds")

tar_plan(

  ## Setup: download + assemble the official HybPiper test dataset --------
  tar_target(
    hp_data,
    run_hybpiper_test_dataset_clean(
      output_dir = "targets_hybpiper_test",
      samples = c("EG30", "EG98", "MWL2"),
      conda_env = "hybpiper_env",
      cpu = 2,
      mapper = "bwa"
    )
  ),
  tar_target(hybpiper_dir, hp_data$hybpiper_dir),
  tar_target(namelist, hp_data$namelist),
  tar_target(targets_file, hp_data$targets_file),

  ## Part 1: SNP assessment ---------------------------------------------
  tar_target(
    hybphaser_output_dir,
    run_consensus_sequences(
      hybpiper_dir = hybpiper_dir,
      output_dir = file.path(hybpiper_dir, "hybphaser_output"),
      namelist = namelist,
      threads = 2
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
      paralog_threshold = "outliers"
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

  ## Clade-association / phasing input CSVs ------------------------------
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
        hybphaser_output_dir, "04_clade_association",
        "clade_references_example.csv"
      )
    ),
    format = "file"
  ),
  tar_target(
    phasing_prep_csv,
    {
      # In a real analysis the phasing table is built by hand from the
      # clade-association results: one row per accession to phase, with its
      # parental references. write_phasing_prep_csv() writes a one-reference
      # template. Here, for illustration, phase the first sample against the
      # other two as if they were its parents.
      path <- file.path(
        hybphaser_output_dir, "05_phasing", "phasing_prep_example.csv"
      )
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(
        data.frame(
          sample = reference_samples[[1]],
          ref1 = reference_samples[[2]], abb1 = "P1",
          ref2 = reference_samples[[3]], abb2 = "P2"
        ),
        path,
        row.names = FALSE
      )
      path
    },
    format = "file"
  ),

  ## Part 2: clade association -----------------------------------------
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
        hybphaser_output_dir, "04_clade_association"
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
        hybphaser_output_dir, "04_clade_association"
      ),
      csv_file_with_clade_reference_names = clade_ref_csv,
      path_to_output_folder = hybphaser_output_dir
    )
  ),

  ## Part 3: phasing --------------------------------------------------
  tar_target(
    phasing_job,
    run_phasing(
      path_to_phasing_folder = file.path(hybphaser_output_dir, "05_phasing"),
      csv_file_with_phasing_prep_info = phasing_prep_csv,
      path_to_read_files_phasing = mapped_reads$output_dir,
      read_type_4phasing = "single-end",
      reference_sequence_folder = reference_dir
    )
  ),
  tar_target(
    phasing_stats,
    collate_phasing_stats(
      path_to_phasing_folder = file.path(hybphaser_output_dir, "05_phasing"),
      csv_file_with_phasing_prep_info = phasing_prep_csv
    )
  ),

  ## Part 1 again, on the phased reads --------------------------------
  tar_target(
    phased_hp_dir,
    run_phased_hybpiper(
      phased_reads_dir = phasing_job$folder_for_phased_reads,
      phased_hp_dir = file.path(hybpiper_dir, "hybpiper_output_phased"),
      targets_file = targets_file,
      conda_env = "hybpiper_env",
      cpu = 2
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
      threads = 2
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

  ## Part 4: merge phased + non-phased sequence lists ---------------
  tar_target(
    merged_sequence_lists,
    {
      sequence_lists        # force both branches before merging
      phased_sequence_lists
      merge_sequence_lists(
        path_to_sequence_lists_normal = file.path(
          hybphaser_output_dir, "03_sequence_lists", "loci_consensus"
        ),
        path_to_sequence_lists_phased = file.path(
          phased_hybphaser_output, "03_sequence_lists", "loci_consensus"
        ),
        path_of_sequence_lists_output = file.path(
          hybphaser_output_dir, "03_sequence_lists_merged", "loci_consensus"
        ),
        path_to_namelist_normal = namelist,
        path_to_namelist_phased = phased_namelist,
        exchange_phased_with_not_phased_samples = "yes",
        include_phased_seqlists_when_non_phased_locus_absent = "no"
      )
    }
  )
)
```

## Running the pipeline

From the project directory:

``` r

targets::tar_make()
```

`targets` works out the order from the dependencies between targets. On
a later run it skips any target that is already up to date.

## Inspecting the results

Any target can be read back with
[`tar_read()`](https://docs.ropensci.org/targets/reference/tar_read.html).
When the pipeline has run, each subsection below also shows its output
for the three samples of the HybPiper test dataset.

### Heterozygosity per locus

`tar_read(snp_results)$tab_snps` is a locus × sample matrix of the
proportion of polymorphic (ambiguity-coded) sites in each consensus
sequence. Higher values mean more within-sample allelic variation.

### Dataset summary

`tar_read(dataset_assessment)$summary_table` collapses that to one row
per sample. `locus_heterozygosity` is the percentage of loci carrying
any polymorphism and `allele_divergence` is the mean per-site divergence
weighted by sequence length. Samples elevated on **both** axes are the
candidates to carry forward to phasing; the useful cut-offs depend on
the dataset.

[`assess_dataset()`](https://joelnitta.github.io/rhybphaser/reference/assess_dataset.md)
also drops loci whose mean SNP proportion across samples is an outlier
(putative paralogs). `loci_removed_paralogs_all` lists them — empty here
because nothing crossed the threshold in this small example.

[`assess_dataset()`](https://joelnitta.github.io/rhybphaser/reference/assess_dataset.md)
also writes this as a plot (`02_assessment/3_LH_vs_AD.png`); the same
view of the summary table:

### Clade association

`tar_read(clade_table)$table_clade_association` gives, per sample, the
percentage of on-target reads that map unambiguously to each clade
reference (`C1`, `C2`, `C3` are the three samples used as references
here). A sample with strong support for two or more divergent references
is a candidate hybrid; the samples below each map mostly to their own
reference.

### Phasing

`tar_read(phasing_stats)$table_phasing_stats` reports, for each phased
accession, the proportion of reads assigned to each parental reference
(here `P1` and `P2`). The ratio is informative about ploidy: a diploid
F1 hybrid is near 1:1, whereas an allotetraploid can be closer to 3:1.

[`run_phasing()`](https://joelnitta.github.io/rhybphaser/reference/run_phasing.md)
also wrote one read file per haplotype (`<sample>_to_<ref>.fastq`),
which the pipeline reassembled with HybPiper and fed back through part
1.

### Merged sequence lists

`tar_read(merged_sequence_lists)` describes the combined per-locus FASTA
files, with the phased haplotypes substituted for their unphased sample,
ready for alignment and phylogenetics.

## Notes

- Set `RHYBPHASER_RUN_TARGETS_VIGNETTE=true` before rendering to execute
  the pipeline; otherwise only the `_targets.R` listing is shown.
- The pipeline requires a `hybpiper_env` conda environment, Docker (with
  the `joelnitta/rhybphaser` image), and network access for the test
  data.
- [`run_consensus_sequences()`](https://joelnitta.github.io/rhybphaser/reference/run_consensus_sequences.md)
  is used rather than
  [`run_generate_consensus_sequences()`](https://joelnitta.github.io/rhybphaser/reference/run_generate_consensus_sequences.md)
  because it returns the output path, which `targets` tracks naturally.
  \`\`\`
