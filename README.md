# rhybphaser

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/joelnitta/rhybphaser/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/joelnitta/rhybphaser/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`rhybphaser` is an R package implementation of the **HybPhaser** method for
detecting and phasing hybrid and polyploid accessions in target-capture
datasets.

It is a reimplementation of the original
[HybPhaser](https://github.com/LarsNauheimer/HybPhaser) collection of R and
bash scripts (Nauheimer et al.) as an installable, tested R package. For the
scientific background, parameter guidance, and citation, see the **original
repository and paper** linked under [Related work](#related-work); this README
covers only how to use the R package.

## Installation

```r
# install.packages("pak")
pak::pak("joelnitta/rhybphaser")
```

### System requirements

- **R** ≥ 4.0
- **Docker** — consensus generation (BWA, SAMtools, BCFtools) and clade
  association / phasing (BBSplit) run inside the `joelnitta/rhybphaser` image.
  These functions pull and run it for you; pass `pull_image = TRUE` on first
  use. See [`DOCKER.md`](DOCKER.md). Not required if you run every step with
  `engine = "local"` and have the tools on `PATH`.
- **HybPiper in a conda environment** (default name `hybpiper_env`) — only
  needed for the functions that *run* HybPiper: `hybpiper_assemble()`,
  `hybpiper_stats()`, `run_hybpiper_test_dataset()`,
  `run_hybpiper_test_dataset_clean()`, `run_phased_hybpiper()`. These use conda
  directly, not Docker. If you already have HybPiper output, they are not
  needed.

## Workflow

`rhybphaser` starts from HybPiper output and produces phased, quality-filtered
sequence lists ready for phylogenetic analysis.

| Step | Function | Runs via |
|------|----------|----------|
| Generate consensus sequences | `run_generate_consensus_sequences()` | Docker |
| Count SNPs / recovered length | `count_snps()` | R |
| Assess dataset, filter samples & loci, flag paralogs | `assess_dataset()` | R |
| Generate per-locus / per-sample FASTA lists | `generate_sequence_lists()` | R |
| Extract on-target reads | `extract_mapped_reads()` | R |
| Clade association (BBSplit) | `run_clade_association()` → `collate_bbsplit_results()` | Docker |
| Phasing (BBSplit) | `run_phasing()` → `collate_phasing_stats()` | Docker |
| Merge phased & non-phased lists | `merge_sequence_lists()` | R |

Every function that shells out to an external tool takes
`engine = c("docker", "local")` and defaults to `"docker"` (the pinned
`rhybphaser` image). Pass `engine = "local"` to use `bwa`/`samtools`/`bcftools`/
`bbsplit.sh` from your `PATH` instead — for example on an HPC system without
Docker.

Each clade-association / phasing / merge function also has a `*_from_config()`
variant that reads its arguments from a `config.txt` file (`read_config()`), for
compatibility with the original script workflow.

An end-to-end example (including running HybPiper on the official test data) is
in the vignettes:

```r
vignette("rhybphaser")          # step-by-step
vignette("rhybphaser_targets")  # the same workflow as a {targets} pipeline
```

## Quick start

Starting from an existing HybPiper output directory:

```r
library(rhybphaser)

out <- "hybphaser_output"

# 1. Remap reads to contigs and call consensus sequences (Docker)
run_generate_consensus_sequences(
  hybpiper_dir = "hybpiper_output",
  output_dir   = out,
  namelist     = "namelist.txt",
  threads      = 4,
  pull_image   = TRUE
)

# 2. Quantify heterozygosity and recovered length
snps <- count_snps(
  path_to_output_folder   = out,
  fasta_file_with_targets = "targets.fasta",
  targets_file_format     = "DNA",
  path_to_namelist        = "namelist.txt"
)

# 3. Assess the dataset and apply quality filters
assessment <- assess_dataset(
  snp_table    = snps$tab_snps,
  length_table = snps$tab_length,
  targets_file = "targets.fasta",
  targets_type = "DNA",
  output_dir   = out,
  min_loci_per_sample_prop   = 0.5,
  min_samples_per_locus_prop = 0.5,
  paralog_threshold          = "outliers"
)

# 4. Write per-locus and per-sample FASTA lists
generate_sequence_lists(
  path_to_output_folder   = out,
  fasta_file_with_targets = "targets.fasta",
  targets_file_format     = "DNA",
  path_to_namelist        = "namelist.txt"
)
```

Continue with `extract_mapped_reads()`, `run_clade_association()`,
`run_phasing()` and `merge_sequence_lists()` as shown in `vignette("rhybphaser")`.

## Getting help

- Function reference: `?count_snps`, `?assess_dataset`,
  `?generate_sequence_lists`, `?run_clade_association`, `?run_phasing`,
  `?merge_sequence_lists`
- Bug reports and questions:
  <https://github.com/joelnitta/rhybphaser/issues>

## Related work

- Original method and code:
  [LarsNauheimer/HybPhaser](https://github.com/LarsNauheimer/HybPhaser)
- HybPhaser paper: Nauheimer et al.,
  [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.10.27.354589v2)
- Assembly pipeline: [HybPiper](https://github.com/mossmatters/HybPiper)

## License

GPL (≥ 3)
