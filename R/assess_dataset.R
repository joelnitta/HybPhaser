#' Assess HybPhaser dataset quality
#'
#' Generate tables and graphs to assess sequence recovery, missing data,
#' and identify putative paralogs based on SNP proportions. This function
#' applies quality filters and produces comprehensive assessment reports.
#'
#' @param snp_table Path to RDS file containing SNP proportions table
#'   (output from [count_snps()]), or a matrix/data.frame with loci as rows
#'   and samples as columns.
#' @param length_table Path to RDS file containing sequence length table
#'   (output from [count_snps()]), or a matrix/data.frame with loci as rows
#'   and samples as columns.
#' @param targets_file Path to FASTA file with target sequences.
#' @param targets_type Character string, either "DNA" or "AA" (amino acid).
#'   If "AA", lengths will be multiplied by 3 to get nucleotide length.
#' @param output_dir Path to output directory where assessment results will
#'   be saved.
#' @param subset_name Optional name for dataset optimization subset. If
#'   provided, creates subfolder "02_assessment_<subset_name>". If empty
#'   string (default), creates "02_assessment" folder.
#' @param min_loci_per_sample_prop Minimum proportion of loci that must be
#'   recovered per sample (default: 0.5). Samples below this threshold are
#'   removed.
#' @param min_target_length_per_sample_prop Minimum proportion of target
#'   sequence length that must be recovered per sample (default: 0.5).
#'   Samples below this threshold are removed.
#' @param min_samples_per_locus_prop Minimum proportion of samples that must
#'   recover a locus (default: 0.5). Loci below this threshold are removed.
#' @param min_target_length_per_locus_prop Minimum proportion of target
#'   sequence length that must be recovered per locus (default: 0.5).
#'   Loci below this threshold are removed.
#' @param paralog_threshold Threshold for removing putative paralog loci
#'   across all samples. Can be:
#'   - Numeric: remove loci with mean SNP proportion above this value
#'   - "outliers": use 1.5*IQR method to identify outliers
#'   - "file": read list of loci to remove from `paralog_file`
#'   - "none" or NA: skip paralog removal (default)
#' @param paralog_file Path to text file with list of paralog loci to remove
#'   (one per line). Only used if `paralog_threshold = "file"`.
#' @param remove_sample_outliers Logical. If TRUE, removes outlier loci for
#'   each sample individually using 1.5*IQR method (default: FALSE).
#'
#' @return Invisibly returns a list with components:
#'   - `snp_table_cleaned`: Cleaned SNP table after filtering
#'   - `length_table_cleaned`: Cleaned length table after filtering
#'   - `samples_removed`: Character vector of removed sample names
#'   - `loci_removed_missing`: Character vector of loci removed for missing
#'     data
#'   - `loci_removed_paralogs_all`: Character vector of loci removed as
#'     paralogs across all samples
#'   - `loci_removed_paralogs_each`: Named list of loci removed per sample
#'   - `summary_table`: Data frame with heterozygosity and allele divergence
#'     metrics per sample
#'
#' @details
#' The function performs three main assessment steps:
#'
#' **Step 1: Missing Data Assessment**
#' - Calculates proportion of loci recovered per sample
#' - Calculates proportion of samples recovering each locus
#' - Calculates proportion of target sequence length recovered
#' - Removes samples and loci failing to meet minimum thresholds
#' - Generates plots and tables showing data recovery
#'
#' **Step 2a: Paralog Detection (All Samples)**
#' - Identifies loci with unusually high SNP proportions across all samples
#' - Applies specified threshold method
#' - Generates barplots showing mean SNP proportions per locus
#'
#' **Step 2b: Paralog Detection (Per Sample)**
#' - Optionally identifies outlier loci for each sample individually
#' - Uses 1.5*IQR method on loci with SNPs
#' - Generates boxplots showing SNP distributions per sample
#'
#' **Step 3: Heterozygosity and Allele Divergence**
#' - Calculates locus heterozygosity (proportion of loci with SNPs)
#' - Calculates allele divergence (weighted mean SNP proportion)
#' - Generates summary table and scatter plots
#'
#' Output files are saved to `<output_dir>/02_assessment[_<subset_name>]/`:
#' - `1_Data_recovered_overview.pdf/png`: Missing data plots
#' - `1_Data_recovered_per_sample.csv`: Recovery metrics per sample
#' - `1_Data_recovered_per_locus.csv`: Recovery metrics per locus
#' - `1_Summary_missing_data.txt`: Text summary of removed samples/loci
#' - `2a_Paralogs_for_all_samples.pdf/png`: Paralog detection plots
#' - `2a_List_of_paralogs_removed_for_all_samples.txt`: List of removed loci
#' - `2b_Paralogs_for_each_sample.pdf/png`: Per-sample outlier plots
#' - `2_Summary_Paralogs.txt`: Text summary of paralog removal
#' - `3_LH_vs_AD.pdf/png`: Locus heterozygosity vs allele divergence plot
#' - `3_varLH_vs_AD.pdf/png`: Multiple heterozygosity threshold plots
#' - `4_Summary_table.csv`: Final metrics per sample
#' - `0_Table_SNPs.csv`: Cleaned SNP table
#' - `0_Table_consensus_length.csv`: Cleaned length table
#' - `0_namelist_included_samples.txt`: List of retained samples
#'
#' R objects are also saved to `<output_dir>/00_R_objects/<subset_name>/`:
#' - `Table_SNPs_cleaned.Rds`
#' - `Table_consensus_length_cleaned.Rds`
#' - `Summary_table.Rds`
#' - Various vectors of removed samples/loci
#'
#' @examples
#' \dontrun{
#' # Run count_snps first
#' snp_results <- count_snps(
#'   consensus_dir = "01_seqs",
#'   output_dir = "output"
#' )
#'
#' # Then assess dataset
#' assess_dataset(
#'   snp_table = file.path("output", "00_R_objects", "Table_SNPs.Rds"),
#'   length_table = file.path("output", "00_R_objects",
#'     "Table_consensus_length.Rds"),
#'   targets_file = "targets.fasta",
#'   targets_type = "DNA",
#'   output_dir = "output",
#'   min_loci_per_sample_prop = 0.5,
#'   min_samples_per_locus_prop = 0.5,
#'   paralog_threshold = "outliers",
#'   remove_sample_outliers = TRUE
#' )
#' }
#'
#' @export
assess_dataset <- function(
  snp_table,
  length_table,
  targets_file,
  targets_type = c("DNA", "AA"),
  output_dir,
  subset_name = "",
  min_loci_per_sample_prop = 0.5,
  min_target_length_per_sample_prop = 0.5,
  min_samples_per_locus_prop = 0.5,
  min_target_length_per_locus_prop = 0.5,
  paralog_threshold = "none",
  paralog_file = NULL,
  remove_sample_outliers = FALSE
) {
  # Input validation
  targets_type <- match.arg(targets_type)

  if (!file.exists(targets_file)) {
    stop("targets_file does not exist: ", targets_file)
  }

  # Load SNP and length tables
  if (is.character(snp_table) && length(snp_table) == 1) {
    if (!file.exists(snp_table)) {
      stop("snp_table file does not exist: ", snp_table)
    }
    tab_snps <- readRDS(snp_table)
  } else {
    tab_snps <- snp_table
  }

  if (is.character(length_table) && length(length_table) == 1) {
    if (!file.exists(length_table)) {
      stop("length_table file does not exist: ", length_table)
    }
    tab_length <- readRDS(length_table)
  } else {
    tab_length <- length_table
  }

  # Convert to numeric matrix
  # Use data.matrix() to ensure proper conversion from data frame to numeric
  # matrix (as.matrix() can produce character matrix if any non-numeric data)
  if (is.data.frame(tab_snps)) {
    tab_snps <- data.matrix(tab_snps)
  } else {
    tab_snps <- as.matrix(tab_snps)
    storage.mode(tab_snps) <- "numeric"
  }

  if (is.data.frame(tab_length)) {
    tab_length <- data.matrix(tab_length)
  } else {
    tab_length <- as.matrix(tab_length)
    storage.mode(tab_length) <- "numeric"
  }

  # Create output directories
  folder_subset_add <- if (subset_name != "") {
    paste0("_", subset_name)
  } else {
    ""
  }

  output_assess <- file.path(
    output_dir,
    paste0("02_assessment", folder_subset_add)
  )
  output_robjects <- file.path(output_dir, "00_R_objects", subset_name)

  dir.create(output_assess, showWarnings = FALSE, recursive = TRUE)
  dir.create(output_robjects, showWarnings = FALSE, recursive = TRUE)

  # Transpose for loci as columns
  loci <- t(tab_snps)

  nloci <- ncol(loci)
  nsamples <- nrow(loci)

  # Identify completely failed loci and samples
  failed_loci <- which(colSums(is.na(loci)) == nsamples)
  failed_samples <- which(colSums(is.na(tab_snps)) == nrow(tab_snps))

  # Step 1: Missing Data Assessment ----

  # Sequences per locus
  seq_per_locus <- colSums(!is.na(loci))
  names(seq_per_locus) <- colnames(loci)
  seq_per_locus_prop <- seq_per_locus / nsamples

  # Sequences per sample
  seq_per_sample <- colSums(!is.na(tab_snps))
  names(seq_per_sample) <- colnames(tab_snps)
  seq_per_sample_prop <- seq_per_sample / nloci

  # Target sequence lengths
  if (targets_type == "AA") {
    targets_length_all <- lengths(ape::read.FASTA(targets_file, type = "AA")) *
      3
  } else {
    targets_length_all <- lengths(ape::read.FASTA(targets_file))
  }

  # Get max target length per gene
  gene_names <- unique(gsub(
    ".*-",
    "",
    gsub(" .*", "", names(targets_length_all))
  ))
  max_target_length <- vapply(
    gene_names,
    function(gene) {
      max(targets_length_all[grep(
        paste0("\\b", gene, "\\b"),
        names(targets_length_all)
      )])
    },
    numeric(1)
  )

  # Calculate proportion of target length recovered
  comb_target_length <- sum(max_target_length)
  comb_seq_length_samples <- colSums(tab_length, na.rm = TRUE)
  prop_target_length_per_sample <- comb_seq_length_samples /
    comb_target_length
  mean_seq_length_loci <- rowMeans(tab_length, na.rm = TRUE)
  # Ensure names match between mean_seq_length_loci and max_target_length
  # Use row names from data as the source of truth
  loci_names <- rownames(tab_length)
  # Create a mapping from loci names to max_target_length
  target_lengths_for_loci <- vapply(
    loci_names,
    function(locus) {
      if (locus %in% names(max_target_length)) {
        max_target_length[locus]
      } else {
        # Try extracting gene name
        gene_part <- gsub(".*-", "", locus)
        if (gene_part %in% names(max_target_length)) {
          max_target_length[gene_part]
        } else {
          NA_real_
        }
      }
    },
    numeric(1)
  )
  prop_target_length_per_locus <- mean_seq_length_loci /
    target_lengths_for_loci
  names(prop_target_length_per_locus) <- loci_names

  # Apply thresholds to identify samples and loci to remove
  outsamples_missing_loci <- seq_per_sample_prop[
    seq_per_sample_prop < min_loci_per_sample_prop
  ]
  outsamples_missing_target <- prop_target_length_per_sample[
    prop_target_length_per_sample < min_target_length_per_sample_prop
  ]
  outsamples_missing <- unique(names(c(
    outsamples_missing_loci,
    outsamples_missing_target
  )))

  outloci_missing_samples <- seq_per_locus_prop[
    seq_per_locus_prop < min_samples_per_locus_prop
  ]
  outloci_missing_target <- prop_target_length_per_locus[
    prop_target_length_per_locus < min_target_length_per_locus_prop
  ]
  outloci_missing <- unique(names(c(
    outloci_missing_samples,
    outloci_missing_target
  )))

  # Remove bad samples and loci
  tab_snps_cl1 <- tab_snps

  if (length(outsamples_missing) > 0) {
    tab_snps_cl1 <- tab_snps_cl1[,
      !(colnames(tab_snps) %in%
        outsamples_missing),
      drop = FALSE
    ]
  }

  if (length(outloci_missing) > 0) {
    tab_snps_cl1 <- tab_snps_cl1[
      !(rownames(tab_snps) %in%
        outloci_missing),
      ,
      drop = FALSE
    ]
  }

  loci_cl1 <- t(tab_snps_cl1)
  nloci_cl1 <- ncol(loci_cl1)
  nsamples_cl1 <- nrow(loci_cl1)

  # Generate missing data plots
  .generate_missing_data_plots(
    seq_per_sample_prop = seq_per_sample_prop,
    prop_target_length_per_sample = prop_target_length_per_sample,
    seq_per_locus_prop = seq_per_locus_prop,
    prop_target_length_per_locus = prop_target_length_per_locus,
    nloci = nloci,
    nsamples = nsamples,
    min_loci_per_sample_prop = min_loci_per_sample_prop,
    min_target_length_per_sample_prop = min_target_length_per_sample_prop,
    min_samples_per_locus_prop = min_samples_per_locus_prop,
    min_target_length_per_locus_prop = min_target_length_per_locus_prop,
    outsamples_missing_loci = outsamples_missing_loci,
    outsamples_missing_target = outsamples_missing_target,
    outloci_missing_samples = outloci_missing_samples,
    outloci_missing_target = outloci_missing_target,
    output_assess = output_assess
  )

  # Save missing data tables
  tab_seq_per_sample <- cbind(
    seq_per_sample,
    round(seq_per_sample_prop, 3),
    round(prop_target_length_per_sample, 3)
  )
  colnames(tab_seq_per_sample) <- c(
    "No. loci",
    "Prop. of loci",
    "Prop. of target length"
  )
  write.csv(
    tab_seq_per_sample,
    file.path(output_assess, "1_Data_recovered_per_sample.csv")
  )

  tab_seq_per_locus <- cbind(
    seq_per_locus,
    round(seq_per_locus_prop, 3),
    round(prop_target_length_per_locus, 3)
  )
  colnames(tab_seq_per_locus) <- c(
    "No. samples",
    "Prop.of samples",
    "Prop. of target length"
  )
  write.csv(
    tab_seq_per_locus,
    file.path(output_assess, "1_Data_recovered_per_locus.csv")
  )

  # Save summary text file
  summary_file <- file.path(output_assess, "1_Summary_missing_data.txt")
  cat(
    file = summary_file,
    append = FALSE,
    "Dataset optimisation: Samples and loci removed to reduce missing ",
    "data\n"
  )

  cat(
    file = summary_file,
    append = TRUE,
    "\n",
    length(failed_samples),
    " samples failed completely:\n",
    paste(names(failed_samples)),
    "\n",
    sep = ""
  )
  cat(
    file = summary_file,
    append = TRUE,
    "\n",
    length(outsamples_missing_loci),
    " samples are below the threshold (",
    min_loci_per_sample_prop,
    ") for proportion of recovered loci:\n",
    paste(
      names(outsamples_missing_loci),
      "\t",
      round(outsamples_missing_loci, 3),
      "\n"
    ),
    sep = ""
  )
  cat(
    file = summary_file,
    append = TRUE,
    "\n",
    length(outsamples_missing_target),
    " samples are below the threshold (",
    min_target_length_per_sample_prop,
    ") for recovered target sequence length\n",
    paste(
      names(outsamples_missing_target),
      "\t",
      round(outsamples_missing_target, 3),
      "\n"
    ),
    sep = ""
  )
  cat(
    file = summary_file,
    append = TRUE,
    "\nIn total ",
    length(outsamples_missing),
    " samples were removed:\n",
    paste(outsamples_missing, "\n"),
    sep = ""
  )

  cat(
    file = summary_file,
    append = TRUE,
    "\n",
    length(failed_loci),
    " loci failed completely:\n",
    paste(names(failed_loci)),
    "\n",
    sep = ""
  )
  cat(
    file = summary_file,
    append = TRUE,
    "\n",
    length(outloci_missing_samples),
    " loci are below the threshold (",
    min_samples_per_locus_prop,
    ") for proportion of recovered samples:\n",
    paste(
      names(outloci_missing_samples),
      "\t",
      round(outloci_missing_samples, 3),
      "\n"
    ),
    sep = ""
  )
  cat(
    file = summary_file,
    append = TRUE,
    "\n",
    length(outloci_missing_target),
    " loci are below the threshold (",
    min_target_length_per_locus_prop,
    ") for proportion of recovered target sequence length:\n",
    paste(
      names(outloci_missing_target),
      "\t",
      round(outloci_missing_target, 3),
      "\n"
    ),
    sep = ""
  )
  cat(
    file = summary_file,
    append = TRUE,
    "\nIn total ",
    length(outloci_missing),
    " loci were removed:\n",
    paste(outloci_missing, "\n"),
    sep = ""
  )

  # Step 2a: Paralog Detection (All Samples) ----

  loci_cl1_colmeans <- colMeans(loci_cl1, na.rm = TRUE)

  # Apply paralog threshold
  paralog_result <- .apply_paralog_threshold(
    loci_cl1_colmeans = loci_cl1_colmeans,
    paralog_threshold = paralog_threshold,
    paralog_file = paralog_file
  )

  outloci_para_all <- paralog_result$outloci_para_all
  outloci_para_all_values <- paralog_result$outloci_para_all_values
  threshold_value <- paralog_result$threshold_value

  # Generate paralog plots (all samples)
  .generate_paralog_all_plots(
    loci_cl1_colmeans = loci_cl1_colmeans,
    outloci_para_all = outloci_para_all,
    nsamples_cl1 = nsamples_cl1,
    nloci_cl1 = nloci_cl1,
    threshold_value = threshold_value,
    paralog_threshold = paralog_threshold,
    output_assess = output_assess
  )

  # Remove paralogs from table
  tab_snps_cl2a <- if (length(outloci_para_all) == 0) {
    tab_snps_cl1
  } else {
    tab_snps_cl1[
      !(rownames(tab_snps_cl1) %in% outloci_para_all),
      ,
      drop = FALSE
    ]
  }

  # Step 2b: Paralog Detection (Per Sample) ----

  # Create table without zeros (only loci with SNPs)
  tab_snps_cl2a_nozero <- tab_snps_cl2a
  tab_snps_cl2a_nozero[tab_snps_cl2a_nozero == 0] <- NA
  tab_snps_cl2b_nozero <- tab_snps_cl2a_nozero
  tab_snps_cl2b <- tab_snps_cl2a

  outloci_para_each <- vector("list", ncol(tab_snps_cl2a))
  names(outloci_para_each) <- colnames(tab_snps_cl2a)
  threshold_para_each <- outloci_para_each

  if (remove_sample_outliers) {
    for (i in seq_along(colnames(tab_snps_cl2a))) {
      threshold_i <- 1.5 *
        IQR(tab_snps_cl2a_nozero[, i], na.rm = TRUE) +
        quantile(tab_snps_cl2a_nozero[, i], na.rm = TRUE)[[4]]
      outlier_loci_i <- tab_snps_cl2a_nozero[
        which(tab_snps_cl2a_nozero[, i] > threshold_i),
        i
      ]
      outloci_para_each[[i]] <- outlier_loci_i
      threshold_para_each[[i]] <- threshold_i

      # Set outliers to NA
      if (length(outlier_loci_i) > 0) {
        outlier_names <- names(outlier_loci_i)
        tab_snps_cl2b[outlier_names, i] <- NA
        tab_snps_cl2b_nozero[outlier_names, i] <- NA
      }
    }
    outliers_color <- "red"
  } else {
    outliers_color <- "black"
  }

  # Get cleaned length table
  tab_length_cl2b <- tab_length[
    rownames(tab_length) %in% rownames(tab_snps_cl2b),
    colnames(tab_length) %in% colnames(tab_snps_cl2b),
    drop = FALSE
  ]

  # Generate per-sample paralog plots
  .generate_paralog_each_plots(
    tab_snps_cl2a_nozero = tab_snps_cl2a_nozero,
    outliers_color = outliers_color,
    output_assess = output_assess
  )

  # Write paralog summary
  .write_paralog_summary(
    output_assess = output_assess,
    paralog_threshold = paralog_threshold,
    paralog_file = paralog_file,
    threshold_value = threshold_value,
    outloci_para_all = outloci_para_all,
    outloci_para_all_values = outloci_para_all_values,
    remove_sample_outliers = remove_sample_outliers,
    outloci_para_each = outloci_para_each,
    threshold_para_each = threshold_para_each
  )

  # Save cleaned tables
  write.csv(tab_snps_cl2b, file = file.path(output_assess, "0_Table_SNPs.csv"))
  write.csv(
    tab_length_cl2b,
    file = file.path(output_assess, "0_Table_consensus_length.csv")
  )

  # Save list of included samples
  included_samples <- rownames(loci)[
    !(rownames(loci) %in%
      outsamples_missing)
  ]
  writeLines(
    included_samples,
    file.path(output_assess, "0_namelist_included_samples.txt")
  )

  # Step 3: Heterozygosity and Allele Divergence ----

  nloci_cl2 <- nrow(tab_snps_cl2b)
  targets_length_cl2b <- sum(max_target_length[
    gsub(".*-", "", names(max_target_length)) %in% rownames(tab_snps_cl2b)
  ])

  # Generate summary table
  summary_table <- .generate_summary_table(
    tab_snps_cl2b = tab_snps_cl2b,
    tab_length_cl2b = tab_length_cl2b,
    targets_length_cl2b = targets_length_cl2b,
    outloci_para_all = outloci_para_all,
    outloci_para_each = outloci_para_each,
    nloci_cl1 = nloci_cl1,
    nloci_cl2 = nloci_cl2
  )

  write.csv(
    summary_table,
    file = file.path(output_assess, "4_Summary_table.csv"),
    row.names = FALSE
  )

  # Generate heterozygosity plots
  .generate_heterozygosity_plots(
    summary_table = summary_table,
    output_assess = output_assess
  )

  # Save R objects
  saveRDS(
    tab_snps_cl2b,
    file = file.path(output_robjects, "Table_SNPs_cleaned.Rds")
  )
  saveRDS(
    tab_length_cl2b,
    file = file.path(output_robjects, "Table_consensus_length_cleaned.Rds")
  )
  saveRDS(summary_table, file = file.path(output_robjects, "Summary_table.Rds"))
  saveRDS(
    outloci_missing,
    file = file.path(output_robjects, "outloci_missing.Rds")
  )
  saveRDS(
    outsamples_missing,
    file = file.path(output_robjects, "outsamples_missing.Rds")
  )
  saveRDS(
    outloci_para_all,
    file = file.path(output_robjects, "outloci_para_all.Rds")
  )
  saveRDS(
    outloci_para_each,
    file = file.path(output_robjects, "outloci_para_each.Rds")
  )

  # Return results
  invisible(list(
    snp_table_cleaned = tab_snps_cl2b,
    length_table_cleaned = tab_length_cl2b,
    samples_removed = outsamples_missing,
    loci_removed_missing = outloci_missing,
    loci_removed_paralogs_all = outloci_para_all,
    loci_removed_paralogs_each = outloci_para_each,
    summary_table = summary_table
  ))
}

# Helper functions ----

#' Generate missing data assessment plots
#' @noRd
.generate_missing_data_plots <- function(
  seq_per_sample_prop,
  prop_target_length_per_sample,
  seq_per_locus_prop,
  prop_target_length_per_locus,
  nloci,
  nsamples,
  min_loci_per_sample_prop,
  min_target_length_per_sample_prop,
  min_samples_per_locus_prop,
  min_target_length_per_locus_prop,
  outsamples_missing_loci,
  outsamples_missing_target,
  outloci_missing_samples,
  outloci_missing_target,
  output_assess
) {
  for (i in 1:2) {
    if (i == 1) {
      pdf(
        file = file.path(output_assess, "1_Data_recovered_overview.pdf"),
        width = 11,
        height = 7
      )
    } else {
      png(
        filename = file.path(output_assess, "1_Data_recovered_overview.png"),
        width = 1400,
        height = 1000
      )
      par(cex.axis = 2, cex.lab = 2, cex.main = 2)
    }
    par(mfrow = c(2, 3))

    # Samples: proportion of loci recovered
    boxplot(
      seq_per_sample_prop,
      main = paste0("Samples: prop. of ", nloci, " loci recovered"),
      xlab = paste0(
        "mean:",
        round(mean(seq_per_sample_prop, na.rm = TRUE), 2),
        " | median:",
        round(median(seq_per_sample_prop, na.rm = TRUE), 2),
        " | threshold:",
        min_loci_per_sample_prop,
        " (",
        length(outsamples_missing_loci),
        " out)"
      )
    )
    abline(h = min_loci_per_sample_prop, lty = 2, col = "red")

    # Samples: proportion of target length recovered
    boxplot(
      prop_target_length_per_sample,
      main = "Samples: prop. of target sequence length recovered",
      xlab = paste0(
        "mean:",
        round(mean(prop_target_length_per_sample, na.rm = TRUE), 2),
        " | median:",
        round(median(prop_target_length_per_sample, na.rm = TRUE), 2),
        " | threshold:",
        min_target_length_per_sample_prop,
        " (",
        length(outsamples_missing_target),
        " out)"
      )
    )
    abline(h = min_target_length_per_sample_prop, lty = 2, col = "red")

    # Scatter: prop loci vs prop target length
    plot(
      prop_target_length_per_sample,
      seq_per_sample_prop,
      main = "Prop. of loci vs\n prop. of target length",
      xlab = "Prop. of target length",
      ylab = "Prop. of loci"
    )
    abline(h = min_loci_per_sample_prop, lty = 2, col = "red")
    abline(v = min_target_length_per_sample_prop, lty = 2, col = "red")

    # Loci: proportion of samples recovered
    boxplot(
      seq_per_locus_prop,
      main = paste0("Loci: prop. of ", nsamples, " samples recovered"),
      xlab = paste0(
        "mean:",
        round(mean(seq_per_locus_prop, na.rm = TRUE), 2),
        " | median:",
        round(median(seq_per_locus_prop, na.rm = TRUE), 2),
        " | threshold:",
        min_samples_per_locus_prop,
        " (",
        length(outloci_missing_samples),
        " out)"
      )
    )
    abline(h = min_samples_per_locus_prop, lty = 2, col = "red")

    # Loci: proportion of target length recovered
    boxplot(
      prop_target_length_per_locus,
      main = "Loci: prop. of target sequence length recovered",
      xlab = paste0(
        "mean:",
        round(mean(prop_target_length_per_locus, na.rm = TRUE), 2),
        " | median:",
        round(median(prop_target_length_per_locus, na.rm = TRUE), 2),
        " | threshold:",
        min_target_length_per_locus_prop,
        " (",
        length(outloci_missing_target),
        " out)"
      )
    )
    abline(h = min_target_length_per_locus_prop, lty = 2, col = "red")

    # Scatter: prop samples vs prop target length
    plot(
      prop_target_length_per_locus,
      seq_per_locus_prop,
      main = "Prop. of samples vs\n prop. of target length",
      xlab = "Prop. of target length",
      ylab = "Prop. of samples"
    )
    abline(h = min_samples_per_locus_prop, lty = 2, col = "red")
    abline(v = min_target_length_per_locus_prop, lty = 2, col = "red")

    dev.off()
  }
}

#' Apply paralog threshold and identify outliers
#' @noRd
.apply_paralog_threshold <- function(
  loci_cl1_colmeans,
  paralog_threshold,
  paralog_file
) {
  if (
    length(paralog_threshold) != 1 ||
      paralog_threshold == "none" ||
      is.na(paralog_threshold)
  ) {
    outloci_para_all <- character(0)
    outloci_para_all_values <- numeric(0)
    threshold_value <- 1
  } else if (paralog_threshold == "outliers") {
    # Remove NA/NaN values before calculating threshold
    valid_values <- loci_cl1_colmeans[
      !is.na(loci_cl1_colmeans) &
        !is.nan(loci_cl1_colmeans)
    ]
    if (length(valid_values) < 2) {
      # Not enough data for IQR
      outloci_para_all <- character(0)
      outloci_para_all_values <- numeric(0)
      threshold_value <- NA_real_
    } else {
      threshold_value <- 1.5 *
        IQR(valid_values, na.rm = TRUE) +
        quantile(valid_values, na.rm = TRUE)[4]
      outloci_para_all_values <- loci_cl1_colmeans[
        !is.na(loci_cl1_colmeans) &
          !is.nan(loci_cl1_colmeans) &
          loci_cl1_colmeans > threshold_value
      ]
      outloci_para_all <- names(outloci_para_all_values)
    }
  } else if (paralog_threshold == "file") {
    if (is.null(paralog_file) || !file.exists(paralog_file)) {
      stop("paralog_threshold is 'file' but paralog_file does not exist")
    }
    threshold_value <- 1
    outloci_para_all <- readLines(paralog_file)
    outloci_para_all_values <- loci_cl1_colmeans[
      names(loci_cl1_colmeans) %in% outloci_para_all
    ]
  } else {
    # Numeric threshold
    threshold_value <- as.numeric(paralog_threshold)
    outloci_para_all_values <- loci_cl1_colmeans[
      loci_cl1_colmeans > threshold_value
    ]
    outloci_para_all <- names(outloci_para_all_values)
  }

  list(
    outloci_para_all = outloci_para_all,
    outloci_para_all_values = outloci_para_all_values,
    threshold_value = threshold_value
  )
}

#' Generate paralog detection plots (all samples)
#' @noRd
.generate_paralog_all_plots <- function(
  loci_cl1_colmeans,
  outloci_para_all,
  nsamples_cl1,
  nloci_cl1,
  threshold_value,
  paralog_threshold,
  output_assess
) {
  # Color outliers red
  colour_outparaall <- rep("black", length(loci_cl1_colmeans))
  sorted_names <- names(sort(loci_cl1_colmeans))
  colour_outparaall[sorted_names %in% outloci_para_all] <- "red"

  for (i in 1:2) {
    if (i == 1) {
      pdf(
        file = file.path(output_assess, "2a_Paralogs_for_all_samples.pdf"),
        width = 11,
        height = 7
      )
    } else {
      png(
        filename = file.path(output_assess, "2a_Paralogs_for_all_samples.png"),
        width = 1400,
        height = 1000
      )
      par(cex.axis = 2, cex.lab = 2, cex.main = 2)
    }

    layout(matrix(c(1, 2), 2, 2, byrow = TRUE), widths = c(5, 1))

    barplot(
      sort(loci_cl1_colmeans),
      col = colour_outparaall,
      border = NA,
      las = 2,
      main = paste0(
        "Mean % SNPs across samples (n=",
        nsamples_cl1,
        ") for each locus (n=",
        nloci_cl1,
        ")"
      )
    )
    if (length(threshold_value) > 0 && paralog_threshold != "file") {
      abline(h = threshold_value, col = "red", lty = 2)
    }

    boxplot(loci_cl1_colmeans, las = 2)
    if (length(threshold_value) > 0 && paralog_threshold != "file") {
      abline(h = threshold_value, col = "red", lty = 2)
    }

    dev.off()
  }
}

#' Generate paralog detection plots (per sample)
#' @noRd
.generate_paralog_each_plots <- function(
  tab_snps_cl2a_nozero,
  outliers_color,
  output_assess
) {
  # Order samples by mean SNP proportion
  sample_order <- order(colMeans(tab_snps_cl2a_nozero, na.rm = TRUE))

  for (i in 1:2) {
    if (i == 1) {
      pdf(
        file = file.path(output_assess, "2b_Paralogs_for_each_sample.pdf"),
        width = 10,
        height = 14
      )
    } else {
      png(
        filename = file.path(output_assess, "2b_Paralogs_for_each_sample.png"),
        width = 1000,
        height = 1400
      )
    }
    par(mfrow = c(1, 1))
    boxplot(
      as.data.frame(tab_snps_cl2a_nozero[, sample_order]),
      horizontal = TRUE,
      las = 1,
      yaxt = "n",
      main = paste0(
        "Proportions of SNPs for all loci per sample\n",
        "(only loci with any SNPs)"
      ),
      ylab = "Samples",
      xlab = "Proportion of SNPs",
      col = "grey",
      pars = list(outcol = outliers_color),
      outpch = 20
    )
    dev.off()
  }
}

#' Write paralog summary text file
#' @noRd
.write_paralog_summary <- function(
  output_assess,
  paralog_threshold,
  paralog_file,
  threshold_value,
  outloci_para_all,
  outloci_para_all_values,
  remove_sample_outliers,
  outloci_para_each,
  threshold_para_each
) {
  cl2b_file <- file.path(output_assess, "2_Summary_Paralogs.txt")
  cat(file = cl2b_file, "Removal of putative paralog loci.\n")
  cat(file = cl2b_file, "Paralogs removed for all samples:\n", append = TRUE)
  cat(
    file = cl2b_file,
    paste0(
      "Variable 'paralog_threshold' set to: ",
      paralog_threshold,
      "\n"
    ),
    append = TRUE
  )

  if (paralog_threshold == "file") {
    cat(
      file = cl2b_file,
      paste0(
        "Loci listed in this file were removed: '",
        paralog_file,
        "'\n"
      ),
      append = TRUE
    )
  } else if (paralog_threshold == "none" || is.na(paralog_threshold)) {
    cat(file = cl2b_file, "None!\n", append = TRUE)
  } else {
    cat(
      file = cl2b_file,
      paste0(
        "Resulting threshold value (mean proportion of SNPs):",
        round(threshold_value, 5),
        "\n"
      ),
      append = TRUE
    )
  }

  if (length(outloci_para_all) > 0) {
    cat(
      file = cl2b_file,
      paste0(length(outloci_para_all), " loci were removed:\n"),
      append = TRUE
    )
    cat(file = cl2b_file, "locus\tmean_prop_SNPs\n", append = TRUE)
    cat(
      file = cl2b_file,
      paste(
        paste(outloci_para_all, round(outloci_para_all_values, 4), sep = "\t"),
        collapse = "\n"
      ),
      append = TRUE
    )
    cat(file = cl2b_file, "\n\n", append = TRUE)
  }

  # Save list of paralogs
  writeLines(
    outloci_para_all,
    file.path(output_assess, "2a_List_of_paralogs_removed_for_all_samples.txt")
  )

  if (remove_sample_outliers) {
    cat(file = cl2b_file, "Paralogs removed for each sample:\n", append = TRUE)
    cat(file = cl2b_file, "Sample\tthreshold\t#removed\tnames\n", append = TRUE)
    for (i in seq_along(names(outloci_para_each))) {
      cat(file = cl2b_file, names(outloci_para_each)[i], "\t", append = TRUE)
      cat(
        file = cl2b_file,
        round(threshold_para_each[[i]], 5),
        length(outloci_para_each[[i]]),
        paste(names(outloci_para_each[[i]]), collapse = ", "),
        sep = "\t",
        append = TRUE
      )
      cat(file = cl2b_file, "\n", append = TRUE)
    }
  } else {
    cat(
      file = cl2b_file,
      "The step for removing paralogs for each samples was skipped.\n",
      append = TRUE
    )
  }
}

#' Generate summary table with heterozygosity metrics
#' @noRd
.generate_summary_table <- function(
  tab_snps_cl2b,
  tab_length_cl2b,
  targets_length_cl2b,
  outloci_para_all,
  outloci_para_each,
  nloci_cl1,
  nloci_cl2
) {
  tab_het_ad <- data.frame(
    sample = colnames(tab_snps_cl2b),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(colnames(tab_snps_cl2b))) {
    tab_het_ad$bp[i] <- sum(tab_length_cl2b[, i], na.rm = TRUE)
    tab_het_ad$bpoftarget[i] <- round(
      sum(tab_length_cl2b[, i], na.rm = TRUE) / targets_length_cl2b,
      3
    ) *
      100
    tab_het_ad$paralogs_all[i] <- length(outloci_para_all)
    tab_het_ad$paralogs_each[i] <- length(outloci_para_each[[i]])
    tab_het_ad$nloci[i] <- nloci_cl1 -
      length(outloci_para_all) -
      sum(is.na(tab_snps_cl2b[, i])) -
      length(outloci_para_each[[i]])
    tab_het_ad$allele_divergence[i] <- 100 *
      round(
        sum(tab_length_cl2b[, i] * tab_snps_cl2b[, i], na.rm = TRUE) /
          sum(tab_length_cl2b[, i], na.rm = TRUE),
        5
      )
    tab_het_ad$locus_heterozygosity[i] <- 100 *
      round(
        1 -
          sum(tab_snps_cl2b[, i] == 0, na.rm = TRUE) /
            (nloci_cl2 - sum(is.na(tab_snps_cl2b[, i]))),
        4
      )
    tab_het_ad$`loci with >0.5% SNPs`[i] <- 100 *
      round(
        1 -
          sum(tab_snps_cl2b[, i] < 0.005, na.rm = TRUE) /
            (nloci_cl2 - sum(is.na(tab_snps_cl2b[, i]))),
        4
      )
    tab_het_ad$`loci with >1% SNPs`[i] <- 100 *
      round(
        1 -
          sum(tab_snps_cl2b[, i] < 0.01, na.rm = TRUE) /
            (nloci_cl2 - sum(is.na(tab_snps_cl2b[, i]))),
        4
      )
    tab_het_ad$`loci with >2% SNPs`[i] <- 100 *
      round(
        1 -
          sum(tab_snps_cl2b[, i] < 0.02, na.rm = TRUE) /
            (nloci_cl2 - sum(is.na(tab_snps_cl2b[, i]))),
        4
      )
  }

  tab_het_ad
}

#' Generate heterozygosity and allele divergence plots
#' @noRd
.generate_heterozygosity_plots <- function(
  summary_table,
  output_assess
) {
  # Simple LH vs AD plot
  for (i in 1:2) {
    if (i == 1) {
      pdf(file.path(output_assess, "3_LH_vs_AD.pdf"), height = 10, width = 10)
    } else {
      png(
        file.path(output_assess, "3_LH_vs_AD.png"),
        height = 1000,
        width = 1000
      )
    }

    plot(
      summary_table$allele_divergence,
      summary_table$locus_heterozygosity,
      xlab = "Allele divergence [%]",
      ylab = "Locus heterozygosity [%]",
      main = "Locus heterozygosity vs allele divergence",
      las = 1
    )
    dev.off()
  }

  # Variable LH thresholds vs AD plots
  for (i in 1:2) {
    if (i == 1) {
      pdf(
        file.path(output_assess, "3_varLH_vs_AD.pdf"),
        height = 10,
        width = 10
      )
    } else {
      png(
        file.path(output_assess, "3_varLH_vs_AD.png"),
        height = 1000,
        width = 1000
      )
    }

    par(mfrow = c(2, 2))
    plot(
      summary_table$allele_divergence,
      summary_table$locus_heterozygosity,
      xlab = "Allele divergence [%]",
      ylab = "Locus heterozygosity (0% SNPs) [%]",
      main = "Locus heterozygosity (0% SNPs) vs allele divergence",
      las = 1
    )
    plot(
      summary_table$allele_divergence,
      summary_table$`loci with >0.5% SNPs`,
      xlab = "Allele divergence [%]",
      ylab = "Locus heterozygosity (>0.5% SNPs) [%]",
      main = "Locus heterozygosity (>0.5% SNPs) vs allele divergence",
      las = 1
    )
    plot(
      summary_table$allele_divergence,
      summary_table$`loci with >1% SNPs`,
      xlab = "Allele divergence [%]",
      ylab = "Locus heterozygosity (>1% SNPs) [%]",
      main = "Locus heterozygosity (>1% SNPs) vs allele divergence",
      las = 1
    )
    plot(
      summary_table$allele_divergence,
      summary_table$`loci with >2% SNPs`,
      xlab = "Allele divergence [%]",
      ylab = "Locus heterozygosity (>2% SNPs) [%]",
      main = "Locus heterozygosity (>2% SNPs) vs allele divergence",
      las = 1
    )
    par(mfrow = c(1, 1))
    dev.off()
  }
}
