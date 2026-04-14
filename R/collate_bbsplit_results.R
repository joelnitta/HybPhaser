#' Collate BBSplit Results
#'
#' Collates BBSplit `refstats` files into clade-association tables,
#' including a normalized table and optional merges with the HybPhaser
#' summary table.
#'
#' @param path_to_clade_association_folder Path to clade association folder
#'   containing `bbsplit_stats/`.
#' @param csv_file_with_clade_reference_names CSV with clade reference sample
#'   names and abbreviations (first two columns).
#' @param path_to_output_folder Optional HybPhaser output folder used to find
#'   `00_R_objects/.../Summary_table.Rds` for merged output tables.
#' @param subset_name Optional subset folder name under `00_R_objects`.
#'
#' @return Invisibly returns a list with generated tables and output paths.
#' @export
#'
#' @examples
#' \dontrun{
#' res <- collate_bbsplit_results(
#'   path_to_clade_association_folder = "04_clade_association",
#'   csv_file_with_clade_reference_names = "clade_references.csv"
#' )
#' }
collate_bbsplit_results <- function(
  path_to_clade_association_folder,
  csv_file_with_clade_reference_names,
  path_to_output_folder = NULL,
  subset_name = ""
) {
  validate_paths(
    paths = list(
      "path_to_clade_association_folder" = path_to_clade_association_folder,
      "csv_file_with_clade_reference_names" = csv_file_with_clade_reference_names
    ),
    must_exist = TRUE,
    type = c("dir", "file")
  )

  folder_bbsplit_stats <- file.path(
    path_to_clade_association_folder,
    "bbsplit_stats"
  )
  if (!dir.exists(folder_bbsplit_stats)) {
    stop("bbsplit_stats folder not found: ", folder_bbsplit_stats)
  }

  stats_files <- list.files(
    folder_bbsplit_stats,
    pattern = "_bbsplit-stats\\.txt$",
    full.names = FALSE
  )
  if (length(stats_files) == 0) {
    stop("No BBSplit stats files found in: ", folder_bbsplit_stats)
  }

  samples <- gsub("(_mapped.*|.fasta.*|_bbsplit-stats.*)", "", stats_files)

  ref_samples <- utils::read.csv(
    csv_file_with_clade_reference_names,
    header = TRUE,
    stringsAsFactors = FALSE
  )

  if (ncol(ref_samples) < 2) {
    stop("csv_file_with_clade_reference_names must have at least 2 columns")
  }

  colnames(ref_samples)[1:2] <- c("sample", "abbreviation")

  if ("abbreviation" %in% colnames(ref_samples)) {
    ref_names <- ref_samples$abbreviation
  } else {
    ref_names <- ref_samples$sample
  }

  tab_clade_assoc <- matrix(
    nrow = length(ref_names),
    ncol = length(samples),
    dimnames = list(ref_names, samples)
  )

  for (i in seq_along(stats_files)) {
    stats_file <- stats_files[[i]]
    p_unamb <- read_bbsplit_stats(
      file.path(folder_bbsplit_stats, stats_file)
    )

    tab_clade_assoc[, i] <- p_unamb$unambiguous_reads[match(
      rownames(tab_clade_assoc),
      p_unamb$name
    )]
  }

  tab_clade_assoc_norm <- normalize_clade_association(
    tab_clade_assoc = tab_clade_assoc,
    ref_samples = ref_samples
  )

  ttab_clade_assoc <- as.data.frame(t(tab_clade_assoc))
  ttab_clade_assoc$sample <- rownames(ttab_clade_assoc)

  ttab_clade_assoc_norm <- as.data.frame(t(tab_clade_assoc_norm))
  ttab_clade_assoc_norm$sample <- rownames(ttab_clade_assoc_norm)

  summary_table <- maybe_read_summary_table(
    path_to_output_folder = path_to_output_folder,
    subset_name = subset_name
  )

  tab_hetad_cladeasso <- NULL
  tab_hetad_cladeasso_norm <- NULL
  if (!is.null(summary_table)) {
    tab_hetad_cladeasso <- merge(
      summary_table,
      ttab_clade_assoc,
      by = "sample",
      incomparables = FALSE
    )

    tab_hetad_cladeasso_norm <- merge(
      summary_table,
      ttab_clade_assoc_norm,
      by = "sample"
    )
  }

  out_raw <- file.path(
    path_to_clade_association_folder,
    "Table_clade_association.csv"
  )
  out_norm <- file.path(
    path_to_clade_association_folder,
    "Table_clade_association_normalised.csv"
  )
  out_raw_summary <- file.path(
    path_to_clade_association_folder,
    "Table_clade_association_and_summary_table.csv"
  )
  out_norm_summary <- file.path(
    path_to_clade_association_folder,
    "Table_clade_association_normalized_and_summary_table.csv"
  )

  utils::write.csv(t(tab_clade_assoc), file = out_raw)
  utils::write.csv(t(tab_clade_assoc_norm), file = out_norm)

  if (!is.null(tab_hetad_cladeasso)) {
    utils::write.csv(tab_hetad_cladeasso, file = out_raw_summary)
    utils::write.csv(tab_hetad_cladeasso_norm, file = out_norm_summary)
  }

  invisible(list(
    table_clade_association = ttab_clade_assoc,
    table_clade_association_normalized = ttab_clade_assoc_norm,
    table_with_summary = tab_hetad_cladeasso,
    table_normalized_with_summary = tab_hetad_cladeasso_norm,
    output_files = c(
      table = out_raw,
      table_normalized = out_norm,
      table_with_summary = out_raw_summary,
      table_normalized_with_summary = out_norm_summary
    )
  ))
}


#' Collate BBSplit Results Using a HybPhaser Config File
#'
#' Reads required settings from `config.txt` and runs
#' `collate_bbsplit_results()`.
#'
#' @param config_file Path to HybPhaser configuration file.
#'
#' @return Invisibly returns the same object as `collate_bbsplit_results()`.
#' @export
collate_bbsplit_results_from_config <- function(config_file = "./config.txt") {
  required <- c(
    "path_to_clade_association_folder",
    "csv_file_with_clade_reference_names",
    "path_to_output_folder",
    "name_for_dataset_optimization_subset"
  )

  cfg <- read_config(config_file, required_vars = required)

  collate_bbsplit_results(
    path_to_clade_association_folder = cfg$path_to_clade_association_folder,
    csv_file_with_clade_reference_names = cfg$csv_file_with_clade_reference_names,
    path_to_output_folder = cfg$path_to_output_folder,
    subset_name = cfg$name_for_dataset_optimization_subset
  )
}


# Read BBSplit refstats with robust first-two-column parsing.
read_bbsplit_stats <- function(stats_file) {
  x <- utils::read.delim(stats_file, stringsAsFactors = FALSE)
  if (ncol(x) < 2) {
    stop("Stats file has fewer than 2 columns: ", stats_file)
  }

  out <- x[, 1:2]
  colnames(out) <- c("name", "unambiguous_reads")
  out
}


# Normalize clade-association values by reference self-match columns.
normalize_clade_association <- function(tab_clade_assoc, ref_samples) {
  tab_clade_assoc_norm <- tab_clade_assoc

  for (i in seq_len(nrow(tab_clade_assoc))) {
    ref_sample <- ref_samples$sample[[i]]
    sample_col <- grep(
      paste0("\\b", ref_sample, "\\b"),
      colnames(tab_clade_assoc)
    )

    if (length(sample_col) > 0) {
      denom <- tab_clade_assoc[i, sample_col[[1]]]
      if (!is.na(denom) && denom != 0) {
        tab_clade_assoc_norm[i, ] <- tab_clade_assoc[i, ] / denom
      }
    }
  }

  tab_clade_assoc_norm
}


# Read summary table if available; otherwise return NULL.
maybe_read_summary_table <- function(path_to_output_folder, subset_name) {
  if (is.null(path_to_output_folder) || path_to_output_folder == "") {
    return(NULL)
  }

  summary_path <- file.path(
    path_to_output_folder,
    "00_R_objects",
    subset_name,
    "Summary_table.Rds"
  )

  if (!file.exists(summary_path)) {
    return(NULL)
  }

  readRDS(summary_path)
}
