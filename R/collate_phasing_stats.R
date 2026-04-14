#' Collate Phasing BBSplit Stats
#'
#' Collates phasing `refstats` files into a sample-by-reference table with
#' values expressed as proportions of unambiguous reads.
#'
#' @param path_to_phasing_folder Path to phasing folder where output table is
#'   written.
#' @param csv_file_with_phasing_prep_info CSV used for phasing preparation.
#'   Expected columns are `sample` or `samples`, then
#'   `ref1, abb1, ref2, abb2, ...`.
#' @param folder_for_phasing_stats Optional folder containing phasing stats
#'   files (default: `path_to_phasing_folder/phasing_stats`).
#'
#' @return Invisibly returns a list with the collated table and output path.
#' @export
#'
#' @examples
#' \dontrun{
#' res <- collate_phasing_stats(
#'   path_to_phasing_folder = "05_phasing",
#'   csv_file_with_phasing_prep_info = "phasing_prep.csv"
#' )
#' }
collate_phasing_stats <- function(
  path_to_phasing_folder,
  csv_file_with_phasing_prep_info,
  folder_for_phasing_stats = ""
) {
  validate_paths(
    paths = list(
      "path_to_phasing_folder" = path_to_phasing_folder,
      "csv_file_with_phasing_prep_info" = csv_file_with_phasing_prep_info
    ),
    must_exist = TRUE,
    type = c("dir", "file")
  )

  if (folder_for_phasing_stats == "") {
    folder_for_phasing_stats <- file.path(
      path_to_phasing_folder,
      "phasing_stats"
    )
  }

  if (!dir.exists(folder_for_phasing_stats)) {
    stop("folder_for_phasing_stats not found: ", folder_for_phasing_stats)
  }

  prep_phasing <- utils::read.csv(
    csv_file_with_phasing_prep_info,
    header = TRUE,
    stringsAsFactors = FALSE
  )
  prep_phasing[is.na(prep_phasing)] <- ""

  sample_col <- detect_phasing_stats_sample_column(prep_phasing)
  used_reference_abbr <- extract_phasing_reference_abbreviations(
    prep_phasing = prep_phasing,
    sample_col = sample_col
  )

  if (length(used_reference_abbr) == 0) {
    stop("No reference abbreviations found in phasing prep CSV")
  }

  stats_files <- list.files(
    folder_for_phasing_stats,
    pattern = "_phasing-stats\\.txt$",
    full.names = FALSE
  )
  if (length(stats_files) == 0) {
    stop("No phasing stats files found in: ", folder_for_phasing_stats)
  }

  sample_names <- gsub("_phasing-stats\\.txt$", "", stats_files)

  tab_phasing_stats <- matrix(
    NA_real_,
    nrow = length(stats_files),
    ncol = length(used_reference_abbr),
    dimnames = list(sample_names, used_reference_abbr)
  )

  for (i in seq_along(stats_files)) {
    stats_file <- file.path(folder_for_phasing_stats, stats_files[[i]])
    p_unamb <- read_phasing_stats(stats_file)

    for (j in seq_len(nrow(p_unamb))) {
      idx <- match(p_unamb$name[[j]], colnames(tab_phasing_stats))
      if (!is.na(idx)) {
        tab_phasing_stats[i, idx] <- round(p_unamb$percent[[j]] / 100, 5)
      }
    }
  }

  table_phasing_stats <- as.data.frame(tab_phasing_stats)
  table_phasing_stats$sample <- rownames(table_phasing_stats)
  table_phasing_stats <- table_phasing_stats[, c("sample", used_reference_abbr)]

  out_file <- file.path(path_to_phasing_folder, "Table_phasing_stats.csv")
  utils::write.csv(
    table_phasing_stats,
    file = out_file,
    row.names = FALSE,
    na = ""
  )

  invisible(list(
    table_phasing_stats = table_phasing_stats,
    output_file = out_file
  ))
}


#' Collate Phasing BBSplit Stats Using a HybPhaser Config File
#'
#' Reads required settings from `config.txt` and runs
#' `collate_phasing_stats()`.
#'
#' @param config_file Path to HybPhaser configuration file.
#'
#' @return Invisibly returns the same object as `collate_phasing_stats()`.
#' @export
collate_phasing_stats_from_config <- function(config_file = "./config.txt") {
  required <- c(
    "path_to_phasing_folder",
    "csv_file_with_phasing_prep_info",
    "folder_for_phasing_stats"
  )

  cfg <- read_config(config_file, required_vars = required)

  collate_phasing_stats(
    path_to_phasing_folder = cfg$path_to_phasing_folder,
    csv_file_with_phasing_prep_info = cfg$csv_file_with_phasing_prep_info,
    folder_for_phasing_stats = cfg$folder_for_phasing_stats
  )
}


# Detect sample column name in phasing prep table.
detect_phasing_stats_sample_column <- function(prep_phasing) {
  cols <- colnames(prep_phasing)
  if ("samples" %in% cols) {
    return("samples")
  }
  if ("sample" %in% cols) {
    return("sample")
  }
  cols[[1]]
}


# Extract non-empty unique reference abbreviations from prep table.
extract_phasing_reference_abbreviations <- function(prep_phasing, sample_col) {
  ref_cols <- setdiff(colnames(prep_phasing), sample_col)
  if (length(ref_cols) < 2) {
    return(character(0))
  }

  abb_cols <- ref_cols[seq(2, length(ref_cols), by = 2)]
  used_reference_abbr <- unique(as.vector(as.matrix(prep_phasing[, abb_cols])))
  used_reference_abbr <- used_reference_abbr[
    !(is.na(used_reference_abbr) | used_reference_abbr == "")
  ]

  as.character(used_reference_abbr)
}


# Read phasing stats table and keep first two columns as name/percent.
read_phasing_stats <- function(stats_file) {
  x <- utils::read.delim(stats_file, stringsAsFactors = FALSE)
  if (ncol(x) < 2) {
    stop("Stats file has fewer than 2 columns: ", stats_file)
  }

  out <- x[, 1:2]
  colnames(out) <- c("name", "percent")
  out
}
