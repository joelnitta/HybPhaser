#' Merge Phased and Non-phased Sequence Lists
#'
#' Combines phased and non-phased sequence lists, optionally applies locus and
#' sample subsets, and can replace non-phased samples with phased derivatives.
#'
#' @param path_to_sequence_lists_normal Folder with non-phased sequence lists.
#' @param path_to_sequence_lists_phased Folder with phased sequence lists.
#' @param path_of_sequence_lists_output Output folder for merged sequence lists.
#' @param path_to_namelist_normal Path to namelist for non-phased samples.
#' @param path_to_namelist_phased Path to namelist for phased samples.
#' @param file_with_samples_included Optional file listing samples to include.
#' @param file_with_samples_excluded Optional file listing samples to exclude.
#' @param file_with_loci_excluded Optional file listing loci to exclude.
#' @param file_with_loci_included Optional file listing loci to include.
#' @param exchange_phased_with_not_phased_samples Whether to remove non-phased
#'   accessions when phased derivatives exist (`"yes"` or `"no"`).
#' @param include_phased_seqlists_when_non_phased_locus_absent Whether to copy
#'   loci present only in phased lists (`"yes"` or `"no"`).
#'
#' @return Invisibly returns a list with output files and retained samples/loci.
#' @export
#'
#' @examples
#' \dontrun{
#' merge_sequence_lists(
#'   path_to_sequence_lists_normal = "03_sequence_lists/loci_consensus",
#'   path_to_sequence_lists_phased = "03_sequence_lists_phased/loci_consensus",
#'   path_of_sequence_lists_output = "03_sequence_lists_merged/loci_consensus",
#'   path_to_namelist_normal = "namelist_normal.txt",
#'   path_to_namelist_phased = "namelist_phased.txt"
#' )
#' }
merge_sequence_lists <- function(
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
) {
  exchange_phased_with_not_phased_samples <- match.arg(
    exchange_phased_with_not_phased_samples
  )
  include_phased_seqlists_when_non_phased_locus_absent <- match.arg(
    include_phased_seqlists_when_non_phased_locus_absent
  )

  tryCatch(
    {
      validate_paths(
        paths = list(
          "path_to_sequence_lists_normal" = path_to_sequence_lists_normal,
          "path_to_sequence_lists_phased" = path_to_sequence_lists_phased,
          "path_to_namelist_normal" = path_to_namelist_normal,
          "path_to_namelist_phased" = path_to_namelist_phased
        ),
        must_exist = TRUE,
        type = c("dir", "dir", "file", "file")
      )
    },
    error = function(e) {
      if (grepl("path_to_sequence_lists_phased not found", e$message)) {
        stop(
          "path_to_sequence_lists_phased not found: ",
          path_to_sequence_lists_phased,
          "\nExpected phased sequence lists from a separate phased run ",
          "(e.g. hybphaser_output_phased/03_sequence_lists/loci_consensus)."
        )
      }
      stop(e)
    }
  )

  validate_optional_list_file(
    file_path = file_with_samples_included,
    arg_name = "file_with_samples_included"
  )
  validate_optional_list_file(
    file_path = file_with_samples_excluded,
    arg_name = "file_with_samples_excluded"
  )
  validate_optional_list_file(
    file_path = file_with_loci_included,
    arg_name = "file_with_loci_included"
  )
  validate_optional_list_file(
    file_path = file_with_loci_excluded,
    arg_name = "file_with_loci_excluded"
  )

  dir.create(
    path_of_sequence_lists_output,
    recursive = TRUE,
    showWarnings = FALSE
  )

  seqlists_normal <- list_sequence_list_files(path_to_sequence_lists_normal)
  seqlists_phased <- list_sequence_list_files(path_to_sequence_lists_phased)

  loci2include <- determine_loci_to_include(
    seqlists_normal = seqlists_normal,
    seqlists_phased = seqlists_phased,
    file_with_loci_included = file_with_loci_included,
    file_with_loci_excluded = file_with_loci_excluded
  )

  copied_normal <- copy_normal_sequence_lists(
    seqlists_normal = seqlists_normal,
    loci2include = loci2include,
    output_dir = path_of_sequence_lists_output
  )

  appended_or_copied_phased <- append_or_copy_phased_sequence_lists(
    seqlists_phased = seqlists_phased,
    seqlists_normal = seqlists_normal,
    loci2include = loci2include,
    output_dir = path_of_sequence_lists_output,
    include_phased_only = include_phased_seqlists_when_non_phased_locus_absent ==
      "yes"
  )

  samples_normal <- read_trimmed_lines(path_to_namelist_normal)
  samples_phased <- read_trimmed_lines(path_to_namelist_phased)

  samples2include <- determine_samples_to_include(
    file_with_samples_included = file_with_samples_included,
    samples_normal = samples_normal,
    samples_phased = samples_phased
  )

  samples2exclude <- character(0)
  if (file_with_samples_excluded != "") {
    samples2exclude <- read_trimmed_lines(file_with_samples_excluded)
  }

  notphased2remove <- character(0)
  if (exchange_phased_with_not_phased_samples == "yes") {
    phased_samples_base <- sub("_to_.*", "", samples_phased)
    notphased2remove <- intersect(samples2include, phased_samples_base)
  }

  excluded_all <- unique(c(samples2exclude, notphased2remove))
  if (length(excluded_all) > 0) {
    samples2include <- setdiff(samples2include, excluded_all)
  }

  output_files <- list_sequence_list_files(path_of_sequence_lists_output)
  for (seqlist_out in output_files) {
    keep_only_selected_samples(
      seqlist_file = seqlist_out,
      samples2include = samples2include
    )
  }

  invisible(list(
    output_dir = path_of_sequence_lists_output,
    output_files = output_files,
    loci_included = loci2include,
    samples_included = samples2include,
    samples_excluded = excluded_all,
    copied_normal_files = copied_normal,
    phased_files_merged = appended_or_copied_phased
  ))
}


#' Merge Sequence Lists Using a HybPhaser Config File
#'
#' Reads required settings from `config.txt` and runs `merge_sequence_lists()`.
#'
#' @param config_file Path to HybPhaser configuration file.
#'
#' @return Invisibly returns the same object as `merge_sequence_lists()`.
#' @export
merge_sequence_lists_from_config <- function(config_file = "./config.txt") {
  required <- c(
    "path_to_sequence_lists_normal",
    "path_to_sequence_lists_phased",
    "path_of_sequence_lists_output",
    "path_to_namelist_normal",
    "path_to_namelist_phased",
    "file_with_samples_included",
    "file_with_samples_excluded",
    "file_with_loci_excluded",
    "file_with_loci_included",
    "exchange_phased_with_not_phased_samples",
    "include_phased_seqlists_when_non_phased_locus_absent"
  )

  cfg <- read_config(config_file, required_vars = required)

  merge_sequence_lists(
    path_to_sequence_lists_normal = cfg$path_to_sequence_lists_normal,
    path_to_sequence_lists_phased = cfg$path_to_sequence_lists_phased,
    path_of_sequence_lists_output = cfg$path_of_sequence_lists_output,
    path_to_namelist_normal = cfg$path_to_namelist_normal,
    path_to_namelist_phased = cfg$path_to_namelist_phased,
    file_with_samples_included = cfg$file_with_samples_included,
    file_with_samples_excluded = cfg$file_with_samples_excluded,
    file_with_loci_excluded = cfg$file_with_loci_excluded,
    file_with_loci_included = cfg$file_with_loci_included,
    exchange_phased_with_not_phased_samples = tolower(
      cfg$exchange_phased_with_not_phased_samples
    ),
    include_phased_seqlists_when_non_phased_locus_absent = tolower(
      cfg$include_phased_seqlists_when_non_phased_locus_absent
    )
  )
}


# Validate optional list file argument.
validate_optional_list_file <- function(file_path, arg_name) {
  if (file_path == "") {
    return(invisible(TRUE))
  }
  if (!file.exists(file_path)) {
    stop(arg_name, " not found: ", file_path)
  }
  invisible(TRUE)
}


# List sequence list files in a directory.
list_sequence_list_files <- function(path) {
  list.files(
    path,
    pattern = "_(consensus|contig)\\.fasta$",
    full.names = TRUE
  )
}


# Derive locus name from sequence list filename.
get_locus_name <- function(path) {
  x <- basename(path)
  x <- sub("_(consensus|contig)\\.fasta$", "", x)
  x
}


# Read and trim non-empty lines.
read_trimmed_lines <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x[nzchar(x)]
}


# Determine loci to include from files or defaults.
determine_loci_to_include <- function(
  seqlists_normal,
  seqlists_phased,
  file_with_loci_included,
  file_with_loci_excluded
) {
  if (file_with_loci_included == "") {
    loci_normal <- vapply(seqlists_normal, get_locus_name, character(1))
    loci_phased <- vapply(seqlists_phased, get_locus_name, character(1))
    loci2include <- unique(c(loci_normal, loci_phased))
  } else {
    loci2include <- read_trimmed_lines(file_with_loci_included)
  }

  if (file_with_loci_excluded != "") {
    loci2exclude <- read_trimmed_lines(file_with_loci_excluded)
    loci2include <- setdiff(loci2include, loci2exclude)
  }

  loci2include
}


# Copy normal sequence lists for included loci.
copy_normal_sequence_lists <- function(
  seqlists_normal,
  loci2include,
  output_dir
) {
  copied <- character(0)

  for (seqln in seqlists_normal) {
    if (get_locus_name(seqln) %in% loci2include) {
      if (file.copy(seqln, to = output_dir, overwrite = TRUE)) {
        copied <- c(copied, file.path(output_dir, basename(seqln)))
      }
    }
  }

  copied
}


# Append phased records to existing loci or copy phased-only loci.
append_or_copy_phased_sequence_lists <- function(
  seqlists_phased,
  seqlists_normal,
  loci2include,
  output_dir,
  include_phased_only
) {
  normal_loci <- vapply(seqlists_normal, get_locus_name, character(1))
  out_files <- list_sequence_list_files(output_dir)
  out_loci <- vapply(out_files, get_locus_name, character(1))

  merged <- character(0)
  for (seqlp in seqlists_phased) {
    loci_p <- get_locus_name(seqlp)
    if (!(loci_p %in% loci2include)) {
      next
    }

    if (loci_p %in% normal_loci) {
      target <- out_files[match(loci_p, out_loci)]
      if (!is.na(target)) {
        file.append(target, seqlp)
        merged <- c(merged, target)
      }
    } else if (isTRUE(include_phased_only)) {
      if (file.copy(seqlp, to = output_dir, overwrite = TRUE)) {
        merged <- c(merged, file.path(output_dir, basename(seqlp)))
      }
    }
  }

  unique(merged)
}


# Determine included samples from include file or namelists.
determine_samples_to_include <- function(
  file_with_samples_included,
  samples_normal,
  samples_phased
) {
  if (file_with_samples_included == "") {
    samples2include <- c(samples_normal, samples_phased)
  } else {
    samples2include <- read_trimmed_lines(file_with_samples_included)
  }

  unique(samples2include)
}


# Keep only FASTA records with headers matching selected samples.
keep_only_selected_samples <- function(seqlist_file, samples2include) {
  raw <- readLines(seqlist_file, warn = FALSE)
  if (length(raw) == 0) {
    return(invisible(NULL))
  }

  header_idx <- which(startsWith(raw, ">"))
  if (length(header_idx) == 0) {
    writeLines(character(0), con = seqlist_file)
    return(invisible(NULL))
  }

  header_idx_end <- c(header_idx[-1] - 1, length(raw))
  keep_ranges <- vector("list", length(header_idx))
  keep_i <- 0L

  for (i in seq_along(header_idx)) {
    h <- raw[header_idx[[i]]]
    sample_id <- trimws(sub("^>", "", h))
    if (sample_id %in% samples2include) {
      keep_i <- keep_i + 1L
      keep_ranges[[keep_i]] <- seq.int(header_idx[[i]], header_idx_end[[i]])
    }
  }

  if (keep_i == 0) {
    writeLines(character(0), con = seqlist_file)
  } else {
    keep_ranges <- keep_ranges[seq_len(keep_i)]
    keep_lines <- raw[unlist(keep_ranges, use.names = FALSE)]
    writeLines(keep_lines, con = seqlist_file)
  }

  invisible(NULL)
}
