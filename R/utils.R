#' Read and Validate HybPhaser Configuration
#'
#' Reads a HybPhaser `config.txt` file and validates required parameters.
#' The file is a list of `key = value` lines, where `value` is a quoted
#' string or a bare number (trailing `# comments` are ignored). Values are
#' taken literally, so Windows paths with backslashes are read correctly.
#'
#' @param config_file Path to configuration file
#' @param required_vars Character vector of required variable names
#'
#' @return Named list of configuration parameters
#' @export
#'
#' @examples
#' \dontrun{
#' config <- read_config("config.txt")
#' }
read_config <- function(config_file, required_vars = NULL) {
  if (!file.exists(config_file)) {
    stop("Configuration file not found: ", config_file)
  }

  lines <- readLines(config_file, warn = FALSE)
  config <- list()

  for (raw in lines) {
    line <- trimws(raw)
    if (line == "" || startsWith(line, "#")) {
      next
    }

    eq <- regexpr("=", line, fixed = TRUE)
    if (eq < 1L) {
      next
    }

    key <- trimws(substr(line, 1L, eq - 1L))
    if (key == "") {
      next
    }

    config[[key]] <- .parse_config_value(
      trimws(substr(line, eq + 1L, nchar(line)))
    )
  }

  if (!is.null(required_vars)) {
    missing_vars <- setdiff(required_vars, names(config))
    if (length(missing_vars) > 0) {
      stop(
        "Missing required config variables: ",
        paste(missing_vars, collapse = ", ")
      )
    }
  }

  config
}


# Parse the right-hand side of a `key = value` config line.
.parse_config_value <- function(rhs) {
  # Quoted string: return the content up to the closing quote, ignoring any
  # trailing comment. Backslashes are kept literally.
  if (startsWith(rhs, '"') || startsWith(rhs, "'")) {
    q <- substr(rhs, 1L, 1L)
    close <- regexpr(q, substr(rhs, 2L, nchar(rhs)), fixed = TRUE)
    if (close < 1L) {
      stop("Unterminated quoted value in config: ", rhs)
    }
    return(substr(rhs, 2L, close))
  }

  # Unquoted: strip a trailing comment, then coerce number / logical.
  hash <- regexpr("#", rhs, fixed = TRUE)
  if (hash > 0L) {
    rhs <- trimws(substr(rhs, 1L, hash - 1L))
  }
  if (rhs %in% c("TRUE", "FALSE")) {
    return(as.logical(rhs))
  }
  num <- suppressWarnings(as.numeric(rhs))
  if (nzchar(rhs) && !is.na(num)) {
    return(num)
  }
  rhs
}


#' Validate File and Directory Paths
#'
#' Checks that specified paths exist and have correct type (file or directory)
#'
#' @param paths Named list where names are path descriptions and values
#'   are paths
#' @param must_exist Logical; if TRUE, paths must exist
#' @param type Character vector matching paths; "file", "dir", or "any"
#'
#' @return Logical; TRUE if all validations pass, error otherwise
#' @keywords internal
validate_paths <- function(paths, must_exist = TRUE, type = "any") {
  type <- rep_len(type, length(paths))

  for (i in seq_along(paths)) {
    path <- paths[[i]]
    path_name <- names(paths)[i]
    path_type <- type[i]

    if (must_exist && !file.exists(path)) {
      stop(path_name, " not found: ", path)
    }

    if (file.exists(path)) {
      is_dir <- dir.exists(path)
      if (path_type == "dir" && !is_dir) {
        stop(path_name, " must be a directory: ", path)
      } else if (path_type == "file" && is_dir) {
        stop(path_name, " must be a file: ", path)
      }
    }
  }

  TRUE
}


#' Create Output Directory Structure
#'
#' Creates output directories with optional subset names
#'
#' @param base_path Base output directory path
#' @param subdirs Character vector of subdirectory names to create
#' @param subset_name Optional subset name to append to directories
#' @param recursive Create parent directories if needed
#'
#' @return Character vector of created directory paths (invisibly)
#' @keywords internal
create_output_dirs <- function(
  base_path,
  subdirs,
  subset_name = "",
  recursive = TRUE
) {
  if (subset_name != "") {
    subset_suffix <- paste0("_", subset_name)
  } else {
    subset_suffix <- ""
  }

  created_dirs <- character(length(subdirs))

  for (i in seq_along(subdirs)) {
    dir_path <- file.path(base_path, paste0(subdirs[i], subset_suffix))
    dir.create(dir_path, showWarnings = FALSE, recursive = recursive)
    created_dirs[i] <- dir_path
  }

  invisible(created_dirs)
}
