#' Get Path to HybPhaser Bash Scripts
#'
#' Returns the full path to bash scripts included with the package
#'
#' @param script_name Name of the script. If NULL, returns the scripts
#'   directory. One of:
#'   - "1_generate_consensus_sequences.sh"
#'   - "1_generate_consensus_sequences_single-core.sh"
#'   - "2_extract_mapped_reads.sh"
#'
#' @return Character path to the script or scripts directory
#' @export
#'
#' @examples
#' # Get path to scripts directory
#' script_dir <- hybphaser_scripts()
#'
#' # Get path to specific script
#' consensus_script <- hybphaser_scripts(
#'   "1_generate_consensus_sequences.sh"
#' )
hybphaser_scripts <- function(script_name = NULL) {
  scripts_dir <- system.file("scripts", package = "rhybphaser")

  if (scripts_dir == "") {
    stop(
      "Cannot find HybPhaser scripts directory. ",
      "Is the package installed correctly?"
    )
  }

  if (is.null(script_name)) {
    return(scripts_dir)
  }

  script_path <- file.path(scripts_dir, script_name)

  if (!file.exists(script_path)) {
    stop(
      "Script not found: ",
      script_name,
      "\n",
      "Available scripts: ",
      paste(list.files(scripts_dir), collapse = ", ")
    )
  }

  script_path
}


#' Get Path to HybPhaser Example Data
#'
#' Returns the full path to example data files included with the package
#'
#' @param file_name Name of the example file. If NULL, returns the extdata
#'   directory. Available files:
#'   - "clade_references.csv"
#'   - "phasing_prep.csv"
#'
#' @return Character path to the example file or extdata directory
#' @export
#'
#' @examples
#' # Get path to extdata directory
#' data_dir <- hybphaser_example()
#'
#' # Get path to specific example file
#' clade_ref <- hybphaser_example("clade_references.csv")
hybphaser_example <- function(file_name = NULL) {
  extdata_dir <- system.file("extdata", package = "rhybphaser")

  if (extdata_dir == "") {
    stop(
      "Cannot find HybPhaser extdata directory. ",
      "Is the package installed correctly?"
    )
  }

  if (is.null(file_name)) {
    return(extdata_dir)
  }

  file_path <- file.path(extdata_dir, file_name)

  if (!file.exists(file_path)) {
    stop(
      "Example file not found: ",
      file_name,
      "\n",
      "Available files: ",
      paste(list.files(extdata_dir), collapse = ", ")
    )
  }

  file_path
}


#' Run HybPhaser Bash Script
#'
#' Execute one of the HybPhaser bash scripts with specified arguments
#'
#' @param script_name Name of the script to run
#' @param args Character vector of arguments to pass to the script
#' @param dry_run Logical; if TRUE, prints the command instead of running it
#'
#' @return If dry_run = FALSE, returns the exit status invisibly. If
#'   dry_run = TRUE, returns the command string.
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate consensus sequences
#' run_hybphaser_script(
#'   "1_generate_consensus_sequences.sh",
#'   args = c("-n", "namelist.txt", "-p", "hybpiper_output")
#' )
#'
#' # Preview command without running
#' run_hybphaser_script(
#'   "1_generate_consensus_sequences.sh",
#'   args = c("-n", "namelist.txt"),
#'   dry_run = TRUE
#' )
#' }
run_hybphaser_script <- function(
  script_name,
  args = character(0),
  dry_run = FALSE
) {
  script_path <- hybphaser_scripts(script_name)

  # Make sure script is executable
  if (!dry_run) {
    Sys.chmod(script_path, mode = "0755")
  }

  cmd <- paste(c(shQuote(script_path), args), collapse = " ")

  if (dry_run) {
    message("Command to run:\n", cmd)
    return(cmd)
  }

  message("Running: ", basename(script_path))
  result <- system(cmd)

  if (result != 0) {
    warning("Script exited with non-zero status: ", result)
  }

  invisible(result)
}
