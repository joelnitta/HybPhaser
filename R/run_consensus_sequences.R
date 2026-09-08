#' Generate Consensus Sequences and Return Output Directory
#'
#' Wrapper around [run_generate_consensus_sequences()] that returns
#' `output_dir` on success. This is useful in `{targets}` pipelines where
#' returning a path can make dependencies easier to track.
#'
#' @inheritParams run_generate_consensus_sequences
#'
#' @return Character scalar. The `output_dir` path.
#' @export
#'
#' @examples
#' \dontrun{
#' out_dir <- run_consensus_sequences(
#'   hybpiper_dir = "path/to/hybpiper_output",
#'   output_dir = "path/to/hybphaser_output",
#'   namelist = "samples.txt",
#'   threads = 4
#' )
#' }
run_consensus_sequences <- function(
  hybpiper_dir,
  output_dir,
  namelist = NULL,
  sample = NULL,
  intronerate = FALSE,
  cleanup = FALSE,
  threads = 1,
  min_depth = 10,
  min_allele_freq = 0.15,
  min_allele_count = 4,
  engine = c("docker", "local"),
  docker_image = "joelnitta/rhybphaser:latest",
  pull_image = FALSE
) {
  engine <- match.arg(engine)

  result <- run_generate_consensus_sequences(
    hybpiper_dir = hybpiper_dir,
    output_dir = output_dir,
    namelist = namelist,
    sample = sample,
    intronerate = intronerate,
    cleanup = cleanup,
    threads = threads,
    min_depth = min_depth,
    min_allele_freq = min_allele_freq,
    min_allele_count = min_allele_count,
    engine = engine,
    docker_image = docker_image,
    pull_image = pull_image
  )

  if (!identical(as.integer(result), 0L)) {
    stop("Consensus sequence generation failed with exit code: ", result)
  }

  # The bash script uses `set +e`, so it can exit 0 even when mkdir/copy
  # operations fail. Verify that at least one consensus FASTA was actually
  # produced.
  consensus_root <- file.path(output_dir, "01_data")
  consensus_files <- list.files(
    consensus_root,
    pattern = "\\.fasta$",
    recursive = TRUE
  )
  if (length(consensus_files) == 0L) {
    msg <- paste0(
      "Consensus sequence generation reported success but no FASTA files ",
      "were produced in '",
      consensus_root,
      "'."
    )
    if (engine == "docker") {
      msg <- paste0(
        msg,
        " This often means the Docker volume mount for the output directory ",
        "failed (a known Docker Desktop issue with freshly-created nested ",
        "directories on macOS). Re-running the pipeline usually resolves it."
      )
    }
    stop(msg)
  }

  normalizePath(output_dir, mustWork = TRUE)
}
