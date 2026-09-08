#' Check if Docker is Available
#'
#' Verifies that Docker is installed and running on the system
#'
#' @param quiet Logical; if TRUE, suppresses messages
#'
#' @return Logical; TRUE if Docker is available, FALSE otherwise
#' @export
#'
#' @examples
#' \dontrun{
#' if (check_docker()) {
#'   message("Docker is available!")
#' }
#' }
check_docker <- function(quiet = FALSE) {
  # Check if docker command exists
  docker_path <- Sys.which("docker")

  if (docker_path == "") {
    if (!quiet) {
      warning("Docker command not found. Please install Docker.")
    }
    return(FALSE)
  }

  # Check if Docker daemon is running
  result <- system2("docker", args = "info", stdout = FALSE, stderr = FALSE)

  if (result != 0) {
    if (!quiet) {
      warning("Docker is installed but not running. ", "Please start Docker.")
    }
    return(FALSE)
  }

  if (!quiet) {
    message("Docker is available")
  }

  TRUE
}


#' Check if the rhybphaser Docker Image Exists
#'
#' Verifies that the rhybphaser Docker image is available locally
#'
#' @param image Name of the Docker image. Default is
#'   "joelnitta/rhybphaser:latest"
#' @param pull Logical; if TRUE and image is not found, attempt to pull it
#'
#' @return Logical; TRUE if image exists,  FALSE otherwise
#' @export
#'
#' @examples
#' \dontrun{
#' check_docker_image("joelnitta/rhybphaser:latest")
#' }
check_docker_image <- function(
  image = "joelnitta/rhybphaser:latest",
  pull = FALSE
) {
  if (!check_docker(quiet = TRUE)) {
    stop("Docker is not available")
  }

  # Check if image exists
  result <- system2(
    "docker",
    args = c("image", "inspect", image),
    stdout = FALSE,
    stderr = FALSE
  )

  if (result == 0) {
    message("Docker image '", image, "' found")
    return(TRUE)
  }

  if (pull) {
    message("Pulling Docker image '", image, "'...")
    pull_result <- system2("docker", args = c("pull", image))

    if (pull_result == 0) {
      message("Successfully pulled image")
      return(TRUE)
    } else {
      stop("Failed to pull Docker image")
    }
  }

  message(
    "Docker image '",
    image,
    "' not found. ",
    "Set pull = TRUE to download it."
  )
  FALSE
}


#' Normalize Path for Docker Volume Mounting
#'
#' Converts a file path to an absolute path suitable for Docker volume
#' mounting
#'
#' @param path Character; file or directory path
#'
#' @return Character; normalized absolute path
#' @keywords internal
.normalize_docker_path <- function(path) {
  # Expand ~ and make absolute
  path <- path.expand(path)
  path <- normalizePath(path, mustWork = FALSE)

  # On Windows, convert to Unix-style path for Docker
  if (.Platform$OS.type == "windows") {
    # Convert C:/path to /c/path format
    path <- gsub("\\\\", "/", path)
    if (grepl("^[A-Za-z]:", path)) {
      drive <- tolower(substr(path, 1, 1))
      path <- paste0("/", drive, substr(path, 3, nchar(path)))
    }
  }

  path
}


#' Compute the relative path from a base to a target directory
#'
#' Returns the relative portion of `path` with respect to `base`, using
#' forward slashes, if `path` is a strict descendant of `base`.
#' Returns `NULL` when `path` is the same as `base` or is not nested
#' inside it.
#'
#' @param path  Character; absolute path of the target.
#' @param base  Character; absolute path of the base directory.
#'
#' @return Character relative path (no leading `/`), or `NULL`.
#' @keywords internal
.relative_path <- function(path, base) {
  path <- gsub("\\\\", "/", path)
  base <- gsub("\\\\", "/", base)
  base <- sub("/$", "", base)
  prefix <- paste0(base, "/")
  if (startsWith(path, prefix) && nchar(path) > nchar(prefix)) {
    substring(path, nchar(prefix) + 1L)
  } else {
    NULL
  }
}


#' Docker `--user` argument for the current host user
#'
#' On Linux the container process must run as the host user so that it can
#' write to (and create sub-directories in) bind-mounted host directories.
#' Docker Desktop on macOS and Windows remaps bind-mount permissions, so
#' `--user` is neither needed nor always safe there.
#'
#' @return `"<uid>:<gid>"`, or `NULL` when it should not be set.
#' @keywords internal
.docker_user <- function() {
  if (Sys.info()[["sysname"]] != "Linux") {
    return(NULL)
  }
  ids <- tryCatch(
    c(
      system2("id", "-u", stdout = TRUE, stderr = FALSE),
      system2("id", "-g", stdout = TRUE, stderr = FALSE)
    ),
    error = function(e) character(0)
  )
  ids <- ids[nzchar(ids)]
  if (length(ids) == 2 && all(grepl("^[0-9]+$", ids))) {
    paste(ids, collapse = ":")
  } else {
    NULL
  }
}


#' Run a Command in the rhybphaser Docker Container
#'
#' Internal function to execute commands in the rhybphaser Docker container
#'
#' @param cmd Character vector; command to run
#' @param volumes Named character vector; host paths to mount (names are
#'   container paths)
#' @param image Docker image name
#' @param ... Additional arguments passed to system2
#'
#' @return Exit code from Docker command
#' @keywords internal
.run_docker <- function(
  cmd,
  volumes = NULL,
  image = "joelnitta/rhybphaser:latest",
  ...
) {
  # Build volume mount arguments
  volume_args <- character(0)
  if (!is.null(volumes) && length(volumes) > 0) {
    for (i in seq_along(volumes)) {
      host_path <- .normalize_docker_path(volumes[i])
      container_path <- names(volumes)[i]
      volume_args <- c(
        volume_args,
        "-v",
        paste0(host_path, ":", container_path)
      )
    }
  }

  user <- .docker_user()
  user_args <- if (is.null(user)) character(0) else c("--user", user)

  # Give the container a writable working directory. The image's WORKDIR
  # (/data) is not writable when running as an arbitrary --user, and some
  # tools (BBSplit) write index files relative to the working directory.
  work_dir <- tempfile("rhybphaser-work-")
  dir.create(work_dir, recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  work_args <- c(
    "-v",
    paste0(.normalize_docker_path(work_dir), ":/work"),
    "--workdir",
    "/work"
  )

  # Build docker run command
  docker_args <- c(
    "run",
    "--rm", # Remove container after exit
    user_args,
    work_args,
    volume_args,
    image,
    cmd
  )

  # Execute
  system2("docker", args = docker_args, ...)
}
