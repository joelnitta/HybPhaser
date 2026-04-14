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


#' Check if HybPhaser Docker Image Exists
#'
#' Verifies that the HybPhaser Docker image is available locally
#'
#' @param image Name of the Docker image. Default is
#'   "joelnitta/hybphaser:latest"
#' @param pull Logical; if TRUE and image is not found, attempt to pull it
#'
#' @return Logical; TRUE if image exists,  FALSE otherwise
#' @export
#'
#' @examples
#' \dontrun{
#' check_docker_image("joelnitta/hybphaser:latest")
#' }
check_docker_image <- function(
  image = "joelnitta/hybphaser:latest",
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

  warning(
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


#' Run Docker Container with HybPhaser
#'
#' Internal function to execute commands in HybPhaser Docker container
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
  image = "joelnitta/hybphaser:latest",
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

  # Build docker run command
  docker_args <- c(
    "run",
    "--rm", # Remove container after exit
    volume_args,
    image,
    cmd
  )

  # Execute
  system2("docker", args = docker_args, ...)
}
