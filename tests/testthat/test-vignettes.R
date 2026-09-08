# The vignette workflow code is almost entirely eval = FALSE (it needs conda,
# Docker and network). These static checks guard against that code silently
# drifting out of sync with the package API. End-to-end execution happens in
# the `vignettes.yaml` workflow.

vignette_source <- function(name) {
  installed <- system.file("doc", name, package = "rhybphaser")
  if (nzchar(installed) && file.exists(installed)) {
    return(installed)
  }
  in_tree <- testthat::test_path("..", "..", "vignettes", name)
  if (file.exists(in_tree)) {
    return(in_tree)
  }
  NA_character_
}

# Pull the contents of every ```{r ...} fenced chunk out of an Rmd.
extract_r_chunks <- function(rmd_path) {
  lines <- readLines(rmd_path, warn = FALSE)
  fence <- grepl("^```", lines)
  starts <- which(fence & grepl("^```\\{r", lines))
  code <- character()
  for (s in starts) {
    close <- which(fence & seq_along(lines) > s)[1]
    if (!is.na(close) && close > s + 1) {
      code <- c(code, lines[(s + 1):(close - 1)])
    }
  }
  code
}

exported_calls_with_bad_args <- function(code) {
  ns <- asNamespace("rhybphaser")
  exported <- getNamespaceExports("rhybphaser")
  problems <- character()

  # TRUE for the empty argument, e.g. the row index in `x[, cols]`
  is_empty_arg <- function(x) {
    tryCatch(
      {
        force(x)
        FALSE
      },
      error = function(e) TRUE
    )
  }

  walk <- function(e) {
    if (is.call(e)) {
      fn <- e[[1]]
      if (is.name(fn) && as.character(fn) %in% exported) {
        fname <- as.character(fn)
        fmls <- names(formals(get(fname, envir = ns)))
        if (!("..." %in% fmls)) {
          nm <- names(e)[-1]
          named <- nm[!is.na(nm) & nzchar(nm)]
          bad <- setdiff(named, fmls)
          if (length(bad) > 0) {
            problems <<- c(
              problems,
              sprintf("%s(): %s", fname, paste(bad, collapse = ", "))
            )
          }
        }
      }
      for (i in seq_along(e)) {
        if (!is_empty_arg(e[[i]])) {
          walk(e[[i]])
        }
      }
    }
  }

  for (e in as.list(code)) {
    walk(e)
  }
  problems
}

for (vig in c("rhybphaser.Rmd", "rhybphaser_targets.Rmd")) {
  test_that(paste(vig, "code parses and calls the package API correctly"), {
    src <- vignette_source(vig)
    skip_if(is.na(src), "vignette source not found")

    chunk_code <- extract_r_chunks(src)
    code <- parse(text = paste(chunk_code, collapse = "\n"))
    expect_type(code, "expression")

    expect_identical(exported_calls_with_bad_args(code), character())
  })
}
