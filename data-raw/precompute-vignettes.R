# Pre-render vignettes whose code is too expensive to run during
# `R CMD check` / pkgdown (they need conda + HybPiper + Docker + network).
#
# `vignettes/rhybphaser_targets.Rmd.orig` is the source (live code);
# this script knits it to `vignettes/rhybphaser_targets.Rmd`, which is the
# committed, pre-rendered vignette with results and the figure baked in.
#
# Usage (from the package root, with rhybphaser installed):
#   R CMD INSTALL .
#   Rscript data-raw/precompute-vignettes.R
#
# CI runs this on a schedule via .github/workflows/precompute-vignettes.yaml.

if (!requireNamespace("knitr", quietly = TRUE)) {
  stop("knitr is required to pre-render the vignettes")
}

Sys.setenv(RHYBPHASER_RUN_TARGETS_VIGNETTE = "true")

old_wd <- setwd("vignettes")
on.exit(setwd(old_wd), add = TRUE)

knitr::opts_chunk$set(comment = "#>", collapse = TRUE)
knitr::knit(
  input = "rhybphaser_targets.Rmd.orig",
  output = "rhybphaser_targets.Rmd"
)

message("Wrote vignettes/rhybphaser_targets.Rmd")
