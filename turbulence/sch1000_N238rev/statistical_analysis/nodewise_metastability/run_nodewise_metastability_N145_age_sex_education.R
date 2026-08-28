#!/usr/bin/env Rscript

# Secondary node-level metastability analysis adding education to the primary
# age/sex nuisance model. The implementation and output format are shared with
# run_nodewise_metastability_N145.R.

options(stringsAsFactors = FALSE)

resolve_script_dir <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) == 1L) {
    return(dirname(normalizePath(
      sub("^--file=", "", script_arg),
      mustWork = TRUE
    )))
  }
  project_candidate <- file.path(
    getwd(),
    "run_nodewise_metastability_N145_age_sex_education.R"
  )
  if (file.exists(project_candidate)) {
    return(dirname(normalizePath(project_candidate, mustWork = TRUE)))
  }
  if (
    interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    nzchar(rstudioapi::getSourceEditorContext()$path)
  ) {
    return(dirname(normalizePath(
      rstudioapi::getSourceEditorContext()$path,
      mustWork = TRUE
    )))
  }
  stop("Could not resolve this script's directory. Run it with Rscript or click Source in RStudio.")
}

core_script <- file.path(resolve_script_dir(), "run_nodewise_metastability_N145.R")
if (!file.exists(core_script)) {
  stop("Missing shared node-wise implementation: ", core_script)
}

old_skip_option <- getOption("nodewise_metastability.skip_main")
options(nodewise_metastability.skip_main = TRUE)
source(core_script, local = .GlobalEnv)
options(nodewise_metastability.skip_main = old_skip_option)

main(c("--include-education", commandArgs(trailingOnly = TRUE)))
