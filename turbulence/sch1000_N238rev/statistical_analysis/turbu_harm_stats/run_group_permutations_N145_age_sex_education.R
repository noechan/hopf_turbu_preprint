#!/usr/bin/env Rscript

# Secondary N145 turbulence group analysis adjusting for age, sex, and
# education. The validated fitting and reporting implementation is shared with
# run_group_permutations_N145.R; this entry point enables education adjustment
# and writes education-specific result filenames so primary outputs are kept.

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
    "run_group_permutations_N145_age_sex_education.R"
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

  stop(
    "Could not resolve this script's directory. Run it with Rscript or open ",
    "turbu_harm_stats.Rproj and click Source."
  )
}

core_script <- file.path(resolve_script_dir(), "run_group_permutations_N145.R")
if (!file.exists(core_script)) {
  stop("Missing shared permutation implementation: ", core_script)
}

old_skip_option <- getOption("turbu_harm_stats.skip_main")
options(turbu_harm_stats.skip_main = TRUE)
source(core_script, local = .GlobalEnv)
options(turbu_harm_stats.skip_main = old_skip_option)

main(c("--include-education", commandArgs(trailingOnly = TRUE)))
