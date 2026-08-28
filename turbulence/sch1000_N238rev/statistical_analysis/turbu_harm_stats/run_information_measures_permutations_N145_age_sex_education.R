#!/usr/bin/env Rscript

# Secondary N145 information-measures analysis adding education to the age- and
# sex-adjusted Freedman-Lane models. The shared driver writes distinct
# education-specific filenames so primary results are preserved.

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
    "run_information_measures_permutations_N145_age_sex_education.R"
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

driver_script <- file.path(
  resolve_script_dir(),
  "run_information_measures_permutations_N145.R"
)
if (!file.exists(driver_script)) {
  stop("Missing information-measures permutation driver: ", driver_script)
}

old_skip_option <- getOption("turbu_harm_stats.skip_information_main")
options(turbu_harm_stats.skip_information_main = TRUE)
source(driver_script, local = .GlobalEnv)
options(turbu_harm_stats.skip_information_main = old_skip_option)

main_information(c("--include-education", commandArgs(trailingOnly = TRUE)))
