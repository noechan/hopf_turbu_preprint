#!/usr/bin/env Rscript

# Recalculate pairwise BH-FDR across the five planned AD-stage contrasts
# separately at each turbulence lambda.
#
# This script does not refit the Freedman-Lane models. It consumes the active
# pairwise result table, preserving its permutation p-values, test statistics,
# effect sizes, and sample sizes, and changes only the multiplicity family.

options(stringsAsFactors = FALSE)

EXPECTED_OUTCOMES <- paste0(
  "Turbu_lam_0_",
  c("27", "24", "21", "18", "15", "12", "09", "06", "03", "01")
)

EXPECTED_COMPARISONS <- c(
  "HC_ABneg_vs_HC_ABpos",
  "HC_ABpos_vs_MCI_ABpos",
  "MCI_ABpos_vs_AD_ABpos",
  "HC_ABneg_vs_MCI_ABpos",
  "HC_ABneg_vs_AD_ABpos"
)

NEW_FDR_COLUMN <- "FDR_BH_across_5_stage_contrasts_within_lambda"
NEW_SIG_COLUMN <- "Significant_FDR_within_lambda_0_05"

usage <- function() {
  cat(paste0(
    "Usage: Rscript recalculate_pairwise_FDR_by_lambda_N145.R [options]\n\n",
    "Options:\n",
    "  --input-file PATH    Pairwise Freedman-Lane result CSV.\n",
    "  --output-file PATH   Destination CSV with lambda-wise FDR.\n",
    "  --alpha NUMBER       Significance threshold (default: 0.05).\n",
    "  --validate-only      Validate without writing an output file.\n",
    "  --help               Show this message.\n"
  ))
}

resolve_script_path <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) == 1L) {
    return(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
  }

  project_candidate <- file.path(
    getwd(),
    "recalculate_pairwise_FDR_by_lambda_N145.R"
  )
  if (file.exists(project_candidate)) {
    return(normalizePath(project_candidate, mustWork = TRUE))
  }

  if (
    interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    nzchar(rstudioapi::getSourceEditorContext()$path)
  ) {
    return(normalizePath(
      rstudioapi::getSourceEditorContext()$path,
      mustWork = TRUE
    ))
  }

  stop(
    "Could not resolve this script's path. Run it with Rscript or open ",
    "turbu_harm_stats.Rproj and click Source."
  )
}

parse_args <- function(args) {
  project_dir <- dirname(resolve_script_path())
  results_dir <- file.path(project_dir, "results", "N145_group_permutations")
  config <- list(
    input_file = file.path(
      results_dir,
      "N145_turbulence_pairwise_FreedmanLane_age_sex.csv"
    ),
    output_file = file.path(
      results_dir,
      "N145_turbulence_pairwise_FreedmanLane_age_sex_FDR_by_lambda.csv"
    ),
    alpha = 0.05,
    validate_only = FALSE
  )

  value_options <- c(
    "--input-file" = "input_file",
    "--output-file" = "output_file",
    "--alpha" = "alpha"
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      usage()
      quit(status = 0L)
    }
    if (arg == "--validate-only") {
      config$validate_only <- TRUE
      i <- i + 1L
      next
    }
    if (!arg %in% names(value_options)) {
      stop("Unknown argument: ", arg, ". Use --help for available options.")
    }
    if (i == length(args)) {
      stop("Missing value after ", arg)
    }
    config[[value_options[[arg]]]] <- args[[i + 1L]]
    i <- i + 2L
  }

  config$alpha <- suppressWarnings(as.numeric(config$alpha))
  if (!is.finite(config$alpha) || config$alpha <= 0 || config$alpha >= 1) {
    stop("--alpha must be a number strictly between 0 and 1.")
  }
  config
}

validate_pairwise_table <- function(results) {
  required <- c(
    "Outcome", "Comparison", "Group_1", "Group_2", "N_group_1",
    "N_group_2", "t", "P_permutation_two_sided", "Hedges_g_adjusted"
  )
  missing <- setdiff(required, names(results))
  if (length(missing)) {
    stop("Input table is missing columns: ", paste(missing, collapse = ", "))
  }

  if (nrow(results) != 50L) {
    stop("Expected 50 rows (10 lambdas x 5 contrasts); found ", nrow(results), ".")
  }
  if (!setequal(results$Outcome, EXPECTED_OUTCOMES)) {
    stop("Input table does not contain the expected ten turbulence lambdas.")
  }
  if (!setequal(results$Comparison, EXPECTED_COMPARISONS)) {
    stop("Input table does not contain the expected five planned contrasts.")
  }

  cell_counts <- table(results$Outcome, results$Comparison)
  if (!all(cell_counts == 1L)) {
    stop("Each lambda-by-comparison combination must occur exactly once.")
  }

  p_values <- suppressWarnings(as.numeric(results$P_permutation_two_sided))
  if (anyNA(p_values) || any(!is.finite(p_values))) {
    stop("Permutation p-values must all be finite numeric values.")
  }
  if (any(p_values < 0 | p_values > 1)) {
    stop("Permutation p-values must lie between 0 and 1.")
  }

  invisible(TRUE)
}

recalculate_fdr <- function(results, alpha) {
  # Remove the former comparison-wise correction from this alternative table
  # so downstream code cannot accidentally select the wrong FDR definition.
  old_fdr_columns <- intersect(
    c("FDR_BH_within_comparison", "Significant_FDR_0_05"),
    names(results)
  )
  if (length(old_fdr_columns)) {
    results[old_fdr_columns] <- NULL
  }

  results[[NEW_FDR_COLUMN]] <- ave(
    results$P_permutation_two_sided,
    results$Outcome,
    FUN = function(values) stats::p.adjust(values, method = "BH")
  )
  results[[NEW_SIG_COLUMN]] <- results[[NEW_FDR_COLUMN]] < alpha
  results
}

print_summary <- function(results) {
  significant <- results[results[[NEW_SIG_COLUMN]], c(
    "Outcome", "Comparison", "P_permutation_two_sided", NEW_FDR_COLUMN,
    "Hedges_g_adjusted"
  )]

  cat("\nFDR-significant tests after correcting five contrasts within each lambda:\n")
  if (!nrow(significant)) {
    cat("  None\n")
  } else {
    print(significant, row.names = FALSE)
  }
}

main <- function() {
  config <- parse_args(commandArgs(trailingOnly = TRUE))
  input_file <- normalizePath(config$input_file, mustWork = TRUE)
  results <- utils::read.csv(input_file, check.names = FALSE)
  validate_pairwise_table(results)

  cat("Validated 50 pairwise tests: 10 lambdas x 5 planned AD-stage contrasts.\n")
  cat("FDR family: five planned contrasts separately within each lambda.\n")

  adjusted <- recalculate_fdr(results, config$alpha)
  print_summary(adjusted)

  if (config$validate_only) {
    cat("\nValidation completed; no output file was written.\n")
    return(invisible(adjusted))
  }

  dir.create(dirname(config$output_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(adjusted, config$output_file, row.names = FALSE)
  cat("\nSaved lambda-wise FDR results to: ", config$output_file, "\n", sep = "")
  invisible(adjusted)
}

main()
