#!/usr/bin/env Rscript

# N145 group comparisons for information flow, information cascade, and
# one-minus information transfer. The validated Freedman-Lane implementation
# is shared with the active turbulence analysis, while outcome validation,
# multiplicity correction, and output names remain measure-specific here.

options(stringsAsFactors = FALSE)

resolve_information_script_dir <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) == 1L) {
    return(dirname(normalizePath(
      sub("^--file=", "", script_arg),
      mustWork = TRUE
    )))
  }

  project_candidate <- file.path(
    getwd(),
    "run_information_measures_permutations_N145.R"
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

information_script_dir <- resolve_information_script_dir()
core_script <- file.path(
  information_script_dir,
  "run_group_permutations_N145.R"
)
if (!file.exists(core_script)) {
  stop("Missing shared permutation implementation: ", core_script)
}

old_skip_option <- getOption("turbu_harm_stats.skip_main")
options(turbu_harm_stats.skip_main = TRUE)
source(core_script, local = .GlobalEnv)
options(turbu_harm_stats.skip_main = old_skip_option)

LAMBDA_ALL <- c(
  "0_27", "0_24", "0_21", "0_18", "0_15",
  "0_12", "0_09", "0_06", "0_03", "0_01"
)
LAMBDA_FLOW <- LAMBDA_ALL[-1L]

EXPECTED_INFORMATION_FLOW <- paste0("InfoFlow_lam_", LAMBDA_FLOW)
EXPECTED_INFORMATION_TRANSFER <- paste0("InfoTransfer_lam_", LAMBDA_ALL)
ONE_MINUS_INFORMATION_TRANSFER <- paste0(
  "OneMinusInfoTransfer_lam_",
  LAMBDA_ALL
)
EXPECTED_INFORMATION_CASCADE <- "InfoCascade"

validate_numeric_features <- function(data, columns, label) {
  require_columns(data, columns, "Harmonized workbook")
  for (column in columns) {
    data[[column]] <- suppressWarnings(as.numeric(data[[column]]))
  }
  values <- as.matrix(data[columns])
  if (anyNA(values)) {
    missing_counts <- colSums(is.na(values))
    missing_counts <- missing_counts[missing_counts > 0L]
    stop(
      label, " contains missing or non-numeric values: ",
      paste(names(missing_counts), missing_counts, sep = "=", collapse = ", ")
    )
  }
  if (any(!is.finite(values))) {
    stop(label, " contains non-finite values.")
  }
  data
}

prepare_information_outcomes <- function(data) {
  data <- validate_numeric_features(
    data,
    EXPECTED_INFORMATION_FLOW,
    "Information flow"
  )
  data <- validate_numeric_features(
    data,
    EXPECTED_INFORMATION_TRANSFER,
    "Information transfer"
  )
  data <- validate_numeric_features(
    data,
    EXPECTED_INFORMATION_CASCADE,
    "Information cascade"
  )

  for (index in seq_along(EXPECTED_INFORMATION_TRANSFER)) {
    raw_column <- EXPECTED_INFORMATION_TRANSFER[[index]]
    transformed_column <- ONE_MINUS_INFORMATION_TRANSFER[[index]]
    data[[transformed_column]] <- 1 - data[[raw_column]]
  }

  transformed_values <- as.matrix(data[ONE_MINUS_INFORMATION_TRANSFER])
  if (anyNA(transformed_values) || any(!is.finite(transformed_values))) {
    stop("One-minus information transfer contains invalid values.")
  }

  list(
    data = data,
    measures = list(
      information_flow = list(
        label = "Information flow",
        outcomes = EXPECTED_INFORMATION_FLOW
      ),
      one_minus_information_transfer = list(
        label = "1 - Information transfer",
        outcomes = ONE_MINUS_INFORMATION_TRANSFER
      ),
      information_cascade = list(
        label = "Information cascade",
        outcomes = EXPECTED_INFORMATION_CASCADE
      )
    )
  )
}

rename_omnibus_fdr <- function(results) {
  old_name <- "FDR_BH_across_10_lambdas"
  new_name <- "FDR_BH_across_scales_within_measure"
  if (!old_name %in% names(results)) {
    stop("Shared omnibus results are missing the expected FDR column.")
  }
  names(results)[names(results) == old_name] <- new_name
  results
}

apply_stage_fdr_within_scale <- function(results) {
  remove_columns <- intersect(
    c("FDR_BH_within_comparison", "Significant_FDR_0_05"),
    names(results)
  )
  if (length(remove_columns)) {
    results[remove_columns] <- NULL
  }

  fdr_column <- "FDR_BH_across_5_stage_contrasts_within_scale"
  significance_column <- "Significant_FDR_within_scale_0_05"
  results[[fdr_column]] <- ave(
    results$P_permutation_two_sided,
    results$Outcome,
    FUN = function(values) stats::p.adjust(values, method = "BH")
  )
  results[[significance_column]] <- results[[fdr_column]] < 0.05

  counts <- table(results$Outcome)
  if (any(counts != length(PAIRWISE_COMPARISONS))) {
    stop("Each outcome must have exactly five planned stage contrasts.")
  }
  results
}

model_suffix <- function(include_education) {
  if (include_education) "age_sex_education" else "age_sex"
}

write_information_results <- function(
  config,
  measure_name,
  omnibus,
  pairwise
) {
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  suffix <- model_suffix(config$include_education)

  utils::write.csv(
    omnibus,
    file.path(
      config$output_dir,
      paste0(
        "N145_", measure_name,
        "_omnibus_FreedmanLane_", suffix, ".csv"
      )
    ),
    row.names = FALSE
  )
  utils::write.csv(
    pairwise,
    file.path(
      config$output_dir,
      paste0(
        "N145_", measure_name,
        "_pairwise_FreedmanLane_", suffix,
        "_FDR_by_scale.csv"
      )
    ),
    row.names = FALSE
  )
}

write_information_run_summary <- function(config, loaded) {
  suffix <- model_suffix(config$include_education)
  model <- if (config$include_education) {
    "Outcome ~ Group + Age + Sex + Education"
  } else {
    "Outcome ~ Group + Age + Sex"
  }

  summary <- data.frame(
    Setting = c(
      "Harmonized_input",
      "Metadata_input",
      "Model",
      "Permutation_method",
      "Permutations",
      "Seed",
      "Information_flow_outcomes",
      "One_minus_information_transfer_outcomes",
      "Information_cascade_outcomes",
      "Information_transfer_transformation",
      "Omnibus_FDR_family",
      "Pairwise_FDR_family",
      "HC_ABneg_N",
      "HC_ABpos_N",
      "MCI_ABpos_N",
      "AD_ABpos_N"
    ),
    Value = c(
      portable_path(config$harmonized_file, config$sch1000_root, "$SCH1000_ROOT"),
      portable_path(config$metadata_file, config$adni3_root, "$ADNI3_ROOT"),
      model,
      "Freedman-Lane",
      as.character(config$n_perm),
      as.character(config$seed),
      as.character(length(EXPECTED_INFORMATION_FLOW)),
      as.character(length(ONE_MINUS_INFORMATION_TRANSFER)),
      "1",
      "1 - harmonized InfoTransfer value",
      "BH across scales separately within each measure",
      "BH across five planned stage contrasts separately within each scale",
      as.character(loaded$counts[["HC_ABneg"]]),
      as.character(loaded$counts[["HC_ABpos"]]),
      as.character(loaded$counts[["MCI_ABpos"]]),
      as.character(loaded$counts[["AD_ABpos"]])
    ),
    stringsAsFactors = FALSE
  )

  utils::write.csv(
    summary,
    file.path(
      config$output_dir,
      paste0(
        "N145_information_measures_permutation_run_summary_",
        suffix,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}

print_information_summary <- function(measure, pairwise) {
  significance_column <- "Significant_FDR_within_scale_0_05"
  significant <- pairwise[pairwise[[significance_column]], c(
    "Outcome", "Comparison", "P_permutation_two_sided",
    "FDR_BH_across_5_stage_contrasts_within_scale",
    "Hedges_g_adjusted"
  )]

  cat("\n", measure$label, " FDR-significant pairwise tests:\n", sep = "")
  if (!nrow(significant)) {
    cat("  None\n")
  } else {
    print(significant, row.names = FALSE)
  }
}

main_information <- function(args = commandArgs(trailingOnly = TRUE)) {
  config <- parse_args(args)
  loaded <- load_analysis_data(
    config$harmonized_file,
    config$metadata_file,
    config$include_education
  )
  prepared <- prepare_information_outcomes(loaded$data)
  loaded$data <- prepared$data

  cat(
    "Validated canonical N145 cohort: ",
    paste(names(loaded$counts), loaded$counts, sep = "=", collapse = " | "),
    "\n",
    sep = ""
  )
  cat(
    "Validated information outcomes: flow=9 | 1-transfer=10 | cascade=1.\n"
  )
  cat(
    "Information transfer analysis uses 1 - harmonized InfoTransfer values.\n"
  )
  covariate_label <- if (config$include_education) {
    "age, sex, and education"
  } else {
    "age and sex"
  }
  cat("Covariates: ", covariate_label, ".\n", sep = "")

  if (config$validate_only) {
    cat(
      "Validation completed; no permutation models were fitted and no files were written.\n"
    )
    return(invisible(NULL))
  }

  require_package("permuco")
  set.seed(config$seed)
  options(contrasts = c("contr.sum", "contr.poly"))

  for (measure_name in names(prepared$measures)) {
    measure <- prepared$measures[[measure_name]]
    cat(
      "\nRunning ", measure$label, ": ",
      length(measure$outcomes), " outcome(s).\n",
      sep = ""
    )

    omnibus <- run_omnibus(
      loaded$data,
      measure$outcomes,
      config$n_perm,
      config$include_education
    )
    omnibus <- rename_omnibus_fdr(omnibus)

    pairwise <- run_pairwise(
      loaded$data,
      measure$outcomes,
      config$n_perm,
      config$include_education
    )
    pairwise <- apply_stage_fdr_within_scale(pairwise)

    write_information_results(
      config,
      measure_name,
      omnibus,
      pairwise
    )
    print_information_summary(measure, pairwise)
  }

  write_information_run_summary(config, loaded)
  cat("\nSaved information-measure results to: ", config$output_dir, "\n", sep = "")
}

if (!isTRUE(getOption("turbu_harm_stats.skip_information_main"))) {
  main_information()
}
