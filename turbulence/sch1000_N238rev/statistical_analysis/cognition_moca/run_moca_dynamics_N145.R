#!/usr/bin/env Rscript

# Cross-sectional associations between harmonized brain dynamics and MOCA in
# the canonical N145 amyloid-status cohort.
#
# Two estimands are kept separate:
#   1. Total association, adjusted for age, sex, and education.
#   2. Association beyond disease stage, additionally adjusted for Group.
#
# Inference for the standardized dynamical predictor uses a two-sided
# Freedman--Lane residual-permutation test. The default confirmatory predictor
# set contains three prespecified global measures. An all-scale predictor set is
# available as an explicitly exploratory option.

options(stringsAsFactors = FALSE)

EXPECTED_GROUP_COUNTS <- c(
  HC_ABneg = 51L,
  HC_ABpos = 37L,
  MCI_ABpos = 31L,
  AD_ABpos = 26L
)
GROUP_LEVELS <- names(EXPECTED_GROUP_COUNTS)

usage <- function(script_name = "run_moca_dynamics_N145.R") {
  cat(paste0(
    "Usage: Rscript ", script_name, " [options]\n\n",
    "Options:\n",
    "  --harmonized-file PATH  Recomputed N145 all-feature ComBat table.\n",
    "  --clinical-file PATH    N145 clinical CSV containing PTID, Group, MOCA.\n",
    "  --metadata-file PATH    Harmonization metadata with age, gender, edu.\n",
    "  --output-dir PATH       Destination for result CSV files.\n",
    "  --n-perm INTEGER        Permutations per model (default: 100000).\n",
    "  --seed INTEGER          Random seed (default: 2025).\n",
    "  --predictor-set VALUE   primary or all-scales (default: primary).\n",
    "  --model-set VALUE       both, total, or group-adjusted (default: both).\n",
    "  --validate-only         Validate inputs without fitting or writing.\n",
    "  --help                  Show this message.\n"
  ))
}

resolve_script_path <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) == 1L) {
    return(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
  }

  project_candidate <- file.path(getwd(), "run_moca_dynamics_N145.R")
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
    "Could not resolve this script's path. Run it with Rscript or open it ",
    "in RStudio and click Source."
  )
}

parse_args <- function(args) {
  script_path <- resolve_script_path()
  analysis_dir <- dirname(script_path)
  statistical_dir <- dirname(analysis_dir)
  sch1000_root <- dirname(statistical_dir)
  adni3_root <- Sys.getenv(
    "ADNI3_ROOT",
    unset = "/path/to/ADNI3"
  )

  config <- list(
    script_path = script_path,
    analysis_dir = analysis_dir,
    sch1000_root = sch1000_root,
    adni3_root = adni3_root,
    harmonized_file = file.path(
      sch1000_root,
      "harmonization_allfeat",
      "recomputed",
      "Turbu_ComBat_ADNI3_allfeatures_N145.xlsx"
    ),
    clinical_file = file.path(
      sch1000_root,
      "visualization",
      "python",
      "data",
      "Turbu_ComBat_clin_ADNI3_allfeatures_N145.csv"
    ),
    metadata_file = file.path(
      sch1000_root,
      "data",
      "covariates",
      "covariates_ADNI3_ABeta_N152.csv"
    ),
    output_dir = file.path(analysis_dir, "results", "N145_MOCA_associations"),
    n_perm = 100000L,
    seed = 2025L,
    predictor_set = "primary",
    model_set = "both",
    validate_only = FALSE
  )

  value_options <- c(
    "--harmonized-file" = "harmonized_file",
    "--clinical-file" = "clinical_file",
    "--metadata-file" = "metadata_file",
    "--output-dir" = "output_dir",
    "--n-perm" = "n_perm",
    "--seed" = "seed",
    "--predictor-set" = "predictor_set",
    "--model-set" = "model_set"
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      usage(basename(script_path))
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

  config$n_perm <- suppressWarnings(as.integer(config$n_perm))
  config$seed <- suppressWarnings(as.integer(config$seed))
  if (is.na(config$n_perm) || config$n_perm < 1L) {
    stop("--n-perm must be a positive integer.")
  }
  if (is.na(config$seed)) {
    stop("--seed must be an integer.")
  }
  if (!config$predictor_set %in% c("primary", "all-scales")) {
    stop("--predictor-set must be primary or all-scales.")
  }
  if (!config$model_set %in% c("both", "total", "group-adjusted")) {
    stop("--model-set must be both, total, or group-adjusted.")
  }
  config
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Missing R package '", package, "'. Install dependencies with:\n",
      "  install.packages(c(\"readxl\", \"permuco\"))"
    )
  }
}

require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

normalize_ptid <- function(values, label) {
  result <- trimws(as.character(values))
  if (anyNA(result) || any(result == "")) {
    stop(label, " contains an empty PTID.")
  }
  duplicates <- unique(result[duplicated(result)])
  if (length(duplicates)) {
    stop(label, " contains duplicate PTIDs: ", paste(duplicates, collapse = ", "))
  }
  result
}

read_harmonized <- function(path) {
  if (!file.exists(path)) {
    stop("Missing harmonized N145 input: ", path)
  }
  extension <- tolower(tools::file_ext(path))
  if (extension == "xlsx") {
    require_package("readxl")
    return(as.data.frame(readxl::read_excel(path)))
  }
  if (extension == "csv") {
    return(utils::read.csv(path, check.names = FALSE))
  }
  stop("Unsupported harmonized input: .", extension, ". Use .xlsx or .csv.")
}

validate_group_counts <- function(groups, label) {
  observed <- table(factor(groups, levels = GROUP_LEVELS))
  observed <- stats::setNames(as.integer(observed), GROUP_LEVELS)
  if (!identical(observed, EXPECTED_GROUP_COUNTS)) {
    stop(
      label, " is not the canonical N145 cohort. Expected ",
      paste(names(EXPECTED_GROUP_COUNTS), EXPECTED_GROUP_COUNTS,
        sep = "=", collapse = ", "
      ),
      "; found ",
      paste(names(observed), observed, sep = "=", collapse = ", "),
      "."
    )
  }
  observed
}

build_predictor_spec <- function(harmonized_names, predictor_set) {
  turbulence <- grep("^Turbu_lam_", harmonized_names, value = TRUE)
  info_flow <- grep("^InfoFlow_lam_", harmonized_names, value = TRUE)
  info_transfer <- grep("^InfoTransfer_lam_", harmonized_names, value = TRUE)
  singleton <- c("InfoCascade", "Metastability", "gKoP")
  missing_singleton <- setdiff(singleton, harmonized_names)
  if (length(missing_singleton)) {
    stop("Harmonized input is missing: ", paste(missing_singleton, collapse = ", "))
  }

  source <- c(turbulence, info_flow, info_transfer, singleton)
  family <- c(
    rep("Turbulence", length(turbulence)),
    rep("InformationFlow", length(info_flow)),
    rep("OneMinusInformationTransfer", length(info_transfer)),
    "InformationCascade",
    "Metastability",
    "GlobalKuramotoOrder"
  )
  transform <- c(
    rep("identity", length(turbulence) + length(info_flow)),
    rep("one_minus", length(info_transfer)),
    rep("identity", length(singleton))
  )
  predictor <- source
  transfer_rows <- family == "OneMinusInformationTransfer"
  predictor[transfer_rows] <- sub(
    "^InfoTransfer_",
    "OneMinus_InfoTransfer_",
    predictor[transfer_rows]
  )

  scale <- rep(NA_real_, length(source))
  scaled_rows <- grepl("_lam_", source)
  scale[scaled_rows] <- suppressWarnings(as.numeric(gsub(
    "_", ".", sub("^.*_lam_", "", source[scaled_rows])
  )))

  primary_sources <- c(
    "Turbu_lam_0_01",
    "InfoFlow_lam_0_01",
    "InfoTransfer_lam_0_01"
  )
  missing_primary <- setdiff(primary_sources, source)
  if (length(missing_primary)) {
    stop("Primary predictor set is incomplete: ", paste(missing_primary, collapse = ", "))
  }

  spec <- data.frame(
    Predictor = predictor,
    Source_variable = source,
    Predictor_family = family,
    Scale_lambda = scale,
    Transform = transform,
    Primary_predictor = source %in% primary_sources,
    stringsAsFactors = FALSE
  )
  if (predictor_set == "primary") {
    spec <- spec[spec$Primary_predictor, , drop = FALSE]
  }
  rownames(spec) <- NULL
  spec
}

load_analysis_data <- function(config) {
  harmonized <- read_harmonized(config$harmonized_file)
  if (!file.exists(config$clinical_file)) {
    stop("Missing N145 clinical CSV: ", config$clinical_file)
  }
  if (!file.exists(config$metadata_file)) {
    stop("Missing harmonization metadata CSV: ", config$metadata_file)
  }
  clinical <- utils::read.csv(
    config$clinical_file,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN")
  )
  metadata <- utils::read.csv(
    config$metadata_file,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN")
  )

  require_columns(harmonized, c("PTID", "Group"), "Harmonized input")
  require_columns(clinical, c("PTID", "Group", "MOCA"), "Clinical input")
  require_columns(
    metadata,
    c("PTID", "Group", "age", "gender", "edu"),
    "Metadata input"
  )

  harmonized$PTID <- normalize_ptid(harmonized$PTID, "Harmonized input")
  clinical$PTID <- normalize_ptid(clinical$PTID, "Clinical input")
  metadata$PTID <- normalize_ptid(metadata$PTID, "Metadata input")
  counts <- validate_group_counts(harmonized$Group, "Harmonized input")

  clinical_index <- match(harmonized$PTID, clinical$PTID)
  metadata_index <- match(harmonized$PTID, metadata$PTID)
  if (anyNA(clinical_index)) {
    stop(
      "PTIDs missing from clinical input: ",
      paste(harmonized$PTID[is.na(clinical_index)], collapse = ", ")
    )
  }
  if (anyNA(metadata_index)) {
    stop(
      "PTIDs missing from metadata input: ",
      paste(harmonized$PTID[is.na(metadata_index)], collapse = ", ")
    )
  }

  harmonized_group <- as.character(harmonized$Group)
  clinical_group <- as.character(clinical$Group[clinical_index])
  metadata_group <- as.character(metadata$Group[metadata_index])
  group_mismatch <- harmonized$PTID[
    harmonized_group != clinical_group | harmonized_group != metadata_group
  ]
  if (length(group_mismatch)) {
    stop(
      "Group labels disagree across inputs for: ",
      paste(group_mismatch, collapse = ", ")
    )
  }

  data <- harmonized
  data$Group <- factor(harmonized_group, levels = GROUP_LEVELS)
  data$MOCA <- suppressWarnings(as.numeric(clinical$MOCA[clinical_index]))
  data$Age <- suppressWarnings(as.numeric(metadata$age[metadata_index]))
  gender <- suppressWarnings(as.numeric(metadata$gender[metadata_index]))
  data$Education <- suppressWarnings(as.numeric(metadata$edu[metadata_index]))
  if (anyNA(data$Age) || anyNA(gender) || anyNA(data$Education)) {
    stop("Age, gender, and education must be complete numeric values for N145.")
  }
  if (!all(gender %in% c(0, 1))) {
    stop("Metadata gender must contain only 0 and 1.")
  }
  data$Sex <- factor(gender, levels = c(0, 1), labels = c("Female", "Male"))

  spec <- build_predictor_spec(names(harmonized), config$predictor_set)
  for (source in spec$Source_variable) {
    data[[source]] <- suppressWarnings(as.numeric(data[[source]]))
  }
  predictor_missing <- colSums(is.na(data[spec$Source_variable]))
  predictor_missing <- predictor_missing[predictor_missing > 0L]
  if (length(predictor_missing)) {
    stop(
      "Dynamical predictors contain missing/non-numeric values: ",
      paste(names(predictor_missing), predictor_missing, sep = "=", collapse = ", ")
    )
  }

  complete_moca <- !is.na(data$MOCA)
  if (sum(complete_moca) < 20L) {
    stop("Fewer than 20 participants have MOCA data.")
  }
  complete_counts <- table(factor(data$Group[complete_moca], levels = GROUP_LEVELS))
  complete_counts <- stats::setNames(as.integer(complete_counts), GROUP_LEVELS)

  list(
    data = data,
    predictor_spec = spec,
    canonical_counts = counts,
    moca_counts = complete_counts,
    n_moca = sum(complete_moca),
    n_moca_missing = sum(!complete_moca),
    moca_missing_ptids = data$PTID[!complete_moca]
  )
}

find_p_column <- function(table) {
  candidates <- c(
    "resampled Pr(>|t|)",
    "resampled P(>|t|)",
    "resampled P(>F)",
    "P(resampled)"
  )
  hit <- candidates[candidates %in% colnames(table)]
  if (!length(hit)) {
    stop(
      "Could not identify the permutation p-value column. Available: ",
      paste(colnames(table), collapse = ", ")
    )
  }
  hit[[1L]]
}

find_term_row <- function(table, term) {
  exact <- which(rownames(table) == term)
  if (length(exact) == 1L) {
    return(exact)
  }
  fallback <- grep(paste0("^", term), rownames(table))
  if (!length(fallback)) {
    stop("Could not find term '", term, "' in permutation results.")
  }
  fallback[[1L]]
}

model_definitions <- function(model_set) {
  definitions <- list(
    Total_association = list(
      formula = "MOCA ~ Predictor_z + Age + Sex + Education",
      reduced = "MOCA ~ Age + Sex + Education",
      estimand = "Total association adjusted for age, sex, and education"
    ),
    Group_adjusted = list(
      formula = "MOCA ~ Predictor_z + Age + Sex + Education + Group",
      reduced = "MOCA ~ Age + Sex + Education + Group",
      estimand = "Association beyond amyloid-status disease stage"
    )
  )
  if (model_set == "total") {
    return(definitions["Total_association"])
  }
  if (model_set == "group-adjusted") {
    return(definitions["Group_adjusted"])
  }
  definitions
}

run_models <- function(loaded, config) {
  require_package("permuco")
  definitions <- model_definitions(config$model_set)
  rows <- list()
  row_index <- 1L

  for (model_name in names(definitions)) {
    definition <- definitions[[model_name]]
    for (i in seq_len(nrow(loaded$predictor_spec))) {
      spec <- loaded$predictor_spec[i, , drop = FALSE]
      source_values <- loaded$data[[spec$Source_variable]]
      predictor_values <- if (spec$Transform == "one_minus") {
        1 - source_values
      } else {
        source_values
      }

      analysis <- data.frame(
        MOCA = loaded$data$MOCA,
        Predictor_raw = predictor_values,
        Age = loaded$data$Age,
        Sex = loaded$data$Sex,
        Education = loaded$data$Education,
        Group = loaded$data$Group
      )
      analysis <- analysis[stats::complete.cases(analysis), , drop = FALSE]
      if (!is.finite(stats::sd(analysis$Predictor_raw)) ||
          stats::sd(analysis$Predictor_raw) == 0) {
        stop("Predictor has zero or invalid variance: ", spec$Predictor)
      }
      analysis$Predictor_z <- as.numeric(scale(analysis$Predictor_raw))
      analysis$Sex <- droplevels(analysis$Sex)
      analysis$Group <- droplevels(analysis$Group)

      full_formula <- stats::as.formula(definition$formula)
      reduced_formula <- stats::as.formula(definition$reduced)
      message(
        model_name, ": ", spec$Predictor,
        " (N=", nrow(analysis), ", permutations=", config$n_perm, ")"
      )

      permutation_fit <- permuco::lmperm(
        formula = full_formula,
        data = analysis,
        np = config$n_perm,
        method = "freedman_lane",
        type = "permutation"
      )
      permutation_table <- permutation_fit$table
      term_row <- find_term_row(permutation_table, "Predictor_z")
      p_column <- find_p_column(permutation_table)

      full_fit <- stats::lm(full_formula, data = analysis)
      reduced_fit <- stats::lm(reduced_formula, data = analysis)
      coefficient_table <- summary(full_fit)$coefficients
      coefficient <- coefficient_table["Predictor_z", "Estimate"]
      standard_error <- coefficient_table["Predictor_z", "Std. Error"]
      t_value <- coefficient_table["Predictor_z", "t value"]
      confidence_interval <- stats::confint(full_fit, "Predictor_z", level = 0.95)
      residual_df <- stats::df.residual(full_fit)
      partial_r2 <- t_value^2 / (t_value^2 + residual_df)
      full_r2 <- summary(full_fit)$r.squared
      reduced_r2 <- summary(reduced_fit)$r.squared
      sample_counts <- table(factor(analysis$Group, levels = GROUP_LEVELS))

      rows[[row_index]] <- data.frame(
        Model = model_name,
        Estimand = definition$estimand,
        Predictor = spec$Predictor,
        Source_variable = spec$Source_variable,
        Predictor_family = spec$Predictor_family,
        Scale_lambda = spec$Scale_lambda,
        Transform = spec$Transform,
        Primary_predictor = spec$Primary_predictor,
        N = nrow(analysis),
        HC_ABneg_N = as.integer(sample_counts[["HC_ABneg"]]),
        HC_ABpos_N = as.integer(sample_counts[["HC_ABpos"]]),
        MCI_ABpos_N = as.integer(sample_counts[["MCI_ABpos"]]),
        AD_ABpos_N = as.integer(sample_counts[["AD_ABpos"]]),
        MOCA_mean = mean(analysis$MOCA),
        MOCA_SD = stats::sd(analysis$MOCA),
        Beta_MOCA_per_predictor_SD = as.numeric(coefficient),
        SE = as.numeric(standard_error),
        CI95_low = as.numeric(confidence_interval[1L]),
        CI95_high = as.numeric(confidence_interval[2L]),
        t = as.numeric(t_value),
        DF_residual = as.integer(residual_df),
        Partial_R2 = as.numeric(partial_r2),
        Full_model_R2 = as.numeric(full_r2),
        Delta_R2_vs_covariates = as.numeric(full_r2 - reduced_r2),
        P_permutation_two_sided = as.numeric(
          permutation_table[term_row, p_column]
        ),
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }

  results <- do.call(rbind, rows)
  results$P_FDR_BH_across_predictors <- ave(
    results$P_permutation_two_sided,
    results$Model,
    FUN = function(values) stats::p.adjust(values, method = "BH")
  )
  family_key <- interaction(
    results$Model,
    results$Predictor_family,
    drop = TRUE
  )
  results$P_FDR_BH_within_family <- ave(
    results$P_permutation_two_sided,
    family_key,
    FUN = function(values) stats::p.adjust(values, method = "BH")
  )
  results$Significant_global_FDR_0_05 <-
    results$P_FDR_BH_across_predictors < 0.05
  results$Significant_within_family_FDR_0_05 <-
    results$P_FDR_BH_within_family < 0.05
  results
}

portable_path <- function(path, root, root_label) {
  resolved <- normalizePath(path, mustWork = TRUE)
  resolved_root <- normalizePath(root, mustWork = TRUE)
  prefix <- paste0(resolved_root, .Platform$file.sep)
  if (startsWith(resolved, prefix)) {
    return(file.path(root_label, substring(resolved, nchar(prefix) + 1L)))
  }
  resolved
}

write_results <- function(results, loaded, config) {
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  set_suffix <- gsub("-", "_", config$predictor_set)
  combined_file <- file.path(
    config$output_dir,
    paste0(
      "N145_MOCA_dynamics_associations_FreedmanLane_age_sex_education_",
      set_suffix,
      ".csv"
    )
  )
  utils::write.csv(results, combined_file, row.names = FALSE)

  for (model_name in unique(results$Model)) {
    model_suffix <- tolower(gsub("_", "_", model_name))
    utils::write.csv(
      results[results$Model == model_name, , drop = FALSE],
      file.path(
        config$output_dir,
        paste0(
          "N145_MOCA_dynamics_", model_suffix,
          "_FreedmanLane_age_sex_education_", set_suffix, ".csv"
        )
      ),
      row.names = FALSE
    )
  }

  summary <- data.frame(
    Setting = c(
      "Harmonized_dynamics_input",
      "Clinical_input",
      "Clinical_columns_used",
      "Metadata_input",
      "Metadata_columns_used",
      "Outcome",
      "Predictor_set",
      "Model_set",
      "Total_association_model",
      "Group_adjusted_model",
      "Permutation_method",
      "Permutations_per_model",
      "Seed",
      "Global_FDR_family",
      "Within_family_FDR",
      "Canonical_N",
      "MOCA_complete_N",
      "MOCA_missing_N",
      "MOCA_missing_PTIDs",
      "MOCA_HC_ABneg_N",
      "MOCA_HC_ABpos_N",
      "MOCA_MCI_ABpos_N",
      "MOCA_AD_ABpos_N"
    ),
    Value = c(
      portable_path(config$harmonized_file, config$sch1000_root, "$SCH1000_ROOT"),
      portable_path(config$clinical_file, config$sch1000_root, "$SCH1000_ROOT"),
      "PTID, Group, MOCA only; dynamical columns in this file are ignored",
      portable_path(config$metadata_file, config$sch1000_root, "$SCH1000_ROOT"),
      "PTID, Group, age, gender, edu",
      "MOCA",
      config$predictor_set,
      config$model_set,
      "MOCA ~ z(dynamics) + Age + Sex + Education",
      "MOCA ~ z(dynamics) + Age + Sex + Education + Group",
      "Freedman-Lane residual permutation; two-sided predictor test",
      as.character(config$n_perm),
      as.character(config$seed),
      paste0(
        "BH across all ", nrow(loaded$predictor_spec),
        " selected predictors, separately within each model"
      ),
      "BH across scales separately within each predictor family and model",
      as.character(sum(loaded$canonical_counts)),
      as.character(loaded$n_moca),
      as.character(loaded$n_moca_missing),
      paste(loaded$moca_missing_ptids, collapse = ";"),
      as.character(loaded$moca_counts[["HC_ABneg"]]),
      as.character(loaded$moca_counts[["HC_ABpos"]]),
      as.character(loaded$moca_counts[["MCI_ABpos"]]),
      as.character(loaded$moca_counts[["AD_ABpos"]])
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    summary,
    file.path(
      config$output_dir,
      paste0("N145_MOCA_dynamics_run_summary_", set_suffix, ".csv")
    ),
    row.names = FALSE
  )

  combined_file
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  config <- parse_args(args)
  loaded <- load_analysis_data(config)

  cat(
    "Validated canonical N145 dynamics and metadata: ",
    paste(
      names(loaded$canonical_counts), loaded$canonical_counts,
      sep = "=", collapse = " | "
    ),
    "\n",
    sep = ""
  )
  cat(
    "MOCA complete sample: N=", loaded$n_moca,
    " (missing=", loaded$n_moca_missing, ") | ",
    paste(names(loaded$moca_counts), loaded$moca_counts, sep = "=", collapse = " | "),
    "\n",
    sep = ""
  )
  cat(
    "Predictor set: ", config$predictor_set,
    " (", nrow(loaded$predictor_spec), " predictors)\n",
    sep = ""
  )

  if (config$validate_only) {
    cat("Validation completed; no models were fitted and no files were written.\n")
    return(invisible(NULL))
  }

  set.seed(config$seed)
  results <- run_models(loaded, config)
  output_file <- write_results(results, loaded, config)

  cat("\nMOCA--dynamics association results:\n")
  print(
    results[c(
      "Model", "Predictor", "N", "Beta_MOCA_per_predictor_SD",
      "P_permutation_two_sided", "P_FDR_BH_across_predictors"
    )],
    row.names = FALSE
  )
  cat("\nSaved combined results to: ", output_file, "\n", sep = "")
}

if (!isTRUE(getOption("cognition_moca.skip_main"))) {
  main()
}
