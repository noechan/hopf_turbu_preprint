#!/usr/bin/env Rscript

# N145 turbulence group comparisons with configurable nuisance covariates.
#
# This is the active replacement for the incomplete legacy
# turbu_R_harm/code/Turbu_harm_ADNI3.R workflow. It reads the current ComBat
# table directly, joins the harmonization metadata by PTID, and runs
# Freedman-Lane permutation models using permuco.

options(stringsAsFactors = FALSE)

EXPECTED_GROUP_COUNTS <- c(
  HC_ABneg = 51L,
  HC_ABpos = 37L,
  MCI_ABpos = 31L,
  AD_ABpos = 26L
)

GROUP_LEVELS <- names(EXPECTED_GROUP_COUNTS)

# These are the five contrasts used by the active N145 MATLAB group analysis.
PAIRWISE_COMPARISONS <- list(
  HC_ABneg_vs_HC_ABpos = c("HC_ABneg", "HC_ABpos"),
  HC_ABpos_vs_MCI_ABpos = c("HC_ABpos", "MCI_ABpos"),
  MCI_ABpos_vs_AD_ABpos = c("MCI_ABpos", "AD_ABpos"),
  HC_ABneg_vs_MCI_ABpos = c("HC_ABneg", "MCI_ABpos"),
  HC_ABneg_vs_AD_ABpos = c("HC_ABneg", "AD_ABpos")
)

usage <- function(script_name = "run_group_permutations_N145.R") {
  cat(paste0(
    "Usage: Rscript ", script_name, " [options]\n\n",
    "Options:\n",
    "  --harmonized-file PATH  Current N145 all-feature ComBat table (.xlsx/.csv).\n",
    "  --metadata-file PATH    CSV with PTID, age, gender, and edu.\n",
    "  --output-dir PATH       Destination for aggregate CSV results.\n",
    "  --n-perm INTEGER        Permutations per model (default: 100000).\n",
    "  --seed INTEGER          Random seed (default: 2025).\n",
    "  --include-education     Also adjust for years of education (edu).\n",
    "  --validate-only         Validate inputs without fitting models.\n",
    "  --help                  Show this message.\n"
  ))
}

parse_args <- function(args) {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) == 1L) {
    script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
  } else {
    # RStudio Source does not provide --file. The project working directory is
    # the directory containing this script, so use that as the first fallback.
    project_candidate <- file.path(getwd(), "run_group_permutations_N145.R")
    if (file.exists(project_candidate)) {
      script_path <- normalizePath(project_candidate, mustWork = TRUE)
    } else if (
      interactive() &&
      requireNamespace("rstudioapi", quietly = TRUE) &&
      nzchar(rstudioapi::getSourceEditorContext()$path)
    ) {
      script_path <- normalizePath(
        rstudioapi::getSourceEditorContext()$path,
        mustWork = TRUE
      )
    } else {
      stop(
        "Could not resolve this script's path. Run it with Rscript or open ",
        "turbu_harm_stats.Rproj and click Source."
      )
    }
  }
  project_dir <- dirname(script_path)
  statistical_dir <- dirname(project_dir)
  sch1000_root <- dirname(statistical_dir)

  adni3_root <- Sys.getenv(
    "ADNI3_ROOT",
    unset = "/path/to/ADNI3"
  )

  config <- list(
    sch1000_root = sch1000_root,
    adni3_root = adni3_root,
    harmonized_file = file.path(
      sch1000_root,
      "harmonization_allfeat",
      "recomputed",
      "Turbu_ComBat_ADNI3_allfeatures_N145.xlsx"
    ),
    metadata_file = file.path(
      sch1000_root,
      "data",
      "covariates",
      "covariates_ADNI3_ABeta_N152.csv"
    ),
    output_dir = file.path(project_dir, "results", "N145_group_permutations"),
    n_perm = 100000L,
    seed = 2025L,
    include_education = FALSE,
    validate_only = FALSE
  )

  value_options <- c(
    "--harmonized-file" = "harmonized_file",
    "--metadata-file" = "metadata_file",
    "--output-dir" = "output_dir",
    "--n-perm" = "n_perm",
    "--seed" = "seed"
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--help") {
      usage(basename(script_path))
      quit(status = 0L)
    }
    if (arg == "--include-education") {
      config$include_education <- TRUE
      i <- i + 1L
      next
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
  config
}

portable_path <- function(path, root, root_label) {
  resolved <- normalizePath(path, mustWork = TRUE)
  resolved_root <- normalizePath(root, mustWork = TRUE)
  prefix <- paste0(resolved_root, .Platform$file.sep)
  if (startsWith(resolved, prefix)) {
    relative <- substring(resolved, nchar(prefix) + 1L)
    return(file.path(root_label, relative))
  }
  resolved
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Missing R package '", package, "'. Install the active R dependencies with:\n",
      "  install.packages(c(\"readxl\", \"permuco\"))"
    )
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

require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

load_analysis_data <- function(
  harmonized_file,
  metadata_file,
  include_education = FALSE
) {
  if (!file.exists(harmonized_file)) {
    stop("Missing harmonized N145 workbook: ", harmonized_file)
  }
  if (!file.exists(metadata_file)) {
    stop("Missing demographic metadata: ", metadata_file)
  }

  extension <- tolower(tools::file_ext(harmonized_file))
  if (extension == "xlsx") {
    require_package("readxl")
    harmonized <- as.data.frame(readxl::read_excel(harmonized_file))
  } else if (extension == "csv") {
    harmonized <- utils::read.csv(harmonized_file, check.names = FALSE)
  } else {
    stop("Unsupported harmonized input extension: .", extension, ". Use .xlsx or .csv.")
  }
  metadata <- utils::read.csv(metadata_file, check.names = FALSE)
  require_columns(harmonized, c("PTID", "Group"), "Harmonized workbook")
  required_metadata <- c("PTID", "age", "gender")
  if (include_education) {
    required_metadata <- c(required_metadata, "edu")
  }
  require_columns(metadata, required_metadata, "Metadata CSV")

  harmonized$PTID <- normalize_ptid(harmonized$PTID, "Harmonized workbook")
  metadata$PTID <- normalize_ptid(metadata$PTID, "Metadata CSV")

  observed_counts <- table(factor(harmonized$Group, levels = GROUP_LEVELS))
  observed_counts <- stats::setNames(as.integer(observed_counts), GROUP_LEVELS)
  if (!identical(observed_counts, EXPECTED_GROUP_COUNTS)) {
    stop(
      "The harmonized workbook is not the canonical N145 cohort. Expected ",
      paste(names(EXPECTED_GROUP_COUNTS), EXPECTED_GROUP_COUNTS, sep = "=", collapse = ", "),
      "; found ",
      paste(names(observed_counts), observed_counts, sep = "=", collapse = ", "),
      "."
    )
  }

  metadata_index <- match(harmonized$PTID, metadata$PTID)
  if (anyNA(metadata_index)) {
    stop(
      "PTIDs missing from demographic metadata: ",
      paste(harmonized$PTID[is.na(metadata_index)], collapse = ", ")
    )
  }

  data <- harmonized
  data$Age <- suppressWarnings(as.numeric(metadata$age[metadata_index]))
  gender <- suppressWarnings(as.numeric(metadata$gender[metadata_index]))
  if (anyNA(data$Age) || anyNA(gender)) {
    stop("Age or gender could not be converted to numeric values.")
  }
  if (!all(gender %in% c(0, 1))) {
    stop("Expected gender to contain only 0 and 1.")
  }
  data$Sex <- factor(gender, levels = c(0, 1), labels = c("Female", "Male"))
  if (include_education) {
    data$Education <- suppressWarnings(as.numeric(metadata$edu[metadata_index]))
    if (anyNA(data$Education)) {
      stop("Education could not be converted to complete numeric values.")
    }
  }
  data$Group <- factor(data$Group, levels = GROUP_LEVELS)

  outcomes <- grep("^Turbu_lam_", names(data), value = TRUE)
  expected_outcomes <- paste0(
    "Turbu_lam_",
    c("0_27", "0_24", "0_21", "0_18", "0_15", "0_12", "0_09", "0_06", "0_03", "0_01")
  )
  if (!identical(outcomes, expected_outcomes)) {
    stop(
      "Unexpected turbulence columns or order. Expected: ",
      paste(expected_outcomes, collapse = ", "),
      "; found: ", paste(outcomes, collapse = ", ")
    )
  }

  for (outcome in outcomes) {
    data[[outcome]] <- suppressWarnings(as.numeric(data[[outcome]]))
  }
  analysis_columns <- c(outcomes, "Age", "Sex")
  if (include_education) {
    analysis_columns <- c(analysis_columns, "Education")
  }
  analysis_columns <- c(analysis_columns, "Group")
  if (anyNA(data[analysis_columns])) {
    missing_counts <- colSums(is.na(data[analysis_columns]))
    missing_counts <- missing_counts[missing_counts > 0]
    stop(
      "Missing analysis values: ",
      paste(names(missing_counts), missing_counts, sep = "=", collapse = ", ")
    )
  }

  list(
    data = data,
    outcomes = outcomes,
    counts = observed_counts,
    include_education = include_education
  )
}

model_formula <- function(outcome, include_education) {
  nuisance_terms <- c("Age", "Sex")
  if (include_education) {
    nuisance_terms <- c(nuisance_terms, "Education")
  }
  stats::as.formula(
    paste(outcome, "~ Group +", paste(nuisance_terms, collapse = " + "))
  )
}

find_p_column <- function(table, candidates) {
  hit <- candidates[candidates %in% colnames(table)]
  if (!length(hit)) {
    stop(
      "Could not find a resampled p-value column. Available columns: ",
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
    stop("Could not locate term '", term, "' in permutation table.")
  }
  fallback[[1L]]
}

run_omnibus <- function(data, outcomes, n_perm, include_education = FALSE) {
  results <- vector("list", length(outcomes))
  for (index in seq_along(outcomes)) {
    outcome <- outcomes[[index]]
    message("Omnibus: ", outcome)
    formula <- model_formula(outcome, include_education)
    # coding_sum=TRUE requests marginal main-effect tests for the factor.
    fit <- permuco::aovperm(
      formula = formula,
      data = data,
      np = n_perm,
      method = "freedman_lane",
      type = "permutation",
      coding_sum = TRUE
    )
    table <- fit$table
    row <- find_term_row(table, "Group")
    p_column <- find_p_column(
      table,
      c("resampled P(>F)", "resampled Pr(>F)", "P(resampled)")
    )
    f_column <- intersect(c("F", "F value"), colnames(table))
    results[[index]] <- data.frame(
      Outcome = outcome,
      N = nrow(data),
      F = if (length(f_column)) as.numeric(table[row, f_column[[1L]]]) else NA_real_,
      P_permutation = as.numeric(table[row, p_column]),
      stringsAsFactors = FALSE
    )
  }
  result <- do.call(rbind, results)
  result$FDR_BH_across_10_lambdas <- stats::p.adjust(
    result$P_permutation,
    method = "BH"
  )
  result$Significant_FDR_0_05 <- result$FDR_BH_across_10_lambdas < 0.05
  result
}

adjusted_hedges_g <- function(
  model,
  data,
  levels_two,
  include_education = FALSE
) {
  mean_age <- mean(data$Age)
  modal_sex <- names(which.max(table(data$Sex)))[[1L]]
  prediction_data <- data.frame(
    Group = factor(levels_two, levels = levels_two),
    Age = rep(mean_age, 2L),
    Sex = factor(rep(modal_sex, 2L), levels = levels(data$Sex))
  )
  if (include_education) {
    prediction_data$Education <- rep(mean(data$Education), 2L)
  }
  adjusted_means <- stats::predict(model, newdata = prediction_data)
  residual_df <- stats::df.residual(model)
  correction <- 1 - 3 / (4 * residual_df - 1)
  correction * ((adjusted_means[[1L]] - adjusted_means[[2L]]) / stats::sigma(model))
}

run_pairwise <- function(data, outcomes, n_perm, include_education = FALSE) {
  rows <- list()
  row_index <- 1L

  for (comparison_name in names(PAIRWISE_COMPARISONS)) {
    levels_two <- PAIRWISE_COMPARISONS[[comparison_name]]
    pair_data <- data[data$Group %in% levels_two, , drop = FALSE]
    pair_data$Group <- factor(as.character(pair_data$Group), levels = levels_two)
    pair_data$Sex <- droplevels(pair_data$Sex)

    for (outcome in outcomes) {
      message("Pairwise ", comparison_name, ": ", outcome)
      formula <- model_formula(outcome, include_education)
      fit <- permuco::lmperm(
        formula = formula,
        data = pair_data,
        np = n_perm,
        method = "freedman_lane",
        type = "permutation"
      )
      table <- fit$table
      row <- find_term_row(table, "Group")
      p_column <- find_p_column(
        table,
        c("resampled Pr(>|t|)", "resampled P(>F)", "P(resampled)")
      )
      t_column <- intersect(c("t", "t value"), colnames(table))

      parametric_model <- stats::lm(formula, data = pair_data)
      rows[[row_index]] <- data.frame(
        Outcome = outcome,
        Comparison = comparison_name,
        Group_1 = levels_two[[1L]],
        Group_2 = levels_two[[2L]],
        N_group_1 = sum(pair_data$Group == levels_two[[1L]]),
        N_group_2 = sum(pair_data$Group == levels_two[[2L]]),
        t = if (length(t_column)) as.numeric(table[row, t_column[[1L]]]) else NA_real_,
        P_permutation_two_sided = as.numeric(table[row, p_column]),
        Hedges_g_adjusted = as.numeric(
          adjusted_hedges_g(
            parametric_model,
            pair_data,
            levels_two,
            include_education
          )
        ),
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }

  result <- do.call(rbind, rows)
  result$FDR_BH_within_comparison <- ave(
    result$P_permutation_two_sided,
    result$Comparison,
    FUN = function(values) stats::p.adjust(values, method = "BH")
  )
  result$Significant_FDR_0_05 <- result$FDR_BH_within_comparison < 0.05
  result
}

write_results <- function(config, loaded, omnibus, pairwise) {
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  model_suffix <- if (config$include_education) {
    "age_sex_education"
  } else {
    "age_sex"
  }
  model_label <- if (config$include_education) {
    "Outcome ~ Group + Age + Sex + Education"
  } else {
    "Outcome ~ Group + Age + Sex"
  }
  summary_filename <- if (config$include_education) {
    "N145_group_permutation_run_summary_age_sex_education.csv"
  } else {
    "N145_group_permutation_run_summary.csv"
  }

  utils::write.csv(
    omnibus,
    file.path(
      config$output_dir,
      paste0("N145_turbulence_omnibus_FreedmanLane_", model_suffix, ".csv")
    ),
    row.names = FALSE
  )
  utils::write.csv(
    pairwise,
    file.path(
      config$output_dir,
      paste0("N145_turbulence_pairwise_FreedmanLane_", model_suffix, ".csv")
    ),
    row.names = FALSE
  )

  run_summary <- data.frame(
    Setting = c(
      "Harmonized_input",
      "Metadata_input",
      "Model",
      "Permutation_method",
      "Permutations",
      "Seed",
      "Omnibus_FDR_family",
      "Pairwise_FDR_family",
      "HC_ABneg_N",
      "HC_ABpos_N",
      "MCI_ABpos_N",
      "AD_ABpos_N"
    ),
    Value = c(
      portable_path(config$harmonized_file, config$sch1000_root, "$SCH1000_ROOT"),
      portable_path(config$metadata_file, config$sch1000_root, "$SCH1000_ROOT"),
      model_label,
      "Freedman-Lane",
      as.character(config$n_perm),
      as.character(config$seed),
      "BH across 10 turbulence lambdas",
      "BH across 10 turbulence lambdas, separately within each comparison",
      as.character(loaded$counts[["HC_ABneg"]]),
      as.character(loaded$counts[["HC_ABpos"]]),
      as.character(loaded$counts[["MCI_ABpos"]]),
      as.character(loaded$counts[["AD_ABpos"]])
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    run_summary,
    file.path(config$output_dir, summary_filename),
    row.names = FALSE
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  config <- parse_args(args)
  loaded <- load_analysis_data(
    config$harmonized_file,
    config$metadata_file,
    config$include_education
  )

  cat(
    "Validated canonical N145 cohort: ",
    paste(names(loaded$counts), loaded$counts, sep = "=", collapse = " | "),
    "\n",
    sep = ""
  )
  covariate_label <- if (config$include_education) {
    "age, sex, and education"
  } else {
    "age and sex"
  }
  cat(
    "Validated 10 turbulence outcomes with complete ",
    covariate_label,
    " covariates.\n",
    sep = ""
  )

  if (config$validate_only) {
    cat("Validation completed; no permutation models were fitted and no files were written.\n")
    return(invisible(NULL))
  }

  require_package("permuco")
  set.seed(config$seed)
  options(contrasts = c("contr.sum", "contr.poly"))

  omnibus <- run_omnibus(
    loaded$data,
    loaded$outcomes,
    config$n_perm,
    config$include_education
  )
  pairwise <- run_pairwise(
    loaded$data,
    loaded$outcomes,
    config$n_perm,
    config$include_education
  )
  write_results(config, loaded, omnibus, pairwise)

  cat("\nOmnibus results:\n")
  print(omnibus, row.names = FALSE)
  cat("\nPairwise FDR-significant tests by comparison:\n")
  print(
    aggregate(
      Significant_FDR_0_05 ~ Comparison,
      data = pairwise,
      FUN = sum
    ),
    row.names = FALSE
  )
  cat("\nSaved aggregate results to: ", config$output_dir, "\n", sep = "")
}

if (!isTRUE(getOption("turbu_harm_stats.skip_main"))) {
  main()
}
