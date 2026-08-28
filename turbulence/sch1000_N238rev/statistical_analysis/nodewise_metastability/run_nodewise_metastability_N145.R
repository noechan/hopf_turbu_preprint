#!/usr/bin/env Rscript

# Covariate-adjusted node-level metastability/turbulence analysis for the
# canonical harmonized N145 cohort.
#
# Each Schaefer-1000 parcel is tested with a two-sided Freedman-Lane
# permutation test. The primary model adjusts for age and sex:
#
#   NodeMetastability ~ Group + Age + Sex
#
# The --include-education option adds Education as a secondary sensitivity
# covariate. Benjamini-Hochberg FDR is applied across the 1,000 parcels
# separately within each planned contrast.

options(stringsAsFactors = FALSE)

EXPECTED_GROUP_COUNTS <- c(
  HC_ABneg = 51L,
  HC_ABpos = 37L,
  MCI_ABpos = 31L,
  AD_ABpos = 26L
)
GROUP_LEVELS <- names(EXPECTED_GROUP_COUNTS)

PLANNED_CONTRASTS <- list(
  HC_ABneg_vs_AD_ABpos = c("HC_ABneg", "AD_ABpos"),
  MCI_ABpos_vs_AD_ABpos = c("MCI_ABpos", "AD_ABpos")
)

usage <- function(script_name = "run_nodewise_metastability_N145.R") {
  cat(paste0(
    "Usage: Rscript ", script_name, " [options]\n\n",
    "Options:\n",
    "  --node-file PATH       Harmonized N145 Schaefer-1000 node table (.xlsx/.csv).\n",
    "  --metadata-file PATH   CSV with PTID, age, gender, and edu.\n",
    "  --output-dir PATH      Destination for aggregate CSV results.\n",
    "  --n-perm INTEGER       Permutations per contrast (default: 10000).\n",
    "  --seed INTEGER         Random seed, reset per contrast (default: 1).\n",
    "  --include-education    Add years of education as a nuisance covariate.\n",
    "  --validate-only        Validate inputs without fitting models.\n",
    "  --help                 Show this message.\n"
  ))
}

resolve_script_path <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) == 1L) {
    return(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
  }

  project_candidate <- file.path(getwd(), "run_nodewise_metastability_N145.R")
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
    "Could not resolve this script's path. Run it with Rscript, or open it ",
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
    node_file = file.path(
      sch1000_root,
      "harmonization_allfeat",
      "recomputed",
      "Turbu_ComBat_ADNI3_HC_MCI_AD_ABeta_lambda_0_01_sch1000_N145.xlsx"
    ),
    metadata_file = file.path(
      sch1000_root,
      "data",
      "covariates",
      "covariates_ADNI3_ABeta_N152.csv"
    ),
    output_dir = file.path(
      analysis_dir,
      "results",
      "N145_nodewise_metastability"
    ),
    n_perm = 10000L,
    seed = 1L,
    include_education = FALSE,
    validate_only = FALSE
  )

  value_options <- c(
    "--node-file" = "node_file",
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

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Missing R package '", package, "'. Install it with:\n",
      "  install.packages(\"", package, "\")"
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

load_node_table <- function(path) {
  if (!file.exists(path)) {
    stop("Missing harmonized node table: ", path)
  }
  extension <- tolower(tools::file_ext(path))
  if (extension == "xlsx") {
    require_package("readxl")
    data <- as.data.frame(readxl::read_excel(path), check.names = FALSE)
  } else if (extension == "csv") {
    data <- utils::read.csv(path, check.names = FALSE)
  } else {
    stop("Unsupported node input extension: .", extension, ". Use .xlsx or .csv.")
  }
  data
}

load_analysis_data <- function(node_file, metadata_file, include_education) {
  nodes <- load_node_table(node_file)
  if (!file.exists(metadata_file)) {
    stop("Missing harmonization metadata: ", metadata_file)
  }
  metadata <- utils::read.csv(metadata_file, check.names = FALSE)

  require_columns(nodes, c("PTID", "Group"), "Node table")
  required_metadata <- c("PTID", "Group", "age", "gender")
  if (include_education) {
    required_metadata <- c(required_metadata, "edu")
  }
  require_columns(metadata, required_metadata, "Metadata CSV")

  nodes$PTID <- normalize_ptid(nodes$PTID, "Node table")
  metadata$PTID <- normalize_ptid(metadata$PTID, "Metadata CSV")

  observed_counts <- table(factor(nodes$Group, levels = GROUP_LEVELS))
  observed_counts <- stats::setNames(as.integer(observed_counts), GROUP_LEVELS)
  if (!identical(observed_counts, EXPECTED_GROUP_COUNTS)) {
    stop(
      "Node table is not the canonical N145 cohort. Expected ",
      paste(names(EXPECTED_GROUP_COUNTS), EXPECTED_GROUP_COUNTS,
        sep = "=", collapse = ", "
      ),
      "; found ",
      paste(names(observed_counts), observed_counts, sep = "=", collapse = ", "),
      "."
    )
  }

  node_columns <- grep("^Schaefer_", names(nodes), value = TRUE)
  if (length(node_columns) != 1000L) {
    stop("Expected 1,000 Schaefer node columns; found ", length(node_columns), ".")
  }
  node_matrix <- as.matrix(nodes[, node_columns, drop = FALSE])
  storage.mode(node_matrix) <- "double"
  if (any(!is.finite(node_matrix))) {
    stop("Node table contains non-finite Schaefer values.")
  }

  metadata_index <- match(nodes$PTID, metadata$PTID)
  if (anyNA(metadata_index)) {
    stop(
      "PTIDs missing from harmonization metadata: ",
      paste(nodes$PTID[is.na(metadata_index)], collapse = ", ")
    )
  }
  metadata_group <- as.character(metadata$Group[metadata_index])
  group_mismatch <- nodes$PTID[metadata_group != nodes$Group]
  if (length(group_mismatch)) {
    stop(
      "Group labels differ between node and metadata inputs for: ",
      paste(group_mismatch, collapse = ", ")
    )
  }

  age <- suppressWarnings(as.numeric(metadata$age[metadata_index]))
  sex <- suppressWarnings(as.numeric(metadata$gender[metadata_index]))
  if (anyNA(age) || anyNA(sex)) {
    stop("Age or sex contains missing/non-numeric values for N145.")
  }
  if (!all(sex %in% c(0, 1))) {
    stop("Expected metadata gender coding 0=Female and 1=Male.")
  }

  covariates <- data.frame(Age = age, Sex = sex)
  if (include_education) {
    education <- suppressWarnings(as.numeric(metadata$edu[metadata_index]))
    if (anyNA(education)) {
      stop("Education contains missing/non-numeric values for N145.")
    }
    covariates$Education <- education
  }

  list(
    PTID = nodes$PTID,
    Group = factor(nodes$Group, levels = GROUP_LEVELS),
    node_names = node_columns,
    Y = node_matrix,
    covariates = covariates,
    observed_counts = observed_counts
  )
}

center_covariates <- function(covariates) {
  centered <- as.data.frame(lapply(covariates, function(values) {
    values - mean(values)
  }))
  names(centered) <- names(covariates)
  centered
}

# Vectorized Freedman-Lane implementation for multiple node outcomes.
# A common row permutation is used for all parcels at each iteration, which
# preserves the observed cross-parcel covariance structure and avoids fitting
# 1,000 separate permutation models.
freedman_lane_nodewise <- function(Y, group_indicator, nuisance, n_perm, seed) {
  Y <- as.matrix(Y)
  nuisance <- as.matrix(nuisance)
  n_subjects <- nrow(Y)
  n_nodes <- ncol(Y)

  X0 <- cbind(Intercept = 1, nuisance)
  X <- cbind(X0, Group = group_indicator)
  if (qr(X0)$rank != ncol(X0) || qr(X)$rank != ncol(X)) {
    stop("Reduced or full design matrix is rank deficient.")
  }

  qr_x0 <- qr(X0)
  coefficients_reduced <- qr.coef(qr_x0, Y)
  fitted_reduced <- X0 %*% coefficients_reduced
  residuals_reduced <- Y - fitted_reduced

  group_residual <- as.numeric(qr.resid(qr_x0, group_indicator))
  group_ss <- sum(group_residual^2)
  if (!is.finite(group_ss) || group_ss <= .Machine$double.eps) {
    stop("Group contrast is not estimable after nuisance adjustment.")
  }

  numerator_observed <- as.numeric(crossprod(group_residual, Y))
  beta_observed <- numerator_observed / group_ss
  sse_reduced_observed <- colSums(residuals_reduced^2)
  ss_group_observed <- numerator_observed^2 / group_ss
  residual_df <- n_subjects - ncol(X)
  sse_full_observed <- pmax(
    sse_reduced_observed - ss_group_observed,
    .Machine$double.eps
  )
  t_observed <- numerator_observed /
    sqrt((sse_full_observed / residual_df) * group_ss)

  xtx0_inverse <- solve(crossprod(X0))
  exceedances <- integer(n_nodes)
  progress_every <- max(1L, floor(n_perm / 10L))
  set.seed(seed)

  for (permutation_index in seq_len(n_perm)) {
    permutation <- sample.int(n_subjects, n_subjects, replace = FALSE)
    permuted_residuals <- residuals_reduced[permutation, , drop = FALSE]

    x0_cross_y <- crossprod(X0, permuted_residuals)
    reduced_coefficients_perm <- xtx0_inverse %*% x0_cross_y
    sse_reduced_perm <- colSums(permuted_residuals^2) -
      colSums(reduced_coefficients_perm * x0_cross_y)

    numerator_perm <- as.numeric(crossprod(group_residual, permuted_residuals))
    ss_group_perm <- numerator_perm^2 / group_ss
    sse_full_perm <- pmax(
      sse_reduced_perm - ss_group_perm,
      .Machine$double.eps
    )
    t_perm <- numerator_perm /
      sqrt((sse_full_perm / residual_df) * group_ss)
    exceedances <- exceedances + (abs(t_perm) >= abs(t_observed))

    if (
      permutation_index %% progress_every == 0L ||
      permutation_index == n_perm
    ) {
      cat("  completed ", permutation_index, "/", n_perm, " permutations\n", sep = "")
    }
  }

  p_values <- (exceedances + 1) / (n_perm + 1)
  partial_eta_squared <- t_observed^2 / (t_observed^2 + residual_df)

  list(
    beta = beta_observed,
    t = t_observed,
    p = p_values,
    partial_eta_squared = partial_eta_squared,
    residual_df = residual_df
  )
}

run_contrast <- function(data, contrast_name, groups, config) {
  group_1 <- groups[[1L]]
  group_2 <- groups[[2L]]
  keep <- data$Group %in% groups

  Y <- data$Y[keep, , drop = FALSE]
  group <- droplevels(data$Group[keep])
  group_indicator <- as.numeric(group == group_1)
  nuisance <- center_covariates(data$covariates[keep, , drop = FALSE])
  n_1 <- sum(group == group_1)
  n_2 <- sum(group == group_2)

  cat(
    "\n", contrast_name, ": ", group_1, " (N=", n_1, ") vs ",
    group_2, " (N=", n_2, ")\n",
    sep = ""
  )
  cat(
    "Model: NodeMetastability ~ Group + ",
    paste(names(nuisance), collapse = " + "),
    "\n",
    sep = ""
  )
  cat("Input: ", nrow(Y), " subjects x ", ncol(Y), " nodes\n", sep = "")

  fit <- freedman_lane_nodewise(
    Y = Y,
    group_indicator = group_indicator,
    nuisance = nuisance,
    n_perm = config$n_perm,
    seed = config$seed
  )

  mean_1 <- colMeans(Y[group == group_1, , drop = FALSE])
  mean_2 <- colMeans(Y[group == group_2, , drop = FALSE])
  adjusted_p <- stats::p.adjust(fit$p, method = "BH")

  results <- data.frame(
    Node = data$node_names,
    Group_1 = group_1,
    Group_2 = group_2,
    N_1 = n_1,
    N_2 = n_2,
    Mean_Group_1 = mean_1,
    Mean_Group_2 = mean_2,
    Raw_Mean_Difference = mean_1 - mean_2,
    Adjusted_Group_Beta = fit$beta,
    T = fit$t,
    Residual_DF = fit$residual_df,
    Partial_Eta_Squared = fit$partial_eta_squared,
    P_Permutation = fit$p,
    P_FDR_BH = adjusted_p,
    Significant_FDR_0_05 = adjusted_p < 0.05,
    stringsAsFactors = FALSE
  )

  suffix <- if (config$include_education) {
    "FreedmanLane_age_sex_education"
  } else {
    "FreedmanLane_age_sex"
  }
  output_file <- file.path(
    config$output_dir,
    paste0(
      "N145_nodewise_metastability_",
      contrast_name,
      "_",
      suffix,
      ".csv"
    )
  )
  utils::write.csv(results, output_file, row.names = FALSE)

  cat(
    "Significant parcels after BH-FDR: ",
    sum(results$Significant_FDR_0_05),
    "/",
    nrow(results),
    "\n",
    sep = ""
  )
  cat("Saved: ", output_file, "\n", sep = "")

  data.frame(
    Contrast = contrast_name,
    Group_1 = group_1,
    Group_2 = group_2,
    N_1 = n_1,
    N_2 = n_2,
    Significant_FDR_0_05 = sum(results$Significant_FDR_0_05),
    Minimum_P_Permutation = min(results$P_Permutation),
    Minimum_P_FDR_BH = min(results$P_FDR_BH),
    Maximum_P_FDR_BH = max(results$P_FDR_BH),
    Output_File = normalizePath(output_file, mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

write_run_summary <- function(config, data, contrast_summaries) {
  suffix <- if (config$include_education) "age_sex_education" else "age_sex"
  summary_file <- file.path(
    config$output_dir,
    paste0("N145_nodewise_metastability_run_summary_", suffix, ".csv")
  )

  model <- if (config$include_education) {
    "NodeMetastability ~ Group + Age + Sex + Education"
  } else {
    "NodeMetastability ~ Group + Age + Sex"
  }

  provenance <- data.frame(
    Item = c(
      "Node_input", "Metadata_input", "Cohort", "Scale", "Parcellation",
      "Model", "Permutation_method", "Permutation_statistic",
      "Permutations", "Random_seed", "Permutation_scope", "FDR_family"
    ),
    Value = c(
      normalizePath(config$node_file, mustWork = TRUE),
      normalizePath(config$metadata_file, mustWork = TRUE),
      paste(names(data$observed_counts), data$observed_counts,
        sep = "=", collapse = " | "
      ),
      "lambda=0.01",
      "Schaefer-1000",
      model,
      "Freedman-Lane residual permutation under the reduced nuisance model",
      "two-sided t statistic for Group",
      as.character(config$n_perm),
      as.character(config$seed),
      "one common participant-row permutation across all 1000 parcels per iteration",
      "Benjamini-Hochberg across 1000 parcels separately within each contrast"
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(provenance, summary_file, row.names = FALSE)

  contrast_file <- file.path(
    config$output_dir,
    paste0("N145_nodewise_metastability_contrast_summary_", suffix, ".csv")
  )
  utils::write.csv(contrast_summaries, contrast_file, row.names = FALSE)
  cat("Run summary: ", summary_file, "\n", sep = "")
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  config <- parse_args(args)
  data <- load_analysis_data(
    config$node_file,
    config$metadata_file,
    config$include_education
  )

  cat(
    "Validated canonical N145 node table: ",
    paste(names(data$observed_counts), data$observed_counts,
      sep = "=", collapse = " | "
    ),
    "\n",
    sep = ""
  )
  cat("Validated ", length(data$node_names), " Schaefer parcels at lambda=0.01.\n", sep = "")
  cat(
    "Covariates: ",
    paste(names(data$covariates), collapse = ", "),
    "\n",
    sep = ""
  )

  if (config$validate_only) {
    cat("Validation completed; no permutation models were fitted.\n")
    return(invisible(NULL))
  }

  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  summaries <- do.call(rbind, lapply(names(PLANNED_CONTRASTS), function(name) {
    run_contrast(data, name, PLANNED_CONTRASTS[[name]], config)
  }))
  write_run_summary(config, data, summaries)
  invisible(summaries)
}

if (!isTRUE(getOption("nodewise_metastability.skip_main"))) {
  main()
}
