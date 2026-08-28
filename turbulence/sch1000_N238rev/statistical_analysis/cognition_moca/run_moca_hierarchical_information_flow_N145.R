#!/usr/bin/env Rscript

# Focused blockwise hierarchical regression for cross-sectional MOCA.
#
# M0: MOCA ~ Age + Sex + Education + Group
# M1: M0 + z(Turbulence at lambda=0.01)
# M2: M1 + z(Information flow at lambda=0.01)
#
# The confirmatory estimand is the incremental contribution of information
# flow in M2 beyond demographics, amyloid-status group, and turbulence.
# Freedman--Lane inference is retained for continuity with the repository;
# HC3-robust inference is also reported because residual heteroscedasticity was
# detected in the broader cognition models.

options(stringsAsFactors = FALSE)

resolve_this_script <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) == 1L) {
    return(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
  }
  if (
    interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    nzchar(rstudioapi::getSourceEditorContext()$path)
  ) {
    return(normalizePath(rstudioapi::getSourceEditorContext()$path, mustWork = TRUE))
  }
  stop("Run this script with Rscript or open it in RStudio and click Source.")
}

SCRIPT_PATH <- resolve_this_script()
ANALYSIS_DIR <- dirname(SCRIPT_PATH)
options(cognition_moca.skip_main = TRUE)
source(file.path(ANALYSIS_DIR, "run_moca_dynamics_N145.R"), local = .GlobalEnv)

usage <- function(script_name = basename(SCRIPT_PATH)) {
  cat(paste0(
    "Usage: Rscript ", script_name, " [options]\n\n",
    "Options:\n",
    "  --harmonized-file PATH  Recomputed N145 all-feature ComBat table.\n",
    "  --clinical-file PATH    N145 clinical CSV containing PTID, Group, MOCA.\n",
    "  --metadata-file PATH    Harmonization metadata with age, gender, edu.\n",
    "  --output-dir PATH       Destination for result CSV files.\n",
    "  --n-perm INTEGER        Permutations per block test (default: 100000).\n",
    "  --seed INTEGER          Random seed (default: 2025).\n",
    "  --validate-only         Validate inputs without fitting or writing.\n",
    "  --help                  Show this message.\n"
  ))
}

find_permutation_p <- function(fit, term) {
  table <- fit$table
  row <- find_term_row(table, term)
  column <- find_p_column(table)
  as.numeric(table[row, column])
}

standardize <- function(values, label) {
  value_sd <- stats::sd(values)
  if (!is.finite(value_sd) || value_sd == 0) {
    stop(label, " has zero or invalid variance.")
  }
  as.numeric(scale(values))
}

hc3_row <- function(model, term) {
  robust <- lmtest::coeftest(model, vcov. = sandwich::vcovHC(model, type = "HC3"))
  robust[term, , drop = TRUE]
}

model_metrics <- function(model, model_name, added_block, previous_r2 = NA_real_) {
  model_summary <- summary(model)
  data.frame(
    Model = model_name,
    Added_block = added_block,
    N = stats::nobs(model),
    Parameters = length(stats::coef(model)),
    Residual_DF = stats::df.residual(model),
    R2 = model_summary$r.squared,
    Adjusted_R2 = model_summary$adj.r.squared,
    Delta_R2 = if (is.na(previous_r2)) NA_real_ else model_summary$r.squared - previous_r2,
    RMSE = sqrt(mean(stats::residuals(model)^2)),
    AIC = stats::AIC(model),
    stringsAsFactors = FALSE
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if ("--help" %in% args) {
    usage()
    return(invisible(NULL))
  }

  output_dir_supplied <- "--output-dir" %in% args
  config <- parse_args(args)
  config$predictor_set <- "primary"
  config$model_set <- "group-adjusted"
  if (!output_dir_supplied) {
    config$output_dir <- file.path(
      config$analysis_dir,
      "results",
      "N145_MOCA_hierarchical_information_flow"
    )
  }

  for (package in c("permuco", "lmtest", "sandwich", "car")) {
    require_package(package)
  }
  loaded <- load_analysis_data(config)

  required_dynamics <- c("Turbu_lam_0_01", "InfoFlow_lam_0_01")
  require_columns(loaded$data, required_dynamics, "Harmonized input")
  analysis <- data.frame(
    PTID = loaded$data$PTID,
    Group = loaded$data$Group,
    MOCA = loaded$data$MOCA,
    Age = loaded$data$Age,
    Sex = loaded$data$Sex,
    Education = loaded$data$Education,
    Turbulence_raw = as.numeric(loaded$data$Turbu_lam_0_01),
    InformationFlow_raw = as.numeric(loaded$data$InfoFlow_lam_0_01),
    stringsAsFactors = FALSE
  )
  analysis <- analysis[stats::complete.cases(analysis), , drop = FALSE]
  analysis$Group <- droplevels(analysis$Group)
  analysis$Sex <- droplevels(analysis$Sex)
  analysis$Turbulence_z <- standardize(analysis$Turbulence_raw, "Turbulence")
  analysis$InformationFlow_z <- standardize(
    analysis$InformationFlow_raw,
    "Information flow"
  )

  group_counts <- table(factor(analysis$Group, levels = GROUP_LEVELS))
  cat(
    "Validated hierarchical MOCA sample: N=", nrow(analysis), " | ",
    paste(names(group_counts), as.integer(group_counts), sep = "=", collapse = " | "),
    "\n", sep = ""
  )
  cat(
    "Predictor correlation: r(Turbulence, InformationFlow)=",
    sprintf("%.4f", stats::cor(analysis$Turbulence_z, analysis$InformationFlow_z)),
    "\n", sep = ""
  )
  if (config$validate_only) {
    cat("Validation completed; no models were fitted and no files were written.\n")
    return(invisible(NULL))
  }

  m0 <- stats::lm(MOCA ~ Age + Sex + Education + Group, data = analysis)
  m1 <- stats::lm(
    MOCA ~ Age + Sex + Education + Group + Turbulence_z,
    data = analysis
  )
  m2 <- stats::lm(
    MOCA ~ Age + Sex + Education + Group + Turbulence_z + InformationFlow_z,
    data = analysis
  )

  set.seed(config$seed)
  message("Testing M0 -> M1 with ", config$n_perm, " Freedman--Lane permutations")
  permutation_m1 <- permuco::lmperm(
    MOCA ~ Age + Sex + Education + Group + Turbulence_z,
    data = analysis,
    np = config$n_perm,
    method = "freedman_lane",
    type = "permutation"
  )
  message("Testing M1 -> M2 with ", config$n_perm, " Freedman--Lane permutations")
  permutation_m2 <- permuco::lmperm(
    MOCA ~ Age + Sex + Education + Group + Turbulence_z + InformationFlow_z,
    data = analysis,
    np = config$n_perm,
    method = "freedman_lane",
    type = "permutation"
  )

  m0_summary <- summary(m0)
  m1_summary <- summary(m1)
  m2_summary <- summary(m2)
  models <- rbind(
    model_metrics(m0, "M0", "Age + Sex + Education + Group"),
    model_metrics(m1, "M1", "Turbulence", m0_summary$r.squared),
    model_metrics(m2, "M2", "Information flow", m1_summary$r.squared)
  )

  nested_m01 <- stats::anova(m0, m1)
  nested_m12 <- stats::anova(m1, m2)
  block_terms <- c("Turbulence_z", "InformationFlow_z")
  block_models <- list(m1, m2)
  block_previous <- list(m0, m1)
  block_permutation <- list(permutation_m1, permutation_m2)
  block_names <- c("M0_to_M1_add_turbulence", "M1_to_M2_add_information_flow")
  nested_tables <- list(nested_m01, nested_m12)

  block_rows <- lapply(seq_along(block_terms), function(i) {
    model <- block_models[[i]]
    previous <- block_previous[[i]]
    term <- block_terms[[i]]
    ordinary <- summary(model)$coefficients[term, ]
    robust <- hc3_row(model, term)
    robust_critical <- stats::qt(0.975, df = stats::df.residual(model))
    data.frame(
      Block_test = block_names[[i]],
      Added_term = term,
      N = stats::nobs(model),
      Beta_MOCA_per_predictor_SD = ordinary[["Estimate"]],
      OLS_SE = ordinary[["Std. Error"]],
      OLS_t = ordinary[["t value"]],
      OLS_p = ordinary[["Pr(>|t|)"]],
      HC3_SE = robust[["Std. Error"]],
      HC3_t = robust[["t value"]],
      HC3_p = robust[["Pr(>|t|)"]],
      HC3_CI95_low = ordinary[["Estimate"]] - robust_critical * robust[["Std. Error"]],
      HC3_CI95_high = ordinary[["Estimate"]] + robust_critical * robust[["Std. Error"]],
      Delta_R2 = summary(model)$r.squared - summary(previous)$r.squared,
      Nested_F = nested_tables[[i]][2L, "F"],
      Nested_F_p = nested_tables[[i]][2L, "Pr(>F)"],
      P_FreedmanLane_two_sided = find_permutation_p(block_permutation[[i]], term),
      Permutations = config$n_perm,
      stringsAsFactors = FALSE
    )
  })
  block_tests <- do.call(rbind, block_rows)

  coefficient_names <- rownames(m2_summary$coefficients)
  robust_m2 <- lmtest::coeftest(m2, vcov. = sandwich::vcovHC(m2, type = "HC3"))
  coefficients <- data.frame(
    Term = coefficient_names,
    Estimate = m2_summary$coefficients[, "Estimate"],
    OLS_SE = m2_summary$coefficients[, "Std. Error"],
    OLS_t = m2_summary$coefficients[, "t value"],
    OLS_p = m2_summary$coefficients[, "Pr(>|t|)"],
    HC3_SE = robust_m2[coefficient_names, "Std. Error"],
    HC3_t = robust_m2[coefficient_names, "t value"],
    HC3_p = robust_m2[coefficient_names, "Pr(>|t|)"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  residuals_m2 <- stats::residuals(m2)
  cook <- stats::cooks.distance(m2)
  leverage <- stats::hatvalues(m2)
  studentized <- stats::rstudent(m2)
  p <- length(stats::coef(m2))
  n <- stats::nobs(m2)
  vif_table <- car::vif(m2)
  if (is.matrix(vif_table)) {
    vif_turbulence <- vif_table["Turbulence_z", "GVIF"]
    vif_flow <- vif_table["InformationFlow_z", "GVIF"]
  } else {
    vif_turbulence <- vif_table[["Turbulence_z"]]
    vif_flow <- vif_table[["InformationFlow_z"]]
  }
  shapiro <- stats::shapiro.test(residuals_m2)
  bp <- lmtest::bptest(m2)
  reset <- lmtest::resettest(m2, power = 2:3, type = "fitted")
  diagnostics <- data.frame(
    Diagnostic = c(
      "N", "Unique_PTIDs", "MOCA_missing_PTID", "Predictor_correlation",
      "Design_rank", "Design_columns", "VIF_Turbulence", "VIF_InformationFlow",
      "RESET_p", "Shapiro_Wilk_W", "Shapiro_Wilk_p", "Breusch_Pagan_p",
      "Cook_threshold_4_over_N", "Cook_flagged_N", "Cook_max",
      "Leverage_threshold_2p_over_N", "Leverage_flagged_N",
      "Externally_studentized_abs_gt_3_N"
    ),
    Value = as.character(c(
      n, length(unique(analysis$PTID)), paste(loaded$moca_missing_ptids, collapse = ";"),
      stats::cor(analysis$Turbulence_z, analysis$InformationFlow_z),
      qr(stats::model.matrix(m2))$rank, ncol(stats::model.matrix(m2)),
      vif_turbulence, vif_flow, reset$p.value, shapiro$statistic,
      shapiro$p.value, bp$p.value, 4 / n, sum(cook > 4 / n), max(cook),
      2 * p / n, sum(leverage > 2 * p / n), sum(abs(studentized) > 3)
    )),
    stringsAsFactors = FALSE
  )
  influence <- data.frame(
    PTID = analysis$PTID,
    Group = analysis$Group,
    Cook_D = cook,
    Leverage = leverage,
    Externally_studentized_residual = studentized,
    Cook_flag = cook > 4 / n,
    Leverage_flag = leverage > 2 * p / n,
    Studentized_abs_gt_3 = abs(studentized) > 3,
    stringsAsFactors = FALSE
  )
  influence <- influence[
    influence$Cook_flag | influence$Leverage_flag | influence$Studentized_abs_gt_3,
    , drop = FALSE
  ]

  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  files <- c(
    models = file.path(config$output_dir, "N145_MOCA_hierarchical_information_flow_models.csv"),
    block_tests = file.path(config$output_dir, "N145_MOCA_hierarchical_information_flow_block_tests.csv"),
    coefficients = file.path(config$output_dir, "N145_MOCA_hierarchical_information_flow_M2_coefficients.csv"),
    diagnostics = file.path(config$output_dir, "N145_MOCA_hierarchical_information_flow_diagnostics.csv"),
    influence = file.path(config$output_dir, "N145_MOCA_hierarchical_information_flow_influence.csv")
  )
  utils::write.csv(models, files[["models"]], row.names = FALSE)
  utils::write.csv(block_tests, files[["block_tests"]], row.names = FALSE)
  utils::write.csv(coefficients, files[["coefficients"]], row.names = FALSE)
  utils::write.csv(diagnostics, files[["diagnostics"]], row.names = FALSE)
  utils::write.csv(influence, files[["influence"]], row.names = FALSE)

  cat("\nHierarchical model summary:\n")
  print(models, row.names = FALSE)
  cat("\nBlock tests:\n")
  print(block_tests, row.names = FALSE)
  cat("\nSaved results under: ", config$output_dir, "\n", sep = "")
  invisible(list(models = models, block_tests = block_tests, diagnostics = diagnostics))
}

main()
