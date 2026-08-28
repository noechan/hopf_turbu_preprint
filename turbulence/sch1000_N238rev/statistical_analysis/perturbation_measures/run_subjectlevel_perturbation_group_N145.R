#!/usr/bin/env Rscript

# Education-adjusted group comparisons for participant-level information
# capacity and susceptibility in the canonical N145 amyloid-status cohort.
#
# Both outcomes are read from the PTID-matched combined export produced by
# data_export/data_for_ML_4Staging_with_InfoCap_Susceptibility.m. Age, sex,
# and education are joined by PTID from the canonical ComBat metadata.

options(stringsAsFactors = FALSE)

EXPECTED_GROUP_COUNTS <- c(
  HC_ABneg = 51L,
  HC_ABpos = 37L,
  MCI_ABpos = 31L,
  AD_ABpos = 26L
)
GROUP_LEVELS <- names(EXPECTED_GROUP_COUNTS)
OUTCOMES <- c("Info_Cap", "Susceptibility")

PAIRWISE_COMPARISONS <- list(
  HC_ABneg_vs_HC_ABpos = c("HC_ABneg", "HC_ABpos"),
  HC_ABpos_vs_MCI_ABpos = c("HC_ABpos", "MCI_ABpos"),
  MCI_ABpos_vs_AD_ABpos = c("MCI_ABpos", "AD_ABpos"),
  HC_ABneg_vs_MCI_ABpos = c("HC_ABneg", "MCI_ABpos"),
  HC_ABneg_vs_AD_ABpos = c("HC_ABneg", "AD_ABpos")
)

usage <- function(script_name = "run_subjectlevel_perturbation_group_N145.R") {
  cat(paste0(
    "Usage: Rscript ", script_name, " [options]\n\n",
    "Options:\n",
    "  --hopf-file PATH      N145 PTID/Group/Info_Cap/Susceptibility table.\n",
    "  --metadata-file PATH  Canonical CSV with PTID, Group, age, gender, edu.\n",
    "  --output-dir PATH     Destination for aggregate result CSV files.\n",
    "  --n-perm INTEGER      Permutations per model (default: 100000).\n",
    "  --seed INTEGER        Random seed (default: 2025).\n",
    "  --validate-only       Validate inputs without fitting or writing.\n",
    "  --help                Show this message.\n"
  ))
}

resolve_script_path <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) == 1L) {
    return(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
  }
  project_candidate <- file.path(
    getwd(), "run_subjectlevel_perturbation_group_N145.R"
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
    hopf_file = file.path(
      sch1000_root,
      "data_export",
      paste0(
        "ML_Input_ADNI3_4STAGINGBYABETA_ComBat_N145_",
        "with_infocap_suscep_sch1000.csv"
      )
    ),
    metadata_file = file.path(
      sch1000_root,
      "data",
      "covariates",
      "covariates_ADNI3_ABeta_N152.csv"
    ),
    output_dir = file.path(
      analysis_dir, "results", "N145_subjectlevel_perturbation"
    ),
    n_perm = 100000L,
    seed = 2025L,
    validate_only = FALSE
  )

  value_options <- c(
    "--hopf-file" = "hopf_file",
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

portable_path <- function(path, root, root_label) {
  resolved <- normalizePath(path, mustWork = TRUE)
  if (!dir.exists(root)) {
    return(resolved)
  }
  resolved_root <- normalizePath(root, mustWork = TRUE)
  prefix <- paste0(resolved_root, .Platform$file.sep)
  if (startsWith(resolved, prefix)) {
    return(file.path(root_label, substring(resolved, nchar(prefix) + 1L)))
  }
  resolved
}

load_analysis_data <- function(config) {
  if (!file.exists(config$hopf_file)) {
    stop(
      "Missing participant-level Hopf input: ", config$hopf_file, "\n",
      "Run data_export/data_for_ML_4Staging_with_InfoCap_Susceptibility.m."
    )
  }
  if (!file.exists(config$metadata_file)) {
    stop(
      "Missing canonical harmonisation metadata: ", config$metadata_file, "\n",
      "Mount ADNI or pass --metadata-file PATH."
    )
  }

  hopf <- utils::read.csv(
    config$hopf_file,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN")
  )
  metadata <- utils::read.csv(
    config$metadata_file,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN")
  )
  require_columns(
    hopf,
    c("PTID", "Group", OUTCOMES),
    "Participant-level Hopf input"
  )
  require_columns(
    metadata,
    c("PTID", "Group", "age", "gender", "edu"),
    "Canonical metadata"
  )

  hopf$PTID <- normalize_ptid(hopf$PTID, "Participant-level Hopf input")
  metadata$PTID <- normalize_ptid(metadata$PTID, "Canonical metadata")
  counts <- table(factor(hopf$Group, levels = GROUP_LEVELS))
  counts <- stats::setNames(as.integer(counts), GROUP_LEVELS)
  if (!identical(counts, EXPECTED_GROUP_COUNTS)) {
    stop(
      "Hopf input is not the canonical N145 cohort. Expected ",
      paste(names(EXPECTED_GROUP_COUNTS), EXPECTED_GROUP_COUNTS,
        sep = "=", collapse = ", "
      ),
      "; found ", paste(names(counts), counts, sep = "=", collapse = ", "), "."
    )
  }

  metadata_index <- match(hopf$PTID, metadata$PTID)
  if (anyNA(metadata_index)) {
    stop(
      "PTIDs missing from canonical metadata: ",
      paste(hopf$PTID[is.na(metadata_index)], collapse = ", ")
    )
  }
  metadata_group <- as.character(metadata$Group[metadata_index])
  group_mismatch <- hopf$PTID[as.character(hopf$Group) != metadata_group]
  if (length(group_mismatch)) {
    stop(
      "Group labels disagree between inputs for: ",
      paste(group_mismatch, collapse = ", ")
    )
  }

  data <- hopf[c("PTID", "Group", OUTCOMES)]
  for (outcome in OUTCOMES) {
    data[[outcome]] <- suppressWarnings(as.numeric(data[[outcome]]))
  }
  data$Age <- suppressWarnings(as.numeric(metadata$age[metadata_index]))
  gender <- suppressWarnings(as.numeric(metadata$gender[metadata_index]))
  data$Education <- suppressWarnings(as.numeric(metadata$edu[metadata_index]))
  if (anyNA(data[c(OUTCOMES, "Age", "Education")]) || anyNA(gender)) {
    missing_counts <- colSums(is.na(data[c(OUTCOMES, "Age", "Education")]))
    stop(
      "Missing or non-numeric analysis values: ",
      paste(names(missing_counts), missing_counts, sep = "=", collapse = ", ")
    )
  }
  if (!all(gender %in% c(0, 1))) {
    stop("Canonical metadata gender must contain only 0 and 1.")
  }
  data$Sex <- factor(gender, levels = c(0, 1), labels = c("Female", "Male"))
  data$Group <- factor(data$Group, levels = GROUP_LEVELS)

  list(data = data, counts = counts)
}

find_p_column <- function(table, candidates) {
  hit <- candidates[candidates %in% colnames(table)]
  if (!length(hit)) {
    stop(
      "Could not identify permutation p-value column. Available: ",
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

model_formula <- function(outcome) {
  stats::as.formula(
    paste(outcome, "~ Group + Age + Sex + Education")
  )
}

run_omnibus <- function(data, n_perm) {
  rows <- vector("list", length(OUTCOMES))
  for (i in seq_along(OUTCOMES)) {
    outcome <- OUTCOMES[[i]]
    message("Omnibus: ", outcome)
    fit <- permuco::aovperm(
      formula = model_formula(outcome),
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
    rows[[i]] <- data.frame(
      Outcome = outcome,
      N = nrow(data),
      F = if (length(f_column)) as.numeric(table[row, f_column[[1L]]]) else NA_real_,
      P_permutation = as.numeric(table[row, p_column]),
      stringsAsFactors = FALSE
    )
  }
  results <- do.call(rbind, rows)
  results$P_FDR_BH_across_2_outcomes <- stats::p.adjust(
    results$P_permutation,
    method = "BH"
  )
  results$Significant_FDR_0_05 <- results$P_FDR_BH_across_2_outcomes < 0.05
  results
}

adjusted_hedges_g <- function(model, data, levels_two) {
  modal_sex <- names(which.max(table(data$Sex)))[[1L]]
  prediction_data <- data.frame(
    Group = factor(levels_two, levels = levels_two),
    Age = rep(mean(data$Age), 2L),
    Sex = factor(rep(modal_sex, 2L), levels = levels(data$Sex)),
    Education = rep(mean(data$Education), 2L)
  )
  adjusted_means <- stats::predict(model, newdata = prediction_data)
  residual_df <- stats::df.residual(model)
  small_sample_correction <- 1 - 3 / (4 * residual_df - 1)
  small_sample_correction * (
    (adjusted_means[[1L]] - adjusted_means[[2L]]) / stats::sigma(model)
  )
}

run_pairwise <- function(data, n_perm) {
  rows <- list()
  row_index <- 1L
  for (comparison_name in names(PAIRWISE_COMPARISONS)) {
    levels_two <- PAIRWISE_COMPARISONS[[comparison_name]]
    pair_data <- data[data$Group %in% levels_two, , drop = FALSE]
    pair_data$Group <- factor(as.character(pair_data$Group), levels = levels_two)
    pair_data$Sex <- droplevels(pair_data$Sex)

    for (outcome in OUTCOMES) {
      message("Pairwise ", comparison_name, ": ", outcome)
      formula <- model_formula(outcome)
      permutation_fit <- permuco::lmperm(
        formula = formula,
        data = pair_data,
        np = n_perm,
        method = "freedman_lane",
        type = "permutation"
      )
      table <- permutation_fit$table
      row <- find_term_row(table, "Group")
      p_column <- find_p_column(
        table,
        c("resampled Pr(>|t|)", "resampled P(>F)", "P(resampled)")
      )
      t_column <- intersect(c("t", "t value"), colnames(table))
      parametric_fit <- stats::lm(formula, data = pair_data)

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
          adjusted_hedges_g(parametric_fit, pair_data, levels_two)
        ),
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }

  results <- do.call(rbind, rows)
  results$P_FDR_BH_across_5_stage_comparisons_within_outcome <- ave(
    results$P_permutation_two_sided,
    results$Outcome,
    FUN = function(values) stats::p.adjust(values, method = "BH")
  )
  results$Significant_FDR_0_05 <-
    results$P_FDR_BH_across_5_stage_comparisons_within_outcome < 0.05
  results
}

write_results <- function(config, loaded, omnibus, pairwise) {
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  omnibus_file <- file.path(
    config$output_dir,
    paste0(
      "N145_subjectlevel_perturbation_omnibus_FreedmanLane_",
      "age_sex_education.csv"
    )
  )
  pairwise_file <- file.path(
    config$output_dir,
    paste0(
      "N145_subjectlevel_perturbation_pairwise_FreedmanLane_",
      "age_sex_education_FDR_by_outcome.csv"
    )
  )
  utils::write.csv(omnibus, omnibus_file, row.names = FALSE)
  utils::write.csv(pairwise, pairwise_file, row.names = FALSE)

  summary <- data.frame(
    Setting = c(
      "Participant_level_Hopf_input",
      "Metadata_input",
      "Outcomes",
      "Outcome_Pearson_correlation",
      "Model",
      "Permutation_method",
      "Permutations_per_model",
      "Seed",
      "Omnibus_FDR_family",
      "Pairwise_FDR_family",
      "HC_ABneg_N",
      "HC_ABpos_N",
      "MCI_ABpos_N",
      "AD_ABpos_N"
    ),
    Value = c(
      portable_path(config$hopf_file, config$sch1000_root, "$SCH1000_ROOT"),
      portable_path(config$metadata_file, config$sch1000_root, "$SCH1000_ROOT"),
      paste(OUTCOMES, collapse = ", "),
      format(
        stats::cor(
          loaded$data$Info_Cap,
          loaded$data$Susceptibility,
          method = "pearson"
        ),
        digits = 16
      ),
      "Outcome ~ Group + Age + Sex + Education",
      "Freedman-Lane residual permutation; two-sided pairwise tests",
      as.character(config$n_perm),
      as.character(config$seed),
      "BH across the 2 participant-level perturbation outcomes",
      paste0(
        "BH across the 5 planned stage comparisons separately within ",
        "each outcome"
      ),
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
      "N145_subjectlevel_perturbation_group_run_summary.csv"
    ),
    row.names = FALSE
  )
  c(omnibus_file, pairwise_file)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  config <- parse_args(args)
  loaded <- load_analysis_data(config)
  cat(
    "Validated canonical N145 cohort: ",
    paste(names(loaded$counts), loaded$counts, sep = "=", collapse = " | "),
    "\nValidated complete information capacity, susceptibility, age, sex, ",
    "and education.\n",
    sep = ""
  )
  cat(
    "Information capacity--susceptibility Pearson r: ",
    sprintf("%.6f", stats::cor(
      loaded$data$Info_Cap,
      loaded$data$Susceptibility
    )),
    "\n",
    sep = ""
  )
  if (config$validate_only) {
    cat("Validation completed; no models were fitted and no files were written.\n")
    return(invisible(NULL))
  }

  require_package("permuco")
  set.seed(config$seed)
  options(contrasts = c("contr.sum", "contr.poly"))
  omnibus <- run_omnibus(loaded$data, config$n_perm)
  pairwise <- run_pairwise(loaded$data, config$n_perm)
  output_files <- write_results(config, loaded, omnibus, pairwise)

  cat("\nOmnibus results:\n")
  print(omnibus, row.names = FALSE)
  cat("\nPairwise results:\n")
  print(pairwise, row.names = FALSE)
  cat("\nSaved aggregate results:\n", paste(output_files, collapse = "\n"), "\n")
}

if (!isTRUE(getOption("subjectlevel_perturbation.skip_main"))) {
  main()
}
