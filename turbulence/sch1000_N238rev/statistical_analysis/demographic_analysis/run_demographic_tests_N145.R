#!/usr/bin/env Rscript

# Recompute the demographic and clinical table for the canonical N145
# amyloid-status cohort.
#
# Legacy audit: the archived turbu_R and turbu_R_harm scripts modelled
# turbulence outcomes with demographic covariates, but did not contain a
# reproducible demographic-table analysis. This active script therefore uses
# the tests documented in the manuscript table: Kruskal-Wallis tests for
# continuous variables, Pearson's chi-square test for sex, and Holm-corrected
# Dunn tests after significant continuous-variable omnibus tests.

options(stringsAsFactors = FALSE)

EXPECTED_GROUP_COUNTS <- c(
  HC_ABneg = 51L,
  HC_ABpos = 37L,
  MCI_ABpos = 31L,
  AD_ABpos = 26L
)
GROUP_LEVELS <- names(EXPECTED_GROUP_COUNTS)
GROUP_CODES <- c(
  HC_ABneg = "a",
  HC_ABpos = "b",
  MCI_ABpos = "c",
  AD_ABpos = "d"
)

usage <- function(script_name = "run_demographic_tests_N145.R") {
  cat(paste0(
    "Usage: Rscript ", script_name, " [options]\n\n",
    "Options:\n",
    "  --clinical-file PATH  N145 clinical CSV used by the active visualisation pipeline.\n",
    "  --metadata-file PATH  Harmonisation metadata CSV with PTID, age, gender, and edu.\n",
    "  --output-dir PATH     Destination for aggregate CSV and LaTeX outputs.\n",
    "  --validate-only       Validate cohort membership and inputs without testing.\n",
    "  --help                Show this message.\n"
  ))
}

resolve_script_path <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) == 1L) {
    return(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
  }

  project_candidate <- file.path(getwd(), "run_demographic_tests_N145.R")
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
    output_dir = file.path(analysis_dir, "results", "N145_demographics"),
    validate_only = FALSE
  )

  value_options <- c(
    "--clinical-file" = "clinical_file",
    "--metadata-file" = "metadata_file",
    "--output-dir" = "output_dir"
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
  config
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

load_analysis_data <- function(clinical_file, metadata_file) {
  if (!file.exists(clinical_file)) {
    stop("Missing active N145 clinical CSV: ", clinical_file)
  }
  if (!file.exists(metadata_file)) {
    stop("Missing harmonisation metadata CSV: ", metadata_file)
  }

  clinical <- utils::read.csv(
    clinical_file,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN")
  )
  metadata <- utils::read.csv(
    metadata_file,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN")
  )

  clinical_columns <- c(
    "PTID", "Group", "Sex", "Age", "EDUYEARS", "MMMEM", "MMSCORE",
    "MOCA", "TOTAL13", "CL_pvc"
  )
  metadata_columns <- c("PTID", "Group", "age", "gender", "edu")
  require_columns(clinical, clinical_columns, "Clinical CSV")
  require_columns(metadata, metadata_columns, "Metadata CSV")

  clinical$PTID <- normalize_ptid(clinical$PTID, "Clinical CSV")
  metadata$PTID <- normalize_ptid(metadata$PTID, "Metadata CSV")

  observed_counts <- table(factor(clinical$Group, levels = GROUP_LEVELS))
  observed_counts <- stats::setNames(as.integer(observed_counts), GROUP_LEVELS)
  if (!identical(observed_counts, EXPECTED_GROUP_COUNTS)) {
    stop(
      "Clinical CSV is not the canonical N145 cohort. Expected ",
      paste(names(EXPECTED_GROUP_COUNTS), EXPECTED_GROUP_COUNTS,
        sep = "=", collapse = ", "
      ),
      "; found ",
      paste(names(observed_counts), observed_counts, sep = "=", collapse = ", "),
      "."
    )
  }

  metadata_index <- match(clinical$PTID, metadata$PTID)
  if (anyNA(metadata_index)) {
    stop(
      "PTIDs missing from harmonisation metadata: ",
      paste(clinical$PTID[is.na(metadata_index)], collapse = ", ")
    )
  }

  metadata_group <- as.character(metadata$Group[metadata_index])
  group_mismatch <- clinical$PTID[metadata_group != clinical$Group]
  if (length(group_mismatch)) {
    stop(
      "Group labels differ between clinical and metadata inputs for: ",
      paste(group_mismatch, collapse = ", ")
    )
  }

  sex <- toupper(trimws(as.character(clinical$Sex)))
  sex <- substr(sex, 1L, 1L)
  if (!all(sex %in% c("F", "M"))) {
    stop("Clinical Sex must contain only F/Female and M/Male values.")
  }
  metadata_gender <- suppressWarnings(as.numeric(metadata$gender[metadata_index]))
  expected_gender <- ifelse(sex == "F", 0, 1)
  if (anyNA(metadata_gender) || !all(metadata_gender %in% c(0, 1))) {
    stop("Metadata gender must contain only numeric 0 and 1 values.")
  }
  if (!all(metadata_gender == expected_gender)) {
    stop("Sex coding differs between the clinical and metadata inputs.")
  }

  numeric_columns <- c("Age", "EDUYEARS", "MMMEM", "MMSCORE", "MOCA", "TOTAL13", "CL_pvc")
  for (column in numeric_columns) {
    clinical[[column]] <- suppressWarnings(as.numeric(clinical[[column]]))
  }

  canonical_education <- suppressWarnings(as.numeric(metadata$edu[metadata_index]))
  if (anyNA(canonical_education)) {
    stop("Metadata education contains missing or non-numeric values for N145.")
  }
  education_mismatch <- clinical$PTID[
    !is.na(clinical$EDUYEARS) & clinical$EDUYEARS != canonical_education
  ]
  if (length(education_mismatch)) {
    message(
      "Education provenance: using the harmonisation metadata value for ",
      paste(education_mismatch, collapse = ", "),
      " because it differs from the clinical export."
    )
  }

  clinical$Education <- canonical_education
  clinical$Sex <- factor(sex, levels = c("M", "F"))
  clinical$Group <- factor(clinical$Group, levels = GROUP_LEVELS)

  list(
    data = clinical,
    observed_counts = observed_counts,
    education_mismatch = education_mismatch
  )
}

dunn_holm <- function(values, groups) {
  keep <- stats::complete.cases(values, groups)
  values <- as.numeric(values[keep])
  groups <- droplevels(factor(groups[keep], levels = GROUP_LEVELS))
  n_total <- length(values)
  ranks <- rank(values, ties.method = "average")

  tie_counts <- table(values)
  tie_sum <- sum(as.numeric(tie_counts)^3 - as.numeric(tie_counts))
  rank_variance <- n_total * (n_total + 1) / 12 -
    tie_sum / (12 * (n_total - 1))

  pairs <- utils::combn(levels(groups), 2L, simplify = FALSE)
  result <- do.call(rbind, lapply(pairs, function(pair) {
    index_1 <- groups == pair[[1L]]
    index_2 <- groups == pair[[2L]]
    z_value <- (mean(ranks[index_1]) - mean(ranks[index_2])) /
      sqrt(rank_variance * (1 / sum(index_1) + 1 / sum(index_2)))
    data.frame(
      Group_1 = pair[[1L]],
      Group_2 = pair[[2L]],
      N_1 = sum(index_1),
      N_2 = sum(index_2),
      Z = z_value,
      P_Raw = 2 * stats::pnorm(-abs(z_value)),
      stringsAsFactors = FALSE
    )
  }))
  result$P_Holm <- stats::p.adjust(result$P_Raw, method = "holm")
  result$Significant_Holm_0_05 <- result$P_Holm < 0.05
  result
}

continuous_spec <- data.frame(
  Variable = c(
    "Age", "Education", "MMSE Memory", "MMSE Overall", "MOCA",
    "ADAS-Cog 13", "Centiloid"
  ),
  Column = c(
    "Age", "Education", "MMMEM", "MMSCORE", "MOCA", "TOTAL13", "CL_pvc"
  ),
  Latex_Label = c(
    "Age (years)", "Education (years)", "MMSE Memory", "MMSE Overall", "MOCA",
    "ADAS-Cog 13", "Centiloid (A$\\beta$)"
  ),
  stringsAsFactors = FALSE
)

analyze_demographics <- function(data) {
  descriptives <- list()
  omnibus <- list()
  posthoc <- list()

  for (i in seq_len(nrow(continuous_spec))) {
    variable <- continuous_spec$Variable[[i]]
    column <- continuous_spec$Column[[i]]
    values <- data[[column]]

    descriptives[[variable]] <- do.call(rbind, lapply(GROUP_LEVELS, function(group) {
      sample <- values[data$Group == group]
      sample <- sample[is.finite(sample)]
      data.frame(
        Variable = variable,
        Group = group,
        N = length(sample),
        Mean = mean(sample),
        SD = stats::sd(sample),
        Median = stats::median(sample),
        Q1 = unname(stats::quantile(sample, 0.25)),
        Q3 = unname(stats::quantile(sample, 0.75)),
        Minimum = min(sample),
        Maximum = max(sample),
        stringsAsFactors = FALSE
      )
    }))

    test_data <- data.frame(values = values, Group = data$Group)
    test_data <- test_data[stats::complete.cases(test_data), , drop = FALSE]
    kw <- stats::kruskal.test(values ~ Group, data = test_data)
    n_test <- nrow(test_data)
    k <- length(unique(test_data$Group))
    epsilon_squared <- max(0, (unname(kw$statistic) - k + 1) / (n_test - k))
    omnibus[[variable]] <- data.frame(
      Variable = variable,
      Test = "Kruskal-Wallis",
      Statistic = unname(kw$statistic),
      DF = unname(kw$parameter),
      N = n_test,
      P_Value = kw$p.value,
      Effect_Size_Name = "epsilon_squared",
      Effect_Size = epsilon_squared,
      stringsAsFactors = FALSE
    )

    pairwise <- dunn_holm(values, data$Group)
    pairwise$Variable <- variable
    posthoc[[variable]] <- pairwise[, c(
      "Variable", "Group_1", "Group_2", "N_1", "N_2", "Z", "P_Raw",
      "P_Holm", "Significant_Holm_0_05"
    )]
  }

  sex_table <- table(data$Group, data$Sex)
  sex_test <- stats::chisq.test(sex_table, correct = FALSE)
  n_sex <- sum(sex_table)
  cramers_v <- sqrt(
    unname(sex_test$statistic) /
      (n_sex * min(nrow(sex_table) - 1, ncol(sex_table) - 1))
  )
  sex_omnibus <- data.frame(
    Variable = "Sex",
    Test = "Pearson chi-square",
    Statistic = unname(sex_test$statistic),
    DF = unname(sex_test$parameter),
    N = n_sex,
    P_Value = sex_test$p.value,
    Effect_Size_Name = "Cramers_V",
    Effect_Size = cramers_v,
    stringsAsFactors = FALSE
  )
  sex_counts <- as.data.frame.matrix(sex_table)
  sex_counts$Group <- rownames(sex_counts)
  rownames(sex_counts) <- NULL
  sex_counts <- sex_counts[, c("Group", "M", "F")]

  list(
    descriptives = do.call(rbind, descriptives),
    omnibus = rbind(do.call(rbind, omnibus), sex_omnibus),
    posthoc = do.call(rbind, posthoc),
    sex_counts = sex_counts
  )
}

format_p_latex <- function(p_value) {
  if (p_value < 0.001) {
    return("$< 0.001$")
  }
  sprintf("%.3f", p_value)
}

significance_codes <- function(variable, omnibus_p, posthoc) {
  codes <- stats::setNames(rep("", length(GROUP_LEVELS)), GROUP_LEVELS)
  if (!is.finite(omnibus_p) || omnibus_p >= 0.05) {
    return(codes)
  }
  rows <- posthoc[
    posthoc$Variable == variable & posthoc$Significant_Holm_0_05,
    ,
    drop = FALSE
  ]
  if (!nrow(rows)) {
    return(codes)
  }
  for (i in seq_len(nrow(rows))) {
    group_1 <- rows$Group_1[[i]]
    group_2 <- rows$Group_2[[i]]
    codes[[group_1]] <- paste(c(codes[[group_1]], GROUP_CODES[[group_2]]), collapse = ",")
    codes[[group_2]] <- paste(c(codes[[group_2]], GROUP_CODES[[group_1]]), collapse = ",")
  }
  vapply(codes, function(value) {
    pieces <- sort(unique(strsplit(value, ",", fixed = TRUE)[[1L]]))
    paste(pieces[nzchar(pieces)], collapse = ",")
  }, character(1L))
}

format_mean_sd <- function(mean_value, sd_value, code = "") {
  result <- sprintf("%.2f $\\pm$ %.2f", mean_value, sd_value)
  if (nzchar(code)) {
    result <- paste0(result, "\\textsuperscript{", code, "}")
  }
  result
}

make_continuous_row <- function(variable, latex_label, results) {
  test <- results$omnibus[results$omnibus$Variable == variable, , drop = FALSE]
  codes <- significance_codes(variable, test$P_Value, results$posthoc)
  values <- vapply(GROUP_LEVELS, function(group) {
    row <- results$descriptives[
      results$descriptives$Variable == variable &
        results$descriptives$Group == group,
      ,
      drop = FALSE
    ]
    format_mean_sd(row$Mean, row$SD, codes[[group]])
  }, character(1L))
  effect <- if (test$P_Value < 0.05) sprintf("%.3f", test$Effect_Size) else "--"
  c(
    paste0("        ", latex_label, " &"),
    paste0("        ", values[[1L]], " &"),
    paste0("        ", values[[2L]], " &"),
    paste0("        ", values[[3L]], " &"),
    paste0("        ", values[[4L]], " &"),
    "        KW &",
    paste0("        ", format_p_latex(test$P_Value), " &"),
    paste0("        ", effect, " \\\\ ")
  )
}

write_latex_table <- function(results, output_file) {
  sex_test <- results$omnibus[results$omnibus$Variable == "Sex", , drop = FALSE]
  sex_values <- vapply(GROUP_LEVELS, function(group) {
    row <- results$sex_counts[results$sex_counts$Group == group, , drop = FALSE]
    paste0(row$M, " / ", row$F)
  }, character(1L))

  lines <- c(
    "% Generated by run_demographic_tests_N145.R",
    "% Requires: \\usepackage{booktabs,makecell,graphicx}",
    "\\begin{table*}[t]",
    "    \\caption{",
    "        \\textbf{Demographic and clinical characteristics of the ADNI3 cohort} (N = 145 after harmonisation and exclusion of motion outliers).",
    "        Values are mean $\\pm$ SD for continuous variables and counts for categorical variables.",
    "        Group comparisons were performed using Kruskal--Wallis tests for continuous variables and Pearson's chi-square test for sex.",
    "        Significant continuous-variable omnibus tests were followed by Dunn pairwise tests with Holm correction.",
    "        Effect sizes for significant continuous-variable omnibus tests are reported as epsilon squared ($\\varepsilon^2$).",
    "        Education was taken from the same PTID-level metadata used for ComBat and the education-adjusted permutation analyses.",
    "        Group sample sizes refer to the full cohort; tests used all available observations for each variable.",
    "        \\textsuperscript{a}$p_{\\mathrm{Holm}} < 0.05$ vs HC$^-$;",
    "        \\textsuperscript{b}$p_{\\mathrm{Holm}} < 0.05$ vs HC$^+$;",
    "        \\textsuperscript{c}$p_{\\mathrm{Holm}} < 0.05$ vs MCI$^+$;",
    "        \\textsuperscript{d}$p_{\\mathrm{Holm}} < 0.05$ vs AD$^+$.",
    "        \\textit{Abbreviations:} MMSE, Mini-Mental State Examination; MOCA, Montreal Cognitive Assessment; ADAS-Cog 13, Alzheimer's Disease Assessment Scale--Cognitive Subscale (13-item version).",
    "    }",
    "    \\label{tab:adni3-demo-clinical}",
    "    \\centering",
    "    \\footnotesize",
    "    \\setlength{\\tabcolsep}{3pt}",
    "    \\renewcommand{\\arraystretch}{1.15}",
    "    \\resizebox{\\textwidth}{!}{%",
    "    \\begin{tabular}{lccccccc}",
    "        \\toprule",
    "        \\textbf{Variable} &",
    "        \\textbf{HC$^-$} &",
    "        \\textbf{HC$^+$} &",
    "        \\textbf{MCI$^+$} &",
    "        \\textbf{AD$^+$} &",
    "        \\makecell{\\textbf{Test}} &",
    "        \\makecell{\\textbf{$p$-value}} &",
    "        \\makecell{\\textbf{Effect} \\\\ \\textbf{size}} \\\\ ",
    "        & (n = 51) & (n = 37) & (n = 31) & (n = 26) & & & \\\\ ",
    "        \\midrule"
  )

  lines <- c(lines, make_continuous_row("Age", "Age (years)", results), "")
  lines <- c(
    lines,
    "        Sex (M/F) &",
    paste0("        ", sex_values[[1L]], " &"),
    paste0("        ", sex_values[[2L]], " &"),
    paste0("        ", sex_values[[3L]], " &"),
    paste0("        ", sex_values[[4L]], " &"),
    "        $\\chi^2$ &",
    paste0("        ", format_p_latex(sex_test$P_Value), " &"),
    paste0(
      "        ",
      if (sex_test$P_Value < 0.05) sprintf("%.3f", sex_test$Effect_Size) else "--",
      " \\\\ "
    ),
    ""
  )
  lines <- c(lines, make_continuous_row("Education", "Education (years)", results), "")

  remaining <- continuous_spec[!continuous_spec$Variable %in% c("Age", "Education"), ]
  for (i in seq_len(nrow(remaining))) {
    lines <- c(
      lines,
      make_continuous_row(
        remaining$Variable[[i]],
        remaining$Latex_Label[[i]],
        results
      )
    )
    if (i < nrow(remaining)) {
      lines <- c(lines, "")
    }
  }

  lines <- c(
    lines,
    "        \\bottomrule",
    "    \\end{tabular}%",
    "    }",
    "\\end{table*}"
  )
  writeLines(lines, output_file, useBytes = TRUE)
}

write_outputs <- function(config, loaded, results) {
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  utils::write.csv(
    results$descriptives,
    file.path(config$output_dir, "N145_demographic_descriptives.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    results$omnibus,
    file.path(config$output_dir, "N145_demographic_omnibus_tests.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    results$posthoc,
    file.path(config$output_dir, "N145_demographic_posthoc_Dunn_Holm.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    results$sex_counts,
    file.path(config$output_dir, "N145_demographic_sex_counts.csv"),
    row.names = FALSE
  )

  latex_file <- file.path(
    config$output_dir,
    "N145_demographic_clinical_table.tex"
  )
  write_latex_table(results, latex_file)

  summary <- data.frame(
    Item = c(
      "Cohort", "Clinical_input", "Education_input", "Education_mismatches",
      "Continuous_test", "Categorical_test", "Posthoc_test",
      "Posthoc_adjustment", "Continuous_effect_size"
    ),
    Value = c(
      paste(names(loaded$observed_counts), loaded$observed_counts,
        sep = "=", collapse = " | "
      ),
      normalizePath(config$clinical_file, mustWork = TRUE),
      normalizePath(config$metadata_file, mustWork = TRUE),
      if (length(loaded$education_mismatch)) {
        paste(loaded$education_mismatch, collapse = ";")
      } else {
        "None"
      },
      "Kruskal-Wallis",
      "Pearson chi-square",
      "Dunn rank-sum test after significant omnibus test",
      "Holm across six pairwise group contrasts per variable",
      "epsilon squared = (H - k + 1) / (N - k)"
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    summary,
    file.path(config$output_dir, "N145_demographic_run_summary.csv"),
    row.names = FALSE
  )

  cat("\nSaved demographic results to: ", config$output_dir, "\n", sep = "")
  cat("LaTeX table: ", latex_file, "\n", sep = "")
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  config <- parse_args(args)
  loaded <- load_analysis_data(config$clinical_file, config$metadata_file)
  cat(
    "Validated canonical N145 cohort: ",
    paste(names(loaded$observed_counts), loaded$observed_counts,
      sep = "=", collapse = " | "
    ),
    "\n",
    sep = ""
  )
  cat("Education source: harmonisation metadata column 'edu'.\n")

  if (config$validate_only) {
    cat("Validation completed; no statistical outputs were written.\n")
    return(invisible(NULL))
  }

  results <- analyze_demographics(loaded$data)
  write_outputs(config, loaded, results)

  cat("\nOmnibus results:\n")
  print(results$omnibus[, c(
    "Variable", "Test", "Statistic", "DF", "N", "P_Value", "Effect_Size"
  )], row.names = FALSE)
  invisible(results)
}

if (!isTRUE(getOption("demographic_analysis.skip_main"))) {
  main()
}
