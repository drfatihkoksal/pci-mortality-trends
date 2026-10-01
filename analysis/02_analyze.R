#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages(library(survival))
source("analysis/config.R")

for (path in unname(study_config[c("results_dir", "tables_dir", "figures_dir", "audit_dir")])) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(study_config$input_csv)) {
  stop("Protected analysis cohort not found. Run 01_extract_deidentify.py first.")
}

d <- read.csv(study_config$input_csv, na.strings = c("", "NA"), check.names = FALSE)
d$index_date <- as.Date(d$index_date)
d$death_date <- as.Date(d$death_date)
d$last_vital_status_date <- as.Date(d$last_vital_status_date)
d$index_year <- as.integer(d$index_year)
d$age <- as.numeric(d$age)
d$sex <- factor(d$sex, levels = c("Male", "Female"))
d$clinical_presentation <- factor(d$clinical_presentation)
d$source <- factor(d$source)
d$time_days <- ifelse(d$dead == 1 & d$followup_days == 0, 0.5, d$followup_days)
d$valid_time <- !is.na(d$time_days) & d$time_days >= 0 & d$death_before_index == 0

binary_columns <- intersect(
  c(
    "primary_pci", "left_main", "graft_pci", "cto", "restenosis", "des", "bms",
    "diabetes", "hypertension", "current_smoker", "heart_failure", "prior_mi",
    "prior_pci", "prior_cabg", "prior_stroke", "peripheral_artery_disease",
    "dialysis", "cardiogenic_shock", "cardiac_arrest", "radial_access",
    "multivessel_disease", "procedural_success", "intravascular_imaging"
  ),
  names(d)
)
for (column in binary_columns) d[[column]] <- as.numeric(d[[column]])
for (column in intersect(c("egfr", "hemoglobin", "lvef"), names(d))) {
  d[[column]] <- as.numeric(d[[column]])
}
if ("ehr_presentation" %in% names(d)) d$ehr_presentation <- factor(d$ehr_presentation)

ehr_signal_columns <- c("diabetes", "egfr", "lvef", "cardiogenic_shock", "cardiac_arrest")
ehr_enriched <- any(ehr_signal_columns %in% names(d))
analysis_mode <- if (ehr_enriched) "EHR-enriched" else "limited-adjustment preliminary"
selected_model_name <- if (ehr_enriched) "available_case_mix" else "age_sex_adjusted"
selected_model_label <- if (ehr_enriched) "Case-mix adjusted" else "Age/sex adjusted"

percent <- function(x, digits = 1) {
  if (length(x) == 0 || all(is.na(x))) return("NA")
  sprintf(paste0("%.", digits, "f"), mean(x, na.rm = TRUE) * 100)
}

median_iqr <- function(x, digits = 0) {
  if (length(x) == 0 || all(is.na(x))) return("NA")
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf(paste0("%.", digits, "f (%.", digits, "f–%.", digits, "f)"), q[2], q[1], q[3])
}

count_percent <- function(x, condition = function(z) z == 1) {
  keep <- !is.na(x)
  if (!any(keep)) return("NA")
  n <- sum(condition(x[keep]))
  sprintf("%s (%.1f)", format(n, big.mark = ","), 100 * n / sum(keep))
}

write_csv_safe <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "")
}

eligible_horizon <- function(data, start_year, end_year, horizon) {
  x <- data[
    data$valid_time & data$index_year >= start_year & data$index_year <= end_year,
    , drop = FALSE
  ]
  x$event <- as.integer(x$dead == 1 & x$time_days <= horizon)
  x$mature <- x$time_days >= horizon | x$event == 1
  x <- x[x$mature, , drop = FALSE]
  rownames(x) <- seq_len(nrow(x))
  x
}

usable_numeric <- function(data, column, threshold = 0.60) {
  column %in% names(data) && mean(!is.na(data[[column]])) >= threshold &&
    length(unique(data[[column]][!is.na(data[[column]])])) > 1
}

available_model_terms <- function(data, enriched = TRUE) {
  terms <- c("splines::ns(age, df = 3)", "sex")
  if (!enriched) return(terms)

  presentation <- if ("ehr_presentation" %in% names(data) &&
    mean(!is.na(data$ehr_presentation)) >= 0.60) "ehr_presentation" else "clinical_presentation"
  if (length(unique(data[[presentation]][!is.na(data[[presentation]])])) > 1) {
    terms <- c(terms, presentation)
  }
  for (column in c(
    "left_main", "graft_pci", "diabetes", "hypertension", "current_smoker",
    "heart_failure", "prior_mi", "prior_pci", "prior_cabg", "prior_stroke",
    "peripheral_artery_disease", "dialysis", "cardiogenic_shock", "cardiac_arrest"
  )) {
    if (usable_numeric(data, column)) terms <- c(terms, column)
  }
  for (column in c("egfr", "hemoglobin", "lvef")) {
    if (usable_numeric(data, column)) {
      terms <- c(terms, sprintf("splines::ns(%s, df = 3)", column))
    }
  }
  unique(terms)
}

formula_variables <- function(formula) unique(all.vars(formula))

year_spline <- function(df = 3) sprintf("splines::ns(index_year, df = %d)", df)
year_levels <- as.character(sort(unique(d$index_year)))

# Patient-clustered sandwich covariance for a logistic model; used when one
# patient can contribute several procedure episodes.
cluster_vcov <- function(model, cluster) {
  X <- model.matrix(model)
  keep <- !is.na(coef(model))
  X <- X[, keep, drop = FALSE]
  scores <- X * as.vector(model$y - fitted(model))
  scores_by_cluster <- rowsum(scores, cluster)
  bread <- vcov(model)[keep, keep, drop = FALSE]
  g <- nrow(scores_by_cluster)
  bread %*% crossprod(scores_by_cluster) %*% bread * g / (g - 1)
}

fit_logistic_model <- function(data, covariate_terms, year_term = year_spline(3), cluster = NULL) {
  rhs <- paste(c(year_term, covariate_terms), collapse = " + ")
  formula <- as.formula(paste("event ~", rhs))
  needed <- formula_variables(formula)
  if (!is.null(cluster)) needed <- unique(c(needed, cluster))
  data$year_f <- factor(data$index_year, levels = year_levels)
  model_data <- data[complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]
  model_data$year_f <- droplevels(model_data$year_f)
  model <- glm(formula, data = model_data, family = binomial(), x = TRUE, model = TRUE)
  keep <- !is.na(coef(model))

  if (is.null(cluster)) {
    reduced_rhs <- if (length(covariate_terms) == 0) "1" else paste(covariate_terms, collapse = " + ")
    reduced_formula <- as.formula(paste("event ~", reduced_rhs))
    reduced <- glm(reduced_formula, data = model_data, family = binomial())
    lrt <- anova(reduced, model, test = "Chisq")
    p_year <- lrt$`Pr(>Chi)`[2]
    V <- vcov(model)[keep, keep, drop = FALSE]
  } else {
    # Robust Wald test of the calendar-year terms.
    V <- cluster_vcov(model, model_data[[cluster]])
    beta <- coef(model)[keep]
    year_index <- grep("index_year|year_f", names(beta))
    b <- beta[year_index]
    stat <- as.numeric(t(b) %*% solve(V[year_index, year_index, drop = FALSE]) %*% b)
    p_year <- pchisq(stat, df = length(year_index), lower.tail = FALSE)
  }
  list(
    model = model, data = model_data, p_year = p_year, V = V,
    formula = paste(deparse(formula), collapse = " "),
    n_clusters = if (is.null(cluster)) nrow(model_data) else length(unique(model_data[[cluster]]))
  )
}

standardized_risk <- function(fit, year, fixed = list()) {
  newdata <- fit$data
  newdata$index_year <- year
  newdata$year_f <- factor(year, levels = levels(fit$data$year_f))
  for (column in names(fixed)) {
    newdata[[column]] <- factor(fixed[[column]], levels = levels(fit$data[[column]]))
  }
  tt <- delete.response(terms(fit$model))
  X <- model.matrix(tt, newdata, contrasts.arg = fit$model$contrasts)
  beta <- coef(fit$model)
  keep <- !is.na(beta)
  X <- X[, keep, drop = FALSE]
  beta <- beta[keep]
  V <- fit$V
  eta <- as.vector(X %*% beta)
  pr <- plogis(eta)
  risk <- mean(pr)
  gradient <- colMeans(X * as.vector(pr * (1 - pr)))
  variance <- as.numeric(t(gradient) %*% V %*% gradient)
  se <- sqrt(max(variance, 0))
  se_logit <- se / max(risk * (1 - risk), 1e-9)
  ci <- plogis(qlogis(risk) + c(-1, 1) * qnorm(0.975) * se_logit)
  list(risk = risk, lower = ci[1], upper = ci[2], gradient = gradient, V = V)
}

risk_difference <- function(fit, start_year, end_year, fixed = list()) {
  a <- standardized_risk(fit, start_year, fixed)
  b <- standardized_risk(fit, end_year, fixed)
  gradient <- b$gradient - a$gradient
  variance <- as.numeric(t(gradient) %*% a$V %*% gradient)
  estimate <- b$risk - a$risk
  se <- sqrt(max(variance, 0))
  c(
    estimate = estimate,
    lower = estimate - qnorm(0.975) * se,
    upper = estimate + qnorm(0.975) * se
  )
}

tidy_glm <- function(model, model_name, outcome_name) {
  s <- summary(model)$coefficients
  data.frame(
    outcome = outcome_name,
    model = model_name,
    term = rownames(s),
    estimate = s[, 1],
    std_error = s[, 2],
    statistic = s[, 3],
    p_value = s[, 4],
    row.names = NULL,
    check.names = FALSE
  )
}

fit_horizon <- function(data, start_year, end_year, horizon, outcome_name) {
  x <- eligible_horizon(data, start_year, end_year, horizon)
  crude <- do.call(rbind, lapply(split(x, x$index_year), function(g) {
    data.frame(
      outcome = outcome_name,
      index_year = unique(g$index_year),
      n = nrow(g),
      deaths = sum(g$event),
      crude_risk = mean(g$event),
      stringsAsFactors = FALSE
    )
  }))
  rownames(crude) <- NULL

  models <- list(
    unadjusted = fit_logistic_model(x, character()),
    age_sex_adjusted = fit_logistic_model(x, available_model_terms(x, enriched = FALSE)),
    available_case_mix = fit_logistic_model(x, available_model_terms(x, enriched = TRUE))
  )

  annual <- do.call(rbind, lapply(names(models), function(model_name) {
    fit <- models[[model_name]]
    do.call(rbind, lapply(start_year:end_year, function(year) {
      z <- standardized_risk(fit, year)
      data.frame(
        outcome = outcome_name,
        model = model_name,
        index_year = year,
        risk = z$risk,
        lower = z$lower,
        upper = z$upper,
        p_year = fit$p_year,
        model_n = nrow(fit$data),
        stringsAsFactors = FALSE
      )
    }))
  }))

  contrasts <- do.call(rbind, lapply(names(models), function(model_name) {
    fit <- models[[model_name]]
    rd <- risk_difference(fit, start_year, end_year)
    data.frame(
      outcome = outcome_name,
      model = model_name,
      start_year = start_year,
      end_year = end_year,
      start_risk = standardized_risk(fit, start_year)$risk,
      end_risk = standardized_risk(fit, end_year)$risk,
      risk_difference = rd["estimate"],
      rd_lower = rd["lower"],
      rd_upper = rd["upper"],
      p_year = fit$p_year,
      model_n = nrow(fit$data),
      formula = paste(fit$formula, collapse = " "),
      stringsAsFactors = FALSE
    )
  }))

  coefficients <- do.call(rbind, lapply(names(models), function(model_name) {
    tidy_glm(models[[model_name]]$model, model_name, outcome_name)
  }))

  # Cox model is a secondary time-to-event estimand and diagnostic; it uses the
  # same covariates as the selected logistic model.
  cox_terms <- available_model_terms(x, enriched = ehr_enriched)
  cox_formula <- as.formula(paste(
    "Surv(pmin(time_days,", horizon, "), event) ~",
    paste(c(year_spline(3), cox_terms), collapse = " + ")
  ))
  cox_needed <- setdiff(formula_variables(cox_formula), c("Surv", "pmin"))
  cox_data <- x[complete.cases(x[, intersect(cox_needed, names(x)), drop = FALSE]), , drop = FALSE]
  cox <- coxph(cox_formula, data = cox_data, x = TRUE, model = TRUE)
  cox_reduced <- coxph(
    as.formula(paste("Surv(pmin(time_days,", horizon, "), event) ~", paste(cox_terms, collapse = " + "))),
    data = cox_data
  )
  cox_p_year <- anova(cox_reduced, cox)$`Pr(>|Chi|)`[2]
  ph <- tryCatch(cox.zph(cox), error = function(e) NULL)
  ph_year_p <- if (is.null(ph)) NA_real_ else ph$table[grep("index_year", rownames(ph$table)), "p"]
  ph_global_p <- if (is.null(ph)) NA_real_ else ph$table["GLOBAL", "p"]
  cox_summary <- data.frame(
    outcome = outcome_name,
    n = cox$n,
    events = cox$nevent,
    p_year_lrt = cox_p_year,
    ph_p_year = ph_year_p,
    ph_p_global = ph_global_p
  )
  ph_table <- if (is.null(ph)) data.frame() else {
    z <- as.data.frame(ph$table)
    z$term <- rownames(z)
    z$outcome <- outcome_name
    rownames(z) <- NULL
    z
  }

  list(
    data = x,
    crude = crude,
    annual = annual,
    contrasts = contrasts,
    coefficients = coefficients,
    cox = cox,
    cox_summary = cox_summary,
    ph = ph_table,
    models = models
  )
}

primary <- fit_horizon(
  d,
  study_config$primary_start_year,
  study_config$primary_end_year,
  study_config$primary_horizon_days,
  "1-year all-cause mortality"
)
early <- fit_horizon(
  d,
  study_config$primary_start_year,
  study_config$early_end_year,
  study_config$early_horizon_days,
  "30-day all-cause mortality"
)
long <- fit_horizon(
  d,
  study_config$primary_start_year,
  study_config$long_term_end_year,
  study_config$long_term_horizon_days,
  "5-year all-cause mortality"
)

all_crude <- rbind(primary$crude, early$crude, long$crude)
all_annual <- rbind(primary$annual, early$annual, long$annual)
all_contrasts <- rbind(primary$contrasts, early$contrasts, long$contrasts)
all_coefficients <- rbind(primary$coefficients, early$coefficients, long$coefficients)
all_ph <- rbind(primary$ph, early$ph, long$ph)
all_cox <- rbind(primary$cox_summary, early$cox_summary, long$cox_summary)
write_csv_safe(all_cox, file.path(study_config$results_dir, "cox_summary.csv"))

write_csv_safe(all_crude, file.path(study_config$tables_dir, "annual_crude_mortality.csv"))
write_csv_safe(all_annual, file.path(study_config$tables_dir, "annual_adjusted_mortality.csv"))
write_csv_safe(all_contrasts, file.path(study_config$tables_dir, "model_contrasts.csv"))
write_csv_safe(all_coefficients, file.path(study_config$tables_dir, "model_coefficients.csv"))
write_csv_safe(all_ph, file.path(study_config$audit_dir, "proportional_hazards_tests.csv"))

# ---------------------------------------------------------------------------
# Sensitivity analyses of the calendar-year pattern.
# ---------------------------------------------------------------------------
selected_terms <- available_model_terms(primary$data, enriched = ehr_enriched)
primary_label <- sprintf("Primary (%s, spline df 3)", tolower(selected_model_label))

summarize_fit <- function(fit, outcome_name, analysis, start_year, end_year, fixed = list(), events = NULL) {
  years <- start_year:end_year
  risks <- sapply(years, function(y) standardized_risk(fit, y, fixed)$risk)
  peak_year <- years[which.max(risks)]
  rd_es <- risk_difference(fit, start_year, end_year, fixed)
  rd_ps <- risk_difference(fit, start_year, peak_year, fixed)
  rd_ep <- risk_difference(fit, peak_year, end_year, fixed)
  data.frame(
    outcome = outcome_name,
    analysis = analysis,
    start_year = start_year,
    end_year = end_year,
    model_n = nrow(fit$data),
    patients = fit$n_clusters,
    events = sum(fit$data$event),
    start_risk = risks[1],
    peak_year = peak_year,
    peak_risk = max(risks),
    end_risk = risks[length(risks)],
    rd_end_start = rd_es["estimate"], rd_end_start_lower = rd_es["lower"], rd_end_start_upper = rd_es["upper"],
    rd_peak_start = rd_ps["estimate"], rd_peak_start_lower = rd_ps["lower"], rd_peak_start_upper = rd_ps["upper"],
    rd_end_peak = rd_ep["estimate"], rd_end_peak_lower = rd_ep["lower"], rd_end_peak_upper = rd_ep["upper"],
    p_year = fit$p_year,
    aic = AIC(fit$model),
    stringsAsFactors = FALSE
  )
}

annual_from_fit <- function(fit, outcome_name, analysis, start_year, end_year, fixed = list()) {
  do.call(rbind, lapply(start_year:end_year, function(y) {
    z <- standardized_risk(fit, y, fixed)
    data.frame(outcome = outcome_name, analysis = analysis, index_year = y,
      risk = z$risk, lower = z$lower, upper = z$upper, stringsAsFactors = FALSE)
  }))
}

sens_rows <- list()
sens_annual <- list()
add_sensitivity <- function(fit, outcome_name, analysis, start_year, end_year, fixed = list()) {
  sens_rows[[length(sens_rows) + 1]] <<- summarize_fit(fit, outcome_name, analysis, start_year, end_year, fixed)
  sens_annual[[length(sens_annual) + 1]] <<- annual_from_fit(fit, outcome_name, analysis, start_year, end_year, fixed)
}

horizon_specs <- list(
  list(result = primary, outcome = "1-year all-cause mortality", end = study_config$primary_end_year,
    horizon = study_config$primary_horizon_days),
  list(result = early, outcome = "30-day all-cause mortality", end = study_config$early_end_year,
    horizon = study_config$early_horizon_days)
)

for (spec in horizon_specs) {
  x <- spec$result$data
  start <- study_config$primary_start_year
  # (a) Selected model and the alternative adjustment tiers.
  add_sensitivity(spec$result$models[[selected_model_name]], spec$outcome, primary_label, start, spec$end)
  add_sensitivity(spec$result$models$unadjusted, spec$outcome, "Unadjusted", start, spec$end)
  if (!ehr_enriched) {
    add_sensitivity(spec$result$models$available_case_mix, spec$outcome,
      "Age/sex + text-derived presentation, left main, graft PCI", start, spec$end)
  }
  # (b) Calendar-year specification.
  for (df in c(2, 4, 5)) {
    add_sensitivity(fit_logistic_model(x, selected_terms, year_spline(df)), spec$outcome,
      sprintf("Spline df %d", df), start, spec$end)
  }
  add_sensitivity(fit_logistic_model(x, selected_terms, "year_f"), spec$outcome,
    "Year as categorical", start, spec$end)
  # (c) Excluding the first, incompletely captured source year.
  add_sensitivity(fit_logistic_model(x[x$index_year >= start + 1, ], selected_terms), spec$outcome,
    sprintf("Index years %d–%d", start + 1, spec$end), start + 1, spec$end)
  # (d) Report-source indicator; risks standardized as if all reports were PDF.
  source_data <- x[x$source %in% c("PDF", "PARS"), , drop = FALSE]
  source_data$source <- factor(as.character(source_data$source), levels = c("PDF", "PARS"))
  add_sensitivity(fit_logistic_model(source_data, c(selected_terms, "source")), spec$outcome,
    "Adjusted for report source (PDF vs PARS)", start, spec$end, fixed = list(source = "PDF"))
  # (e) Fixed 2-year lookback: every procedure episode with no observed PCI
  # during the preceding 730 days, from 2017 so all episodes share the same
  # observable lookback; patient-clustered robust variance.
  days_path <- file.path(dirname(study_config$input_csv), "procedure_days.csv")
  pdays <- read.csv(days_path, na.strings = c("", "NA"))
  pdays$sex <- factor(pdays$sex, levels = c("Male", "Female"))
  pdays$time_days <- ifelse(pdays$dead == 1 & pdays$followup_days == 0, 0.5, pdays$followup_days)
  pdays$valid_time <- !is.na(pdays$time_days) & pdays$time_days >= 0 & pdays$death_before_index == 0
  pdays <- pdays[is.na(pdays$days_since_prior_procedure) | pdays$days_since_prior_procedure > 730, ]
  lookback_start <- start + 2
  episodes <- eligible_horizon(pdays, lookback_start, spec$end, spec$horizon)
  add_sensitivity(fit_logistic_model(episodes, selected_terms, cluster = "study_id"), spec$outcome,
    sprintf("Fixed 2-year lookback episodes, %d–%d", lookback_start, spec$end), lookback_start, spec$end)
}

sensitivity <- do.call(rbind, sens_rows)
rownames(sensitivity) <- NULL
sensitivity_annual <- do.call(rbind, sens_annual)
write_csv_safe(sensitivity, file.path(study_config$tables_dir, "sensitivity_analyses.csv"))
write_csv_safe(sensitivity_annual, file.path(study_config$tables_dir, "sensitivity_annual.csv"))

# PDF-to-PARS bridge year (2022): both systems were in use, sequentially.
bridge <- do.call(rbind, lapply(horizon_specs, function(spec) {
  x <- spec$result$data
  x <- x[x$index_year == 2022 & x$source %in% c("PDF", "PARS"), , drop = FALSE]
  x$source <- factor(as.character(x$source), levels = c("PDF", "PARS"))
  fit <- glm(as.formula(paste("event ~ source +", paste(selected_terms, collapse = " + "))),
    data = x, family = binomial())
  est <- coef(summary(fit))["sourcePARS", ]
  do.call(rbind, lapply(levels(x$source), function(src) {
    g <- x[x$source == src, ]
    data.frame(
      outcome = spec$outcome,
      source = src,
      first_month = format(min(as.Date(g$index_date)), "%Y-%m"),
      last_month = format(max(as.Date(g$index_date)), "%Y-%m"),
      n = nrow(g),
      median_age = median(g$age, na.rm = TRUE),
      female_percent = 100 * mean(g$sex == "Female", na.rm = TRUE),
      stemi_primary_percent = 100 * mean(g$clinical_presentation == "STEMI/primary PCI", na.rm = TRUE),
      deaths = sum(g$event),
      crude_percent = 100 * mean(g$event),
      adjusted_or_pars = exp(est[1]),
      or_lower = exp(est[1] - qnorm(0.975) * est[2]),
      or_upper = exp(est[1] + qnorm(0.975) * est[2]),
      or_p = est[4],
      stringsAsFactors = FALSE
    )
  }))
}))
write_csv_safe(bridge, file.path(study_config$tables_dir, "source_bridge_2022.csv"))

# Calendar-time (period) mortality among earlier cohorts. Patients treated in
# 2015–2018 who survived 1 year are followed through calendar years 2016–2025;
# a rise in their death rate during specific calendar years indicates
# background (period) mortality independent of the index procedure.
period_base <- d[d$valid_time & d$index_year >= 2015 & d$index_year <= 2018 & !is.na(d$age) & !is.na(d$sex), ]
period_base$fu_start <- period_base$index_date + 365
period_base$fu_end <- pmin(period_base$last_vital_status_date, as.Date("2026-01-01"))
period_base <- period_base[period_base$fu_end > period_base$fu_start, ]
period_rows <- do.call(rbind, lapply(2016:2025, function(y) {
  y0 <- as.Date(sprintf("%d-01-01", y)); y1 <- as.Date(sprintf("%d-01-01", y + 1))
  a <- pmax(period_base$fu_start, y0); b <- pmin(period_base$fu_end, y1)
  keep <- b > a
  g <- period_base[keep, ]
  a <- a[keep]; b <- b[keep]
  died <- g$dead == 1 & !is.na(g$death_date) & g$death_date > a & g$death_date <= b
  mid <- a + (b - a) / 2
  data.frame(
    calendar_year = y,
    person_years = as.numeric(b - a) / 365.25,
    death = as.integer(died),
    attained_age = g$age + as.numeric(mid - g$index_date) / 365.25,
    years_since_pci = as.numeric(mid - g$index_date) / 365.25,
    sex = g$sex
  )
}))
period_rows$calendar_f <- relevel(factor(period_rows$calendar_year), ref = "2019")
# Time since PCI is not included: in a closed 2015–2018 cohort it is nearly
# collinear with calendar year (age-period-cohort identification problem).
period_fit <- glm(death ~ calendar_f + splines::ns(attained_age, df = 3) + sex +
  offset(log(person_years)), family = poisson(), data = period_rows)
period_coef <- coef(summary(period_fit))
period_table <- do.call(rbind, lapply(2016:2025, function(y) {
  g <- period_rows[period_rows$calendar_year == y, ]
  term <- paste0("calendar_f", y)
  est <- if (term %in% rownames(period_coef)) period_coef[term, 1:2] else c(0, 0)
  data.frame(
    calendar_year = y,
    person_years = sum(g$person_years),
    deaths = sum(g$death),
    crude_rate_per_100py = 100 * sum(g$death) / sum(g$person_years),
    adjusted_rate_ratio_vs_2019 = exp(est[1]),
    rr_lower = if (y == 2019) NA else exp(est[1] - qnorm(0.975) * est[2]),
    rr_upper = if (y == 2019) NA else exp(est[1] + qnorm(0.975) * est[2])
  )
}))
period_lrt <- anova(update(period_fit, . ~ . - calendar_f), period_fit, test = "Chisq")$`Pr(>Chi)`[2]
period_table$p_calendar_year <- period_lrt
write_csv_safe(period_table, file.path(study_config$tables_dir, "period_mortality_2015_2018_cohorts.csv"))

# Figure S1: calendar-year pattern across sensitivity analyses.
png(file.path(study_config$figures_dir, "figureS1_sensitivity_trends.png"), width = 2600, height = 1300, res = 220)
par(mfrow = c(1, 2), mar = c(4.2, 4.5, 2.5, 1))
sens_show <- c(primary_label, "Year as categorical",
  "Adjusted for report source (PDF vs PARS)")
sens_cols <- c("#1F4E79", "#7F7F7F", "#B9770E", "#117864")
for (outcome_name in c("30-day all-cause mortality", "1-year all-cause mortality")) {
  z <- sensitivity_annual[sensitivity_annual$outcome == outcome_name, ]
  lookback_label <- grep("^Fixed 2-year lookback", unique(z$analysis), value = TRUE)
  labels <- c(sens_show, lookback_label)
  z <- z[z$analysis %in% labels, ]
  ylim <- c(0, 1.18 * max(z$risk) * 100)
  plot(NA, xlim = range(z$index_year), ylim = ylim, xlab = "Index year", ylab = "Standardized mortality, %",
    main = outcome_name, las = 1)
  for (i in seq_along(labels)) {
    zi <- z[z$analysis == labels[i], ]
    lines(zi$index_year, 100 * zi$risk, type = "b", pch = c(19, 1, 17, 15)[i], col = sens_cols[i],
      lwd = c(2.5, 1.5, 1.5, 1.5)[i])
  }
  legend("topleft", c("Primary (spline df 3)", "Year categorical", "Source-adjusted", "2-year lookback episodes"),
    col = sens_cols, pch = c(19, 1, 17, 15), lty = 1, bty = "n", cex = 0.8)
}
dev.off()

# Cohort maturity by year.
maturity <- do.call(rbind, lapply(sort(unique(d$index_year[d$index_year >= 2015 & d$index_year <= 2026])), function(year) {
  g <- d[d$valid_time & d$index_year == year, , drop = FALSE]
  result <- data.frame(index_year = year, n = nrow(g))
  for (horizon in c(30, 365, 1825)) {
    event <- g$dead == 1 & g$time_days <= horizon
    mature <- g$time_days >= horizon | event
    result[[paste0("mature_", horizon, "d")]] <- sum(mature, na.rm = TRUE)
    result[[paste0("deaths_", horizon, "d")]] <- sum(event, na.rm = TRUE)
  }
  result
}))
write_csv_safe(maturity, file.path(study_config$audit_dir, "followup_maturity_by_year.csv"))

# Table 1 in equal two-year periods.
table1_data <- primary$data
table1_data$period <- cut(
  table1_data$index_year,
  breaks = c(2014, 2016, 2018, 2020, 2022, 2024),
  labels = c("2015–16", "2017–18", "2019–20", "2021–22", "2023–24")
)
summarize_period <- function(g, label) {
  data.frame(
    period = label,
    n = nrow(g),
    age = median_iqr(g$age),
    female = count_percent(g$sex, function(z) z == "Female"),
    primary_pci = count_percent(g$primary_pci),
    stemi_primary = count_percent(g$clinical_presentation, function(z) z == "STEMI/primary PCI"),
    left_main = count_percent(g$left_main),
    graft_pci = count_percent(g$graft_pci),
    cto = count_percent(g$cto),
    restenosis = count_percent(g$restenosis),
    des = count_percent(g$des),
    bms = count_percent(g$bms),
    one_year_death = count_percent(g$event),
    stringsAsFactors = FALSE
  )
}
table1 <- do.call(rbind, lapply(levels(table1_data$period), function(level) {
  summarize_period(table1_data[table1_data$period == level, , drop = FALSE], level)
}))
table1 <- rbind(summarize_period(table1_data, "Overall"), table1)
write_csv_safe(table1, file.path(study_config$tables_dir, "table1_baseline_by_period.csv"))

# Annual case-mix trends.
case_mix <- do.call(rbind, lapply(split(primary$data, primary$data$index_year), function(g) {
  data.frame(
    index_year = unique(g$index_year),
    n = nrow(g),
    median_age = median(g$age, na.rm = TRUE),
    female_percent = 100 * mean(g$sex == "Female", na.rm = TRUE),
    primary_pci_percent = 100 * mean(g$primary_pci == 1, na.rm = TRUE),
    left_main_percent = 100 * mean(g$left_main == 1, na.rm = TRUE),
    graft_pci_percent = 100 * mean(g$graft_pci == 1, na.rm = TRUE),
    des_percent = 100 * mean(g$des == 1, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
rownames(case_mix) <- NULL
write_csv_safe(case_mix, file.path(study_config$tables_dir, "annual_case_mix.csv"))

# Figure 1: cohort flow.
in_primary_period <- d$index_year >= study_config$primary_start_year &
  d$index_year <= study_config$primary_end_year
outside_primary_period <- sum(!in_primary_period, na.rm = TRUE)
before_primary <- sum(d$index_year < study_config$primary_start_year, na.rm = TRUE)
after_primary_early <- sum(d$index_year == study_config$early_end_year, na.rm = TRUE)
after_primary_partial <- sum(d$index_year > study_config$early_end_year, na.rm = TRUE)
death_before_index_primary <- sum(
  in_primary_period & d$death_before_index == 1,
  na.rm = TRUE
)
audit_json <- function(path, key) {
  line <- grep(sprintf('"%s"', key), readLines(path, warn = FALSE), value = TRUE)[1]
  as.numeric(gsub("[^0-9]", "", sub(".*:", "", line)))
}
procedure_days_n <- audit_json(file.path(study_config$audit_dir, "data_quality_summary.json"), "procedure_days")
source_rows <- sum(d$source_record_count, na.rm = TRUE)
selected_model_n <- unique(primary$annual$model_n[primary$annual$model == selected_model_name])
fmt <- function(x) format(x, big.mark = ",")

png(file.path(study_config$figures_dir, "figure1_cohort_flow.png"), width = 2400, height = 2200, res = 220)
par(mar = c(0, 0, 0, 0))
plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
flow_box <- function(x, y, text, w = 0.40, h = 0.075, fill = "#EAF2F8", cex = 0.82) {
  rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2, col = fill, border = "#1F4E79", lwd = 1.6)
  text(x, y, text, cex = cex)
}
down <- function(x, y1, y2) arrows(x, y1, x, y2, length = 0.07, lwd = 1.6, col = "#1F4E79")
side <- function(y, x1, x2) arrows(x1, y, x2, y, length = 0.07, lwd = 1.6, col = "#1F4E79")
mx <- 0.30; sx <- 0.76
flow_box(mx, 0.93, sprintf("PCI reports in source workbook\nn = %s reports", fmt(source_rows)))
down(mx, 0.892, 0.818)
flow_box(mx, 0.78, sprintf("Procedure days after merging\nsame-day reports\nn = %s", fmt(procedure_days_n)))
side(0.855, mx, sx - 0.20)
flow_box(sx, 0.855, sprintf("Same-day reports merged\nn = %s reports", fmt(source_rows - procedure_days_n)), fill = "#FDEDEC")
down(mx, 0.742, 0.668)
flow_box(mx, 0.63, sprintf("First observed PCI per patient\nn = %s patients", fmt(nrow(d))))
side(0.705, mx, sx - 0.20)
flow_box(sx, 0.705, sprintf("Repeat procedure days\nnot used as index\nn = %s", fmt(procedure_days_n - nrow(d))), fill = "#FDEDEC")
down(mx, 0.592, 0.518)
flow_box(mx, 0.48, sprintf("Primary 1-year cohort, 2015–2024\nn = %s patients; %s deaths", fmt(nrow(primary$data)), fmt(sum(primary$data$event))), fill = "#E8F8F5")
side(0.555, mx, sx - 0.20)
flow_box(sx, 0.555, sprintf(
  "Excluded from primary cohort\nIndex before 2015: %s\nIndex 2025: %s (30-day cohort only)\nIndex 2026, partial year: %s\nDeath recorded before index: %s",
  fmt(before_primary), fmt(after_primary_early), fmt(after_primary_partial), fmt(death_before_index_primary)
), h = 0.13, fill = "#FDEDEC")
down(mx, 0.442, 0.368)
flow_box(mx, 0.33, sprintf("Age/sex-adjusted model\nn = %s", fmt(selected_model_n)), fill = "#E8F8F5")
side(0.405, mx, sx - 0.20)
flow_box(sx, 0.405, sprintf("Missing or implausible age\nor missing sex\nn = %s", fmt(nrow(primary$data) - selected_model_n)), fill = "#FDEDEC")
text(0.5, 0.20, "Secondary fixed-horizon cohorts", font = 2, cex = 0.85)
flow_box(0.27, 0.12, sprintf("30-day cohort, 2015–2025\nn = %s; %s deaths", fmt(nrow(early$data)), fmt(sum(early$data$event))), w = 0.38, fill = "#F4F6F7")
flow_box(0.73, 0.12, sprintf("5-year cohort, 2015–2020\nn = %s; %s deaths", fmt(nrow(long$data)), fmt(sum(long$data$event))), w = 0.38, fill = "#F4F6F7")
dev.off()

# Figure 2: annual case-mix trends.
png(file.path(study_config$figures_dir, "figure2_case_mix_trends.png"), width = 2600, height = 1800, res = 220)
par(mfrow = c(2, 2), mar = c(4, 4.5, 2.5, 1))
plot(case_mix$index_year, case_mix$n, type = "b", pch = 19, col = "#1F4E79", xlab = "Index year", ylab = "Patients", main = "Annual index PCI cohort")
plot(case_mix$index_year, case_mix$median_age, type = "b", pch = 19, col = "#7D3C98", xlab = "Index year", ylab = "Median age, years", main = "Age")
plot(case_mix$index_year, case_mix$female_percent, type = "b", pch = 19, col = "#B03A2E", xlab = "Index year", ylab = "Percent", main = "Women")
plot(case_mix$index_year, case_mix$primary_pci_percent, type = "b", pch = 19, col = "#117864", xlab = "Index year", ylab = "Percent", main = "Text-classified primary PCI")
dev.off()

# Figure 3: crude and adjusted fixed-horizon mortality.
png(file.path(study_config$figures_dir, "figure3_adjusted_mortality_trends.png"), width = 2600, height = 1300, res = 220)
par(mfrow = c(1, 2), mar = c(4.2, 4.5, 2.5, 1))
for (outcome_name in c("30-day all-cause mortality", "1-year all-cause mortality")) {
  adj <- all_annual[all_annual$outcome == outcome_name & all_annual$model == selected_model_name, ]
  cr <- all_crude[all_crude$outcome == outcome_name, ]
  ylim <- c(0, 1.18 * max(c(adj$upper, cr$crude_risk), na.rm = TRUE))
  plot(adj$index_year, 100 * adj$risk, type = "n", ylim = 100 * ylim, las = 1,
    xlab = "Index year", ylab = "Mortality, %", main = outcome_name)
  polygon(c(adj$index_year, rev(adj$index_year)), c(100 * adj$lower, rev(100 * adj$upper)),
    col = adjustcolor("#2E86C1", alpha.f = 0.20), border = NA)
  lines(adj$index_year, 100 * adj$risk, type = "b", pch = 19, col = "#1F4E79", lwd = 2)
  lines(cr$index_year, 100 * cr$crude_risk, type = "b", pch = 1, col = "#922B21", lwd = 1.5)
  legend("topleft", c(paste(selected_model_label, "(95% CI)"), "Crude"), col = c("#1F4E79", "#922B21"),
    pch = c(19, 1), lty = 1, bty = "n", cex = 0.85)
}
dev.off()

# Figure 4: 5-year Kaplan-Meier survival for mature eras.
km <- long$data
km$era <- factor(ifelse(km$index_year <= 2017, "2015–2017", "2018–2020"))
km_fit <- survfit(Surv(pmin(time_days, 1825), event) ~ era, data = km)
png(file.path(study_config$figures_dir, "figure4_five_year_survival.png"), width = 1800, height = 1500, res = 220)
par(mar = c(4.5, 4.7, 2.5, 1))
plot(km_fit, col = c("#1F4E79", "#B03A2E"), lwd = 2, conf.int = TRUE,
  xlab = "Days after index PCI", ylab = "Survival probability", xlim = c(0, 1825), ylim = c(0.70, 1.00), mark.time = FALSE)
legend("bottomleft", levels(km$era), col = c("#1F4E79", "#B03A2E"), lwd = 2, bty = "n")
title("Five-year all-cause survival")
dev.off()

# Scalar result contract for manuscript generation.
get_contrast <- function(result, outcome, model = selected_model_name) {
  result[result$outcome == outcome & result$model == model, , drop = FALSE][1, ]
}
p1 <- get_contrast(all_contrasts, "1-year all-cause mortality")
p30 <- get_contrast(all_contrasts, "30-day all-cause mortality")
p5 <- get_contrast(all_contrasts, "5-year all-cause mortality")

key_values <- data.frame(
  key = c(
    "analysis_mode", "selected_model", "source_records", "unique_patients",
    "outside_primary_period", "invalid_death_before_index", "invalid_death_before_index_primary",
    "primary_n", "primary_events", "primary_crude_percent", "primary_model_n",
    "primary_start_adjusted_percent", "primary_end_adjusted_percent",
    "primary_risk_difference_points", "primary_rd_lower_points", "primary_rd_upper_points",
    "primary_p_year", "early_n", "early_events", "early_crude_percent",
    "early_start_adjusted_percent", "early_end_adjusted_percent", "early_p_year",
    "long_n", "long_events", "long_crude_percent", "long_start_adjusted_percent",
    "long_end_adjusted_percent", "long_p_year", "median_age", "female_percent",
    "study_start_date", "study_end_date"
  ),
  value = c(
    analysis_mode,
    selected_model_name,
    source_rows,
    nrow(d),
    outside_primary_period,
    sum(d$death_before_index == 1, na.rm = TRUE),
    death_before_index_primary,
    nrow(primary$data),
    sum(primary$data$event),
    100 * mean(primary$data$event),
    p1$model_n,
    100 * p1$start_risk,
    100 * p1$end_risk,
    100 * p1$risk_difference,
    100 * p1$rd_lower,
    100 * p1$rd_upper,
    p1$p_year,
    nrow(early$data),
    sum(early$data$event),
    100 * mean(early$data$event),
    100 * p30$start_risk,
    100 * p30$end_risk,
    p30$p_year,
    nrow(long$data),
    sum(long$data$event),
    100 * mean(long$data$event),
    100 * p5$start_risk,
    100 * p5$end_risk,
    p5$p_year,
    median(primary$data$age, na.rm = TRUE),
    100 * mean(primary$data$sex == "Female", na.rm = TRUE),
    format(min(primary$data$index_date), "%Y-%m-%d"),
    format(max(primary$data$index_date), "%Y-%m-%d")
  ),
  stringsAsFactors = FALSE
)
write_csv_safe(key_values, file.path(study_config$results_dir, "key_results.csv"))

capture.output(summary(primary$cox), file = file.path(study_config$results_dir, "primary_cox_summary.txt"))
capture.output(sessionInfo(), file = file.path(study_config$audit_dir, "r_session_info.txt"))

cat(sprintf("Analysis mode: %s\n", analysis_mode))
cat(sprintf("Primary cohort: %s patients, %s one-year deaths\n", format(nrow(primary$data), big.mark = ","), format(sum(primary$data$event), big.mark = ",")))
cat(sprintf("Adjusted one-year mortality: %.2f%% in 2015 vs %.2f%% in 2024; p(year)=%.4g\n", 100 * p1$start_risk, 100 * p1$end_risk, p1$p_year))
