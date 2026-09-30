#!/usr/bin/env Rscript
# efficacy_analysis.R
#
# Automated efficacy analysis across all 8 trial programs.
# Primary analysis: ANCOVA of change-from-baseline endpoint on treatment
# arm, adjusting for baseline covariate (standard approach for a
# continuous primary endpoint in a randomized trial). Population: the
# efficacy-evaluable / modified-intent-to-treat set (excludes patients with
# a missing primary endpoint value per the pre-specified handling of
# missing data — see docs/statistical_analysis_plan_sections.md).
#
# Fitted with base-R lm() (QR-decomposition least squares) — a standard,
# well-validated numerical method, used here as the "pipeline" method.

suppressWarnings(suppressMessages({
  patients <- read.csv("data/patient_level_data.csv", stringsAsFactors = FALSE)
}))

dir.create("output", showWarnings = FALSE)

programs <- sort(unique(patients$program_id))
results <- data.frame()

for (pid in programs) {
  d <- patients[patients$program_id == pid, ]
  d_eval <- d[!is.na(d$endpoint_chg), ]  # mITT / evaluable set: observed cases

  fit <- lm(endpoint_chg ~ trt + baseline, data = d_eval)
  s <- summary(fit)

  est   <- unname(coef(fit)["trt"])
  se    <- s$coefficients["trt", "Std. Error"]
  pval  <- s$coefficients["trt", "Pr(>|t|)"]
  df_res <- fit$df.residual
  tcrit <- qt(0.975, df_res)
  ci_lo <- est - tcrit * se
  ci_hi <- est + tcrit * se

  n_total <- nrow(d)
  n_eval  <- nrow(d_eval)

  results <- rbind(results, data.frame(
    program_id       = pid,
    n_randomized     = n_total,
    n_evaluable      = n_eval,
    pct_missing      = round(100 * (n_total - n_eval) / n_total, 2),
    treatment_effect = round(est, 4),
    se               = round(se, 4),
    ci_lower_95      = round(ci_lo, 4),
    ci_upper_95      = round(ci_hi, 4),
    p_value          = signif(pval, 6),
    df_residual      = df_res,
    method           = "ANCOVA (lm, baseline-adjusted)",
    stringsAsFactors = FALSE
  ))
}

write.csv(results, "output/efficacy_results.csv", row.names = FALSE)
cat("Efficacy analysis complete for", length(programs), "programs.\n")
print(results[, c("program_id", "treatment_effect", "ci_lower_95", "ci_upper_95", "p_value")])
