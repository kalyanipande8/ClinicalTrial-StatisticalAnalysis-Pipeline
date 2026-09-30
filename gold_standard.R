#!/usr/bin/env Rscript
# gold_standard.R
#
# Independent "trusted reference" computation, coded separately from
# efficacy_analysis.R / safety_analysis.R, used to validate the pipeline.
#
# Efficacy gold standard: solves the same ANCOVA model (endpoint ~ trt +
# baseline) directly via the normal equations (X'X)^-1 X'y using base
# matrix algebra (solve()), rather than lm()'s QR decomposition. Same
# statistical method (OLS ANCOVA), independent numerical path.
#
# Safety gold standard: AE incidence proportions computed directly from
# table() counts (same point estimate as the pipeline), but the 95% CI is
# computed with the Wilson score interval — a different, generally more
# accurate small-sample method than the pipeline's Wald normal-approximation
# CI. This is a deliberate methodological difference (not a bug) and is
# exactly the kind of real-world discrepancy a validation step is meant to
# catch and characterize.

suppressWarnings(suppressMessages({
  patients <- read.csv("data/patient_level_data.csv", stringsAsFactors = FALSE)
}))

dir.create("output", showWarnings = FALSE)

programs <- sort(unique(patients$program_id))

# ---- efficacy gold standard: manual OLS via normal equations -----------
gold_efficacy <- data.frame()
for (pid in programs) {
  d <- patients[patients$program_id == pid, ]
  d_eval <- d[!is.na(d$endpoint_chg), ]

  y <- d_eval$endpoint_chg
  X <- cbind(intercept = 1, trt = d_eval$trt, baseline = d_eval$baseline)
  n <- nrow(X); k <- ncol(X)

  XtX_inv <- solve(t(X) %*% X)
  beta <- XtX_inv %*% t(X) %*% y
  resid <- y - X %*% beta
  sigma2 <- sum(resid^2) / (n - k)
  vcov <- sigma2 * XtX_inv
  se_trt <- sqrt(vcov["trt", "trt"])
  est_trt <- beta["trt", 1]
  df_res <- n - k
  tstat <- est_trt / se_trt
  pval <- 2 * pt(-abs(tstat), df_res)
  tcrit <- qt(0.975, df_res)

  gold_efficacy <- rbind(gold_efficacy, data.frame(
    program_id = pid,
    treatment_effect = round(est_trt, 4),
    se = round(se_trt, 4),
    ci_lower_95 = round(est_trt - tcrit * se_trt, 4),
    ci_upper_95 = round(est_trt + tcrit * se_trt, 4),
    p_value = signif(pval, 6),
    df_residual = df_res,
    method = "OLS via normal equations (solve)",
    stringsAsFactors = FALSE
  ))
}
write.csv(gold_efficacy, "output/gold_efficacy_results.csv", row.names = FALSE)

# ---- safety gold standard: table() counts + Wilson score CI ------------
wilson_ci <- function(x, n, z = 1.96) {
  p <- x / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half <- (z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / denom
  c(lower = max(0, center - half), upper = min(1, center + half))
}

gold_safety <- data.frame()
for (pid in programs) {
  d <- patients[patients$program_id == pid, ]
  for (arm_name in c("Control", "Treatment")) {
    da <- d[d$arm == arm_name, ]
    tab <- table(factor(da$any_ae, levels = c(0, 1)))
    n <- sum(tab); x <- tab[["1"]]
    p <- x / n
    ci <- wilson_ci(x, n)
    gold_safety <- rbind(gold_safety, data.frame(
      program_id = pid, arm = arm_name, n = n, n_with_ae = x,
      ae_rate = round(p, 4),
      ci_lower_95 = round(ci["lower"], 4),
      ci_upper_95 = round(ci["upper"], 4),
      stringsAsFactors = FALSE
    ))
  }
}
write.csv(gold_safety, "output/gold_safety_ae_incidence.csv", row.names = FALSE)

cat("Gold-standard reference computation complete for", length(programs), "programs.\n")
