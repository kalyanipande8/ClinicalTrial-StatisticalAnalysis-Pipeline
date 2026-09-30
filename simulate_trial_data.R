#!/usr/bin/env Rscript
# simulate_trial_data.R
#
# Simulates patient-level data for 8 distinct clinical trial programs
# (Phase II/III style, 2-arm randomized, placebo/active-control vs
# investigational drug). For each program we simulate:
#   - demographics / randomization (treatment arm, baseline covariate)
#   - one continuous primary efficacy endpoint (e.g. change from baseline
#     in a disease severity score), generated with a KNOWN, programmed-in
#     true treatment effect so we have ground truth to validate against
#   - a binary safety endpoint per patient: at least one treatment-emergent
#     adverse event (AE), plus an AE count, with a KNOWN true AE rate per
#     arm and per program
#
# All data is synthetic. Fixed seed -> fully reproducible.

set.seed(42)

dir.create("data", showWarnings = FALSE)

# ---- 8 trial program definitions ---------------------------------------
# true_effect: true mean treatment - control difference on the primary
#              continuous endpoint (negative = improvement, e.g. score
#              reduction), used as ground truth for the sanity check
# true_ae_rate_ctrl / true_ae_rate_trt: true underlying AE incidence
#              probabilities used to generate the safety data
programs <- data.frame(
  program_id      = sprintf("PROG-%02d", 1:8),
  indication      = c("Rheumatoid Arthritis", "Type 2 Diabetes", "Major Depressive Disorder",
                       "Chronic Migraine", "Atopic Dermatitis", "COPD",
                       "Osteoarthritis", "Ulcerative Colitis"),
  n_per_arm       = c(120, 150, 200, 90, 110, 175, 130, 95),
  true_effect     = c(-4.5, -0.8, -3.2, -2.1, -5.8, -0.15, -3.9, -2.6),
  sd_endpoint     = c(9.0, 1.6, 7.0, 5.0, 10.5, 0.45, 8.2, 6.5),
  true_ae_rate_ctrl = c(0.22, 0.18, 0.30, 0.25, 0.35, 0.28, 0.20, 0.32),
  true_ae_rate_trt  = c(0.31, 0.24, 0.34, 0.29, 0.48, 0.33, 0.27, 0.44),
  stringsAsFactors = FALSE
)
write.csv(programs, "data/program_truth.csv", row.names = FALSE)

all_patients <- list()

for (i in seq_len(nrow(programs))) {
  p <- programs[i, ]
  n <- p$n_per_arm

  # randomize 1:1, control = 0, treatment = 1
  arm <- rep(c("Control", "Treatment"), each = n)
  trt <- rep(c(0, 1), each = n)

  # baseline covariate (disease severity at baseline), balanced by design
  baseline <- round(rnorm(2 * n, mean = 50, sd = 8), 1)

  # primary continuous efficacy endpoint = change from baseline
  # true model: endpoint = beta0 + beta_baseline*baseline_centered + true_effect*trt + noise
  baseline_c <- baseline - mean(baseline)
  beta_baseline <- 0.15  # small baseline-dependent drift, ANCOVA will adjust for it
  noise <- rnorm(2 * n, mean = 0, sd = p$sd_endpoint)
  endpoint <- round(beta_baseline * baseline_c + p$true_effect * trt + noise, 2)

  # safety: per-patient AE flag (>=1 TEAE) and AE count (Poisson given flag)
  ae_prob <- ifelse(trt == 1, p$true_ae_rate_trt, p$true_ae_rate_ctrl)
  any_ae <- rbinom(2 * n, 1, ae_prob)
  ae_count <- ifelse(any_ae == 1, 1 + rpois(2 * n, lambda = 0.8), 0)

  # a couple of realistic missing-data / discontinuation flags (MCAR-ish, small %)
  missing_flag <- rbinom(2 * n, 1, 0.04)
  endpoint_obs <- endpoint
  endpoint_obs[missing_flag == 1] <- NA

  df <- data.frame(
    program_id   = p$program_id,
    subject_id   = sprintf("%s-%04d", p$program_id, seq_len(2 * n)),
    arm          = arm,
    trt          = trt,
    baseline     = baseline,
    endpoint_chg = endpoint_obs,
    completer    = 1 - missing_flag,
    any_ae       = any_ae,
    ae_count     = ae_count,
    stringsAsFactors = FALSE
  )
  all_patients[[i]] <- df
}

patients <- do.call(rbind, all_patients)
write.csv(patients, "data/patient_level_data.csv", row.names = FALSE)

cat(sprintf("Simulated %d patients across %d trial programs.\n",
            nrow(patients), nrow(programs)))
cat("Wrote data/patient_level_data.csv and data/program_truth.csv\n")
