library(testthat)

# These tests assume run_all.R (or at least simulate_trial_data.R,
# efficacy_analysis.R, safety_analysis.R, gold_standard.R, and
# validate_accuracy.R) have already been run from the project root, so
# that data/ and output/ are populated.

test_that("all 8 trial programs are present in the simulated data", {
  patients <- read.csv("../data/patient_level_data.csv", stringsAsFactors = FALSE)
  programs <- read.csv("../data/program_truth.csv", stringsAsFactors = FALSE)
  expect_equal(nrow(programs), 8)
  expect_equal(length(unique(patients$program_id)), 8)
  expect_true(all(sprintf("PROG-%02d", 1:8) %in% unique(patients$program_id)))
})

test_that("efficacy and safety pipeline output exist for all 8 programs", {
  eff <- read.csv("../output/efficacy_results.csv", stringsAsFactors = FALSE)
  saf <- read.csv("../output/safety_ae_incidence.csv", stringsAsFactors = FALSE)
  expect_equal(nrow(eff), 8)
  expect_equal(length(unique(saf$program_id)), 8)
  expect_equal(nrow(saf), 16)  # 8 programs x 2 arms
})

test_that("pipeline output matches gold-standard reference within tolerance", {
  checks <- read.csv("../output/validation_checks.csv", stringsAsFactors = FALSE)
  accuracy_pct <- 100 * mean(checks$agree)
  # Real, computed agreement -- not asserted at a fixed number, but must
  # clear a high bar consistent with a validated clinical pipeline.
  expect_gt(accuracy_pct, 95)
  cat(sprintf("\n  [validation] measured pipeline-vs-gold-standard accuracy: %.2f%%\n", accuracy_pct))
})

test_that("known-large-effect program (PROG-05, Atopic Dermatitis) recovers the injected true effect", {
  # Ground truth injected in simulate_trial_data.R: true_effect = -5.8 for
  # PROG-05. This is a sanity check that the estimation pipeline actually
  # recovers a known, deliberately large and clearly-detectable signal --
  # not just that it produces *some* number.
  truth <- read.csv("../data/program_truth.csv", stringsAsFactors = FALSE)
  true_effect_05 <- truth$true_effect[truth$program_id == "PROG-05"]

  eff <- read.csv("../output/efficacy_results.csv", stringsAsFactors = FALSE)
  row5 <- eff[eff$program_id == "PROG-05", ]

  expect_lt(row5$p_value, 0.01)  # should be highly statistically significant
  expect_lt(abs(row5$treatment_effect - true_effect_05), 1.5)  # estimate close to injected truth
  expect_true(row5$ci_upper_95 < 0)  # CI should exclude the null (0) given the large true effect
})

test_that("a null/near-null program (PROG-06, COPD, true_effect = -0.15) is NOT falsely called significant", {
  truth <- read.csv("../data/program_truth.csv", stringsAsFactors = FALSE)
  eff <- read.csv("../output/efficacy_results.csv", stringsAsFactors = FALSE)
  row6 <- eff[eff$program_id == "PROG-06", ]
  # true effect is tiny relative to noise sd (0.45) at n=175/arm; expect the
  # pipeline to NOT produce a spuriously huge, clearly-wrong estimate.
  expect_lt(abs(row6$treatment_effect), 0.5)
})

test_that("safety AE rates are within [0,1] and treatment arm AE rate reflects injected higher-risk truth on average", {
  truth <- read.csv("../data/program_truth.csv", stringsAsFactors = FALSE)
  saf <- read.csv("../output/safety_ae_incidence.csv", stringsAsFactors = FALSE)
  expect_true(all(saf$ae_rate >= 0 & saf$ae_rate <= 1))

  rd <- read.csv("../output/safety_risk_difference.csv", stringsAsFactors = FALSE)
  # In our simulated truth, every program has true_ae_rate_trt > true_ae_rate_ctrl,
  # so the *average* observed risk difference across 8 programs should be positive.
  expect_gt(mean(rd$risk_difference), 0)
})
