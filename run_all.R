#!/usr/bin/env Rscript
# run_all.R -- end-to-end pipeline runner
#
# 1. simulate patient-level data for 8 trial programs
# 2. run the automated efficacy + safety analysis pipeline
# 3. compute the independent gold-standard reference
# 4. validate pipeline output against gold standard, report accuracy
# 5. run the test suite

cat("=== Step 1/5: Simulating patient-level trial data ===\n")
source("simulate_trial_data.R")

cat("\n=== Step 2/5: Running efficacy analysis pipeline ===\n")
source("efficacy_analysis.R")

cat("\n=== Step 3/5: Running safety analysis pipeline ===\n")
source("safety_analysis.R")

cat("\n=== Step 4/5: Computing gold-standard reference ===\n")
source("gold_standard.R")

cat("\n=== Step 5/5: Validating pipeline vs gold standard ===\n")
source("validate_accuracy.R")

cat("\n=== Running test suite ===\n")
test_result <- testthat::test_dir("tests", reporter = "summary")

cat("\nPipeline run complete. See output/ for all result files.\n")
