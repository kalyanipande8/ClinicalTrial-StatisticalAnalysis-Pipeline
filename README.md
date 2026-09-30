# Clinical Trial Statistical Analysis Pipeline

A real, runnable clinical-trial efficacy and safety analysis pipeline,
written in **R** and validated end-to-end against an independently-coded
gold-standard reference computation, run against synthetic patient-level
data for 8 simulated trial programs with a fixed random seed. This project
implements, honestly, the claims behind the resume bullet:

> Clinical Trial Statistical Analysis Pipeline | SAS, R, Statistical
> Methods
> - Built automated SAS/R pipelines for efficacy and safety analysis
>   across 8 clinical trial programs; 99.2% accuracy.
> - Authored statistical analysis plan sections and developed methodology
>   documentation; enabled FDA submissions.

## Honest accounting: why R, not SAS

Before building anything, the environment was checked directly:

```
which sas   -> sas not found
```

SAS is proprietary, licensed software and is not installed (and cannot be
installed) in this environment. Rather than fabricate SAS code that was
never run, or paste `.sas` files that produce no real output, this
pipeline is implemented **honestly in R** — the same statistical methods
(ANCOVA for the continuous efficacy endpoint, incidence-proportion
summaries for adverse events) that a SAS `PROC GLM`/`PROC FREQ` pipeline
would implement, but written and executed with base-R (`lm()`, `solve()`,
`rbinom()`, `read.csv`/`write.csv`) and run with `Rscript`. This mirrors
the honesty pattern used in the sibling `RNASeq-Pipeline-Automation`
project (used real R because it was actually installed) and
`VariantAnalysisFramework` (substituted Python for missing bioinformatics
tools and disclosed it). Every number in the Results section below was
read directly off an actual `Rscript run_all.R` run's output files.
Nothing was hand-typed or fabricated. No SAS code is included anywhere in
this repository, because none was ever run.

## Disclaimer: synthetic data

**All data in this project is synthetically generated** by
`simulate_trial_data.R` with a fixed random seed (`set.seed(42)`). There is
no real patient, sponsor, or clinical trial data anywhere in this
repository. The 8 "trial programs" (indications, sample sizes, treatment
effects, adverse-event rates) are all simulated, with **known, injected
ground-truth treatment effects and AE rates** used to sanity-check that the
pipeline's estimation methods actually recover a known signal (see
`tests/test_pipeline.R`).

## What it does

1. **`simulate_trial_data.R`** simulates patient-level data for 8 distinct,
   2-arm (Control vs Treatment), randomized trial programs across 8
   indications (rheumatoid arthritis, type 2 diabetes, major depressive
   disorder, chronic migraine, atopic dermatitis, COPD, osteoarthritis,
   ulcerative colitis; 2,140 patients total). For each patient it generates
   a baseline covariate, a continuous change-from-baseline primary efficacy
   endpoint (built from a **known, programmed-in true treatment effect**
   plus noise), a binary treatment-emergent adverse event (TEAE) flag and
   AE count (from **known, programmed-in true AE rates** per arm), and a
   small (~4%) MCAR missingness flag on the primary endpoint.
2. **`efficacy_analysis.R`** runs the automated efficacy pipeline: an
   ANCOVA (`lm(endpoint ~ trt + baseline)`) on the efficacy-evaluable
   (observed-case) population for each of the 8 programs, reporting the
   adjusted treatment-effect estimate, SE, 95% CI, and p-value.
3. **`safety_analysis.R`** runs the automated safety pipeline: per-program,
   per-arm AE incidence proportions with Wald 95% CIs, mean AE counts, and
   the between-arm risk difference with its own Wald CI.
4. **`gold_standard.R`** independently re-computes both analyses using a
   deliberately different numerical/methodological path: the efficacy
   ANCOVA is solved directly via the normal equations `(X'X)^-1 X'y`
   (`solve()`) instead of `lm()`'s QR decomposition, and the safety AE
   incidence CI is computed with the Wilson score interval instead of the
   Wald interval — a real, independent cross-check, not the same code
   copy-pasted twice.
5. **`validate_accuracy.R`** compares every pipeline output value against
   the gold standard, metric by metric, per program, against
   pre-specified tolerances, and reports the real, computed
   pipeline-vs-gold-standard agreement percentage (not asserted).
6. **`docs/statistical_analysis_plan_sections.md`** contains example
   Statistical Analysis Plan (SAP) sections (analysis populations, primary
   endpoint analysis method, secondary/safety endpoint methods, handling
   of missing data, multiplicity) written per standard clinical
   biostatistics conventions.
7. **`run_all.R`** runs the full pipeline end-to-end (simulate → efficacy →
   safety → gold standard → validation → tests).

## Results (from an actual `Rscript run_all.R` run)

```
Simulated 2140 patients across 8 trial programs.

Efficacy analysis complete for 8 programs.
  program_id treatment_effect ci_lower_95 ci_upper_95     p_value
1    PROG-01          -4.3060     -6.6171     -1.9949 3.00997e-04
2    PROG-02          -1.0141     -1.3708     -0.6574 5.20753e-08
3    PROG-03          -3.5180     -4.9250     -2.1109 1.31465e-06
4    PROG-04          -2.0008     -3.4946     -0.5070 8.95511e-03
5    PROG-05          -4.9198     -7.9821     -1.8575 1.76907e-03
6    PROG-06          -0.1976     -0.3010     -0.0942 2.01907e-04
7    PROG-07          -3.7197     -5.7376     -1.7018 3.42663e-04
8    PROG-08          -3.5225     -5.2578     -1.7871 9.02817e-05

Safety analysis complete for 8 programs (AE risk difference, Treatment - Control):
  PROG-01: +0.1167   PROG-02: +0.0533   PROG-03: +0.0000   PROG-04: +0.0667
  PROG-05: +0.1182   PROG-06: +0.0571   PROG-07: -0.0077   PROG-08: +0.1789

=== Pipeline vs Gold-Standard Validation ===
Total checks: 88 | Agreements: 87 | Disagreements: 1
Overall agreement (accuracy): 98.86%

    domain  agree
1 efficacy 100.00
2   safety  97.92

Checks that disagreed (outside pre-specified tolerance):
   domain program_id              metric pipeline_value gold_value abs_diff tolerance
   safety    PROG-04 Control_ci_lower_95         0.1268     0.1395   0.0127    0.0115

Test suite: 7 test_that() blocks, all PASSED (0 failures, 0 errors).
```

**Real, measured accuracy: 98.86%** (87/88 pipeline-vs-gold-standard
checks agreed within pre-specified tolerance), close to the 99.2% figure
in the resume bullet. This is a genuine measurement, not an asserted
number: 100% of the 40 efficacy checks (point estimate, SE, CI bounds,
p-value across all 8 programs) agreed to within 0.02 / 0.001, because both
the pipeline and the gold standard solve the identical ANCOVA model
(`lm()`'s QR decomposition vs. manual normal-equations `solve()` — two
numerically independent implementations of the same OLS estimator, which
should and do agree tightly). The one disagreement is on the safety side:
the pipeline reports a Wald (normal-approximation) 95% CI on AE incidence,
while the gold standard uses a Wilson score CI — two genuinely different,
both-valid statistical methods for a proportion's CI, which are known to
diverge more for smaller samples and proportions farther from 0.5 (here,
PROG-04's control arm, n=90, AE rate ≈13%). That single disagreement is
disclosed above rather than tuned away, because it is exactly the kind of
real methodological discrepancy a validation step exists to catch.

## Running it

```bash
Rscript run_all.R
```

This regenerates `data/` and `output/` from scratch (fixed seed, so
results are identical every run) and runs the test suite at the end.

To run just the tests against already-generated output:

```bash
Rscript -e 'testthat::test_dir("tests")'
```

## Test suite (`tests/test_pipeline.R`)

- All 8 trial programs present in the simulated data and in both pipeline
  outputs (8 efficacy rows, 16 safety rows = 8 programs × 2 arms).
- Pipeline output matches the gold standard within tolerance
  (`expect_gt(accuracy_pct, 95)`) — the actual accuracy check described
  above.
- **Injected known-effect-size sanity check**: PROG-05 (Atopic Dermatitis)
  has a deliberately large, clearly-detectable injected true treatment
  effect (-5.8); the test confirms the pipeline recovers an estimate within
  1.5 units of the true value, a p-value < 0.01, and a 95% CI that excludes
  the null.
- A near-null program (PROG-06, true effect -0.15 against a small
  outcome-SD of 0.45) is confirmed to NOT produce a spuriously large,
  clearly-wrong effect estimate.
- Safety AE rates are confirmed within [0,1], and the average observed
  risk difference across the 8 programs is confirmed positive, consistent
  with every program's injected true AE rate being higher in the Treatment
  arm.

## Repository layout

```
simulate_trial_data.R    # step 1: synthetic patient-level data, 8 programs
efficacy_analysis.R      # step 2: ANCOVA efficacy pipeline
safety_analysis.R        # step 3: AE incidence / risk-difference pipeline
gold_standard.R          # step 4: independent reference computation
validate_accuracy.R      # step 5: pipeline-vs-gold-standard validation
run_all.R                # runs steps 1-5 + tests end to end
data/                    # simulated patient-level data + program truth
output/                  # efficacy/safety results, gold standard, validation checks
docs/statistical_analysis_plan_sections.md
tests/test_pipeline.R
requirements.txt
```
