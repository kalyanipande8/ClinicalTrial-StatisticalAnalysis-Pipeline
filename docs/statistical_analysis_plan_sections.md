# Statistical Analysis Plan — Example Sections

These sections are written in the style and structure of a real clinical
trial Statistical Analysis Plan (SAP), for a hypothetical Phase III,
randomized, double-blind, placebo/active-controlled, 2-arm superiority
trial (the design mirrored by the 8 simulated programs in this repository).
They are illustrative/template sections authored for this project, not
sections from any real, sponsor-specific, or submitted SAP.

## 1. Analysis Populations

- **Intent-to-Treat (ITT) Population**: All randomized subjects, analyzed
  according to their randomized treatment assignment regardless of
  treatment actually received. Used for the primary efficacy analysis
  unless otherwise specified, consistent with ICH E9 guidance on
  preserving the benefits of randomization.
- **Modified Intent-to-Treat (mITT) / Efficacy-Evaluable Population**: All
  randomized subjects who received at least one dose of study drug and
  have at least one post-baseline assessment of the primary endpoint. In
  this repository's pipeline (`efficacy_analysis.R`), the evaluable set is
  operationalized as randomized subjects with a non-missing primary
  endpoint value (`endpoint_chg`), and the ANCOVA is fit on that set.
- **Safety Population**: All randomized subjects who received at least one
  dose of study drug, analyzed according to treatment actually received.
  This is the population used for all adverse event summaries
  (`safety_analysis.R`); it is not restricted by efficacy-endpoint
  completeness, which is why the safety population can differ in size from
  the efficacy-evaluable population within the same program.
- **Per-Protocol (PP) Population**: All mITT subjects without major
  protocol deviations judged to materially impact the assessment of the
  primary endpoint. Used only as a sensitivity analysis for the primary
  endpoint, not the primary analysis. (Not separately implemented in this
  repository's simulated example, since no protocol-deviation variable was
  simulated; noted here for completeness of a real SAP.)

## 2. Primary Endpoint Analysis

**Primary endpoint**: Change from baseline to the primary assessment
timepoint in the disease-specific severity score (continuous).

**Primary analysis method**: Analysis of Covariance (ANCOVA), with the
change-from-baseline endpoint as the dependent variable, treatment group
as a fixed factor, and the baseline value of the endpoint as a covariate.
Baseline adjustment is pre-specified because it (a) reflects the
randomized-but-imbalanced-by-chance baseline severity that commonly
remains after randomization in moderate sample sizes, and (b) increases
statistical power relative to an unadjusted comparison of change scores,
per standard practice (e.g., FDA/EMA guidance on covariate adjustment in
clinical trials).

The treatment effect is reported as the ANCOVA-adjusted between-group
least-squares mean difference (Treatment − Control), together with its
two-sided 95% confidence interval and the associated p-value from the
t-test on the treatment coefficient. A two-sided alpha of 0.05 is used for
the primary hypothesis test unless a program's protocol specifies a
different alpha (e.g., after adjustment for multiplicity across co-primary
endpoints).

**Estimand framework**: The primary estimand targets the treatment effect
that would have been observed had all randomized subjects adhered to
assigned treatment for the full study duration (a "treatment policy"-like
strategy for intercurrent events other than the handling of missing data
itself, which is addressed separately in Section 4). Sensitivity analyses
under alternative intercurrent-event strategies (e.g., a hypothetical
strategy censoring data after rescue medication use) are pre-specified
where applicable to a given program's protocol but are not implemented in
this repository's simplified simulation.

## 3. Secondary and Safety Endpoint Analyses

- **Key secondary continuous endpoints**: Analyzed with the same ANCOVA
  approach as the primary endpoint, adjusted for multiplicity using a
  pre-specified fixed-sequence or Hochberg procedure when hierarchically
  ordered with the primary endpoint (not implemented numerically in this
  repository, which focuses on one primary endpoint per program).
- **Adverse events (safety)**: Treatment-emergent adverse events (TEAEs)
  are summarized by treatment arm as the number and percentage of subjects
  in the safety population reporting at least one TEAE, with an exact or
  approximate 95% confidence interval on the incidence proportion. This
  repository's pipeline (`safety_analysis.R`) reports the Wald
  (normal-approximation) CI as the primary reporting method, and the risk
  difference (Treatment − Control) in AE incidence with its own Wald CI —
  a standard, widely used approach for AE incidence tables in clinical
  study reports. `gold_standard.R` independently validates the AE rate
  point estimates and cross-checks the CI against the Wilson score
  interval, an alternative method with better small-sample coverage
  properties; see `README.md` for the resulting measured agreement and an
  honest discussion of where the two CI methods diverge.
- AE severity grading, seriousness, relatedness, and MedDRA
  System-Organ-Class-level tabulation are standard components of a real
  clinical study report safety analysis but are out of scope for this
  repository's simplified binary any-AE simulation.

## 4. Handling of Missing Data

- **Primary endpoint**: The primary analysis uses an "observed cases" /
  complete-case ANCOVA on the mITT-evaluable population (subjects with a
  non-missing primary endpoint value), which is a valid analysis under a
  Missing-At-Random (MAR) assumption conditional on the covariates in the
  model (here, baseline severity and treatment arm).
- **Sensitivity analyses** (specified in a real SAP, not separately
  re-implemented in this simplified repository): multiple imputation under
  MAR (e.g., using a baseline- and treatment-arm-conditional imputation
  model), and a tipping-point analysis to assess robustness of the primary
  conclusion to plausible Missing-Not-At-Random (MNAR) departures — for
  example, imputing worse outcomes for treatment-arm dropouts than for
  control-arm dropouts, and confirming the treatment effect remains
  directionally and (where feasible) statistically significant.
- **Safety data**: AE reporting is generally complete for subjects who
  remain on study (AEs are actively solicited/recorded at each visit up to
  discontinuation), so no imputation is applied to the safety population;
  subjects who discontinue early simply contribute a shorter at-risk
  period, which is a limitation noted narratively in the clinical study
  report rather than statistically imputed in the incidence-proportion
  analysis presented here.
- In this repository's simulation, missingness on the primary endpoint is
  generated at a small (~4%) Missing Completely At Random (MCAR) rate per
  `simulate_trial_data.R`, which is simpler than a realistic MAR/MNAR
  dropout mechanism but is sufficient to exercise the "evaluable
  population != randomized population" logic in the pipeline and tests.

## 5. Multiplicity

Where a program has multiple co-primary or key secondary endpoints, or
multiple dose arms compared against a common control, a pre-specified
multiplicity-control procedure (e.g., Hochberg, fixed-sequence gatekeeping,
or a Bonferroni-type correction) is applied to preserve the family-wise
type I error rate at the nominal two-sided alpha (typically 0.05). Each of
the 8 simulated programs in this repository has a single primary endpoint
and a single active-vs-control comparison, so no multiplicity adjustment
was required for the example analyses shown in `output/`.
