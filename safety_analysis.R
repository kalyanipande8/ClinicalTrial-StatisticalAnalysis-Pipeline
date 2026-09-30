#!/usr/bin/env Rscript
# safety_analysis.R
#
# Automated safety analysis across all 8 trial programs.
# Population: the safety set = all randomized/dosed patients (no exclusion
# for missing efficacy data, per standard convention that the safety
# population is broader than the efficacy population).
#
# For each program x arm: AE incidence (n with >=1 TEAE / N), AE rate,
# Wald (normal-approximation) 95% CI on the incidence proportion, and the
# risk difference (Treatment - Control) with its own Wald CI. Also reports
# a simple overall AE rate summary table.

suppressWarnings(suppressMessages({
  patients <- read.csv("data/patient_level_data.csv", stringsAsFactors = FALSE)
}))

dir.create("output", showWarnings = FALSE)

wald_ci <- function(x, n, z = 1.96) {
  p <- x / n
  se <- sqrt(p * (1 - p) / n)
  c(p = p, lower = max(0, p - z * se), upper = min(1, p + z * se), se = se)
}

programs <- sort(unique(patients$program_id))
ae_table <- data.frame()
rd_table <- data.frame()

for (pid in programs) {
  d <- patients[patients$program_id == pid, ]

  for (arm_name in c("Control", "Treatment")) {
    da <- d[d$arm == arm_name, ]
    n <- nrow(da)
    x <- sum(da$any_ae)
    ci <- wald_ci(x, n)
    ae_table <- rbind(ae_table, data.frame(
      program_id = pid, arm = arm_name, n = n, n_with_ae = x,
      ae_rate = round(ci["p"], 4),
      ci_lower_95 = round(ci["lower"], 4),
      ci_upper_95 = round(ci["upper"], 4),
      mean_ae_count = round(mean(da$ae_count), 3),
      stringsAsFactors = FALSE
    ))
  }

  dc <- d[d$arm == "Control", ]
  dt <- d[d$arm == "Treatment", ]
  pc <- mean(dc$any_ae); pt <- mean(dt$any_ae)
  nc <- nrow(dc); nt <- nrow(dt)
  rd <- pt - pc
  se_rd <- sqrt(pc * (1 - pc) / nc + pt * (1 - pt) / nt)
  rd_table <- rbind(rd_table, data.frame(
    program_id = pid,
    ae_rate_control = round(pc, 4),
    ae_rate_treatment = round(pt, 4),
    risk_difference = round(rd, 4),
    rd_ci_lower_95 = round(rd - 1.96 * se_rd, 4),
    rd_ci_upper_95 = round(rd + 1.96 * se_rd, 4),
    stringsAsFactors = FALSE
  ))
}

write.csv(ae_table, "output/safety_ae_incidence.csv", row.names = FALSE)
write.csv(rd_table, "output/safety_risk_difference.csv", row.names = FALSE)

cat("Safety analysis complete for", length(programs), "programs.\n")
print(rd_table)
