# Sensitivity of the EWMA VaR models to the decay parameter lambda.
# The five main models are untouched; this re-runs only the EWMA volatility
# with other lambdas, for both the normal and the Student-t (df 5) multiplier.
#
# Needs R/var_models.R, R/backtest.R and R/charts.R to be sourced first.

source(file.path("R", "config.R"))

LAMBDAS <- c(0.90, 0.92, 0.94, 0.95, 0.96, 0.966, 0.97, 0.977, 0.98, 0.99)
MULTIPLIERS <- c(normal = "normal", student_t = "student_t")

IN_SAMPLE  <- as.Date(c("2015-01-01", "2019-12-31"))
OUT_SAMPLE <- as.Date(c("2020-01-01", "2025-12-31"))

# Number of days for a return's weight to halve: lambda^h = 1/2.
half_life <- function(lambda) log(0.5) / log(lambda)

multiplier_value <- function(m) switch(m, normal = z_normal(), student_t = z_student_t())

# Breach flags at the forecast dates for one lambda and one multiplier.
# Uses the same EWMA recursion and seed as the main model (ewma_variance()).
ewma_breaches <- function(ret, fc, lambda, multiplier) {
  idx <- match(fc$date, ret$date)
  v <- ewma_variance(ret$ret_port, lambda = lambda)[idx]
  -fc$ret_port > multiplier_value(multiplier) * sqrt(v)
}

lambda_sensitivity <- function(ret, fc) {
  n  <- nrow(fc)
  yr <- format(fc$date, "%Y")
  rows <- list()
  for (m in MULTIPLIERS) for (lam in LAMBDAS) {
    b  <- ewma_breaches(ret, fc, lam, m)
    x  <- sum(b)
    kp <- kupiec_pof(x, n)
    per_year <- tapply(b, yr, sum)
    rows[[length(rows) + 1]] <- data.frame(
      multiplier        = m,
      lambda            = lam,
      half_life_days    = half_life(lam),
      breaches          = x,
      expected_breaches = n * ALPHA,
      breach_rate       = x / n,
      kupiec_lr         = kp$statistic,
      kupiec_p_value    = kp$p_value,
      reject_at_5pct    = kp$p_value < 0.05,
      min_year_breaches = min(per_year),
      max_year_breaches = max(per_year)
    )
  }
  do.call(rbind, rows)
}

# Choose lambda on 2015-2019 only, then test it on 2020-2025.
# Selection rule: the smallest |breaches - expected| in-sample. Ties (one
# lambda above, one below by the same amount) go to the higher in-sample
# Kupiec p-value, then to the lambda closest to 0.94.
out_of_sample_check <- function(ret, fc) {
  ins  <- fc$date >= IN_SAMPLE[1]  & fc$date <= IN_SAMPLE[2]
  outs <- fc$date >= OUT_SAMPLE[1] & fc$date <= OUT_SAMPLE[2]
  rows <- list()
  for (m in MULTIPLIERS) {
    # Out-of-sample results are computed for every candidate only so that ties
    # can be reported; they play no part in the choice.
    cand <- do.call(rbind, lapply(LAMBDAS, function(lam) {
      b <- ewma_breaches(ret, fc, lam, m)
      x <- sum(b[ins])
      data.frame(multiplier = m, lambda = lam,
                 in_breaches = x, in_expected = sum(ins) * ALPHA,
                 in_gap = abs(x - sum(ins) * ALPHA),
                 in_p = kupiec_pof(x, sum(ins))$p_value,
                 out_breaches = sum(b[outs]),
                 out_p = kupiec_pof(sum(b[outs]), sum(outs))$p_value)
    }))
    cand <- cand[order(cand$in_gap, -cand$in_p, abs(cand$lambda - 0.94)), ]
    cand$tied_with_chosen <- cand$in_gap == cand$in_gap[1] & cand$in_p == cand$in_p[1]
    cand$chosen <- seq_len(nrow(cand)) == 1
    best <- cand[1, ]
    b  <- ewma_breaches(ret, fc, best$lambda, m)
    x  <- sum(b[outs])
    kp <- kupiec_pof(x, sum(outs))
    rows[[m]] <- data.frame(
      multiplier         = m,
      chosen_lambda      = best$lambda,
      half_life_days     = half_life(best$lambda),
      in_sample_days     = sum(ins),
      in_sample_breaches = best$in_breaches,
      in_sample_expected = best$in_expected,
      out_sample_days     = sum(outs),
      out_sample_breaches = x,
      out_sample_expected = sum(outs) * ALPHA,
      out_sample_rate     = x / sum(outs),
      out_kupiec_lr       = kp$statistic,
      out_kupiec_p_value  = kp$p_value,
      out_reject_at_5pct  = kp$p_value < 0.05,
      n_tied_lambdas      = sum(cand$tied_with_chosen)
    )
    attr(rows[[m]], "candidates") <- cand
  }
  rows
}

plot_lambda_sensitivity <- function(sens, fc, file = "ewma_lambda_sensitivity.png") {
  # The two EWMA curves keep the colours of the EWMA normal / EWMA t models
  # used in the other charts, plus distinct point shapes. The rolling 250-day
  # results are grey dotted reference lines with direct labels, so no colour
  # is reused for a different model.
  col <- c(normal = MODEL_COL[["ewma_normal"]], student_t = MODEL_COL[["ewma_t"]])
  pch <- c(normal = 16, student_t = 17)
  rolling <- c(normal    = sum(fc$breach_normal),
               student_t = sum(fc$breach_student_t),
               historical = sum(fc$breach_historical))
  expected <- nrow(fc) * ALPHA
  hl <- half_life(LAMBDAS)

  open_png(file, width = 2400, height = 1400)
  par(mar = c(7.6, 5, 7.4, 13))
  xlim <- range(hl) * c(0.9, 1.1)
  ylim <- c(0, max(sens$breaches, rolling) * 1.08)
  plot(NA, xlim = xlim, ylim = ylim, log = "x", axes = FALSE, bty = "n",
       xlab = "", ylab = "Total breaches, 2015-2025")
  mtext("EWMA half-life in trading days (log scale)", side = 1, line = 2.7, col = COL$text)
  at <- pretty(ylim)
  abline(h = at, col = COL$grid, lwd = 0.8)
  axis(2, at = at, lwd = 0)
  axis(1, at = c(5, 7, 10, 15, 20, 30, 50, 70), lwd = 0, lwd.ticks = 1, col.ticks = COL$axis)
  # Top axis: the lambda behind each half-life.
  # Labels drawn by hand: axis() silently drops labels it thinks are too close.
  axis(3, at = hl, labels = FALSE, lwd = 0, lwd.ticks = 1, col.ticks = COL$axis, line = -0.2)
  mtext(format(LAMBDAS), side = 3, at = hl, line = 0.5, cex = 0.72, col = COL$text2)
  mtext("lambda", side = 3, line = 1.6, at = xlim[1], adj = 0, cex = 0.8, col = COL$text2)

  # Reference lines with labels at the right edge.
  par(xpd = NA)
  right <- 10^(par("usr")[2])
  ref <- function(y, lty, lcol, label) {
    segments(10^par("usr")[1], y, right, y, lty = lty, lwd = 1.6, col = lcol)
    text(right * 1.03, y, label, adj = 0, cex = 0.78, col = COL$text)
  }
  ref(expected, 2, COL$text, sprintf("Expected at 1%%: %.2f", expected))
  ref(rolling[["normal"]], 3, COL$text2,
      sprintf("Rolling 250-day, normal: %d", rolling[["normal"]]))
  ref(rolling[["student_t"]], 3, COL$text2,
      sprintf("Rolling 250-day, Student-t: %d", rolling[["student_t"]]))
  ref(rolling[["historical"]], 3, COL$text2,
      sprintf("Rolling 250-day, historical: %d", rolling[["historical"]]))
  par(xpd = FALSE)

  for (m in names(col)) {
    s <- sens[sens$multiplier == m, ]
    lines(s$half_life_days, s$breaches, col = col[[m]], lwd = 2.2)
    points(s$half_life_days, s$breaches, pch = pch[[m]], col = col[[m]], cex = 1.3)
  }
  # Ring the default lambda 0.94 used by the main EWMA models.
  d <- sens[sens$lambda == 0.94, ]
  points(d$half_life_days, d$breaches, pch = 1, cex = 2.6, lwd = 1.5, col = COL$text)

  mtext("EWMA VaR breaches against decay speed, 60/40 SPY/AGG portfolio, 2015-2025",
        side = 3, line = 5.6, adj = 0, cex = 1.2, font = 2, col = COL$text)
  mtext(paste("1-day 99% VaR. Each point is one lambda; shorter half-life = faster-reacting volatility.",
              "Ringed points: lambda 0.94, used by the main EWMA models."),
        side = 3, line = 4.2, adj = 0, cex = 0.82, col = COL$text2)
  bottom_legend(c("EWMA, normal multiplier", "EWMA, Student-t (df 5) multiplier",
                  "Rolling 250-day (dotted)", "Expected (dashed)"),
                cex = 0.85, inset = -0.31,
                col = c(col[["normal"]], col[["student_t"]], COL$text2, COL$text),
                lty = c(1, 1, 3, 2), lwd = c(2.2, 2.2, 1.6, 1.6),
                pch = c(pch[["normal"]], pch[["student_t"]], NA, NA), seg.len = 2)
  dev.off()
}

run_ewma_sensitivity <- function(ret, fc) {
  sens <- lambda_sensitivity(ret, fc)

  # Consistency check: lambda 0.94 must reproduce the main EWMA models.
  stopifnot(
    sens$breaches[sens$multiplier == "normal"    & sens$lambda == 0.94] == sum(fc$breach_ewma_normal),
    sens$breaches[sens$multiplier == "student_t" & sens$lambda == 0.94] == sum(fc$breach_ewma_t)
  )

  write.csv(sens, file.path(OUTPUT_DIR, "ewma_lambda_sensitivity.csv"), row.names = FALSE)
  cat("\n===== EWMA lambda sensitivity (2015-2025, n =", nrow(fc), ") =====\n")
  print(transform(sens,
                  half_life_days = round(half_life_days, 2),
                  breach_rate = sprintf("%.3f%%", 100 * breach_rate),
                  kupiec_lr = round(kupiec_lr, 3),
                  kupiec_p_value = signif(kupiec_p_value, 3)), row.names = FALSE)

  oos <- out_of_sample_check(ret, fc)
  oos_tab <- do.call(rbind, oos)
  cand_tab <- do.call(rbind, lapply(oos, function(o) {
    c <- attr(o, "candidates"); c[order(c$lambda), ]
  }))
  write.csv(oos_tab, file.path(OUTPUT_DIR, "ewma_lambda_out_of_sample.csv"), row.names = FALSE)
  write.csv(cand_tab, file.path(OUTPUT_DIR, "ewma_lambda_out_of_sample_candidates.csv"),
            row.names = FALSE)
  cat("\n===== Out-of-sample check: choose lambda on 2015-2019, test on 2020-2025 =====\n")
  for (m in names(oos)) {
    o <- oos[[m]]
    cat(sprintf("\n[%s multiplier] in-sample candidates (2015-2019, %d days, expected %.2f):\n",
                m, o$in_sample_days, o$in_sample_expected))
    cand <- attr(o, "candidates")
    cat("(out_* columns are for reference only and are not used in the choice)\n")
    print(transform(cand[order(cand$lambda), names(cand) != "multiplier"],
                    in_p = signif(in_p, 3), out_p = signif(out_p, 3)), row.names = FALSE)
    cat(sprintf("chosen lambda %.3f (half-life %.1f days): %d breaches in-sample vs %.2f expected\n",
                o$chosen_lambda, o$half_life_days, o$in_sample_breaches, o$in_sample_expected))
    cat(sprintf("2020-2025: %d breaches vs %.2f expected (%d days, rate %.3f%%) | Kupiec LR %.3f, p = %.3g -> %s at 5%%\n",
                o$out_sample_breaches, o$out_sample_expected, o$out_sample_days,
                100 * o$out_sample_rate, o$out_kupiec_lr, o$out_kupiec_p_value,
                if (o$out_reject_at_5pct) "REJECT" else "do not reject"))
    if (o$n_tied_lambdas > 1) {
      tied <- cand[cand$tied_with_chosen, ]
      cat(sprintf("NOTE: %d lambdas tie on the in-sample rule (%s); the tie-break picked %.3f.\n",
                  nrow(tied), paste(sort(tied$lambda), collapse = ", "), o$chosen_lambda))
      cat(sprintf("      2020-2025 results across the tied lambdas: %s\n",
                  paste(sprintf("%.3f -> %d breaches, p = %.3g (%s)", tied$lambda,
                                tied$out_breaches, tied$out_p,
                                ifelse(tied$out_p < 0.05, "reject", "not rejected")),
                        collapse = "; ")))
    }
  }

  plot_lambda_sensitivity(sens, fc)
  invisible(list(sensitivity = sens, out_of_sample = oos_tab))
}
