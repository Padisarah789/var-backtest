# Backtest the VaR forecasts.
#
# A breach (also called an exception) is a day where the realised loss is
# larger than the VaR forecast for that day:
#     loss_t = -ret_port_t  >  VaR_t
# A loss exactly equal to VaR is not a breach.
#
# Under the null hypothesis that the model is correct at 99%, each day is a
# breach with probability p = 0.01, independently, so the breach count over n
# days follows Binomial(n, p): mean n*p, standard deviation sqrt(n*p*(1-p)).

source(file.path("R", "config.R"))
# MODELS (model name -> VaR column) is defined in R/var_models.R.

# Kupiec (1995) proportion-of-failures likelihood-ratio test.
#   LR = -2 * [ log L(p) - log L(x/n) ]
#   log L(q) = x * log(q) + (n - x) * log(1 - q)
# LR is asymptotically chi-squared with 1 degree of freedom.
# 0 * log(0) is taken as 0, so x = 0 or x = n are handled.
kupiec_pof <- function(x, n, p = ALPHA) {
  xlogy <- function(a, b) if (a == 0) 0 else a * log(b)
  loglik <- function(q) xlogy(x, q) + xlogy(n - x, 1 - q)
  lr <- -2 * (loglik(p) - loglik(x / n))
  list(statistic = lr, p_value = pchisq(lr, df = 1, lower.tail = FALSE))
}

add_breach_flags <- function(fc) {
  for (m in names(MODELS)) {
    fc[[paste0("breach_", m)]] <- -fc$ret_port > fc[[MODELS[[m]]]]
  }
  fc
}

backtest_summary <- function(fc, p = ALPHA, level = 0.05) {
  n <- nrow(fc)
  rows <- lapply(names(MODELS), function(m) {
    x <- sum(fc[[paste0("breach_", m)]])
    k <- kupiec_pof(x, n, p)
    data.frame(
      model             = m,
      n_days            = n,
      breaches          = x,
      breach_rate       = x / n,
      expected_breaches = n * p,
      binom_sd          = sqrt(n * p * (1 - p)),
      z_score           = (x - n * p) / sqrt(n * p * (1 - p)),
      kupiec_lr         = k$statistic,
      kupiec_p_value    = k$p_value,
      reject_at_5pct    = k$p_value < level
    )
  })
  do.call(rbind, rows)
}

breaches_by_year <- function(fc, p = ALPHA) {
  yr <- format(fc$date, "%Y")
  out <- data.frame(year = sort(unique(yr)))
  out$n_days   <- as.integer(table(yr)[out$year])
  out$expected <- out$n_days * p
  for (m in names(MODELS)) {
    out[[paste0(m, "_breaches")]] <-
      as.integer(tapply(fc[[paste0("breach_", m)]], yr, sum)[out$year])
  }
  out
}

# Excess = realised loss minus VaR: how far past the forecast the loss went.
largest_breaches <- function(fc, k = 10) {
  rows <- lapply(names(MODELS), function(m) {
    b <- fc[fc[[paste0("breach_", m)]], ]
    b$excess <- -b$ret_port - b[[MODELS[[m]]]]
    b <- head(b[order(b$excess, decreasing = TRUE), ], k)
    data.frame(model = m, rank = seq_len(nrow(b)), date = b$date,
               ret_port = b$ret_port, var = b[[MODELS[[m]]]], excess = b$excess)
  })
  do.call(rbind, rows)
}

# Severity of breaches: how far past VaR the losses went, per model.
# Excess = loss - VaR on breach days, in percentage points of portfolio
# value (a loss of 3.0% against a VaR of 2.5% is an excess of 0.5).
# Average VaR is over all forecast days, also in % of portfolio value.
breach_severity <- function(fc) {
  rows <- lapply(names(MODELS), function(m) {
    v <- fc[[MODELS[[m]]]]
    b <- fc[[paste0("breach_", m)]]
    ex <- 100 * (-fc$ret_port[b] - v[b])
    has <- length(ex) > 0
    data.frame(
      model               = m,
      avg_var_pct         = 100 * mean(v),
      breaches            = sum(b),
      mean_excess_pp      = if (has) mean(ex) else NA_real_,
      median_excess_pp    = if (has) median(ex) else NA_real_,
      largest_excess_pp   = if (has) max(ex) else NA_real_,
      total_excess_pp     = sum(ex)
    )
  })
  do.call(rbind, rows)
}

# FHS data-implied multiplier (minus the 3rd worst standardised return in
# each window) over the forecast period, next to the fixed multipliers.
fhs_multiplier_summary <- function(fc) {
  k <- fc$fhs_mult
  data.frame(
    statistic = c("fhs_min", "fhs_median", "fhs_max",
                  "normal_fixed", "student_t_df5_fixed"),
    multiplier = c(min(k), median(k), max(k), z_normal(), z_student_t()),
    date = c(format(fc$date[which.min(k)]), NA, format(fc$date[which.max(k)]), NA, NA),
    share_of_days_above = c(NA, NA, NA, mean(k > z_normal()), mean(k > z_student_t()))
  )
}

# Sensitivity of the two Student-t models to the degrees of freedom.
# Re-uses the stored rolling mean/sd and EWMA sd, swapping only the multiplier.
t_df_sensitivity <- function(fc, dfs = c(4, T_DF, 8), p = ALPHA) {
  n <- nrow(fc)
  rows <- list()
  for (df in dfs) {
    k <- z_student_t(df)
    v <- list(student_t = var_parametric(fc$roll_mean, fc$roll_sd, k),
              ewma_t    = var_parametric(0, fc$ewma_sd, k))
    for (m in names(v)) {
      x  <- sum(-fc$ret_port > v[[m]])
      kp <- kupiec_pof(x, n, p)
      rows[[length(rows) + 1]] <- data.frame(
        model = m, df = df, multiplier = k, breaches = x, breach_rate = x / n,
        expected_breaches = n * p, kupiec_lr = kp$statistic,
        kupiec_p_value = kp$p_value, reject_at_5pct = kp$p_value < 0.05)
    }
  }
  out <- do.call(rbind, rows)
  out[order(out$model, out$df), ]
}

run_backtest <- function(fc) {
  fc <- add_breach_flags(fc)
  summ <- backtest_summary(fc)
  year <- breaches_by_year(fc)
  big  <- largest_breaches(fc)
  sens <- t_df_sensitivity(fc)
  sev  <- breach_severity(fc)
  fhsm <- fhs_multiplier_summary(fc)

  write.csv(summ, file.path(OUTPUT_DIR, "backtest_summary.csv"), row.names = FALSE)
  write.csv(year, file.path(OUTPUT_DIR, "backtest_breaches_by_year.csv"), row.names = FALSE)
  write.csv(big,  file.path(OUTPUT_DIR, "backtest_largest_breaches.csv"), row.names = FALSE)
  write.csv(sens, file.path(OUTPUT_DIR, "backtest_t_df_sensitivity.csv"), row.names = FALSE)
  write.csv(sev,  file.path(OUTPUT_DIR, "backtest_severity.csv"), row.names = FALSE)
  write.csv(fhsm, file.path(OUTPUT_DIR, "backtest_fhs_multiplier.csv"), row.names = FALSE)

  op <- options(width = 150)   # keep the seven-model tables on one line
  on.exit(options(op))

  cat("\n===== BACKTEST: hand-check numbers =====\n")
  n <- nrow(fc)
  cat(sprintf("n = %d forecast days, p = %.2f\n", n, ALPHA))
  cat(sprintf("binomial null: mean n*p = %.2f, sd sqrt(n*p*(1-p)) = %.4f\n",
              n * ALPHA, sqrt(n * ALPHA * (1 - ALPHA))))
  for (i in seq_len(nrow(summ))) {
    s <- summ[i, ]
    cat(sprintf("%-11s breaches = %d | rate = %.4f%% | z = %.2f | Kupiec LR = %.4f, p = %.4g -> %s at 5%%\n",
                s$model, s$breaches, 100 * s$breach_rate, s$z_score, s$kupiec_lr,
                s$kupiec_p_value, if (s$reject_at_5pct) "REJECT" else "do not reject"))
  }

  cat("\nbreaches per calendar year (all seven models):\n")
  yp <- year
  names(yp) <- sub("_breaches$", "", names(yp))
  yp$expected <- round(yp$expected, 2)
  tot <- yp[1, ]
  tot$year <- "total"
  tot[, -1] <- lapply(yp[, -1], sum)
  print(rbind(yp, tot), row.names = FALSE)

  cat("\nbreach severity (VaR in % of portfolio, excess = loss - VaR in percentage points):\n")
  print(transform(sev, avg_var_pct = round(avg_var_pct, 4),
                  mean_excess_pp = round(mean_excess_pp, 4),
                  median_excess_pp = round(median_excess_pp, 4),
                  largest_excess_pp = round(largest_excess_pp, 4),
                  total_excess_pp = round(total_excess_pp, 4)), row.names = FALSE)

  cat("\nFHS data-implied multiplier over the forecast period:\n")
  print(transform(fhsm, multiplier = round(multiplier, 6),
                  share_of_days_above = round(share_of_days_above, 4)), row.names = FALSE)

  cat("\n10 largest breaches per model (returns, VaR and excess in %):\n")
  big_print <- transform(big, ret_port = round(100 * ret_port, 3),
                         var = round(100 * var, 3), excess = round(100 * excess, 3))
  print(big_print, row.names = FALSE)

  cat("\nSensitivity: Student-t models with different degrees of freedom\n")
  print(transform(sens, multiplier = round(multiplier, 6),
                  breach_rate = sprintf("%.3f%%", 100 * breach_rate),
                  kupiec_lr = round(kupiec_lr, 3),
                  kupiec_p_value = signif(kupiec_p_value, 3)), row.names = FALSE)

  invisible(list(fc = fc, summary = summ, by_year = year, largest = big,
                 t_sensitivity = sens, severity = sev, fhs_multiplier = fhsm))
}
