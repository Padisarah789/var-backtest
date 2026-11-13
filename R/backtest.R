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

MODELS <- c(historical = "var_hist", normal = "var_normal")

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

run_backtest <- function(fc) {
  fc <- add_breach_flags(fc)
  summ <- backtest_summary(fc)
  year <- breaches_by_year(fc)
  big  <- largest_breaches(fc)

  write.csv(summ, file.path(OUTPUT_DIR, "backtest_summary.csv"), row.names = FALSE)
  write.csv(year, file.path(OUTPUT_DIR, "backtest_breaches_by_year.csv"), row.names = FALSE)
  write.csv(big,  file.path(OUTPUT_DIR, "backtest_largest_breaches.csv"), row.names = FALSE)

  cat("\n===== BACKTEST: hand-check numbers =====\n")
  n <- nrow(fc)
  cat(sprintf("n = %d forecast days, p = %.2f\n", n, ALPHA))
  cat(sprintf("binomial null: mean n*p = %.2f, sd sqrt(n*p*(1-p)) = %.4f\n",
              n * ALPHA, sqrt(n * ALPHA * (1 - ALPHA))))
  for (i in seq_len(nrow(summ))) {
    s <- summ[i, ]
    cat(sprintf("%-10s breaches = %d | rate = %.4f%% | z = %.2f | Kupiec LR = %.4f, p = %.4g -> %s at 5%%\n",
                s$model, s$breaches, 100 * s$breach_rate, s$z_score, s$kupiec_lr,
                s$kupiec_p_value, if (s$reject_at_5pct) "REJECT" else "do not reject"))
  }

  cat("\nbreaches per calendar year:\n")
  print(transform(year, expected = round(expected, 2)), row.names = FALSE)

  cat("\n10 largest breaches per model (returns, VaR and excess in %):\n")
  big_print <- transform(big, ret_port = round(100 * ret_port, 3),
                         var = round(100 * var, 3), excess = round(100 * excess, 3))
  print(big_print, row.names = FALSE)

  invisible(list(fc = fc, summary = summ, by_year = year, largest = big))
}
