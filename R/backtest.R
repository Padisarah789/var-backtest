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

# a * log(b), with 0 * log(0) taken as 0 (the limit as a probability -> 0).
xlogy <- function(a, b) if (a == 0) 0 else a * log(b)

# Kupiec (1995) proportion-of-failures likelihood-ratio test.
#   LR = -2 * [ log L(p) - log L(x/n) ]
#   log L(q) = x * log(q) + (n - x) * log(1 - q)
# LR is asymptotically chi-squared with 1 degree of freedom.
# 0 * log(0) is taken as 0, so x = 0 or x = n are handled.
kupiec_pof <- function(x, n, p = ALPHA) {
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

# ---- Clustering of breaches ---------------------------------------------

# Transition counts of the daily breach indicator I_t over consecutive
# forecast days t-1 -> t (n - 1 pairs for n days):
#   n00: no breach -> no breach    n01: no breach -> breach
#   n10: breach -> no breach       n11: breach -> breach
transition_counts <- function(b) {
  prev <- head(b, -1)
  cur  <- tail(b, -1)
  c(n00 = sum(!prev & !cur), n01 = sum(!prev & cur),
    n10 = sum(prev & !cur),  n11 = sum(prev & cur))
}

# Christoffersen (1998) independence test (first-order Markov alternative).
#   pi01 = n01 / (n00 + n01)   P(breach | no breach the day before)
#   pi11 = n11 / (n10 + n11)   P(breach | breach the day before)
#   pi   = (n01 + n11) / (n00 + n01 + n10 + n11)
#   LR_ind = -2 [ log L(pi) - log L(pi01, pi11) ],  chi-squared, 1 df
#   log L(pi)        = (n00 + n10) log(1 - pi) + (n01 + n11) log(pi)
#   log L(pi01,pi11) = n00 log(1 - pi01) + n01 log(pi01)
#                    + n10 log(1 - pi11) + n11 log(pi11)
# n11 = 0: pi11 is estimated as 0 and the n11 * log(pi11) term is taken as
# 0 * log(0) = 0 (its limit), so log L(pi01, pi11) stays finite and the test
# is computed as usual, with no adjustment to the counts. The statistic is
# then driven by pi01 vs pi only. If there is no day after a breach at all
# (n10 + n11 = 0), pi11 is undefined; its terms are both 0 and drop out.
christoffersen_ind <- function(n00, n01, n10, n11) {
  pi01 <- if (n00 + n01 > 0) n01 / (n00 + n01) else NA_real_
  pi11 <- if (n10 + n11 > 0) n11 / (n10 + n11) else NA_real_
  pi   <- (n01 + n11) / (n00 + n01 + n10 + n11)
  ll0 <- xlogy(n00 + n10, 1 - pi) + xlogy(n01 + n11, pi)
  ll1 <- xlogy(n00, 1 - pi01) + xlogy(n01, pi01) +
         xlogy(n10, 1 - pi11) + xlogy(n11, pi11)
  lr <- max(0, -2 * (ll0 - ll1))   # clip tiny negative rounding error
  list(statistic = lr, p_value = pchisq(lr, df = 1, lower.tail = FALSE),
       pi01 = pi01, pi11 = pi11, pi = pi)
}

# Breaches that came within `window` trading days after the previous breach
# (gap of 1 to `window` forecast days). The share's denominator is all
# breaches, including each model's first breach, which can never count.
clustered_breaches <- function(b, window = 5) {
  pos <- which(b)
  sum(diff(pos) <= window)
}

# Per model: transition counts, conditional breach probabilities, the
# independence test, conditional coverage, and the 5-day clustering share.
# Conditional coverage = Kupiec LR (all n days) + independence LR,
# chi-squared with 2 df. (Christoffersen's original LR_cc computes the
# unconditional part on the n - 1 transitions instead; the two differ only
# through the first day.)
clustering_tests <- function(fc, p = ALPHA, window = 5) {
  n <- nrow(fc)
  rows <- lapply(names(MODELS), function(m) {
    b  <- fc[[paste0("breach_", m)]]
    tc <- transition_counts(b)
    ci <- christoffersen_ind(tc[["n00"]], tc[["n01"]], tc[["n10"]], tc[["n11"]])
    kp <- kupiec_pof(sum(b), n, p)
    lr_cc <- kp$statistic + ci$statistic
    k <- clustered_breaches(b, window)
    data.frame(
      model             = m,
      breaches          = sum(b),
      n00 = tc[["n00"]], n01 = tc[["n01"]], n10 = tc[["n10"]], n11 = tc[["n11"]],
      p_breach_after_no_breach = ci$pi01,
      p_breach_after_breach    = ci$pi11,
      p_breach_pooled          = ci$pi,
      independence_lr   = ci$statistic,
      independence_p    = ci$p_value,
      kupiec_lr         = kp$statistic,
      cond_coverage_lr  = lr_cc,
      cond_coverage_p   = pchisq(lr_cc, df = 2, lower.tail = FALSE),
      clustered_breaches = k,
      share_clustered    = if (sum(b) > 0) k / sum(b) else NA_real_
    )
  })
  do.call(rbind, rows)
}

# One row per model, pulling together the backtest tables.
final_summary <- function(summ, clus, sev) {
  stopifnot(identical(summ$model, clus$model), identical(summ$model, sev$model))
  data.frame(
    model           = summ$model,
    breaches        = summ$breaches,
    kupiec_p        = summ$kupiec_p_value,
    independence_p  = clus$independence_p,
    cond_coverage_p = clus$cond_coverage_p,
    share_clustered = clus$share_clustered,
    avg_var_pct     = sev$avg_var_pct,
    mean_excess_pp  = sev$mean_excess_pp,
    total_excess_pp = sev$total_excess_pp
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
  clus <- clustering_tests(fc)
  fin  <- final_summary(summ, clus, sev)

  write.csv(summ, file.path(OUTPUT_DIR, "backtest_summary.csv"), row.names = FALSE)
  write.csv(year, file.path(OUTPUT_DIR, "backtest_breaches_by_year.csv"), row.names = FALSE)
  write.csv(big,  file.path(OUTPUT_DIR, "backtest_largest_breaches.csv"), row.names = FALSE)
  write.csv(sens, file.path(OUTPUT_DIR, "backtest_t_df_sensitivity.csv"), row.names = FALSE)
  write.csv(sev,  file.path(OUTPUT_DIR, "backtest_severity.csv"), row.names = FALSE)
  write.csv(fhsm, file.path(OUTPUT_DIR, "backtest_fhs_multiplier.csv"), row.names = FALSE)
  write.csv(clus, file.path(OUTPUT_DIR, "backtest_clustering.csv"), row.names = FALSE)
  write.csv(fin,  file.path(OUTPUT_DIR, "final_summary.csv"), row.names = FALSE)

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

  cat("\nbreach clustering: transition counts over consecutive forecast days (n - 1 = %d pairs)\n" |>
        sprintf(n - 1))
  cat("  p(breach | no breach yesterday) = n01 / (n00 + n01); p(breach | breach yesterday) = n11 / (n10 + n11)\n")
  print(data.frame(model = clus$model, n00 = clus$n00, n01 = clus$n01, n10 = clus$n10, n11 = clus$n11,
                   n00_n01 = clus$n00 + clus$n01, n10_n11 = clus$n10 + clus$n11,
                   p_after_no_breach = sprintf("%.6f", clus$p_breach_after_no_breach),
                   p_after_breach = sprintf("%.6f", clus$p_breach_after_breach),
                   ind_lr = round(clus$independence_lr, 4),
                   ind_p = signif(clus$independence_p, 4),
                   cc_lr = round(clus$cond_coverage_lr, 4),
                   cc_p = signif(clus$cond_coverage_p, 4)), row.names = FALSE)
  cat("  models with n11 = 0 use 0 * log(0) = 0 in the likelihood (no change to the counts):",
      if (any(clus$n11 == 0)) paste(clus$model[clus$n11 == 0], collapse = ", ") else "none", "\n")
  cat(sprintf("\nbreaches within 5 trading days after a previous breach (reference: 1 - 0.99^5 = %.4f if breaches are independent at 1%%):\n",
              1 - (1 - ALPHA)^5))
  print(data.frame(model = clus$model, breaches = clus$breaches,
                   clustered = clus$clustered_breaches,
                   share = sprintf("%.4f", clus$share_clustered)), row.names = FALSE)

  cat("\n===== FINAL SUMMARY (output/final_summary.csv) =====\n")
  print(transform(fin, kupiec_p = signif(kupiec_p, 4), independence_p = signif(independence_p, 4),
                  cond_coverage_p = signif(cond_coverage_p, 4),
                  share_clustered = round(share_clustered, 4), avg_var_pct = round(avg_var_pct, 4),
                  mean_excess_pp = round(mean_excess_pp, 4), total_excess_pp = round(total_excess_pp, 4)),
        row.names = FALSE)

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
                 t_sensitivity = sens, severity = sev, fhs_multiplier = fhsm,
                 clustering = clus, final = fin))
}
