# Printouts for checking individual forecasts by hand.

source(file.path("R", "config.R"))

# Rolling-window check: the window's first and last dates, its mean and sd,
# its five worst returns, and the historical and normal VaRs.
hand_check <- function(ret, fc, date) {
  date <- as.Date(date)
  i <- match(date, ret$date)
  if (is.na(i)) stop(date, " is not a trading day in the data")
  w_rows <- window_rows(i)
  w <- ret$ret_port[w_rows]
  f <- fc[fc$date == date, ]
  worst <- order(w)[1:5]

  cat(sprintf("\n===== Hand check: forecast for %s =====\n", date))
  cat(sprintf("window: %d returns, %s to %s\n",
              length(w), ret$date[w_rows[1]], ret$date[w_rows[length(w_rows)]]))
  cat(sprintf("window mean: %.6f (%.4f%%)\n", mean(w), 100 * mean(w)))
  cat(sprintf("window sd:   %.6f (%.4f%%)\n", sd(w), 100 * sd(w)))
  cat("five worst portfolio returns in the window:\n")
  print(data.frame(rank = 1:5, date = ret$date[w_rows][worst],
                   ret_port = sprintf("%.6f", w[worst]),
                   ret_spy = sprintf("%.6f", ret$ret_spy[w_rows][worst]),
                   ret_agg = sprintf("%.6f", ret$ret_agg[w_rows][worst])),
        row.names = FALSE)
  cat(sprintf("historical VaR (3rd worst, sign flipped): %.6f (%.4f%%)\n",
              f$var_hist, 100 * f$var_hist))
  cat(sprintf("normal VaR (%.6f * sd - mean):          %.6f (%.4f%%)\n",
              z_normal(), f$var_normal, 100 * f$var_normal))
  cat(sprintf("realised portfolio return that day: %.6f (%.4f%%)\n",
              f$ret_port, 100 * f$ret_port))
}

# EWMA check: yesterday's EWMA variance and return, today's EWMA variance
# (variance_t = lambda * variance_{t-1} + (1 - lambda) * return_{t-1}^2),
# and all five VaR forecasts with the inputs needed to rebuild them.
hand_check_ewma <- function(ret, fc, date) {
  date <- as.Date(date)
  i <- match(date, ret$date)
  if (is.na(i)) stop(date, " is not a trading day in the data")
  v <- ewma_variance(ret$ret_port)
  f <- fc[fc$date == date, ]
  stopifnot(isTRUE(all.equal(v[i], f$ewma_var)))

  cat(sprintf("\n===== EWMA hand check: forecast for %s =====\n", date))
  cat(sprintf("yesterday (%s) EWMA variance: %.10f  (sd %.6f)\n",
              ret$date[i - 1], v[i - 1], sqrt(v[i - 1])))
  cat(sprintf("yesterday (%s) return:        %.6f\n", ret$date[i - 1], ret$ret_port[i - 1]))
  cat(sprintf("today EWMA variance = %.2f * %.10f + %.2f * (%.6f)^2 = %.10f  (sd %.6f)\n",
              EWMA_LAMBDA, v[i - 1], 1 - EWMA_LAMBDA, ret$ret_port[i - 1], v[i], sqrt(v[i])))
  cat(sprintf("multipliers: normal %.6f, Student-t (df %d, unit variance) %.6f\n",
              z_normal(), T_DF, z_student_t()))
  cat(sprintf("rolling 250-day window: mean %.6f, sd %.6f\n", f$roll_mean, f$roll_sd))
  cat("VaR forecasts:\n")
  cat(sprintf("  historical  %.6f  (%.4f%%)  3rd worst of window\n", f$var_hist, 100 * f$var_hist))
  cat(sprintf("  normal      %.6f  (%.4f%%)  %.6f * rolling sd - rolling mean\n",
              f$var_normal, 100 * f$var_normal, z_normal()))
  cat(sprintf("  student_t   %.6f  (%.4f%%)  %.6f * rolling sd - rolling mean\n",
              f$var_t, 100 * f$var_t, z_student_t()))
  cat(sprintf("  ewma_normal %.6f  (%.4f%%)  %.6f * EWMA sd\n",
              f$var_ewma_normal, 100 * f$var_ewma_normal, z_normal()))
  cat(sprintf("  ewma_t      %.6f  (%.4f%%)  %.6f * EWMA sd\n",
              f$var_ewma_t, 100 * f$var_ewma_t, z_student_t()))
  cat(sprintf("realised portfolio return that day: %.6f (%.4f%%)\n",
              f$ret_port, 100 * f$ret_port))
}
