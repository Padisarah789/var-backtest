# Print everything needed to verify one day's forecast by hand: the window's
# first and last dates, its mean and sd, its five worst returns, and both VaRs.

source(file.path("R", "config.R"))

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
              qnorm(0.99), f$var_normal, 100 * f$var_normal))
  cat(sprintf("realised portfolio return that day: %.6f (%.4f%%)\n",
              f$ret_port, 100 * f$ret_port))
}
