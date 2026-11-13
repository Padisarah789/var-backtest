# One-day 99% Value at Risk forecasts from a rolling 250-day window.
#
# Sign convention: VaR is a positive number giving a loss as a fraction of
# portfolio value. VaR = 0.02 means "on 99% of days we expect to lose less
# than 2%". A day breaches VaR when its return is below -VaR.
#
# No look-ahead: the forecast for day t uses returns from days t-250 to t-1
# only. Day t's own return is stored next to the forecast for the later
# backtest but is never used to make it.

source(file.path("R", "config.R"))

# Historical simulation VaR.
# Quantile definition: R's type = 1, the inverse of the empirical CDF:
# the smallest observed return r with (share of returns <= r) >= 1%.
# With 250 returns, 250 * 0.01 = 2.5, so it rounds up to the 3rd worst return.
# Why type 1 rather than R's default type 7:
#   * the VaR is a loss that actually happened, with no interpolation
#     between two days' returns;
#   * it is easy to check by hand: sort the window, take the 3rd worst;
#   * it follows the textbook definition VaR = -inf{x : F(x) >= alpha}.
# Type 7 would interpolate between the 3rd and 4th worst (position 3.49),
# giving a slightly smaller VaR. quality_report.R shows the size of that gap.
var_historical <- function(window_returns, alpha = ALPHA, type = 1) {
  -unname(quantile(window_returns, probs = alpha, type = type))
}

# Parametric (variance-covariance) VaR, assuming normal returns:
#   VaR = z * sd - mean,  z = qnorm(0.99) = 2.3263...
# sd is the sample standard deviation (divides by n - 1).
# We use the exact qnorm(0.99) rather than the rounded 2.326; they differ by
# 0.00035 * sd, about 0.015% of the VaR.
var_normal <- function(window_returns, alpha = ALPHA) {
  z <- qnorm(1 - alpha)
  z * sd(window_returns) - mean(window_returns)
}

# Rows of `ret` (returned by build_returns()) used for the forecast at row i.
window_rows <- function(i, window = WINDOW) (i - window):(i - 1)

forecast_var <- function(ret, window = WINDOW,
                         from = FORECAST_FROM, to = FORECAST_TO) {
  stopifnot(!is.unsorted(ret$date), !anyNA(ret$ret_port))

  idx <- which(ret$date >= from & ret$date <= to)
  if (min(idx) <= window) {
    stop("Not enough history: the first forecast date needs ", window,
         " prior returns but only has ", min(idx) - 1)
  }

  out <- data.frame(
    date         = ret$date[idx],
    ret_port     = ret$ret_port[idx],
    var_hist     = NA_real_,
    var_normal   = NA_real_,
    window_start = as.Date(NA),
    window_end   = as.Date(NA)
  )
  for (k in seq_along(idx)) {
    rows <- window_rows(idx[k], window)
    w <- ret$ret_port[rows]
    out$var_hist[k]     <- var_historical(w)
    out$var_normal[k]   <- var_normal(w)
    out$window_start[k] <- ret$date[rows[1]]
    out$window_end[k]   <- ret$date[rows[length(rows)]]
  }
  # Guard against look-ahead: every window must end before its forecast date.
  stopifnot(all(out$window_end < out$date))
  out
}
