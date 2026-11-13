# One-day 99% Value at Risk forecasts: five models.
#
# Sign convention: VaR is a positive number giving a loss as a fraction of
# portfolio value. VaR = 0.02 means "on 99% of days we expect to lose less
# than 2%". A day breaches VaR when its return is below -VaR.
#
# No look-ahead: the forecast for day t uses returns up to day t-1 only.
# Day t's own return is stored next to the forecast for the backtest but is
# never used to make it.
#
# Models (name -> column in the forecast table):
#   historical   var_hist         3rd worst of the last 250 returns
#   normal       var_normal       z_normal * rolling sd - rolling mean
#   student_t    var_t            z_t      * rolling sd - rolling mean
#   ewma_normal  var_ewma_normal  z_normal * EWMA sd   (mean taken as zero)
#   ewma_t       var_ewma_t       z_t      * EWMA sd   (mean taken as zero)

source(file.path("R", "config.R"))

MODELS <- c(historical  = "var_hist",
            normal      = "var_normal",
            student_t   = "var_t",
            ewma_normal = "var_ewma_normal",
            ewma_t      = "var_ewma_t")

T_DF        <- 5      # degrees of freedom for the Student-t models
EWMA_LAMBDA <- 0.94   # RiskMetrics daily decay factor

# ---- Multipliers ---------------------------------------------------------

z_normal <- function(alpha = ALPHA) qnorm(1 - alpha)          # 2.326348

# A t distribution with df degrees of freedom has variance df / (df - 2).
# Multiplying by sqrt((df - 2) / df) rescales it to variance 1, so that
# "multiplier * sd" uses the same sd as the normal model and only the shape
# of the tail differs. For df = 5: -qt(0.01, 5) * sqrt(3/5) = 2.606464.
z_student_t <- function(df = T_DF, alpha = ALPHA) {
  stopifnot(df > 2)
  -qt(alpha, df) * sqrt((df - 2) / df)
}

# ---- Historical and rolling-moment models --------------------------------

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

# Parametric VaR from a mean and standard deviation: multiplier * sd - mean.
# sd is the sample standard deviation (divides by n - 1).
# The normal model uses qnorm(0.99) exactly rather than the rounded 2.326;
# they differ by 0.00035 * sd, about 0.015% of the VaR.
var_parametric <- function(mu, sigma, multiplier) multiplier * sigma - mu

var_normal <- function(window_returns, alpha = ALPHA) {
  var_parametric(mean(window_returns), sd(window_returns), z_normal(alpha))
}

# Rows of `ret` (returned by build_returns()) used for the forecast at row i.
window_rows <- function(i, window = WINDOW) (i - window):(i - 1)

# ---- EWMA volatility -----------------------------------------------------

# Exponentially weighted variance forecast, RiskMetrics style:
#   v[i] = lambda * v[i-1] + (1 - lambda) * r[i-1]^2
# v[i] is the variance forecast for row i and uses returns up to row i-1 only.
# Seed: v[init_n + 1] = sample variance (n - 1 divisor) of r[1..init_n], i.e.
# the forecast for the first day after the first 250 returns is their sample
# variance. Rows 1..init_n have no forecast (NA).
# The mean is taken as zero in the recursion (r^2, not (r - mean)^2).
ewma_variance <- function(r, lambda = EWMA_LAMBDA, init_n = WINDOW) {
  stopifnot(length(r) > init_n, !anyNA(r))
  v <- rep(NA_real_, length(r))
  v[init_n + 1] <- var(r[1:init_n])
  if (length(r) > init_n + 1) {
    for (i in (init_n + 2):length(r)) {
      v[i] <- lambda * v[i - 1] + (1 - lambda) * r[i - 1]^2
    }
  }
  v
}

# ---- Forecast table ------------------------------------------------------

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
    roll_mean    = NA_real_,
    roll_sd      = NA_real_,
    window_start = as.Date(NA),
    window_end   = as.Date(NA)
  )
  for (k in seq_along(idx)) {
    rows <- window_rows(idx[k], window)
    w <- ret$ret_port[rows]
    out$var_hist[k]     <- var_historical(w)
    out$roll_mean[k]    <- mean(w)
    out$roll_sd[k]      <- sd(w)
    out$window_start[k] <- ret$date[rows[1]]
    out$window_end[k]   <- ret$date[rows[length(rows)]]
  }

  # EWMA runs over the whole return history (seeded on its first 250 returns,
  # which end well before 2015) and is then read off at the forecast dates.
  ewma_v <- ewma_variance(ret$ret_port)
  out$ewma_var <- ewma_v[idx]
  out$ewma_sd  <- sqrt(out$ewma_var)

  out$var_normal      <- var_parametric(out$roll_mean, out$roll_sd, z_normal())
  out$var_t           <- var_parametric(out$roll_mean, out$roll_sd, z_student_t())
  out$var_ewma_normal <- var_parametric(0, out$ewma_sd, z_normal())
  out$var_ewma_t      <- var_parametric(0, out$ewma_sd, z_student_t())

  # Guards against look-ahead and gaps.
  stopifnot(all(out$window_end < out$date), !anyNA(out[, MODELS]))
  out
}
