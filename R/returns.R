# Turn raw prices into daily simple returns for each fund and for the
# 60/40 portfolio. Defines functions; run_all.R calls them.

source(file.path("R", "config.R"))

read_prices <- function(ticker) {
  df <- read.csv(file.path(RAW_DIR, paste0(ticker, ".csv")))
  df$date <- as.Date(df$date)
  df[order(df$date), ]
}

# Simple return on day t: P_t / P_{t-1} - 1, using dividend-adjusted closes,
# so the return includes dividends (a total return).
simple_returns <- function(p) c(NA, p[-1] / p[-length(p)] - 1)

build_returns <- function() {
  spy <- read_prices("SPY")
  agg <- read_prices("AGG")

  # Keep only dates where both funds have a price. Any date one fund is
  # missing is reported by the quality report rather than silently filled.
  px <- merge(spy[, c("date", "adjusted")], agg[, c("date", "adjusted")],
              by = "date", suffixes = c("_spy", "_agg"))
  px <- px[order(px$date), ]

  ret <- data.frame(
    date    = px$date,
    ret_spy = simple_returns(px$adjusted_spy),
    ret_agg = simple_returns(px$adjusted_agg)
  )
  # Daily rebalancing to fixed weights means the portfolio return each day is
  # just the weighted average of the two fund returns.
  ret$ret_port <- WEIGHTS["SPY"] * ret$ret_spy + WEIGHTS["AGG"] * ret$ret_agg

  ret[-1, ]   # first row has no previous price
}
