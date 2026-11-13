# Shared settings for the whole project. Every other script sources this file,
# so a change here (e.g. a different confidence level) flows through everywhere.

TICKERS <- c("SPY", "AGG")
WEIGHTS <- c(SPY = 0.60, AGG = 0.40)   # rebalanced back to these weights daily

# Download start: about 18 months before the first forecast, so the 250-day
# window for 2015-01-02 is full with plenty of margin.
DOWNLOAD_FROM <- as.Date("2013-06-01")
DOWNLOAD_TO   <- as.Date("2025-12-31")  # inclusive

FORECAST_FROM <- as.Date("2015-01-01")
FORECAST_TO   <- as.Date("2025-12-31")

WINDOW <- 250     # trading days in the estimation window
ALPHA  <- 0.01    # tail probability: 1-day 99% VaR

RAW_DIR    <- file.path("data", "raw")
OUTPUT_DIR <- file.path("output")
