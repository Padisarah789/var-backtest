# Data quality and summary report. Prints to the console and saves a copy to
# output/quality_report.txt.

source(file.path("R", "config.R"))

pct <- function(x, d = 3) sprintf(paste0("%.", d, "f%%"), 100 * x)

report_raw <- function(ticker) {
  df <- read.csv(file.path(RAW_DIR, paste0(ticker, ".csv")))
  df$date <- as.Date(df$date)
  cat(sprintf("\n--- %s raw prices ---\n", ticker))
  cat(sprintf("rows: %d | dates: %s to %s\n", nrow(df), min(df$date), max(df$date)))
  cat("missing values per column:\n")
  print(colSums(is.na(df)))
  cat(sprintf("duplicate dates: %d | non-positive adjusted prices: %d | zero-volume days: %d\n",
              sum(duplicated(df$date)), sum(df$adjusted <= 0, na.rm = TRUE),
              sum(df$volume == 0, na.rm = TRUE)))
  gaps <- as.numeric(diff(sort(df$date)))
  cat(sprintf("largest calendar gap between rows: %d days (ending %s)\n",
              max(gaps), sort(df$date)[which.max(gaps) + 1]))
  # Price above the day's high or below its low suggests a bad print.
  bad_hl <- with(df, close > high + 1e-8 | close < low - 1e-8 | low > high)
  cat(sprintf("rows where close is outside [low, high]: %d\n", sum(bad_hl, na.rm = TRUE)))
  # The adjusted/close ratio should only step at ex-dividend dates (fewer than
  # about 12 per year). Many more steps suggests noisy adjustment data.
  ratio <- df$adjusted / df$close
  steps <- sum(abs(diff(ratio)) > 1e-6, na.rm = TRUE)
  cat(sprintf("days where the adjusted/close ratio changes: %d (about %.1f per year)\n",
              steps, steps / (as.numeric(diff(range(df$date))) / 365.25)))
  invisible(df)
}

report_returns <- function(ret) {
  cat("\n--- Daily returns (after aligning both funds on common dates) ---\n")
  cat(sprintf("rows: %d | dates: %s to %s\n", nrow(ret), min(ret$date), max(ret$date)))
  cat("missing values per column:\n")
  print(colSums(is.na(ret)))
  cat(sprintf("exactly-zero returns: SPY %d, AGG %d (many in a row would suggest stale prices)\n",
              sum(ret$ret_spy == 0), sum(ret$ret_agg == 0)))
  for (col in c("ret_spy", "ret_agg", "ret_port")) {
    o <- order(abs(ret[[col]]), decreasing = TRUE)[1:3]
    cat(sprintf("largest absolute %s: %s\n", col,
                paste(sprintf("%s (%s)", ret$date[o], pct(ret[[col]][o], 2)), collapse = ", ")))
  }
  cat("\nreturn summary (daily):\n")
  print(sapply(ret[, c("ret_spy", "ret_agg", "ret_port")], function(x)
    c(mean = mean(x), sd = sd(x), min = min(x), max = max(x))), digits = 4)
}

report_forecasts <- function(fc, ret) {
  cat("\n--- VaR forecasts ---\n")
  cat(sprintf("forecast days: %d | first: %s | last: %s\n",
              nrow(fc), min(fc$date), max(fc$date)))
  cat("missing values per column:\n")
  print(colSums(is.na(fc)))
  cat("forecast days per year:\n")
  print(table(format(fc$date, "%Y")))

  cat("\nsummary of VaR series (positive = loss, % of portfolio value):\n")
  s <- sapply(fc[, MODELS], function(x)
    100 * c(min = min(x), q25 = unname(quantile(x, 0.25)), median = median(x),
            mean = mean(x), q75 = unname(quantile(x, 0.75)), max = max(x), sd = sd(x)))
  colnames(s) <- names(MODELS)
  print(round(s, 3))
  cat("date of max VaR:\n")
  print(sapply(MODELS, function(v) format(fc$date[which.max(fc[[v]])])))
  cat("\ncorrelation between the VaR series:\n")
  cm <- cor(fc[, MODELS]); dimnames(cm) <- list(names(MODELS), names(MODELS))
  print(round(cm, 3))
  cat(sprintf("EWMA burn-in: seeded %d returns before the first forecast date (seed weight left: %.1e)\n",
              match(fc$date[1], ret$date) - WINDOW - 1,
              EWMA_LAMBDA^(match(fc$date[1], ret$date) - WINDOW - 1)))
  i1 <- match(fc$date, ret$date)
  gap <- abs(sqrt(ewma_variance(ret$ret_port, init_n = FHS_SEED_N)[i1]) / fc$ewma_sd - 1)
  cat(sprintf("FHS standardisation EWMA (seeded on the first %d returns) vs main EWMA at forecast dates: max relative sd gap %.1e\n",
              FHS_SEED_N, max(gap)))
  cat(sprintf("floor VaR comes from EWMA normal on %d days and from historical on %d days\n",
              sum(fc$var_ewma_normal >= fc$var_hist), sum(fc$var_ewma_normal < fc$var_hist)))
  cat(sprintf("days where historical VaR > normal VaR: %d of %d (%s)\n",
              sum(fc$var_hist > fc$var_normal), nrow(fc),
              pct(mean(fc$var_hist > fc$var_normal), 1)))

  # Sensitivity of historical VaR to the quantile definition.
  idx <- match(fc$date, ret$date)
  v7 <- vapply(idx, function(i) var_historical(ret$ret_port[window_rows(i)], type = 7), 0)
  d <- fc$var_hist - v7
  cat(sprintf("historical VaR, type 1 minus type 7: mean %s, max %s (type 1 is never smaller: %s)\n",
              pct(mean(d)), pct(max(d)), all(d >= -1e-12)))
}

run_quality_report <- function(ret, fc) {
  path <- file.path(OUTPUT_DIR, "quality_report.txt")
  lines <- capture.output({
    cat("DATA QUALITY REPORT  (generated", format(Sys.time(), "%Y-%m-%d %H:%M"), ")\n")
    info <- file.path(RAW_DIR, "download_info.txt")
    if (file.exists(info)) cat(readLines(info), sep = "\n")
    for (tk in TICKERS) report_raw(tk)
    report_returns(ret)
    report_forecasts(fc, ret)
  })
  writeLines(lines, path)
  cat(lines, sep = "\n")
  invisible(path)
}
