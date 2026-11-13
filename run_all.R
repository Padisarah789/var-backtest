# Run the whole pipeline from the project root:
#   Rscript run_all.R            # downloads data only if it is not there yet
#   REFRESH=1 Rscript run_all.R  # force a fresh download from Yahoo

source(file.path("R", "config.R"))
source(file.path("R", "returns.R"))
source(file.path("R", "var_models.R"))
source(file.path("R", "quality_report.R"))
source(file.path("R", "hand_check.R"))
source(file.path("R", "backtest.R"))
source(file.path("R", "charts.R"))

raw_files <- file.path(RAW_DIR, paste0(TICKERS, ".csv"))
if (Sys.getenv("REFRESH") == "1" || !all(file.exists(raw_files))) {
  source(file.path("R", "download.R"))
}

dir.create(OUTPUT_DIR, showWarnings = FALSE)

ret <- build_returns()
fc  <- forecast_var(ret)

write.csv(fc[, c("date", "ret_port", MODELS)],
          file.path(OUTPUT_DIR, "var_forecasts.csv"), row.names = FALSE)

run_quality_report(ret, fc)

for (d in c("2019-06-03", "2020-03-16")) hand_check(ret, fc, d)
for (d in c("2020-03-12", "2020-03-16")) hand_check_ewma(ret, fc, d)

bt <- run_backtest(fc)
run_charts(bt$fc, bt$by_year)
