# Download daily prices for SPY and AGG from Yahoo Finance and save them as
# raw CSVs in data/raw/. These files are never committed to git.
#
# Run from the project root:  Rscript R/download.R

source(file.path("R", "config.R"))
suppressPackageStartupMessages(library(quantmod))

dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)

download_one <- function(ticker) {
  # getSymbols treats `to` as exclusive, so ask for one day past the end.
  x <- getSymbols(ticker, src = "yahoo",
                  from = DOWNLOAD_FROM, to = DOWNLOAD_TO + 1,
                  auto.assign = FALSE)
  out <- data.frame(
    date     = index(x),
    open     = as.numeric(Op(x)),
    high     = as.numeric(Hi(x)),
    low      = as.numeric(Lo(x)),
    close    = as.numeric(Cl(x)),   # split-adjusted only
    volume   = as.numeric(Vo(x)),
    adjusted = as.numeric(Ad(x))    # split- and dividend-adjusted
  )
  path <- file.path(RAW_DIR, paste0(ticker, ".csv"))
  write.csv(out, path, row.names = FALSE)
  message(sprintf("%s: %d rows, %s to %s -> %s",
                  ticker, nrow(out), min(out$date), max(out$date), path))
  invisible(out)
}

for (tk in TICKERS) download_one(tk)

# Yahoo revises its adjusted prices retroactively (every new dividend changes
# the whole history), so record when this snapshot was taken.
writeLines(c(
  paste("downloaded_at_utc:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("source: Yahoo Finance via quantmod", as.character(packageVersion("quantmod"))),
  paste("tickers:", paste(TICKERS, collapse = ", ")),
  paste("requested_range:", DOWNLOAD_FROM, "to", DOWNLOAD_TO)
), file.path(RAW_DIR, "download_info.txt"))
