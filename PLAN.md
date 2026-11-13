# Plan

## Done (core)

- Project setup: folders, `.gitignore`, package versions locked with renv.
- `R/download.R`: daily SPY and AGG prices from Yahoo, 2013-06-03 to 2025-12-31.
- `R/returns.R`: daily simple total returns; 60/40 portfolio rebalanced daily.
- `R/var_models.R`: 1-day 99% VaR from a rolling 250-day window, for
  2,766 trading days from 2015-01-02 to 2025-12-31.
  - Historical: empirical 1% quantile, R `type = 1` (3rd worst of 250).
  - Parametric normal: `qnorm(0.99) * sd - mean`.
- `R/quality_report.R` and `R/hand_check.R`: data checks and a printout for
  checking individual dates by hand.

## Next

1. **Backtest.** Mark each day as an exception when `ret_port < -VaR`. Count
   exceptions per model, in total and per year. With 2,766 days, about 27.7
   are expected at 99%.
2. **Kupiec test (unconditional coverage).** Likelihood-ratio test that the
   exception rate equals 1%: `LR_uc ~ chi-squared(1)`. Report the statistic
   and p-value per model, over the full sample and per year.
   - Possible extension: Christoffersen tests for independence and
     conditional coverage, to see whether exceptions cluster (as they would
     in March 2020).
   - Possible extension: Basel traffic-light zones on rolling 250-day
     exception counts.
3. **Plot.** Realised portfolio returns against `-VaR` for both models, with
   exception days marked. Add a zoomed panel for February to April 2020.

## Later / open questions

- Unit tests (`testthat`) for the VaR functions on small known inputs, and a
  test that no window includes its own forecast date.
- Decide whether to commit a frozen copy of `output/var_forecasts.csv`.
  Yahoo revises adjusted prices, so results can drift slightly between
  downloads.
- Check Yahoo's adjusted prices against a second source (e.g. the issuers'
  published NAV and distribution history).
