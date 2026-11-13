# Plan

## Done

### Block 1: core
- Project setup: folders, `.gitignore`, package versions locked with renv.
- `R/download.R`: daily SPY and AGG prices from Yahoo, 2013-06-03 to 2025-12-31.
- `R/returns.R`: daily simple total returns; 60/40 portfolio rebalanced daily.
- `R/var_models.R`: 1-day 99% VaR from a rolling 250-day window, for
  2,766 trading days from 2015-01-02 to 2025-12-31.
  - Historical: empirical 1% quantile, R `type = 1` (3rd worst of 250).
  - Parametric normal: `qnorm(0.99) * sd - mean`.
- `R/quality_report.R` and `R/hand_check.R`: data checks and a printout for
  checking individual dates by hand.

### Block 3: backtest
- `R/backtest.R`: a breach is a day with `-ret_port > VaR` (strictly
  greater). Per model: total breaches, breach rate, expected count, binomial
  standard deviation, and the Kupiec proportion-of-failures test (likelihood
  ratio, p-value, reject at 5%). Also breaches per calendar year and the 10
  largest breaches by excess loss. Tables are in `output/backtest_*.csv`.
- `R/charts.R`: three PNGs in `figures/`:
  - full-period returns vs both VaR thresholds, with breaches marked;
  - the same, zoomed to January 2020 to April 2021;
  - breaches per year, with expected breaches drawn as a reference line.
- `.gitignore`: figures and backtest tables are committed. Raw prices, the
  daily forecast series and the quality report stay local.

### Block 4: three more VaR models
- `R/var_models.R` now has five models, all with no look-ahead:
  - historical;
  - normal;
  - rolling Student-t: same rolling mean and sd as the normal model, with
    multiplier `-qt(0.01, 5) * sqrt(3/5)` = 2.606464 (t rescaled to unit
    variance);
  - EWMA normal: `v_t = 0.94 v_{t-1} + 0.06 r_{t-1}^2`, mean taken as zero,
    multiplier `qnorm(0.99)`;
  - EWMA Student-t: the same EWMA volatility with the t multiplier.
  - The EWMA seed is the sample variance of the first 250 returns
    (2013-06-04 onwards), used as the forecast for the 251st day. That is
    149 days before the first forecast, so the seed's remaining weight is
    about 1e-4.
- Backtest tables, the Kupiec test and per-year counts cover all five models.
- `output/backtest_t_df_sensitivity.csv`: both t models with df 4, 5 and 8.
- `R/hand_check.R`: `hand_check_ewma()` prints yesterday's EWMA variance and
  return, today's EWMA variance, and all five VaRs.
- Charts: the 2020-2021 zoom adds the EWMA normal line and its breaches; the
  per-year chart shows all five models.

## Log

- 2026-11-13: Five models backtested; only historical passes Kupiec; EWMA reduces crisis breaches but not the total; next: breach size and clustering.
- 2026-11-13: Per-year breaches: rolling models cluster in 2018, 2020 and 2022; EWMA is evenly spread at about twice the expected rate.

## Next

1. **Independence of breaches.** Christoffersen (1998) test for independence
   and conditional coverage, for all five models. Kupiec checks only how many
   breaches there are, not whether they bunch together. This matters most for
   comparing the rolling-window models with the EWMA ones.
2. **Kupiec per year.** The per-year counts are already computed; add the
   test statistic and p-value per year, noting that about 2.5 expected
   breaches a year gives the test very little power.
3. **Basel traffic light.** Rolling 250-day breach counts mapped to the
   green, yellow and red zones.
4. **Unit tests** (`testthat`):
   - VaR functions on small known inputs;
   - Kupiec against a hand-computed value, including the edge cases of 0
     breaches and n breaches;
   - no window includes its own forecast date.
5. **Write-up.** Interpret the results in the README (not done yet, on
   purpose).

## Later / open questions

- More models (not started on purpose): filtered historical simulation,
  GARCH(1,1), a t distribution with fitted rather than fixed df.
- EWMA choices to revisit: the seed (sample variance vs mean of squared
  returns) and lambda (0.94 is the RiskMetrics default, not fitted).
- Decide whether to commit a frozen copy of `output/var_forecasts.csv`.
  Yahoo revises adjusted prices, so results can drift slightly between
  downloads.
- Check Yahoo's adjusted prices against a second source (e.g. the issuers'
  published NAV and distribution history).
