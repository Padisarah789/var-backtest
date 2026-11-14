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

### Block 5: EWMA lambda sensitivity
- `R/ewma_sensitivity.R` covers lambda 0.90, 0.92, 0.94, 0.95, 0.96, 0.966,
  0.97, 0.977, 0.98 and 0.99, each with the normal and the Student-t (df 5)
  multipliers. The five main models are unchanged; lambda 0.94 is checked
  to reproduce them.
- `output/ewma_lambda_sensitivity.csv`: for each lambda and multiplier,
  - half-life `log(0.5) / log(lambda)`;
  - total breaches and breach rate;
  - the Kupiec test;
  - minimum and maximum breaches in any calendar year.
- Out-of-sample check (`output/ewma_lambda_out_of_sample.csv`):
  - lambda is chosen on 2015-2019 only, as the one with breaches closest to
    expected, then tested on 2020-2025;
  - ties go to the higher in-sample Kupiec p-value, then to the lambda
    closest to 0.94;
  - both multipliers had a three-way tie, so the tie-break decided the
    choice. `output/ewma_lambda_out_of_sample_candidates.csv` lists every
    candidate's 2020-2025 result (for reference only, not used in the
    choice).
- `figures/ewma_lambda_sensitivity.png`: total breaches against half-life
  (log scale) for both multipliers, with lambda on the top axis. Reference
  lines show the expected 27.66 and the rolling 250-day results (normal,
  Student-t, historical).

### Day 5, block 1: FHS, floor and breach severity
- `R/var_models.R` adds two models; the five existing ones are unchanged
  (checked: identical `var_forecasts.csv` columns and backtest rows).
  - Filtered historical simulation (`fhs`): each past return is divided by
    the EWMA (lambda 0.94) volatility forecast made for its day; the
    multiplier is minus the type-1 1% quantile (3rd worst) of the previous
    250 standardised returns; VaR = multiplier * today's EWMA volatility.
  - Floor (`floor`): the larger of the EWMA normal VaR and the historical VaR.
  - Standardisation needs EWMA volatilities from 2014-01-06 (row 150), but
    the main EWMA starts at row 251. So the standardised returns use the
    same lambda-0.94 recursion seeded on the sample variance of the first 50
    returns (`FHS_SEED_N`). From the first forecast date on, its volatility
    is within 1.7e-5 (relative) of the main EWMA. Today's volatility is the
    main EWMA's, so the FHS multiplier compares directly with 2.326 and 2.606.
- Backtest tables (summary, Kupiec, per-year, largest breaches) cover all
  seven models; the console prints the full per-year table with totals.
- `output/backtest_severity.csv`: per model, average VaR, breaches, and the
  mean, median, largest and total excess (loss minus VaR on breach days, in
  percentage points of portfolio value).
- `output/backtest_fhs_multiplier.csv`: minimum, median and maximum of the
  FHS multiplier (with dates), next to 2.326 and 2.606, and the share of
  days the FHS multiplier is above each.
- `R/hand_check.R`: `hand_check_fhs()` for 2020-03-12 and 2020-03-16 prints
  today's EWMA volatility, the three worst standardised returns in the
  window with their dates, the FHS multiplier, and the FHS and floor VaRs.
- Charts: the per-year chart shows all seven models (green FHS, violet
  floor); new `figures/mean_excess_per_model.png`.

## Log

- 2026-11-13: Five models backtested; only historical passes Kupiec; EWMA reduces crisis breaches but not the total; next: breach size and clustering.
- 2026-11-13: Per-year breaches: rolling models cluster in 2018, 2020 and 2022; EWMA is evenly spread at about twice the expected rate.
- 2026-11-13: EWMA half-life sensitivity: normal multiplier flat at 62 to 67; t multiplier best near a 30-day half-life (40 to 41) but still rejected; out-of-sample result depends on tie-break.

## Next

1. **Independence of breaches.** Christoffersen (1998) test for independence
   and conditional coverage, for all seven models. Kupiec checks only how many
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

- More models (not started on purpose): GARCH(1,1) (also as the filter for
  FHS), a t distribution with fitted rather than fixed df.
- FHS standardisation seed: the 50-return seed is a choice. An alternative
  is to download earlier history so the main EWMA covers the first FHS
  window, but that would change the EWMA seed of the existing models.
- EWMA choices to revisit: the seed (sample variance vs mean of squared
  returns). Lambda selection by breach count is coarse: counts move by only
  a few breaches across lambdas, so ties are common. A smoother criterion
  (e.g. a likelihood or quantile loss on 2015-2019) would avoid relying on a
  tie-break.
- Decide whether to commit a frozen copy of `output/var_forecasts.csv`.
  Yahoo revises adjusted prices, so results can drift slightly between
  downloads.
- Check Yahoo's adjusted prices against a second source (e.g. the issuers'
  published NAV and distribution history).
