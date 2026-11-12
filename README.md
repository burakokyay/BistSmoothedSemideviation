# Gaussian-Smoothed Mean–Semideviation Portfolios on Borsa İstanbul

Replication code for:

Okyay, B. "Smoothing Small-Sample Risk Estimates: Gaussian-Smoothed Mean–Semideviation Portfolios on Borsa İstanbul." Manuscript submitted to *The Journal of Portfolio Management*.

## Files

| File | What it does |
|---|---|
| `download_prices.R` | Downloads adjusted daily closing prices for the BIST-30 stocks from Yahoo Finance and saves `bist30_prices.csv`. |
| `smoothed_msd_backtest_console.R` | Runs the full backtest and writes every table and figure in the paper to a `results/` folder. |
| `bist30_prices.csv` | The price file used in the paper (2015-01-01 to 2026-11-06), if included. |

## How to run

1. Install R (version 4.0 or later).
2. Put the files in one folder and set it as the working directory in R:
   ```r
   setwd("path/to/this/folder")
   ```
3. If `bist30_prices.csv` is not included, create it:
   ```r
   source("download_prices.R")
   ```
4. Run the backtest (about 30 seconds):
   ```r
   source("smoothed_msd_backtest_console.R")
   ```

The script installs the `nloptr` package if it is missing.

## Output (in `results/`)

- `performance_table.csv`: out-of-sample performance of all strategies (Exhibit 1)
- `bootstrap_tests.csv`: block-bootstrap tests, smoothed minus unsmoothed (Exhibit 5)
- `subperiod_table.csv`: results for 2015–2020 and 2021–2026 (Exhibit 8)
- `closed_form_check.csv`: closed-form formula against numerical integration
- `data_fixes.csv`: daily returns above 25% set to zero (unadjusted corporate actions)
- Figures: `exhibit_risk_diversification.png`, `exhibit_effn_time.png`, `exhibit_optimism.png`, `exhibit_wealth_N20.png`, `exhibit_drawdowns_N20.png`, `exhibit_rolling_semidev_N20.png`

## Data note

Prices are Yahoo Finance adjusted closes. Yahoo revises its data from time to time, so a fresh download may differ slightly from the file used in the paper.

## Contact

Burak Okyay, Department of Business Administration, Faculty of Art and Social Sciences, Mudanya University.
