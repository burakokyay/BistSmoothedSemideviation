# =============================================================================
# Gaussian-smoothed mean-semideviation portfolios on BIST-30 stocks
# Rolling-window out-of-sample backtest  --  R console version
#
# HOW TO RUN
#   1. Put your price file (e.g. bist30_prices.csv) in your working directory.
#      Check the folder with getwd(); change it with setwd("C:/path/to/folder").
#   2. Edit DATA_FILE below if your file has a different name.
#   3. Paste this whole script into the R console, or run:
#        source("smoothed_msd_backtest_console.R")
#   Results (tables, figures) are written to the "results" sub-folder.
#
# INPUT FORMAT (CSV)
#   First column  = dates (YYYY-MM-DD), any header name (e.g. "Index", "Date")
#   Other columns = adjusted closing prices, one column per stock
# =============================================================================

# ---- 0. Settings ------------------------------------------------------------
DATA_FILE <- "bist30_prices.csv"   # file name inside your working directory
OUT_DIR   <- "results"             # output folder (created if missing)

# ---- 0b. Packages -----------------------------------------------------------
if (!requireNamespace("nloptr", quietly = TRUE)) install.packages("nloptr")
suppressPackageStartupMessages(library(nloptr))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
outdir <- OUT_DIR
set.seed(20261108)

# ----------------------------------------------------------------------------
# 1. Data import
# ----------------------------------------------------------------------------
cat("Working directory:", getwd(), "\n")
infile <- file.path(getwd(), DATA_FILE)
if (!file.exists(infile)) {
  csvs <- list.files(getwd(), pattern = "\\.csv$", ignore.case = TRUE)
  cat("File '", DATA_FILE, "' not found. CSV files here: ",
      if (length(csvs)) paste(csvs, collapse = ", ") else "none", "\n", sep = "")
  if (interactive()) {
    cat("Please pick the price file in the dialog...\n")
    infile <- file.choose()
  } else stop("Price file not found; set DATA_FILE or setwd().")
}

# Read: auto-detect comma or semicolon separator (Excel in Turkish locale
# often saves with ';' and ',' as decimal mark)
first_line <- readLines(infile, n = 1)
if (grepl(";", first_line)) {
  px <- read.csv2(infile, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  px <- read.csv(infile, check.names = FALSE, stringsAsFactors = FALSE)
}

# Dates: try common formats
parse_dates <- function(x) {
  x <- as.character(x)
  for (f in c("%Y-%m-%d", "%d.%m.%Y", "%d/%m/%Y", "%m/%d/%Y", "%Y/%m/%d")) {
    d <- as.Date(x, format = f)
    if (mean(!is.na(d)) > 0.95) return(d)
  }
  stop("Could not read the date column; use YYYY-MM-DD format.")
}
dates <- parse_dates(px[[1]])
P <- as.matrix(sapply(px[, -1, drop = FALSE], function(v) as.numeric(as.character(v))))
colnames(P) <- sub("\\.IS$", "", colnames(px)[-1])     # drop Yahoo ".IS" suffix
keep_rows <- !is.na(dates); dates <- dates[keep_rows]; P <- P[keep_rows, , drop = FALSE]
ord <- order(dates); dates <- dates[ord]; P <- P[ord, , drop = FALSE]

# Drop stocks with no data at all
empty <- colSums(!is.na(P)) == 0
if (any(empty)) {
  cat("Dropping stocks with no data:", paste(colnames(P)[empty], collapse = ", "), "\n")
  P <- P[, !empty, drop = FALSE]
}
cat(sprintf("Loaded %d price rows x %d stocks, %s to %s\n",
            nrow(P), ncol(P), min(dates), max(dates)))

for (j in seq_len(ncol(P))) {                # carry last price over short gaps
  for (i in 2:nrow(P)) if (is.na(P[i, j])) P[i, j] <- P[i - 1, j]
}
R <- P[-1, , drop = FALSE] / P[-nrow(P), , drop = FALSE] - 1   # simple daily returns
rdates <- dates[-1]
R[is.na(R)] <- 0

# Exchange holidays appear as rows where no stock moves: drop them
holiday <- apply(abs(R) < 1e-12, 1, all)
R <- R[!holiday, , drop = FALSE]; rdates <- rdates[!holiday]

# Single-stock moves above 25% exceed Borsa Istanbul daily price limits;
# they are unadjusted corporate actions (bonus issues) in the source data.
bad <- which(abs(R) > 0.25, arr.ind = TRUE)
data_fixes <- data.frame(date  = rdates[bad[, 1]],
                         stock = colnames(R)[bad[, 2]],
                         raw_return = R[bad])
R[abs(R) > 0.25] <- 0
write.csv(data_fixes, file.path(outdir, "data_fixes.csv"), row.names = FALSE)

n_assets <- ncol(R)
cat(sprintf("Assets: %d, return days: %d (%s to %s), holidays dropped: %d, fixes: %d\n",
            n_assets, nrow(R), min(rdates), max(rdates), sum(holiday), nrow(data_fixes)))

# ----------------------------------------------------------------------------
# 2. Smoothed mean-semideviation (closed form) and gradient
#    Loss of day j: m_j = -w'D_j.  Smoothing adds N(0, sigma^2 I) to each
#    return vector, so each loss becomes N(m_j, s^2) with s = sigma*||w||.
#    msd^sigma[w] = mu + sqrt( mean_j E[(a_j + s Z)_+^2] ),  a_j = m_j - mu
#    E[(a + sZ)_+^2] = (a^2 + s^2) Phi(a/s) + a s phi(a/s)
# ----------------------------------------------------------------------------
msd_obj <- function(w, D, sigma) {
  Dbar <- colMeans(D)
  C    <- sweep(D, 2, Dbar)                 # centred returns
  a    <- -as.vector(C %*% w)               # centred losses
  mu   <- -sum(Dbar * w)
  if (sigma <= 0) {
    ap <- pmax(a, 0)
    V  <- mean(ap^2)
    dV <- -2 * colMeans(C * ap)
  } else {
    s  <- sigma * sqrt(sum(w^2))
    d  <- a / s
    Ph <- pnorm(d); ph <- dnorm(d)
    V  <- mean((a^2 + s^2) * Ph + a * s * ph)
    dV <- -2 * colMeans(C * (a * Ph + s * ph)) + 2 * sigma^2 * mean(Ph) * w
  }
  V  <- max(V, 1e-16)
  list(objective = mu + sqrt(V), gradient = -Dbar + dV / (2 * sqrt(V)))
}

solve_msd <- function(D, sigma, w0 = rep(1 / ncol(D), ncol(D))) {
  k <- ncol(D)
  res <- nloptr(
    x0 = w0,
    eval_f = function(w) msd_obj(w, D, sigma),
    lb = rep(0, k), ub = rep(1, k),
    eval_g_eq = function(w) list(constraints = sum(w) - 1, jacobian = rep(1, k)),
    opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-10,
                ftol_rel = 1e-12, maxeval = 2000))
  w <- pmax(res$solution, 0); w / sum(w)
}

# ----------------------------------------------------------------------------
# 3. Check closed form against numerical integration
# ----------------------------------------------------------------------------
# Independent check: integrate each Gaussian component numerically
mc_check <- function(D, w, sigma) {
  m  <- -as.vector(D %*% w); mu <- mean(m); a <- m - mu
  s  <- sigma * sqrt(sum(w^2))
  V  <- mean(sapply(a, function(aj)
          integrate(function(z) pmax(aj + s * z, 0)^2 * dnorm(z),
                    -Inf, Inf, rel.tol = 1e-12)$value))
  mu + sqrt(V)
}
chk <- do.call(rbind, lapply(1:20, function(i) {
  Nw <- sample(c(20, 40, 60), 1); st <- sample.int(nrow(R) - Nw, 1)
  D  <- R[st:(st + Nw - 1), ]
  w  <- runif(n_assets); w <- w / sum(w)
  sg <- sample(c(0.25, 0.5, 1), 1) * mean(apply(D, 2, sd))
  cf <- msd_obj(w, D, sg)$objective; mc <- mc_check(D, w, sg)
  data.frame(N = Nw, sigma = sg, closed_form = cf, quadrature = mc,
             rel_diff_pct = 100 * abs(cf - mc) / abs(mc))
}))
write.csv(chk, file.path(outdir, "closed_form_check.csv"), row.names = FALSE)
cat(sprintf("Closed form vs numerical integration: max rel. diff %.3f%%, median %.3f%%\n",
            max(chk$rel_diff_pct), median(chk$rel_diff_pct)))

# ----------------------------------------------------------------------------
# 4. Rolling-window backtest
# ----------------------------------------------------------------------------
windows <- c(20, 40, 60)
mults   <- c(0, 0.25, 0.5, 1)
tc      <- 0.0010                           # 10 bp one-way cost for net results

ym <- format(rdates, "%Y-%m")
reb_idx <- which(c(TRUE, ym[-1] != ym[-length(ym)]))   # first day of each month
reb_idx <- reb_idx[reb_idx > max(windows)]              # common start for all N
cat(sprintf("Out-of-sample: %s to %s, %d monthly rebalances\n",
            rdates[reb_idx[1]], max(rdates), length(reb_idx)))

run_strategy <- function(weight_fun) {
  T  <- nrow(R); start <- reb_idx[1]
  ret <- numeric(0); W <- NULL; turn <- numeric(0); est <- numeric(0)
  w_drift <- rep(0, n_assets)
  ends <- c(reb_idx[-1] - 1, T)
  for (k in seq_along(reb_idx)) {
    t0 <- reb_idx[k]
    wr <- weight_fun(t0)
    w  <- wr$w
    turn <- c(turn, if (k == 1) 1 else sum(abs(w - w_drift)))
    W  <- rbind(W, w); est <- c(est, wr$est)
    h  <- t0:ends[k]
    v  <- w
    for (t in h) {                          # buy-and-hold within the month
      rp <- sum(v * R[t, ]); ret <- c(ret, rp)
      v  <- v * (1 + R[t, ]) / (1 + rp)
    }
    w_drift <- v
  }
  list(ret = ret, dates = rdates[reb_idx[1]:T], W = W, turnover = turn,
       est = est, reb = reb_idx)
}

semidev <- function(x) sqrt(mean(pmax(mean(x) - x, 0)^2))     # downside, returns
results <- list()
for (N in windows) for (m in mults) {
  key <- sprintf("N%d_s%.2f", N, m)
  results[[key]] <- run_strategy(function(t0) {
    D  <- R[(t0 - N):(t0 - 1), ]
    sg <- m * mean(apply(D, 2, sd))
    w  <- solve_msd(D, sg)
    list(w = w, est = semidev(as.vector(D %*% w)))  # in-sample semideviation
  })
  cat("done", key, "\n")
}
results[["EW"]] <- run_strategy(function(t0)
  list(w = rep(1 / n_assets, n_assets), est = NA))

# ----------------------------------------------------------------------------
# 5. Performance statistics
# ----------------------------------------------------------------------------
semidev <- function(x) sqrt(mean(pmax(mean(x) - x, 0)^2))     # downside, returns
maxdd   <- function(x) { v <- cumprod(1 + x); min(v / cummax(v) - 1) }

net_returns <- function(res) {
  r <- res$ret
  pos <- match(res$reb, res$reb[1]:(res$reb[1] + length(r) - 1))
  r[pos] <- r[pos] - tc * res$turnover
  r
}

stats <- function(res) {
  r <- res$ret; rn <- net_returns(res)
  hhi <- rowSums(res$W^2)
  # realised next-month semideviation vs in-sample estimate (optimism ratio)
  ends <- c(res$reb[-1] - 1, res$reb[1] + length(r) - 1)
  real_msd <- sapply(seq_along(res$reb), function(k) {
    x <- r[(res$reb[k]:ends[k]) - res$reb[1] + 1]; semidev(x)
  })
  c(ann_return      = 252 * mean(r),
    ann_vol         = sqrt(252) * sd(r),
    ret_vol_ratio   = sqrt(252) * mean(r) / sd(r),
    ann_semidev     = sqrt(252) * semidev(r),
    max_drawdown    = maxdd(r),
    ann_return_net  = 252 * mean(rn),
    ret_vol_net     = sqrt(252) * mean(rn) / sd(rn),
    avg_hhi         = mean(hhi),
    eff_n           = mean(1 / hhi),
    avg_max_weight  = mean(apply(res$W, 1, max)),
    avg_turnover    = mean(res$turnover[-1]),
    optimism_ratio  = if (all(is.na(res$est))) NA else
                        median(real_msd / res$est, na.rm = TRUE))
}

tab <- do.call(rbind, lapply(results, stats))
tab <- data.frame(strategy = rownames(tab), round(tab, 5), row.names = NULL)
write.csv(tab, file.path(outdir, "performance_table.csv"), row.names = FALSE)
print(tab)

# ----------------------------------------------------------------------------
# 6. Paired stationary block bootstrap: smoothed (sigma>0) minus unsmoothed
# ----------------------------------------------------------------------------
boot_idx <- function(T, L = 21) {
  idx <- integer(0)
  while (length(idx) < T) {
    s <- sample.int(T, 1); len <- rgeom(1, 1 / L) + 1
    idx <- c(idx, ((s - 1 + 0:(len - 1)) %% T) + 1)
  }
  idx[1:T]
}
B <- 2000
Tn <- length(results[["EW"]]$ret)
BI <- replicate(B, boot_idx(Tn))
boot_rows <- list()
for (N in windows) {
  base <- results[[sprintf("N%d_s%.2f", N, 0)]]$ret
  for (m in mults[-1]) {
    alt <- results[[sprintf("N%d_s%.2f", N, m)]]$ret
    f_sd <- function(i) sqrt(252) * (semidev(alt[i]) - semidev(base[i]))
    f_rv <- function(i) sqrt(252) * (mean(alt[i]) / sd(alt[i]) - mean(base[i]) / sd(base[i]))
    d_sd <- apply(BI, 2, f_sd); d_rv <- apply(BI, 2, f_rv)
    boot_rows[[length(boot_rows) + 1]] <- data.frame(
      N = N, mult = m,
      d_semidev = f_sd(1:Tn), p_semidev = mean(d_sd >= 0),
      d_ret_vol = f_rv(1:Tn), p_ret_vol_two_sided = 2 * min(mean(d_rv <= 0), mean(d_rv >= 0)))
  }
}
boot_tab <- do.call(rbind, boot_rows)
write.csv(boot_tab, file.path(outdir, "bootstrap_tests.csv"), row.names = FALSE)
print(boot_tab)

# ----------------------------------------------------------------------------
# 7. Exhibit data: cumulative wealth (N = 20) and weights concentration
# ----------------------------------------------------------------------------
wealth <- data.frame(date = results[["EW"]]$dates)
for (k in c("N20_s0.00", "N20_s0.50", "N20_s1.00", "EW"))
  wealth[[k]] <- cumprod(1 + results[[k]]$ret)
write.csv(wealth, file.path(outdir, "wealth_N20.csv"), row.names = FALSE)

png(file.path(outdir, "exhibit_wealth_N20.png"), width = 1800, height = 1000, res = 200)
matplot(wealth$date, log(as.matrix(wealth[, -1])), type = "l", lty = 1, lwd = 1.6,
        col = c("#B8442E", "#2E6FB8", "#1F8A5B", "#555555"),
        xlab = "", ylab = "Log cumulative wealth (TRY)", xaxt = "n",
        main = "Out-of-sample growth, 20-day estimation window")
axis.Date(1, wealth$date, format = "%Y")
legend("topleft", bty = "n", lwd = 2,
       col = c("#B8442E", "#2E6FB8", "#1F8A5B", "#555555"),
       legend = c("Unsmoothed (sigma = 0)", "Smoothed, 0.5x", "Smoothed, 1x", "Equal weight"))
dev.off()

# ----------------------------------------------------------------------------
# 8. Sub-period robustness (split at 2021-01-01) and summary exhibits
# ----------------------------------------------------------------------------
split_date <- as.Date("2021-01-01")
sub_rows <- list()
for (k in names(results)) {
  r <- results[[k]]$ret; d <- results[[k]]$dates
  for (per in c("2015-2020", "2021-2026")) {
    idx <- if (per == "2015-2020") d < split_date else d >= split_date
    x <- r[idx]
    sub_rows[[length(sub_rows) + 1]] <- data.frame(
      strategy = k, period = per,
      ann_return = 252 * mean(x), ann_vol = sqrt(252) * sd(x),
      ret_vol_ratio = sqrt(252) * mean(x) / sd(x),
      ann_semidev = sqrt(252) * semidev(x), max_drawdown = maxdd(x))
  }
}
sub_tab <- do.call(rbind, sub_rows)
write.csv(sub_tab, file.path(outdir, "subperiod_table.csv"), row.names = FALSE)

cols <- c("#B8442E", "#2E6FB8", "#1F8A5B")
png(file.path(outdir, "exhibit_risk_diversification.png"), width = 2000, height = 900, res = 200)
par(mfrow = c(1, 2), mar = c(4.2, 4.5, 2.5, 1))
plot(NA, xlim = c(0, 1), ylim = range(tab$ann_semidev[1:12] * 100, tab$ann_semidev[13] * 100),
     xlab = "Smoothing level (multiple of average stock volatility)",
     ylab = "Realised annual semideviation (%)", main = "Downside risk")
for (i in 1:3) lines(mults, tab$ann_semidev[(i - 1) * 4 + 1:4] * 100, type = "b",
                     pch = 19, col = cols[i], lwd = 2)
abline(h = tab$ann_semidev[13] * 100, lty = 2, col = "#555555")
legend("topright", bty = "n", lwd = 2, pch = c(19, 19, 19, NA), lty = c(1, 1, 1, 2),
       col = c(cols, "#555555"), legend = c("20-day window", "40-day window", "60-day window", "Equal weight"))
plot(NA, xlim = c(0, 1), ylim = c(0, 11),
     xlab = "Smoothing level (multiple of average stock volatility)",
     ylab = "Effective number of stocks (1/HHI)", main = "Diversification")
for (i in 1:3) lines(mults, tab$eff_n[(i - 1) * 4 + 1:4], type = "b", pch = 19, col = cols[i], lwd = 2)
dev.off()

png(file.path(outdir, "exhibit_optimism.png"), width = 1400, height = 900, res = 200)
par(mar = c(4.2, 4.5, 2.5, 1))
plot(NA, xlim = c(0, 1), ylim = c(0.95, 1.45),
     xlab = "Smoothing level (multiple of average stock volatility)",
     ylab = "Realised / estimated semideviation", main = "How much risk the optimiser hides")
for (i in 1:3) lines(mults, tab$optimism_ratio[(i - 1) * 4 + 1:4], type = "b", pch = 19, col = cols[i], lwd = 2)
abline(h = 1, lty = 3)
legend("topright", bty = "n", lwd = 2, pch = 19, col = cols,
       legend = c("20-day window", "40-day window", "60-day window"))
dev.off()

# ----------------------------------------------------------------------------
# 9. Additional exhibits
# ----------------------------------------------------------------------------
dd_path <- function(x) { v <- cumprod(1 + x); v / cummax(v) - 1 }
d_all <- results[["EW"]]$dates
pal <- c(unsm = "#B8442E", sm = "#1F8A5B", ew = "#555555")

# 9a. Drawdowns through time (20-day window)
png(file.path(outdir, "exhibit_drawdowns_N20.png"), width = 1800, height = 900, res = 200)
par(mar = c(3, 4.5, 2.5, 1))
dd <- cbind(dd_path(results[["N20_s0.00"]]$ret), dd_path(results[["N20_s1.00"]]$ret),
            dd_path(results[["EW"]]$ret)) * 100
matplot(d_all, dd, type = "l", lty = 1, lwd = 1.4, col = pal, xaxt = "n",
        xlab = "", ylab = "Drawdown from previous peak (%)",
        main = "Drawdowns, 20-day estimation window")
axis.Date(1, d_all, format = "%Y"); abline(h = 0, col = "grey70")
legend("bottomleft", bty = "n", lwd = 2, col = pal,
       legend = c("Unsmoothed (sigma = 0)", "Smoothed, 1x", "Equal weight"))
dev.off()

# 9b. Concentration through time: effective number of stocks at each rebalance
png(file.path(outdir, "exhibit_effn_time.png"), width = 1800, height = 900, res = 200)
par(mar = c(3, 4.5, 2.5, 1))
reb_dates <- rdates[reb_idx]
effn <- sapply(c("N20_s0.00", "N20_s1.00", "N60_s0.00", "N60_s1.00"),
               function(k) 1 / rowSums(results[[k]]$W^2))
matplot(reb_dates, effn, type = "l", lwd = c(1.6, 1.6, 1.2, 1.2), lty = c(1, 1, 2, 2),
        col = pal[c(1, 2, 1, 2)], xaxt = "n", xlab = "", ylim = c(0, max(effn) * 1.05),
        ylab = "Effective number of stocks (1/HHI)",
        main = "Diversification at each monthly rebalance")
axis.Date(1, reb_dates, format = "%Y")
legend("topleft", bty = "n", lwd = 2, lty = c(1, 1, 2, 2), col = pal[c(1, 2, 1, 2)], ncol = 2,
       legend = c("20-day, unsmoothed", "20-day, smoothed 1x",
                  "60-day, unsmoothed", "60-day, smoothed 1x"))
dev.off()

# 9d. Rolling one-year semideviation (20-day window)
roll_sd <- function(x, L = 252) sapply(L:length(x), function(i) sqrt(252) * semidev(x[(i - L + 1):i]))
rs <- cbind(roll_sd(results[["N20_s0.00"]]$ret), roll_sd(results[["N20_s1.00"]]$ret),
            roll_sd(results[["EW"]]$ret)) * 100
png(file.path(outdir, "exhibit_rolling_semidev_N20.png"), width = 1800, height = 900, res = 200)
par(mar = c(3, 4.5, 2.5, 1))
matplot(d_all[252:length(d_all)], rs, type = "l", lty = 1, lwd = 1.4, col = pal, xaxt = "n",
        xlab = "", ylab = "Annualised semideviation (%)",
        main = "Rolling one-year downside risk, 20-day estimation window")
axis.Date(1, d_all, format = "%Y")
legend("topleft", bty = "n", lwd = 2, col = pal,
       legend = c("Unsmoothed (sigma = 0)", "Smoothed, 1x", "Equal weight"))
dev.off()
roll_share <- mean(rs[, 2] < rs[, 1])
cat(sprintf("Rolling 1-yr semideviation lower with smoothing on %.1f%% of days\n", 100 * roll_share))
max_dd_tab <- data.frame(strategy = c("N20_s0.00", "N20_s1.00", "EW"),
                         worst_dd = apply(dd, 2, min),
                         days_below_minus20 = colSums(dd < -20))
print(max_dd_tab)

save(results, tab, boot_tab, chk, data_fixes, sub_tab, file = file.path(outdir, "backtest.RData"))
cat("Finished.\n")
