# =============================================================================
# Download BIST-30 adjusted closing prices from Yahoo Finance
# Produces bist30_prices.csv in the working directory (first column = date).
# Requires: install.packages("quantmod")
# Note: Yahoo Finance data are revised from time to time, so a fresh download
# may differ slightly from the file used in the paper.
# =============================================================================

if (!requireNamespace("quantmod", quietly = TRUE)) install.packages("quantmod")
suppressPackageStartupMessages(library(quantmod))

tickers <- c("AKBNK", "ALARK", "ARCLK", "ASELS", "BIMAS", "DOHOL", "EKGYO",
             "ENKAI", "EREGL", "FROTO", "GARAN", "HALKB", "HEKTS", "ISCTR",
             "KCHOL", "KOZAL", "KRDMD", "PETKM", "PGSUS", "SAHOL", "SASA",
             "SISE", "TAVHL", "TCELL", "THYAO", "TOASO", "TTKOM", "TUPRS",
             "VAKBN", "YKBNK")
start_date <- "2015-01-01"
end_date   <- "2026-11-07"

prices <- list()
for (tk in tickers) {
  x <- tryCatch(
    getSymbols(paste0(tk, ".IS"), src = "yahoo", from = start_date,
               to = end_date, auto.assign = FALSE),
    error = function(e) NULL)
  if (is.null(x)) { cat("No data for", tk, "\n"); next }
  prices[[tk]] <- Ad(x)
  cat("Downloaded", tk, "\n")
}

px <- do.call(merge, prices)
colnames(px) <- names(prices)
out <- data.frame(Index = index(px), coredata(px), check.names = FALSE)
write.csv(out, "bist30_prices.csv", row.names = FALSE)
cat("Saved bist30_prices.csv:", nrow(out), "rows,", ncol(out) - 1, "stocks\n")
