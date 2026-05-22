getwd()
setwd("D:/MS ECON/Fall 2025/Capstone/R/Midterm")
#Installing Package

install.packages("quantmod")
install.packages("PerformanceAnalytics")
install.packages("tidyverse")
install.packages("quadprog")
install.packages("flextable")
install.packages("officer")

# Load libraries
library(readxl)
library(ggplot2)
library(quantmod)
library(PerformanceAnalytics)
library(xts)
library(zoo)
library(tidyverse)
library(PerformanceAnalytics)
library(PortfolioAnalytics)
library(ROI)
library(ROI.plugin.quadprog)
library(ROI.plugin.glpk)
library(quadprog)
library(dplyr)
library(flextable)
library(officer)

## -----------------------------
## 0. Parameters & tickers
## -----------------------------
sample_start <- as.Date("2023-01-01")
sample_end   <- as.Date("2025-09-29")

from <- format(sample_start)     # << tighten the pull
to   <- format(sample_end)

tickers <- c(
  BTC="BTC-USD", GOLD="GC=F", SPX="^GSPC",
  DOGE="DOGE-USD", XRP="XRP-USD", ETH="ETH-USD",
  FVX="^FVX", TNX="^TNX",
  EURUSD="EURUSD=X", USDJPY="JPY=X", GBPUSD="GBPUSD=X", AUDUSD="AUDUSD=X"
)

rf_daily <- 0      # risk-free set to zero per methodology
alpha_cvar <- 0.95
roll_win <- 252    # 252-day estimation window
gamma_ceq <- 1     # CEQ parameter
f_speed   <- 0.05  # LIBRO trading speed (tune)
M_notional <- 100    # LIBRO notional

## Helper to extract Adjusted (or Close when Adjusted missing)
AdSafe <- function(x) if (NCOL(x) >= 6) Ad(x) else Cl(x)

## -----------------------------
## 1. Asset and data preparation
## -----------------------------

# 1.1 Load prices (handles uneven calendars by design)
getSymbols(tickers, src="yahoo", from=from, to=to, auto.assign=TRUE, warnings=FALSE)

assets <- list(
  BTC    = AdSafe(`BTC-USD`),
  GOLD   = AdSafe(`GC=F`),
  SPX    = AdSafe(`GSPC`),
  DOGE   = AdSafe(`DOGE-USD`),
  XRP    = AdSafe(`XRP-USD`),
  ETH    = AdSafe(`ETH-USD`),
  FVX    = Cl(`FVX`),   # yields (close)
  TNX    = Cl(`TNX`),   # yields (close)
  EURUSD = AdSafe(`EURUSD=X`),
  USDJPY = AdSafe(`JPY=X`),
  GBPUSD = AdSafe(`GBPUSD=X`),
  AUDUSD = AdSafe(`AUDUSD=X`)
)  
# ===============================
# 1.2 Align to NYSE calendar (NO back-fill)
# ===============================
biz_idx <- index(assets$SPX)

# just align to the S&P 500 business-day index; don't fill here
align_to_biz <- function(x, idx) x[idx]

prices_biz <- do.call(merge, lapply(assets, align_to_biz, idx = biz_idx))
colnames(prices_biz) <- names(assets)

# forward-fill only (carry last observation forward); DO NOT back-fill
prices_biz <- na.locf(prices_biz, na.rm = FALSE)

# start the panel when all assets have started (drop pre-start rows)
prices_biz <- na.omit(prices_biz)
# --- ROBUSTNESS WINDOW: keep only 2023–2025
prices_biz <- window(prices_biz, start = sample_start, end = sample_end)
# ===============================
# 1.3 Returns
#   - price-like assets: log and simple % returns
#   - yields (FVX, TNX): use Δy (NOT %)
# ===============================
yield_cols <- c("FVX","TNX")
price_cols <- setdiff(colnames(prices_biz), yield_cols)

# LOG returns for price-like assets
ret_log_prices <- diff(log(prices_biz[, price_cols]))

# Δy for yields
dyields <- diff(prices_biz[, yield_cols])

# align and combine
dyields   <- dyields[index(ret_log_prices)]
ret_log   <- merge(ret_log_prices, dyields)
colnames(ret_log) <- c(price_cols, yield_cols)
ret_log   <- na.omit(ret_log)

# SIMPLE returns (for inflation-adjusted step)
ret_simple_prices <- ROC(prices_biz[, price_cols], type = "discrete")
dyields2          <- diff(prices_biz[, yield_cols])
dyields2          <- dyields2[index(ret_simple_prices)]
ret_simple        <- merge(ret_simple_prices, dyields2)
colnames(ret_simple) <- c(price_cols, yield_cols)
ret_simple        <- na.omit(ret_simple)

# Yields (FVX, TNX) are levels; use yield changes (NOT % returns). Convert to log-return slots:
# columns
yield_cols  <- c("FVX","TNX")
price_cols  <- setdiff(colnames(prices_biz), yield_cols)

# (A) log returns for price-like assets
ret_log_prices <- diff(log(prices_biz[, price_cols]))

# (B) yield changes for FVX/TNX (NOT % returns)
dyields <- diff(prices_biz[, yield_cols])

# (C) align and combine
common_idx <- index(ret_log_prices)                 # same as prices_biz[-1]
dyields     <- dyields[common_idx]                  # align to ret_log index
ret_log     <- merge(ret_log_prices, dyields)       # order preserved
colnames(ret_log) <- c(price_cols, yield_cols)
ret_log     <- na.omit(ret_log)

# Also keep simple returns for inflation adjustment step
ret_simple_prices <- ROC(prices_biz[, price_cols], type = "discrete")
dyields           <- diff(prices_biz[, yield_cols]) # still Δy
dyields           <- dyields[index(ret_simple_prices)]
ret_simple        <- merge(ret_simple_prices, dyields)
colnames(ret_simple) <- c(price_cols, yield_cols)
ret_simple        <- na.omit(ret_simple)

# ==========================================================
# 1.4 Rolling estimation (252d) + Monthly rebalancing (12 assets)
#     - Min-variance (long-only, fully invested)
#     - Trading-speed smoothing (f_speed)
# ==========================================================

assets_all <- colnames(ret_log)    # should be your 12 assets
R <- ret_log[, assets_all]         # daily "log-style" panel (prices: log ret; FVX/TNX: Δy)

# --- endpoints for monthly rebalancing on the NYSE calendar ---
ep <- endpoints(R, on = "months")
ep <- ep[ep > roll_win]            # start after we have 252d of history

# --- helper: min-variance weights (long-only, sum=1) ---
minvar_weights <- function(R_win) {
  n <- NCOL(R_win)
  S <- cov(R_win, use = "pairwise.complete.obs")
  # Regularize if ill-conditioned
  if (any(!is.finite(S)) || det(S) <= 1e-10) {
    S <- S + diag(1e-6, n)
  }
  Dmat <- as.matrix(2 * S)                 # quadprog uses 1/2 * x' D x
  dvec <- rep(0, n)
  Amat <- cbind(rep(1, n), diag(n))        # sum(w)=1 (equality), w>=0
  bvec <- c(1, rep(0, n))
  sol <- tryCatch(
    quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1),
    error = function(e) NULL
  )
  if (is.null(sol)) {
    w <- rep(1/n, n)                        # fallback: equal-weight
  } else {
    w <- pmax(sol$solution, 0)
    w <- w / sum(w)                         # normalize
  }
  setNames(w, colnames(R_win))
}

# --- containers for results ---
p_ret <- xts(rep(NA_real_, NROW(R)), order.by = index(R))   # daily portfolio log returns
W_hist <- xts(matrix(NA_real_, nrow = NROW(R), ncol = NCOL(R)),
              order.by = index(R))
colnames(W_hist) <- colnames(R)

# --- initialize previous weights (equal-weight) ---
w_prev <- rep(1/NCOL(R), NCOL(R))
names(w_prev) <- colnames(R)

for (k in seq_along(ep)) {
  i <- ep[k]
  if (!(i > roll_win)) next
  i_next <- if (k < length(ep)) ep[k + 1] else NROW(R)
  
  # estimation window: last 252 trading days ending at i
  win_idx <- (i - roll_win + 1):i
  R_win <- R[win_idx, , drop = FALSE]
  
  # target weights from optimizer
  w_tgt <- minvar_weights(R_win)
  
  # trading-speed smoothing at the rebalance date
  w_new <- (1 - f_speed) * w_prev + f_speed * w_tgt
  w_new <- pmax(w_new, 0)
  w_new <- w_new / sum(w_new)
  
  # apply fixed weights over the next holding period: (i+1) ... i_next
  if (i < i_next) {
    idx_hold <- (i + 1):i_next
    # daily portfolio log return = R %*% w_new
    p_ret[idx_hold] <- R[idx_hold, ] %*% w_new
    # record weights for each day in the holding period
    W_hist[idx_hold, ] <- matrix(rep(w_new, each = length(idx_hold)),
                                 nrow = length(idx_hold), byrow = FALSE,
                                 dimnames = list(index(R)[idx_hold], names(w_new)))
  }
  
  # update previous weights
  w_prev <- w_new
}

p_ret   <- na.omit(p_ret)
W_hist  <- na.omit(W_hist)

# positions in "percent units" (sum to 100) for the latest weights
positions_percent_latest <- M_notional * w_prev
names(positions_percent_latest) <- colnames(R)

# ==========================================================
# 1.5 Quick sanity checks + one clean cumulative plot
# ==========================================================

# Inspect weights/exposures
tail(W_hist, 1)                         # last weights (fractions, sum ~ 1)
round(positions_percent_latest, 2)      # last positions in percent units (sum ~ 100)

# Convert daily log -> simple once
p_ret_simple <- exp(p_ret) - 1

# Plot cumulative return (choose ONE path automatically)
if (requireNamespace("PerformanceAnalytics", quietly = TRUE)) {
  PerformanceAnalytics::chart.CumReturns(
    p_ret_simple, main = "Updated Cumulative Return", ylab = "Return"
  )
} else {
  eq <- exp(cumsum(p_ret)) - 1          # equity curve from log returns
  plot(eq, type = "l", lwd = 2,
       main = "Updated Cumulative Return", ylab = "Return", xlab = "")
  abline(h = 0, lty = 3)
}


## -----------------------------------------------------
## 2. Inflation Adjustment of Returns (CPI → daily π_t)
## -----------------------------------------------------
# CPI monthly (headline, SA)
getSymbols("CPIAUCSL", src="FRED", auto.assign=TRUE)
cpi_log_m <- log(CPIAUCSL)

# Interpolate log CPI to business days (only within CPI range), fill edges
idx_all   <- index(prices_biz)
idx_use   <- idx_all[idx_all >= first(index(cpi_log_m)) & idx_all <= last(index(cpi_log_m))]
cpi_log_use <- zoo::na.approx(zoo::zoo(as.numeric(cpi_log_m), order.by=index(cpi_log_m)),
                              x=index(cpi_log_m), xout=idx_use)
cpi_log_d <- xts(as.numeric(cpi_log_use), order.by=idx_use)
cpi_log_d <- cpi_log_d[index(prices_biz)]
cpi_log_d <- na.locf(cpi_log_d); cpi_log_d <- na.locf(cpi_log_d, fromLast=TRUE)

# Daily inflation
pi_log_d  <- diff(cpi_log_d)                              # log daily inflation
pi_daily  <- xts(exp(as.numeric(pi_log_d)) - 1, index(pi_log_d))  # simple daily inflation

# Real daily returns (exact): (1+R)/(1+π) - 1
pi_use <- pi_daily[index(ret_simple)]
ret_real_prices <- sweep(1 + ret_simple[, price_cols], 1, 1 + pi_use, `/`) - 1
ret_real <- merge(ret_real_prices, ret_simple[, yield_cols])  # keep Δy for yields
colnames(ret_real) <- c(price_cols, yield_cols)
ret_real <- na.omit(ret_real)
## -----------------------------------------------------
## 2.a Inflation-hedge diagnostics (before optimization)
## -----------------------------------------------------
inflation_xts <- pi_daily

# Work only with price-like assets (exclude yields FVX/TNX)
yield_cols  <- c("FVX","TNX")
price_cols  <- setdiff(colnames(ret_simple), yield_cols)

infl_hedge_stats <- function(r_nom, infl){
  common <- merge(r_nom, infl, join = "inner") %>% na.omit()
  r  <- as.numeric(common[,1])
  pi <- as.numeric(common[,2])
  
  if (length(r) < 30) {
    return(tibble(
      corr = NA_real_, beta = NA_real_, beta_p = NA_real_,
      mean_real = NA_real_, sd_real = NA_real_
    ))
  }
  
  fit <- lm(r ~ pi)
  tibble(
    corr      = cor(r, pi),
    beta      = unname(coef(fit)[2]),
    beta_p    = summary(fit)$coefficients[2,4],
    mean_real = mean((1 + r)/(1 + pi) - 1),
    sd_real   = sd(  (1 + r)/(1 + pi) - 1)
  )
}

infl_table <- map_dfr(price_cols, function(nm) {
  stats <- infl_hedge_stats(ret_simple[, nm], inflation_xts)
  mutate(stats, asset = nm, .before = 1)
})

print(infl_table)

## -----------------------------------------------------
## 3. Asset allocation models 
## -----------------------------------------------------

# ==========================================================
# 3.1) Equally Weighted (EW) Portfolio — monthly rebalance
#     - w_i = 1/N at each rebalance
#     - same trading-speed smoothing: f_speed
# ==========================================================

assets_all_ew <- colnames(ret_log)      # all 12 assets
R_ew <- ret_log[, assets_all_ew]

# monthly endpoints on NYSE calendar (same as before)
ep_ew <- endpoints(R_ew, on = "months")
ep_ew <- ep_ew[ep_ew > roll_win]        # start after 252d history

# containers (separate from min-var)
p_ret_ew  <- xts(rep(NA_real_, NROW(R_ew)), order.by = index(R_ew))  # daily portfolio log returns
W_hist_ew <- xts(matrix(NA_real_, nrow = NROW(R_ew), ncol = NCOL(R_ew)),
                 order.by = index(R_ew))
colnames(W_hist_ew) <- colnames(R_ew)

# initialize previous weights (equal-weight to start)
n_ew   <- NCOL(R_ew)
w_prev_ew <- rep(1 / n_ew, n_ew)
names(w_prev_ew) <- colnames(R_ew)

for (k in seq_along(ep_ew)) {
  i <- ep_ew[k]
  if (!(i > roll_win)) next
  i_next <- if (k < length(ep_ew)) ep_ew[k + 1] else NROW(R_ew)
  
  # target weights: equal-weight every rebalance
  w_tgt_ew <- rep(1 / n_ew, n_ew)
  names(w_tgt_ew) <- colnames(R_ew)
  
  # trading-speed smoothing at the rebalance date
  w_new_ew <- (1 - f_speed) * w_prev_ew + f_speed * w_tgt_ew
  w_new_ew <- pmax(w_new_ew, 0)
  w_new_ew <- w_new_ew / sum(w_new_ew)
  
  # hold weights over next period: (i+1) .. i_next
  if (i < i_next) {
    idx_hold <- (i + 1):i_next
    p_ret_ew[idx_hold] <- R_ew[idx_hold, ] %*% w_new_ew
    W_hist_ew[idx_hold, ] <- matrix(rep(w_new_ew, each = length(idx_hold)),
                                    nrow = length(idx_hold), byrow = FALSE,
                                    dimnames = list(index(R_ew)[idx_hold], names(w_new_ew)))
  }
  
  w_prev_ew <- w_new_ew
}

p_ret_ew  <- na.omit(p_ret_ew)
W_hist_ew <- na.omit(W_hist_ew)

# positions in percent units for latest EW portfolio
positions_percent_latest_ew <- M_notional * w_prev_ew
names(positions_percent_latest_ew) <- colnames(R_ew)

# Quick checks + plot (EW)

# weights & positions
print(tail(W_hist_ew, 1))                     # last weights (fractions, ~equal)
print(round(positions_percent_latest_ew, 2))  # last positions in percent (sum ~ 100)

# cumulative return (convert log -> simple )
p_ret_simple_ew <- exp(p_ret_ew) - 1

if (requireNamespace("PerformanceAnalytics", quietly = TRUE)) {
  PerformanceAnalytics::chart.CumReturns(
    p_ret_simple_ew, main = "Updated EW Portfolio — Cumulative Return", ylab = "Return"
  )
} else {
  eq_ew <- exp(cumsum(p_ret_ew)) - 1
  plot(eq_ew, type = "l", lwd = 2,
       main = "Updated EW Portfolio — Cumulative Return", ylab = "Return", xlab = "")
  abline(h = 0, lty = 3)
}

# ==========================================================
# Mean–Variance (Target Return) Portfolio — monthly rebalance
#   min w'Σw  s.t. 1'w=1, μ'w=r_target, w>=0, caps, FX group cap
#   + converts FVX/TNX Δy to bond returns via duration
# ==========================================================
suppressPackageStartupMessages(library(quadprog))

# --- Helper: convert Yahoo yield levels to DECIMAL if needed
to_decimal_yield <- function(x_close_xts) {
  medv <- suppressWarnings(median(na.omit(as.numeric(x_close_xts))))
  if (is.finite(medv) && medv > 5)      x_close_xts / 1000  # 42.5 -> 0.0425
  else if (is.finite(medv) && medv > 0.2) x_close_xts / 100  # 4.25 -> 0.0425
  else                                    x_close_xts
}

# --- Build investable return panel: price assets + bond returns
yield_cols <- c("FVX","TNX")
price_cols <- setdiff(colnames(ret_log), yield_cols)

# (1) Δy in DECIMAL units from original price levels
yields_dec <- merge(
  to_decimal_yield(Cl(`FVX`)),
  to_decimal_yield(Cl(`TNX`))
)
colnames(yields_dec) <- yield_cols
yields_dec <- yields_dec[index(prices_biz)]              # align to NYSE days
yields_dec <- na.locf(yields_dec, na.rm = FALSE)
yields_dec <- na.omit(yields_dec)
dyields_dec <- diff(yields_dec)
dyields_dec <- dyields_dec[index(ret_log)]               # align to return dates

# (2) Approximate bond price returns:  r_bond ≈ -D * Δy
D_FVX <- 4.5   # ~5y note
D_TNX <- 8.5   # ~10y note
bond_ret <- dyields_dec
bond_ret[,"FVX"] <- -D_FVX * dyields_dec[,"FVX"]
bond_ret[,"TNX"] <- -D_TNX * dyields_dec[,"TNX"]

# (3) Investable panel (all price assets' log returns + bond returns)
R_inv <- merge(ret_log[, price_cols], bond_ret)
R_inv <- na.omit(R_inv)

# --- Settings
roll_win   <- 252
r_target_ann <- 0.08
r_target_day_fixed <- r_target_ann / 252      # used if not adapting per window
asset_cap  <- 0.20                            # per-asset cap (e.g., 20%)
fx_group_cap <- 0.40                          # total FX cap (EURUSD, USDJPY, GBPUSD, AUDUSD)

# --- Light covariance shrinkage (stabilize Σ)
shrink_cov <- function(S, alpha = 0.15) {
  v <- diag(S); Tgt <- diag(v, nrow = length(v))
  (1 - alpha) * S + alpha * Tgt
}

# --- Optimizer: MV with target r, caps, group cap
mv_weights <- function(R_win, r_target_day, cap = asset_cap, fx_cap = fx_group_cap) {
  n  <- NCOL(R_win)
  mu <- colMeans(R_win, na.rm = TRUE)
  S  <- cov(R_win, use = "pairwise.complete.obs")
  S  <- shrink_cov(S, 0.15)
  if (any(!is.finite(S)) || det(S) <= 1e-10) S <- S + diag(1e-6, n)
  
  # Clip target into feasible band under long-only (helps solver)
  r_use <- min(max(r_target_day, min(mu, na.rm=TRUE) + 1e-8),
               max(mu, na.rm=TRUE) - 1e-8)
  
  Dmat <- as.matrix(2 * S)
  dvec <- rep(0, n)
  
  # Identify FX names present (may include bonds now)
  fx_names <- intersect(colnames(R_win), c("EURUSD","USDJPY","GBPUSD","AUDUSD"))
  G <- rep(0, n); if (length(fx_names)) G[colnames(R_win) %in% fx_names] <- 1
  
  # Constraints for quadprog: t(Amat) %*% w >= bvec
  # 1) sum(w)=1 (eq), 2) mu'w=r_use (eq), 3) w>=0, 4) w<=cap, 5) FX group <= fx_cap
  Amat <- cbind(
    rep(1, n),       # sum = 1
    mu,              # target return
    diag(n),         # w >= 0
    -diag(n),        # w <= cap  ->  -w >= -cap
    -G               # FX group: -G w >= -fx_cap
  )
  bvec <- c(1, r_use, rep(0, n), rep(-cap, n), -fx_cap)
  meq  <- 2
  
  sol <- try(quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = meq), silent = TRUE)
  if (inherits(sol, "try-error") || anyNA(sol$solution)) {
    # fallback: min-var with caps/group cap
    Amat2 <- cbind(rep(1,n), diag(n), -diag(n), -G)
    bvec2 <- c(1,      rep(0,n), rep(-cap,n), -fx_cap)
    sol2 <- try(quadprog::solve.QP(Dmat, dvec, Amat2, bvec2, meq = 1), silent = TRUE)
    w <- if (inherits(sol2, "try-error")) rep(1/n, n) else pmax(sol2$solution, 0)
  } else {
    w <- pmax(sol$solution, 0)
  }
  w / sum(w)
}

# --- Rebalancing scaffold (monthly), with trading-speed smoothing
ep_mv <- endpoints(R_inv, on = "months"); ep_mv <- ep_mv[ep_mv > roll_win]

p_ret_mv  <- xts(rep(NA_real_, NROW(R_inv)), order.by = index(R_inv))
W_hist_mv <- xts(matrix(NA_real_, nrow = NROW(R_inv), ncol = NCOL(R_inv)),
                 order.by = index(R_inv))
colnames(W_hist_mv) <- colnames(R_inv)

w_prev_mv <- rep(1 / NCOL(R_inv), NCOL(R_inv)); names(w_prev_mv) <- colnames(R_inv)

for (k in seq_along(ep_mv)) {
  i <- ep_mv[k]; if (!(i > roll_win)) next
  i_next <- if (k < length(ep_mv)) ep_mv[k + 1] else NROW(R_inv)
  
  win_idx <- (i - roll_win + 1):i
  R_win <- R_inv[win_idx, , drop = FALSE]
  
  # Option: adaptive target (use 60th percentile of window means)
  mu_win <- colMeans(R_win, na.rm = TRUE)
  r_target_day <- quantile(mu_win, 0.60, na.rm = TRUE)
  
  w_tgt_mv <- mv_weights(R_win, r_target_day)
  
  # trading-speed smoothing
  w_new_mv <- (1 - f_speed) * w_prev_mv + f_speed * w_tgt_mv
  w_new_mv <- pmax(w_new_mv, 0); w_new_mv <- w_new_mv / sum(w_new_mv)
  
  if (i < i_next) {
    idx_hold <- (i + 1):i_next
    p_ret_mv[idx_hold]   <- R_inv[idx_hold, ] %*% w_new_mv
    W_hist_mv[idx_hold,] <- matrix(rep(w_new_mv, each = length(idx_hold)),
                                   nrow = length(idx_hold), byrow = FALSE,
                                   dimnames = list(index(R_inv)[idx_hold], names(w_new_mv)))
  }
  w_prev_mv <- w_new_mv
}

p_ret_mv  <- na.omit(p_ret_mv)
W_hist_mv <- na.omit(W_hist_mv)

# --- Latest positions in percent units
positions_percent_latest_mv <- M_notional * w_prev_mv
names(positions_percent_latest_mv) <- colnames(R_inv)

# --- Plot cumulative return
p_ret_simple_mv <- exp(p_ret_mv) - 1
if (requireNamespace("PerformanceAnalytics", quietly = TRUE)) {
  PerformanceAnalytics::chart.CumReturns(
    p_ret_simple_mv, main = "Updated Mean–Variance (Target Return) — Cumulative Return", ylab = "Return"
  )
} else {
  eq_mv <- exp(cumsum(p_ret_mv)) - 1
  plot(eq_mv, type = "l", lwd = 2,
       main = "Updated Mean–Variance (Target Return) — Cumulative Return", ylab = "Return", xlab = "")
  abline(h = 0, lty = 3)
}

# --- Quick printouts
print(tail(W_hist_mv, 1))                    # last weights (fractions)
print(round(positions_percent_latest_mv, 2)) # last positions in percent (sum ~ 100)

# ==========================================================
# 3.3) Maximum Sharpe (MV-S) Portfolio — monthly rebalance
#      max (mu'w) / sqrt(w'Σw), rf=0, with same constraints as MV
#      Implementation: search over target returns -> pick highest Sharpe
# ==========================================================
# --- Helpers (defined here if not already present) ---------------------------
if (!exists("to_decimal_yield")) {
  to_decimal_yield <- function(x_close_xts) {
    medv <- suppressWarnings(median(na.omit(as.numeric(x_close_xts))))
    if (is.finite(medv) && medv > 5)      x_close_xts / 1000  # 42.5 -> 0.0425
    else if (is.finite(medv) && medv > 0.2) x_close_xts / 100  # 4.25 -> 0.0425
    else                                    x_close_xts
  }
}
shrink_cov <- function(S, alpha = 0.15) {
  v <- diag(S); Tgt <- diag(v, nrow = length(v))
  (1 - alpha) * S + alpha * Tgt
}

# --- Build investable return panel (price assets + bond returns) if missing ---
if (!exists("R_inv")) {
  yield_cols <- c("FVX","TNX")
  price_cols <- setdiff(colnames(ret_log), yield_cols)
  
  yields_dec <- merge(to_decimal_yield(Cl(`FVX`)), to_decimal_yield(Cl(`TNX`)))
  colnames(yields_dec) <- yield_cols
  yields_dec <- yields_dec[index(prices_biz)]
  yields_dec <- na.locf(yields_dec, na.rm = FALSE)
  yields_dec <- na.omit(yields_dec)
  dyields_dec <- diff(yields_dec)
  dyields_dec <- dyields_dec[index(ret_log)]
  
  D_FVX <- 4.5; D_TNX <- 8.5
  bond_ret <- dyields_dec
  bond_ret[,"FVX"] <- -D_FVX * dyields_dec[,"FVX"]
  bond_ret[,"TNX"] <- -D_TNX * dyields_dec[,"TNX"]
  
  R_inv <- merge(ret_log[, price_cols], bond_ret)
  R_inv <- na.omit(R_inv)
}

# --- Constraints & knobs (same spirit as MV) ---------------------------------
asset_cap    <- 0.20    # per-asset cap (e.g., 20%)
fx_group_cap <- 0.40    # cap on EURUSD+USDJPY+GBPUSD+AUDUSD total
roll_win     <- 252     # estimation window
rf_daily     <- 0       # Sharpe uses Rf = 0 (per your spec)

# Solve min-variance for a given daily target return r_use with caps + FX group cap
solve_minvar_at_r <- function(mu, S, r_use, cnames, cap = asset_cap, fx_cap = fx_group_cap) {
  n <- length(mu)
  S <- shrink_cov(S, 0.15)
  if (any(!is.finite(S)) || det(S) <= 1e-10) S <- S + diag(1e-6, n)
  
  Dmat <- as.matrix(2 * S)
  dvec <- rep(0, n)
  
  fx_names <- intersect(cnames, c("EURUSD","USDJPY","GBPUSD","AUDUSD"))
  G <- rep(0, n); if (length(fx_names)) G[cnames %in% fx_names] <- 1
  
  # t(Amat) %*% w >= bvec
  Amat <- cbind(
    rep(1, n),    # sum=1 (eq)
    mu,           # mu'w=r (eq)
    diag(n),      # w >= 0
    -diag(n),     # w <= cap  -> -w >= -cap
    -G            # FX group: -G w >= -fx_cap  (i.e., G w <= fx_cap)
  )
  bvec <- c(1, r_use, rep(0, n), rep(-cap, n), -fx_cap)
  meq  <- 2
  
  sol <- try(quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = meq), silent = TRUE)
  if (inherits(sol, "try-error") || anyNA(sol$solution)) return(NULL)
  w <- pmax(sol$solution, 0); w <- w / sum(w)
  w
}

# Choose w that maximizes in-window Sharpe among a grid of target returns
maxsharpe_weights <- function(R_win) {
  cn <- colnames(R_win)
  mu <- colMeans(R_win, na.rm = TRUE)
  S  <- cov(R_win, use = "pairwise.complete.obs")
  
  # Clip feasible band (avoid impossible targets under long-only)
  rmin <- min(mu, na.rm = TRUE) + 1e-8
  rmax <- max(mu, na.rm = TRUE) - 1e-8
  if (!is.finite(rmin) || !is.finite(rmax) || rmin >= rmax) {
    return(rep(1/NCOL(R_win), NCOL(R_win)))  # fallback equal-weight
  }
  
  # Target grid (inside feasible band): mid-to-upper portion of mu spectrum
  r_grid <- seq(rmin + 0.15 * (rmax - rmin), rmin + 0.95 * (rmax - rmin), length.out = 15)
  
  best_w <- NULL; best_sr <- -Inf
  for (rt in r_grid) {
    w <- solve_minvar_at_r(mu, S, rt, cn)
    if (is.null(w)) next
    varp <- as.numeric(t(w) %*% S %*% w)
    if (!is.finite(varp) || varp <= 0) next
    retp <- as.numeric(sum(mu * w))
    sr   <- (retp - rf_daily) / sqrt(varp)
    if (is.finite(sr) && sr > best_sr) { best_sr <- sr; best_w <- w }
  }
  if (is.null(best_w)) best_w <- rep(1/NCOL(R_win), NCOL(R_win))
  names(best_w) <- cn
  best_w
}

# --- Rebalancing (monthly) with trading-speed smoothing ----------------------
ep_mvs <- endpoints(R_inv, on = "months"); ep_mvs <- ep_mvs[ep_mvs > roll_win]

p_ret_mvs  <- xts(rep(NA_real_, NROW(R_inv)), order.by = index(R_inv))
W_hist_mvs <- xts(matrix(NA_real_, nrow = NROW(R_inv), ncol = NCOL(R_inv)),
                  order.by = index(R_inv))
colnames(W_hist_mvs) <- colnames(R_inv)

w_prev_mvs <- rep(1 / NCOL(R_inv), NCOL(R_inv)); names(w_prev_mvs) <- colnames(R_inv)

for (k in seq_along(ep_mvs)) {
  i <- ep_mvs[k]; if (!(i > roll_win)) next
  i_next <- if (k < length(ep_mvs)) ep_mvs[k + 1] else NROW(R_inv)
  
  win_idx <- (i - roll_win + 1):i
  R_win <- R_inv[win_idx, , drop = FALSE]
  
  w_tgt_mvs <- maxsharpe_weights(R_win)
  
  # trading-speed smoothing
  w_new_mvs <- (1 - f_speed) * w_prev_mvs + f_speed * w_tgt_mvs
  w_new_mvs <- pmax(w_new_mvs, 0); w_new_mvs <- w_new_mvs / sum(w_new_mvs)
  
  if (i < i_next) {
    idx_hold <- (i + 1):i_next
    p_ret_mvs[idx_hold]   <- R_inv[idx_hold, ] %*% w_new_mvs
    W_hist_mvs[idx_hold,] <- matrix(rep(w_new_mvs, each = length(idx_hold)),
                                    nrow = length(idx_hold), byrow = FALSE,
                                    dimnames = list(index(R_inv)[idx_hold], names(w_new_mvs)))
  }
  w_prev_mvs <- w_new_mvs
}

p_ret_mvs  <- na.omit(p_ret_mvs)
W_hist_mvs <- na.omit(W_hist_mvs)

# --- Latest positions in percent units & plot --------------------------------
positions_percent_latest_mvs <- M_notional * w_prev_mvs
names(positions_percent_latest_mvs) <- colnames(R_inv)
print(tail(W_hist_mvs, 1))
print(round(positions_percent_latest_mvs, 2))

p_ret_simple_mvs <- exp(p_ret_mvs) - 1
if (requireNamespace("PerformanceAnalytics", quietly = TRUE)) {
  PerformanceAnalytics::chart.CumReturns(
    p_ret_simple_mvs, main = "Updated Maximum Sharpe (MV-S) — Cumulative Return", ylab = "Return"
  )
} else {
  eq_mvs <- exp(cumsum(p_ret_mvs)) - 1
  plot(eq_mvs, type = "l", lwd = 2,
       main = "Updated Maximum Sharpe (MV-S) — Cumulative Return", ylab = "Return", xlab = "")
  abline(h = 0, lty = 3)
}
# ==========================================================
# 3.4) Mean–CVaR Portfolio (Rockafellar–Uryasev LP)
#      min_zeta,u  zeta + (1/((1-alpha)*T)) * sum(u_t)
#      s.t. u_t >= -w'R_t - zeta, u_t >= 0
#           1'w = 1,  w >= 0,  μ'w = r_target_day
#      Rolling 252d window, monthly rebalancing, smoothing f_speed
# ==========================================================
# --- choose returns panel ---
R_cvar <- if (exists("R_inv")) R_inv else ret_log[, setdiff(colnames(ret_log), c("FVX","TNX"))]
R_cvar <- na.omit(R_cvar)

alpha_cvar <- 0.95           # tail probability
roll_win   <- 252            # estimation window
rf_daily   <- 0              # kept at zero (as in your spec)

# --- CVaR LP solver (long-only, fully invested, target daily return) ----------
mean_cvar_weights <- function(R_win, alpha = 0.95, r_target_day = NULL) {
  # R_win: T x N (daily returns in the window)
  Tn <- NROW(R_win); n <- NCOL(R_win)
  mu <- colMeans(R_win, na.rm = TRUE)
  
  # If no target given, set a feasible adaptive target (e.g., 60th pct of means)
  if (is.null(r_target_day)) r_target_day <- quantile(mu, 0.60, na.rm = TRUE)
  
  # Variables: [w (n), zeta (1), u (Tn)]  -> total n + 1 + Tn
  # Objective: minimize c' x  with  c = [0_n, 1, rep(1/((1-alpha)*Tn), Tn)]
  cvec <- c(rep(0, n), 1, rep(1/((1 - alpha) * Tn), Tn))
  
  # Constraints list we’ll rbind:
  A <- NULL; dir <- character(); rhs <- numeric()
  
  # (1) Full investment: sum(w) = 1
  A <- rbind(A, c(rep(1, n), 0, rep(0, Tn))); dir <- c(dir, "=="); rhs <- c(rhs, 1)
  
  # (2) Target return: mu' w = r_target_day
  A <- rbind(A, c(mu, 0, rep(0, Tn)));       dir <- c(dir, "=="); rhs <- c(rhs, r_target_day)
  
  # (3) Scenario constraints: u_t >= -w'R_t - zeta  =>  u_t + zeta + w'R_t >= 0
  Rt <- as.matrix(R_win)
  for (t in 1:Tn) {
    row <- c(Rt[t, ], 1, rep(0, Tn))     # w'R_t + zeta
    row[n + 1 + t] <- row[n + 1 + t] + 1 # + u_t
    A <- rbind(A, row); dir <- c(dir, ">="); rhs <- c(rhs, 0)
  }
  
  # (4) u_t >= 0  (Tn rows)
  for (t in 1:Tn) {
    row <- rep(0, n + 1 + Tn); row[n + 1 + t] <- 1
    A <- rbind(A, row); dir <- c(dir, ">="); rhs <- c(rhs, 0)
  }
  
  # (5) Long-only: w_i >= 0  (n rows)
  for (i in 1:n) {
    row <- rep(0, n + 1 + Tn); row[i] <- 1
    A <- rbind(A, row); dir <- c(dir, ">="); rhs <- c(rhs, 0)
  }
  
  op <- OP(objective = L_objective(cvec),
           constraints = L_constraint(L = A, dir = dir, rhs = rhs),
           maximum = FALSE)
  
  sol <- tryCatch(ROI_solve(op, solver = "glpk"), error = function(e) e)
  if (inherits(sol, "error") || sol$status$code != 0) {
    # Fallback: equal-weight if LP infeasible
    w <- rep(1/n, n)
  } else {
    x <- solution(sol)
    w <- x[1:n]
    # normalize defensively
    w <- pmax(w, 0); w <- w / sum(w)
  }
  names(w) <- colnames(R_win)
  w
}

# --- Monthly rebalancing with trading-speed smoothing -------------------------
ep_cvar <- endpoints(R_cvar, on = "months"); ep_cvar <- ep_cvar[ep_cvar > roll_win]

p_ret_cvar  <- xts(rep(NA_real_, NROW(R_cvar)), order.by = index(R_cvar))
W_hist_cvar <- xts(matrix(NA_real_, nrow = NROW(R_cvar), ncol = NCOL(R_cvar)),
                   order.by = index(R_cvar))
colnames(W_hist_cvar) <- colnames(R_cvar)

w_prev_cvar <- rep(1 / NCOL(R_cvar), NCOL(R_cvar)); names(w_prev_cvar) <- colnames(R_cvar)

for (k in seq_along(ep_cvar)) {
  i <- ep_cvar[k]; if (!(i > roll_win)) next
  i_next <- if (k < length(ep_cvar)) ep_cvar[k + 1] else NROW(R_cvar)
  
  win_idx <- (i - roll_win + 1):i
  R_win <- R_cvar[win_idx, , drop = FALSE]
  
  # adaptive target (same as MV/MV-S style): 60th pct of window means
  mu_win <- colMeans(R_win, na.rm = TRUE)
  r_target_day <- quantile(mu_win, 0.60, na.rm = TRUE)
  
  w_tgt <- mean_cvar_weights(R_win, alpha = alpha_cvar, r_target_day = r_target_day)
  
  # trading-speed smoothing
  w_new <- (1 - f_speed) * w_prev_cvar + f_speed * w_tgt
  w_new <- pmax(w_new, 0); w_new <- w_new / sum(w_new)
  
  if (i < i_next) {
    idx_hold <- (i + 1):i_next
    p_ret_cvar[idx_hold]   <- R_cvar[idx_hold, ] %*% w_new
    W_hist_cvar[idx_hold,] <- matrix(rep(w_new, each = length(idx_hold)),
                                     nrow = length(idx_hold), byrow = FALSE,
                                     dimnames = list(index(R_cvar)[idx_hold], names(w_new)))
  }
  w_prev_cvar <- w_new
}

p_ret_cvar  <- na.omit(p_ret_cvar)
W_hist_cvar <- na.omit(W_hist_cvar)

# --- Outputs & plot -----------------------------------------------------------
positions_percent_latest_cvar <- M_notional * w_prev_cvar
names(positions_percent_latest_cvar) <- colnames(R_cvar)
print(tail(W_hist_cvar, 1))
print(round(positions_percent_latest_cvar, 2))

p_ret_simple_cvar <- exp(p_ret_cvar) - 1
if (requireNamespace("PerformanceAnalytics", quietly = TRUE)) {
  PerformanceAnalytics::chart.CumReturns(
    p_ret_simple_cvar,
    main = sprintf("Updated Mean–CVaR (α=%.2f) — Cumulative Return", alpha_cvar),
    ylab = "Return"
  )
} else {
  eq_cvar <- exp(cumsum(p_ret_cvar)) - 1
  plot(eq_cvar, type = "l", lwd = 2,
       main = sprintf("Updated Mean–CVaR (α=%.2f) — Cumulative Return", alpha_cvar),
       ylab = "Return", xlab = "")
  abline(h = 0, lty = 3)
}


## =========================
## ## =========================
## 4. Performance evaluation
## =========================
## =========================

## =========================
## 4.1 Cumulative Wealth (CW)  — fixes included
## =========================

# Convert strategy log-returns -> simple returns once
r_ew_s   <- exp(p_ret_ew)   - 1
r_mv_s   <- exp(p_ret_mv)   - 1
r_mvs_s  <- exp(p_ret_mvs)  - 1
r_cvar_s <- exp(p_ret_cvar) - 1

# Helper: cumulative wealth from simple returns (W0 = 1)
cw_from_simple <- function(r_simple_xts, W0 = 1) cumprod(1 + r_simple_xts) * W0

# Wealth indices
CW_EW   <- cw_from_simple(r_ew_s);   colnames(CW_EW)   <- "EW"
CW_MV   <- cw_from_simple(r_mv_s);   colnames(CW_MV)   <- "MV"
CW_MVS  <- cw_from_simple(r_mvs_s);  colnames(CW_MVS)  <- "MV-S"
CW_CVaR <- cw_from_simple(r_cvar_s); colnames(CW_CVaR) <- "Mean-CVaR"

# Safe inner-merge for multiple xts
merge_xts_inner <- function(...) Reduce(function(a, b) merge(a, b, join = "inner"),
                                        list(...))
CW_ALL <- merge_xts_inner(CW_EW, CW_MV, CW_MVS, CW_CVaR)

# Plot
if (requireNamespace("PerformanceAnalytics", quietly = TRUE)) {
  PerformanceAnalytics::chart.TimeSeries(
    CW_ALL, main = "Updated Cumulative Wealth (W0 = 1)", ylab = "Wealth",
    legend.loc = "topleft", colorset = 1:NCOL(CW_ALL)
  )
} else {
  plot(CW_ALL[,1], type="l", lwd=2, ylim=range(CW_ALL, na.rm=TRUE),
       main="Updated Cumulative Wealth (W0 = 1)", ylab="Wealth", xlab="")
  for (j in 2:NCOL(CW_ALL)) lines(CW_ALL[,j], lwd=2)
  abline(h = 1, lty = 3)
  legend("topleft", legend = colnames(CW_ALL), lwd = 2, col = 1:NCOL(CW_ALL), bty="n")
}

## =========================
## 4.2 Sharpe Ratio (annualized, Rf = 0) — fixes included
## =========================

ann_factor <- 252

sum_stats <- function(r_simple_xts, cw_xts, rf_daily = 0) {
  n      <- NROW(r_simple_xts)
  mu_d   <- mean(r_simple_xts, na.rm = TRUE) - rf_daily
  sd_d   <- sd(r_simple_xts,   na.rm = TRUE)
  sr_ann <- if (sd_d > 0) (mu_d / sd_d) * sqrt(ann_factor) else NA_real_
  
  last_w <- as.numeric(xts::last(cw_xts))   # <— use xts::last
  
  tibble::tibble(
    Final_Wealth = last_w,
    CAGR         = (last_w^(ann_factor / n)) - 1,
    Vol_Ann      = sd_d * sqrt(ann_factor),
    Sharpe       = sr_ann,
    MaxDD        = as.numeric(PerformanceAnalytics::maxDrawdown(r_simple_xts))
  )
}

stats_tbl <- dplyr::bind_rows(
  EW          = sum_stats(r_ew_s,   CW_EW),
  MV          = sum_stats(r_mv_s,   CW_MV),
  `MV-S`      = sum_stats(r_mvs_s,  CW_MVS),
  `Mean-CVaR` = sum_stats(r_cvar_s, CW_CVaR),
  .id = "Strategy"
) |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~round(.x, 4)))

print(stats_tbl)

## =========================
## 4.3 Certainty Equivalent (CEQ)
## =========================
# CEQ_d = mu_d − (gamma/2)*var_d ; CEQ_ann = 252 * CEQ_d  (iid approx)
ann_factor <- 252
gamma_ceq  <- 1   # per your spec

ceq_ann <- function(r_simple_xts, gamma = gamma_ceq, scale = ann_factor) {
  mu_d  <- mean(r_simple_xts, na.rm = TRUE)
  var_d <- var(r_simple_xts,  na.rm = TRUE)
  ceq_d <- mu_d - 0.5 * gamma * var_d
  scale * ceq_d
}

ceq_tbl <- tibble::tibble(
  Strategy = c("EW","MV","MV-S","Mean-CVaR"),
  CEQ_Ann  = c(
    ceq_ann(r_ew_s),
    ceq_ann(r_mv_s),
    ceq_ann(r_mvs_s),
    ceq_ann(r_cvar_s)
  )
) |> dplyr::mutate(CEQ_Ann = round(CEQ_Ann, 4))

print(ceq_tbl)

# (optional) CEQ bar chart
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::ggplot(ceq_tbl, ggplot2::aes(x = Strategy, y = CEQ_Ann)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = round(CEQ_Ann, 3)),
                       vjust = -0.4, size = 3.5) +
    ggplot2::labs(title = "Certainty-Equivalent Return (annualized, γ = 1)",
                  x = "", y = "CEQ (per year)") +
    ggplot2::theme_minimal()
}

## =========================
## 4.4 Adjusted Sharpe Ratio (ASR)
## =========================
# ASR_d = SR_d * (1 + (S/6)*SR_d − (K/24)*SR_d^2), where
#   SR_d uses daily mean/sd with Rf=0, S is skewness, K is excess kurtosis.
# We scale to annual units as ASR_ann ≈ ASR_d * sqrt(252) to match 4.2.
adj_sharpe_ann <- function(r_simple_xts, scale = ann_factor) {
  r <- as.numeric(r_simple_xts)
  r <- r[is.finite(r)]
  if (length(r) < 3) return(list(ASR_Ann = NA_real_, SR_Ann = NA_real_, S = NA_real_, K = NA_real_))
  mu  <- mean(r); sdv <- stats::sd(r)
  if (!is.finite(sdv) || sdv == 0) return(list(ASR_Ann = NA_real_, SR_Ann = NA_real_, S = NA_real_, K = NA_real_))
  
  # daily Sharpe (Rf=0)
  sr_d <- mu / sdv
  
  # sample skewness and EXCESS kurtosis (unbiasedness not critical for large n)
  z   <- (r - mu) / sdv
  S   <- mean(z^3)                          # skewness
  Kex <- mean(z^4) - 3                      # excess kurtosis
  
  # adjusted Sharpe at daily freq
  asr_d <- sr_d * (1 + (S/6) * sr_d - (Kex/24) * sr_d^2)
  
  # scale to annual to be comparable with your 4.2 Sharpe
  list(
    ASR_Ann = asr_d * sqrt(scale),
    SR_Ann  = (sr_d * sqrt(scale)),
    S = S,
    K = Kex
  )
}

asr_list <- list(
  EW        = adj_sharpe_ann(r_ew_s),
  MV        = adj_sharpe_ann(r_mv_s),
  `MV-S`    = adj_sharpe_ann(r_mvs_s),
  `Mean-CVaR` = adj_sharpe_ann(r_cvar_s)
)

asr_tbl <- dplyr::bind_rows(lapply(names(asr_list), function(nm) {
  tibble::tibble(
    Strategy = nm,
    ASR_Ann  = asr_list[[nm]]$ASR_Ann,
    Sharpe_Ann = asr_list[[nm]]$SR_Ann,
    Skewness = asr_list[[nm]]$S,
    KurtosisEx = asr_list[[nm]]$K
  )
})) |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~round(.x, 4)))

print(asr_tbl)

# (optional) ASR bar chart
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::ggplot(asr_tbl, ggplot2::aes(x = Strategy, y = ASR_Ann)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = round(ASR_Ann, 2)),
                       vjust = -0.4, size = 3.5) +
    ggplot2::labs(title = "Adjusted Sharpe Ratio (annualized)",
                  x = "", y = "ASR (per year)") +
    ggplot2::theme_minimal()
}

## =========================
## 4.5 Diversification Measures (Neff, DR)
## =========================

# Helper: compute Neff and DR from a return panel and a weight history
div_measures <- function(R_panel, W_hist, win = 252) {
  stopifnot(NCOL(R_panel) == NCOL(W_hist))
  last_d <- xts::last(index(W_hist))
  w_last <- as.numeric(W_hist[last_d, ])
  names(w_last) <- colnames(W_hist)
  
  # Align return panel to weights & take last 'win' days up to last_d
  R_use <- R_panel[, colnames(W_hist)]
  R_use <- R_use[index(R_use) <= last_d]
  R_win <- tail(na.omit(R_use), win)
  
  # Individual vols and covariance
  sig_i <- apply(R_win, 2, sd, na.rm = TRUE)
  S     <- cov(R_win, use = "pairwise.complete.obs")
  
  # Portfolio vol
  sig_p <- as.numeric(sqrt(t(w_last) %*% S %*% w_last))
  
  # Metrics
  neff <- 1 / sum(w_last^2)
  dr   <- sum(w_last * sig_i) / sig_p
  
  tibble::tibble(
    Neff = as.numeric(neff),
    DR   = as.numeric(dr)
  )
}

# Compute for each strategy
div_ew   <- div_measures(R_ew,   W_hist_ew,  win = roll_win) |> dplyr::mutate(Strategy = "EW", .before = 1)
div_mv   <- div_measures(R_inv,  W_hist_mv,  win = roll_win) |> dplyr::mutate(Strategy = "MV", .before = 1)
div_mvs  <- div_measures(R_inv,  W_hist_mvs, win = roll_win) |> dplyr::mutate(Strategy = "MV-S", .before = 1)
div_cvar <- div_measures(R_cvar, W_hist_cvar,win = roll_win) |> dplyr::mutate(Strategy = "Mean-CVaR", .before = 1)

div_tbl <- dplyr::bind_rows(div_ew, div_mv, div_mvs, div_cvar) |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 3)))

print(div_tbl)

# (optional) quick bars
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::ggplot(div_tbl, ggplot2::aes(x = Strategy, y = Neff)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = Neff), vjust = -0.4) +
    ggplot2::labs(title = "Effective Number of Assets (latest weights)", x = "", y = "N_eff") +
    ggplot2::theme_minimal()
  
  ggplot2::ggplot(div_tbl, ggplot2::aes(x = Strategy, y = DR)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = DR), vjust = -0.4) +
    ggplot2::labs(title = "Diversification Ratio (252d window)", x = "", y = "DR") +
    ggplot2::theme_minimal()
}

## =========================
## 4.6 Turnover & Transaction-cost estimates
## =========================
# Turnover at rebalance t: TO_t = sum_i | w_{i,t+1}^{post} - w_{i,t+1}^{pre} |
# where w^{pre} is the *drifted* weight vector right before trading at t+1,
# and w^{post} is the new target weight applied for the next period.

# Helper: compute drifted weights over a holding period [s_idx, e_idx]
# given the post-trade weights at s_idx and asset returns in that interval.
drift_to_end <- function(R_panel, s_idx, e_idx, w_start, cols){
  if (e_idx < s_idx) return(w_start)
  gross <- apply(1 + R_panel[s_idx:e_idx, cols, drop = FALSE], 2, prod, na.rm = TRUE)
  w_pre <- as.numeric(w_start * gross)
  w_pre / sum(w_pre)
}

# Core function: returns per-rebalance turnover vector and a summary tibble
turnover_for_strategy <- function(R_panel, W_hist, ep_vec){
  stopifnot(all(colnames(R_panel) %in% colnames(W_hist)))
  cols <- colnames(W_hist)
  datesR <- index(R_panel); datesW <- index(W_hist)
  
  TO_vec <- c(); when <- c()
  
  # We use the same monthly endpoints you used to rebalance (after the 252d burn-in).
  for (k in seq_len(length(ep_vec) - 1)) {
    i      <- ep_vec[k]
    i_next <- ep_vec[k + 1]
    # holding period is (i+1) ... i_next
    s_idx <- i + 1
    e_idx <- i_next
    
    if (s_idx > NROW(R_panel) || e_idx > NROW(R_panel)) next
    if (e_idx >= NROW(W_hist)) next                      # need a *next* period to compare
    
    # post-trade weights at start of (i+1) (these are constant in your W_hist for the period)
    w_post_t   <- as.numeric(W_hist[s_idx, cols])
    # drift them to just before the next trade (end of holding period)
    w_pre_next <- drift_to_end(R_panel, s_idx, e_idx, w_post_t, cols)
    
    # next period's post-trade weights (applied at i_next + 1)
    if (e_idx + 1 > NROW(W_hist)) next
    w_post_next <- as.numeric(W_hist[e_idx + 1, cols])
    
    # turnover at this rebalance
    TO_t <- sum(abs(w_post_next - w_pre_next))
    TO_vec <- c(TO_vec, TO_t)
    when   <- c(when, datesW[e_idx + 1])
  }
  
  tibble::tibble(RebalanceDate = when, Turnover = TO_vec)
}

# Compute per-rebalance turnover series for each portfolio
to_ew   <- turnover_for_strategy(R_ew,   W_hist_ew,  ep_ew)   |> dplyr::mutate(Strategy = "EW")
to_mv   <- turnover_for_strategy(R_inv,  W_hist_mv,  ep_mv)   |> dplyr::mutate(Strategy = "MV")
to_mvs  <- turnover_for_strategy(R_inv,  W_hist_mvs, ep_mvs)  |> dplyr::mutate(Strategy = "MV-S")
to_cvar <- turnover_for_strategy(R_cvar, W_hist_cvar,ep_cvar) |> dplyr::mutate(Strategy = "Mean-CVaR")

turnover_all <- dplyr::bind_rows(to_ew, to_mv, to_mvs, to_cvar)

# Summary (your formula average over rebalances; T-K equals number of rebalances used)
turnover_summary <- turnover_all |>
  dplyr::group_by(Strategy) |>
  dplyr::summarise(
    Rebalances = dplyr::n(),
    TO_mean    = mean(Turnover, na.rm = TRUE),
    TO_median  = median(Turnover, na.rm = TRUE),
    TO_p90     = quantile(Turnover, 0.90, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~round(.x, 4)))

print(turnover_summary)

# ---- Optional: translate turnover into transaction costs ---------------------
# Assume a one-way trading cost of 'tc_bps' (in basis points of notional traded).
# Cost per rebalance = tc * sum_i |Δw_i| (since |Δw_i| is fraction of portfolio traded).
# With monthly rebal, annualized cost ≈ 12 * tc * TO_mean.
tc_bps <- 5   # <-- set your assumed one-way cost here (e.g., 5 bps = 0.05%)
tc      <- tc_bps / 1e4

tc_summary <- turnover_summary |>
  dplyr::mutate(
    Cost_per_rebalance = round(tc * TO_mean, 6),
    Cost_per_year      = round(12 * tc * TO_mean, 6)
  )

print(tc_summary)

# (optional) bar chart of average turnover per rebalance
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::ggplot(turnover_summary, ggplot2::aes(x = Strategy, y = TO_mean)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = round(TO_mean, 3)), vjust = -0.4) +
    ggplot2::labs(title = "Average Turnover per Rebalance (L1, drift-adjusted)",
                  x = "", y = "TO (sum |Δw|)") +
    ggplot2::theme_minimal()
}
