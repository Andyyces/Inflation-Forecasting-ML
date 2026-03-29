# Create the output folder first, otherwise the save step will fail
# dir.create(here("Output", "Forecasts"), recursive = TRUE)

# Test on first 3 origins only before committing to the full run
n_origins_test <- 3

# =============================================================================
# 01_ucsv.R  [FIXED - numerical stability + correct yoy scaling]
# UCSV Model: Estimation and Recursive Forecasting
# Translation of Chan (2018), Econometric Reviews, 37(8), 807-823
# Non-centred parameterisation of Stock & Watson (2007)
# As used in Banbura & Bobeica (2023)
# =============================================================================

library(here)
library(readxl)
library(dplyr)
library(lubridate)

# =============================================================================
# SECTION 1: LOAD AND PREPARE DATA
# =============================================================================

panel <- read_excel(here("Data", "EA-MD-QD", "EA-MD-QD-2026-02",
                         "data_TR2", "eadataM_NA_TR2.xlsx"))
raw   <- read_excel(here("Data", "EA-MD-QD", "EA-MD-QD-2026-02", "EAdata.xlsx"),
                    sheet = "data")

panel <- panel %>% mutate(Time = as.Date(Time))
raw   <- raw   %>% mutate(Time = as.Date(Time))

start_date <- as.Date("2000-01-01")
panel <- panel %>% filter(Time >= start_date)
raw   <- raw   %>% filter(Time >= start_date)

dates  <- panel$Time
T_full <- nrow(panel)

cat("Panel dimensions:", T_full, "x", ncol(panel) - 1, "\n")
cat("Date range:", format(min(dates)), "to", format(max(dates)), "\n")

# --- Construct inflation series from raw HICP levels ---
core_level     <- as.numeric(raw$HICPNEF_EA)
headline_level <- as.numeric(raw$HICPOV_EA)

# Month-on-month annualised x1200: INPUT to all models (UCSV, DFM, ML)
pi_core_mom <- c(NA, 1200 * diff(log(core_level)))
pi_hl_mom   <- c(NA, 1200 * diff(log(headline_level)))

# Year-on-year x100: EVALUATION TARGET for RMSE tables
# Note: x100 only, NOT x1200 -- yoy diff already covers 12 months
pi_core_yoy <- c(rep(NA, 12), 100 * diff(log(core_level),     lag = 12))
pi_hl_yoy   <- c(rep(NA, 12), 100 * diff(log(headline_level), lag = 12))

# Sanity check: recent values should be ~1-3% for core HICP
cat("\nRecent core HICP (m/m annualised and y/y):\n")
print(tail(data.frame(date = dates,
                      mom  = round(pi_core_mom, 2),
                      yoy  = round(pi_core_yoy, 2)), 12))

# =============================================================================
# SECTION 2: UCSV MODEL FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# svrw_gam(): Draw log-volatilities via auxiliary mixture sampler
# Translation of SVRW_gam.m (Chan 2018)
# -----------------------------------------------------------------------------
svrw_gam <- function(Ystar, htilde, h0, omegah, b0, Vh0, Vh, Vomegah) {
  
  T <- length(htilde)
  
  # Kim, Shephard & Chib (1998) 7-component normal mixture constants
  pj       <- c(0.0073, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.2575)
  mj       <- c(-10.12999, -3.97281, -8.56686, 2.77786,
                0.61942,   1.79518, -1.08819) - 1.2704
  sigj     <- c(5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261)
  sqrtsigj <- sqrt(sigj)
  
  # --- Step 1: Sample mixture component S ---
  fitted <- h0 + omegah * htilde
  q_mat  <- matrix(NA, T, 7)
  for (j in 1:7) {
    q_mat[, j] <- pj[j] * dnorm(Ystar, mean = fitted + mj[j], sd = sqrtsigj[j])
  }
  row_sums <- rowSums(q_mat)
  row_sums[row_sums <= 0] <- 1e-300
  q_mat <- q_mat / row_sums
  
  cumq     <- t(apply(q_mat, 1, cumsum))
  temprand <- runif(T)
  S        <- rowSums(cumq < temprand) + 1L
  S        <- pmin(pmax(S, 1L), 7L)
  
  # --- Step 2: Sample htilde ---
  dconst   <- mj[S]
  invOmega <- 1 / sigj[S]
  
  invSh_diag       <- c(1/Vh, rep(1, T - 1))
  HtinvShH_diag    <- invSh_diag + c(invSh_diag[-1], 0)
  HtinvShH_diag[T] <- invSh_diag[T]
  HtinvShH_off     <- -invSh_diag[-1]
  
  Kh_diag    <- HtinvShH_diag + omegah^2 * invOmega
  Kh_offdiag <- HtinvShH_off
  
  Kh_mat <- diag(Kh_diag)
  if (T > 1) {
    Kh_mat[cbind(2:T, 1:(T-1))] <- Kh_offdiag
    Kh_mat[cbind(1:(T-1), 2:T)] <- Kh_offdiag
  }
  Kh_mat <- 0.5 * (Kh_mat + t(Kh_mat)) + diag(1e-6, T)
  
  rhs <- omegah * invOmega * (Ystar - dconst - h0)
  
  htilde_draw <- tryCatch({
    Kh_chol   <- chol(Kh_mat)
    htildehat <- backsolve(Kh_chol, forwardsolve(t(Kh_chol), rhs))
    htildehat + backsolve(Kh_chol, rnorm(T))
  }, error = function(e) htilde)
  
  # --- Step 3: Sample h0 and omegah jointly ---
  Xbeta      <- cbind(1, htilde_draw)
  invVbeta   <- diag(c(1/Vh0, 1/Vomegah))
  XtinvOmega <- t(Xbeta) * invOmega
  invDbeta   <- invVbeta + XtinvOmega %*% Xbeta
  invDbeta   <- 0.5 * (invDbeta + t(invDbeta)) + diag(1e-6, 2)
  
  rhs_beta <- invVbeta %*% c(b0, 0) + XtinvOmega %*% (Ystar - dconst)
  
  betahat <- tryCatch(
    as.numeric(solve(invDbeta, rhs_beta)),
    error = function(e) c(h0, omegah)
  )
  
  Dbeta <- tryCatch({
    D <- solve(invDbeta)
    D <- 0.5 * (D + t(D)) + diag(1e-6, 2)
    D
  }, error = function(e) diag(c(Vh0, Vomegah)))
  
  beta_draw <- tryCatch(
    as.numeric(betahat + t(chol(Dbeta)) %*% rnorm(2)),
    error = function(e) c(rnorm(1, b0, sqrt(Vh0)),
                          rnorm(1, 0,  sqrt(Vomegah)))
  )
  
  h0_new     <- beta_draw[1]
  omegah_new <- beta_draw[2]
  
  # Sign normalisation
  U          <- ifelse(runif(1) > 0.5, 1, -1)
  htilde_out <- U * htilde_draw
  omegah_out <- U * omegah_new
  
  omegahhat <- betahat[2]
  Domegah   <- Dbeta[2, 2]
  
  return(list(htilde    = as.numeric(htilde_out),
              h0        = as.numeric(h0_new),
              omegah    = as.numeric(omegah_out),
              omegahhat = as.numeric(omegahhat),
              Domegah   = as.numeric(Domegah)))
}


# -----------------------------------------------------------------------------
# ucsv_gam(): Main Gibbs sampler
# Translation of UCSV_gam.m (Chan 2018)
# -----------------------------------------------------------------------------
ucsv_gam <- function(y, nloop = 11000, burnin = 1000, verbose = FALSE) {
  
  T <- length(y)
  
  # Priors (identical to Chan 2018)
  tau0    <- 0;   Vtau    <- 10
  Vh      <- 10;  Vg      <- 10
  b0      <- 0;   Vh0     <- 10
  c0      <- 0;   Vg0     <- 10
  Vomegah <- 0.2; Vomegag <- 0.2
  
  # Initialise
  h0     <- log(max(var(y, na.rm = TRUE), 1e-4)) / 2
  g0     <- log(max(var(y, na.rm = TRUE), 1e-4)) / 2
  omegah <- sqrt(0.2)
  omegag <- sqrt(0.2)
  htilde <- rep(0, T)
  gtilde <- rep(0, T)
  h      <- h0 + omegah * htilde
  g      <- g0 + omegag * gtilde
  tau    <- rep(mean(y, na.rm = TRUE), T)
  
  # Storage
  n_save      <- nloop - burnin
  store_tau   <- matrix(NA, n_save, T)
  store_h     <- matrix(NA, n_save, T)
  store_g     <- matrix(NA, n_save, T)
  store_theta <- matrix(NA, n_save, 4,
                        dimnames = list(NULL,
                                        c("omegah2", "omegag2", "h0", "g0")))
  
  for (loop in 1:nloop) {
    
    # === Block 1: Sample tau ===
    exp_neg_h <- exp(-h)
    exp_neg_g <- exp(-g)
    
    invStau_diag         <- c((1/Vtau) * exp_neg_g[1], exp_neg_g[-1])
    HtinvStauH_main      <- invStau_diag + c(invStau_diag[-1], 0)
    HtinvStauH_main[T]   <- invStau_diag[T]
    HtinvStauH_off       <- -invStau_diag[-1]
    
    invDtau_main <- HtinvStauH_main + exp_neg_h
    invDtau_off  <- HtinvStauH_off
    
    HtinvStauH_x_tau0    <- numeric(T)
    HtinvStauH_x_tau0[1] <- HtinvStauH_main[1] * tau0 +
      HtinvStauH_off[1]  * tau0
    for (t in 2:(T - 1)) {
      HtinvStauH_x_tau0[t] <- HtinvStauH_off[t-1] * tau0 +
        HtinvStauH_main[t]  * tau0 +
        HtinvStauH_off[t]   * tau0
    }
    HtinvStauH_x_tau0[T] <- HtinvStauH_off[T-1] * tau0 +
      HtinvStauH_main[T]  * tau0
    
    rhs_tau <- HtinvStauH_x_tau0 + y * exp_neg_h
    
    invDtau_mat <- diag(invDtau_main)
    invDtau_mat[cbind(2:T, 1:(T-1))] <- invDtau_off
    invDtau_mat[cbind(1:(T-1), 2:T)] <- invDtau_off
    invDtau_mat <- 0.5 * (invDtau_mat + t(invDtau_mat)) + diag(1e-6, T)
    
    tau <- tryCatch({
      chol_invDtau <- chol(invDtau_mat)
      tauhat       <- backsolve(chol_invDtau,
                                forwardsolve(t(chol_invDtau), rhs_tau))
      tauhat + backsolve(chol_invDtau, rnorm(T))
    }, error = function(e) tau)
    
    # === Block 2: Sample htilde (SV on transitory residuals) ===
    eps_pi  <- y - tau
    Ystar_h <- log(eps_pi^2 + 1e-4)
    
    res_h  <- svrw_gam(Ystar_h, htilde, h0, omegah, b0, Vh0, Vh, Vomegah)
    htilde <- res_h$htilde
    h0     <- res_h$h0
    omegah <- res_h$omegah
    h      <- h0 + omegah * htilde
    
    # === Block 3: Sample gtilde (SV on trend innovations) ===
    eps_tau <- c((tau[1] - tau0) / sqrt(Vtau), diff(tau))
    Ystar_g <- log(eps_tau^2 + 1e-4)
    
    res_g  <- svrw_gam(Ystar_g, gtilde, g0, omegag, c0, Vg0, Vg, Vomegag)
    gtilde <- res_g$htilde
    g0     <- res_g$h0
    omegag <- res_g$omegah
    g      <- g0 + omegag * gtilde
    
    # Clamp log-variances to prevent numerical explosions
    h <- pmax(pmin(h, 10), -10)
    g <- pmax(pmin(g, 10), -10)
    
    # Store after burn-in
    if (loop > burnin) {
      i                <- loop - burnin
      store_tau[i, ]   <- tau
      store_h[i, ]     <- h
      store_g[i, ]     <- g
      store_theta[i, ] <- c(omegah^2, omegag^2, h0, g0)
    }
    
    if (verbose && loop %% 2000 == 0) cat(loop, "/", nloop, "loops\n")
  }
  
  list(tau   = store_tau,
       h     = store_h,
       g     = store_g,
       theta = store_theta)
}


# -----------------------------------------------------------------------------
# ucsv_forecast(): Forward simulation to generate forecasts
# Path-average formula from Banbura & Bobeica (2023)
# Output is in same units as input (m/m annualised, x1200)
# RMSE evaluation converts to y/y using x100 series
# -----------------------------------------------------------------------------
ucsv_forecast <- function(mcmc_output, h_max = 12, n_sim = 1000) {
  
  n_draws <- nrow(mcmc_output$theta)
  T       <- ncol(mcmc_output$tau)
  
  idx <- sample(1:n_draws, min(n_sim, n_draws))
  S   <- length(idx)
  
  pi_paths <- matrix(NA, S, h_max)
  
  for (i in seq_along(idx)) {
    s <- idx[i]
    
    tau_cur <- mcmc_output$tau[s, T]
    h_cur   <- mcmc_output$h[s, T]
    g_cur   <- mcmc_output$g[s, T]
    omegah  <- sqrt(max(mcmc_output$theta[s, "omegah2"], 0))
    omegag  <- sqrt(max(mcmc_output$theta[s, "omegag2"], 0))
    
    pi_sim <- numeric(h_max)
    
    for (j in 1:h_max) {
      h_cur    <- pmax(pmin(h_cur + omegah * rnorm(1), 10), -10)
      g_cur    <- pmax(pmin(g_cur + omegag * rnorm(1), 10), -10)
      sd_trend <- exp(0.5 * g_cur)
      sd_trans <- exp(0.5 * h_cur)
      tau_cur  <- tau_cur + rnorm(1, 0, sd_trend)
      pi_sim[j] <- tau_cur + rnorm(1, 0, sd_trans)
    }
    
    pi_paths[i, ] <- pi_sim
  }
  
  # Path-average point forecasts (m/m annualised units, x1200)
  # For h=12: this approximates the year-on-year rate / 12 * 12 = yoy
  # when compared against pi_core_yoy (x100), scale by 1/12
  forecasts        <- numeric(h_max)
  names(forecasts) <- paste0("h", 1:h_max)
  
  for (h in 1:h_max) {
    path_avgs    <- rowMeans(pi_paths[, 1:h, drop = FALSE])
    forecasts[h] <- median(path_avgs)
  }
  
  list(forecasts = forecasts,
       paths     = pi_paths)
}


# =============================================================================
# SECTION 3: RECURSIVE FORECASTING LOOP
# =============================================================================

horizons   <- c(1, 3, 12)
h_max      <- max(horizons)
nloop      <- 11000
burnin     <- 1000
n_sim      <- 1000

# Change n_origins to 3 for a test run, then set back to full for production
eval_start <- which(dates == as.Date("2010-01-01"))
# n_origins  <- T_full - eval_start + 1 - h_max   # full run (~182)
 n_origins <- 3                                 # test run

cat("\nRecursive forecasting setup:\n")
cat("Estimation window start: Jan 2000\n")
cat("Evaluation period start:", format(dates[eval_start]), "\n")
cat("Number of forecast origins:", n_origins, "\n")
cat("Horizons:", horizons, "\n\n")

fc_ucsv <- matrix(NA, n_origins, length(horizons),
                  dimnames = list(
                    format(dates[eval_start:(eval_start + n_origins - 1)]),
                    paste0("h", horizons)
                  ))

set.seed(2024)
start_time <- Sys.time()

for (i in 1:n_origins) {
  
  t_now <- eval_start - 1 + i
  
  # Expanding window: all data up to forecast origin, drop leading NA
  y_est <- pi_core_mom[1:t_now]
  y_est <- y_est[!is.na(y_est)]
  
  cat("Origin", i, "/", n_origins,
      "| Date:", format(dates[t_now]),
      "| T:", length(y_est), "\n")
  
  mcmc_out <- ucsv_gam(y_est, nloop = nloop, burnin = burnin, verbose = FALSE)
  fc_out   <- ucsv_forecast(mcmc_out, h_max = h_max, n_sim = n_sim)
  
  fc_ucsv[i, ] <- fc_out$forecasts[horizons]
}

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat("\nRecursive loop completed in", round(elapsed, 1), "minutes\n")

# =============================================================================
# SECTION 4: CONVERT FORECASTS TO y/y SCALE FOR EVALUATION
# =============================================================================
# UCSV forecasts are in m/m annualised (x1200) units
# Path-average at h=12 approximates the year-on-year rate
# But x1200 vs x100 means we need to divide by 12 for comparison
# Alternatively: scale the forecast to match evaluation units

# Convert path-average m/m annualised to y/y percentage points
# path_avg_mom (x1200) / 12 * 12 months = cumulative, but path average
# already averages over h months, so:
# E[sum_{j=1}^{12} pi_{t+j}^{mom}] / 12 * 12 = sum = approx yoy in x1200 units
# To get x100 units: divide by 12
fc_ucsv_yoy        <- fc_ucsv / 12
colnames(fc_ucsv_yoy) <- paste0("h", horizons, "_yoy")

# =============================================================================
# SECTION 5: SAVE OUTPUT
# =============================================================================

dir.create(here("Output", "Forecasts"), recursive = TRUE, showWarnings = FALSE)

# Save raw forecasts (m/m annualised, x1200) -- use for model comparison
saveRDS(fc_ucsv,     here("Output", "Forecasts", "fc_ucsv_core_mom.rds"))
write.csv(fc_ucsv,   here("Output", "Forecasts", "fc_ucsv_core_mom.csv"))

# Save y/y scaled forecasts -- use for RMSE tables against pi_core_yoy
saveRDS(fc_ucsv_yoy, here("Output", "Forecasts", "fc_ucsv_core_yoy.rds"))
write.csv(fc_ucsv_yoy, here("Output", "Forecasts", "fc_ucsv_core_yoy.csv"))

cat("\nForecasts saved.\n")
cat("\nPreview (m/m annualised):\n"); print(head(fc_ucsv))
cat("\nPreview (y/y):\n");            print(head(fc_ucsv_yoy))

# =============================================================================
# SECTION 6: DIAGNOSTIC PLOT
# =============================================================================

fc_dates   <- dates[eval_start:(eval_start + n_origins - 1)]

# Recompute realised with correct scaling
pi_core_yoy <- c(rep(NA, 12), 100 * diff(log(core_level), lag = 12))

# Redo diagnostic
target_idx <- (eval_start + 12 - 1):(eval_start + n_origins + 12 - 2)
realised   <- pi_core_yoy[target_idx]
errors_h12 <- fc_ucsv[, "h12"] - realised
rmse_h12   <- sqrt(mean(errors_h12^2, na.rm = TRUE))
cat("RMSE at h=12:", round(rmse_h12, 4), "\n")

plot(fc_dates, fc_ucsv[, "h12"],
     type = "l", col = "blue", lwd = 1.5,
     ylim = range(c(fc_ucsv[, "h12"], realised), na.rm = TRUE),
     main = "UCSV: h=12 forecasts vs realised core HICP",
     ylab = "% (year-on-year)", xlab = "")
lines(fc_dates, realised, col = "black", lwd = 1.5)
legend("topleft", legend = c("UCSV forecast", "Realised"),
       col = c("blue", "black"), lty = 1, bty = "n")

