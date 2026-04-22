# =============================================================================
# 01_ucsv.R
# UCSV Model: Estimation and Recursive Forecasting
# Translation of Chan (2018), Econometric Reviews, 37(8), 807-823
# Non-centred parameterisation of Stock & Watson (2007)
# As used in Banbura & Bobeica (2023)
#
# EVALUATION NOTE:
#   - Models are estimated on pi_core_mom (m/m annualised, x1200)
#   - Forecasts are path-averages of simulated m/m annualised paths
#   - Evaluation uses the h-period average of realised m/m annualised,
#     which is the natural target consistent with the path-average formula
#   - This is consistent across all horizons h=1,2,3,4,5,6,12
# =============================================================================

Sys.setlocale("LC_TIME", "English")
library(pacman)
p_load(here,
       readxl,
       dplyr,
       lubridate)

# =============================================================================
# SECTION 1: LOAD AND PREPARE DATA
# =============================================================================

panel_a <- readRDS(here("Data", "panel_a.rds"))
raw   <- read_excel(here("Data", "EA-MD-QD", "EA-MD-QD-2026-02", "EAdata.xlsx"),
                    sheet = "data")

raw   <- raw   %>% mutate(Time = as.Date(Time))

start_date <- as.Date("2000-01-01")
raw   <- raw   %>% filter(Time >= start_date)

dates  <- panel_a$Time
T_full <- nrow(panel_a)

cat("Panel dimensions:", T_full, "x", ncol(panel_a) - 1, "\n")
cat("Date range:", format(min(dates)), "to", format(max(dates)), "\n")

saveRDS(raw, here("Data", "raw.rds"))

# --- Construct inflation series from raw HICP levels ---
core_level     <- as.numeric(raw$HICPNEF_EA)
headline_level <- as.numeric(raw$HICPOV_EA)

# Month-on-month annualised (x1200): input to all models
pi_core_mom <- c(NA, 1200 * diff(log(core_level)))
pi_hl_mom   <- c(NA, 1200 * diff(log(headline_level)))

# Year-on-year (x100): for reference and robustness plots only
pi_core_yoy <- c(rep(NA, 12), 100 * diff(log(core_level),     lag = 12))
pi_hl_yoy   <- c(rep(NA, 12), 100 * diff(log(headline_level), lag = 12))

# Sanity check
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

  # Step 1: Sample mixture component S
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

  # Step 2: Sample htilde
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

  # Step 3: Sample h0 and omegah jointly
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

  return(list(htilde    = as.numeric(htilde_out),
              h0        = as.numeric(h0_new),
              omegah    = as.numeric(omegah_out),
              omegahhat = as.numeric(betahat[2]),
              Domegah   = as.numeric(Dbeta[2, 2])))
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

    # Block 1: Sample tau
    exp_neg_h <- exp(-h)
    exp_neg_g <- exp(-g)

    invStau_diag         <- c((1/Vtau) * exp_neg_g[1], exp_neg_g[-1])
    HtinvStauH_main      <- invStau_diag + c(invStau_diag[-1], 0)
    HtinvStauH_main[T]   <- invStau_diag[T]
    HtinvStauH_off       <- -invStau_diag[-1]

    invDtau_main <- HtinvStauH_main + exp_neg_h
    invDtau_off  <- HtinvStauH_off

    HtinvStauH_x_tau0    <- numeric(T)
    HtinvStauH_x_tau0[1] <- HtinvStauH_main[1] * tau0 + HtinvStauH_off[1]  * tau0
    for (t in 2:(T - 1)) {
      HtinvStauH_x_tau0[t] <- HtinvStauH_off[t-1] * tau0 +
                                HtinvStauH_main[t]  * tau0 +
                                HtinvStauH_off[t]   * tau0
    }
    HtinvStauH_x_tau0[T] <- HtinvStauH_off[T-1] * tau0 + HtinvStauH_main[T] * tau0

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

    # Block 2: Sample htilde (SV on transitory residuals)
    eps_pi  <- y - tau
    Ystar_h <- log(eps_pi^2 + 1e-4)

    res_h  <- svrw_gam(Ystar_h, htilde, h0, omegah, b0, Vh0, Vh, Vomegah)
    htilde <- res_h$htilde
    h0     <- res_h$h0
    omegah <- res_h$omegah
    h      <- h0 + omegah * htilde

    # Block 3: Sample gtilde (SV on trend innovations)
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
# ucsv_forecast(): Forward simulation
# Path-average formula from Banbura & Bobeica (2023)
# Output units: same as input (m/m annualised, x1200)
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

  # Path-average: median over draws of (1/h * sum pi_{t+1:t+h})
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

horizons   <- c(1, 2, 3, 4, 5, 6, 12)
h_max      <- max(horizons)
nloop      <- 11000
burnin     <- 1000
n_sim      <- 1000

eval_start <- which(dates == as.Date("2010-01-01"))
n_origins  <- T_full - eval_start + 1 - h_max   # full run
# n_origins <- 3                                 # test run

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
  y_est <- pi_core_mom[1:t_now]
  y_est <- y_est[!is.na(y_est)]

  cat("Origin", i, "/", n_origins,
      "| Date:", format(dates[t_now]),
      "| T:", length(y_est), "\n")

  mcmc_out <- ucsv_gam(y_est, nloop = nloop, burnin = burnin, verbose = FALSE)
  fc_out   <- ucsv_forecast(mcmc_out, h_max = h_max, n_sim = n_sim)

  fc_ucsv[i, ] <- fc_out$forecasts[paste0("h", horizons)]
}

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat("\nRecursive loop completed in", round(elapsed, 1), "minutes\n")

# =============================================================================
# SECTION 4: SAVE OUTPUT
# =============================================================================

dir.create(here("Output", "Forecasts"), recursive = TRUE, showWarnings = FALSE)
saveRDS(fc_ucsv, here("Output", "Forecasts", "fc_ucsv_core.rds"))
write.csv(fc_ucsv, here("Output", "Forecasts", "fc_ucsv_core.csv"))
cat("Saved:", nrow(fc_ucsv), "origins x", ncol(fc_ucsv), "horizons\n")

# =============================================================================
# SECTION 5: RANDOM WALK BENCHMARKS
# =============================================================================
# RW2: forecast = 3-month trailing average of m/m annualised
# RW3: forecast = 12-month trailing average of m/m annualised (~current y/y)

make_rw_matrix <- function(n_origins, horizons, dates, eval_start) {
  matrix(NA, n_origins, length(horizons),
         dimnames = list(
           format(dates[eval_start:(eval_start + n_origins - 1)]),
           paste0("h", horizons)
         ))
}

fc_rw2 <- make_rw_matrix(n_origins, horizons, dates, eval_start)
fc_rw3 <- make_rw_matrix(n_origins, horizons, dates, eval_start)

for (i in 1:n_origins) {
  t_now <- eval_start - 1 + i

  idx3         <- max(1, t_now - 2)
  fc_rw2[i, ] <- rep(mean(pi_core_mom[idx3:t_now],  na.rm = TRUE), length(horizons))

  idx12        <- max(1, t_now - 11)
  fc_rw3[i, ] <- rep(mean(pi_core_mom[idx12:t_now], na.rm = TRUE), length(horizons))
}

saveRDS(fc_rw2, here("Output", "Forecasts", "fc_rw2_core.rds"))
write.csv(fc_rw2, here("Output", "Forecasts", "fc_rw2_core.csv"))
saveRDS(fc_rw3, here("Output", "Forecasts", "fc_rw3_core.rds"))
write.csv(fc_rw3, here("Output", "Forecasts", "fc_rw3_core.csv"))

# =============================================================================
# SECTION 6: RMSE BY HORIZON
# =============================================================================
# Realised at horizon h = average of pi_core_mom[t+1 : t+h]
# This is the natural evaluation target consistent with path-average forecasts

compute_rmse_by_horizon <- function(fc_matrix, horizons, pi_mom, eval_start) {
  rmse_vec <- numeric(length(horizons))
  names(rmse_vec) <- paste0("h", horizons)

  for (k in seq_along(horizons)) {
    h         <- horizons[k]
    n_fc      <- nrow(fc_matrix)
    realised  <- numeric(n_fc)

    for (i in 1:n_fc) {
      t_now       <- eval_start - 1 + i
      future_vals <- pi_mom[(t_now + 1):(t_now + h)]
      realised[i] <- mean(future_vals, na.rm = TRUE)
    }

    errors       <- fc_matrix[, paste0("h", h)] - realised
    rmse_vec[k]  <- sqrt(mean(errors^2, na.rm = TRUE))
  }

  rmse_vec
}

rmse_ucsv <- compute_rmse_by_horizon(fc_ucsv, horizons, pi_core_mom, eval_start)
rmse_rw2  <- compute_rmse_by_horizon(fc_rw2,  horizons, pi_core_mom, eval_start)
rmse_rw3  <- compute_rmse_by_horizon(fc_rw3,  horizons, pi_core_mom, eval_start)

# Print results table
cat("\nRMSE by horizon:\n")
cat(sprintf("%-6s  %8s  %8s  %8s  %8s  %8s\n",
            "h", "RW(3m)", "RW(12m)", "UCSV", "UCSV/RW2", "UCSV/RW3"))
cat(strrep("-", 58), "\n")
for (k in seq_along(horizons)) {
  cat(sprintf("h=%-4s  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n",
              horizons[k],
              rmse_rw2[k], rmse_rw3[k], rmse_ucsv[k],
              rmse_ucsv[k] / rmse_rw2[k],
              rmse_ucsv[k] / rmse_rw3[k]))
}

# Save RMSE table
dir.create(here("Output", "Tables"), recursive = TRUE, showWarnings = FALSE)
rmse_table <- data.frame(
  horizon        = paste0("h", horizons),
  rmse_rw2       = rmse_rw2,
  rmse_rw3       = rmse_rw3,
  rmse_ucsv      = rmse_ucsv,
  ratio_ucsv_rw2 = rmse_ucsv / rmse_rw2,
  ratio_ucsv_rw3 = rmse_ucsv / rmse_rw3
)
write.csv(rmse_table, here("Output", "Tables", "rmse_ucsv_rw.csv"), row.names = FALSE)

# =============================================================================
# SECTION 7: PLOTS
# =============================================================================

dir.create(here("Output", "Figures"), recursive = TRUE, showWarnings = FALSE)
fc_ucsv  <- readRDS(here("Output", "Forecasts", "fc_ucsv_core.rds"))
fc_rw2   <- readRDS(here("Output", "Forecasts", "fc_rw2_core.rds"))
fc_rw3   <- readRDS(here("Output", "Forecasts", "fc_rw3_core.rds"))
horizons <- c(1, 2, 3, 4, 5, 6, 12)

ratio_rw2 <- rmse_ucsv / rmse_rw2
ratio_rw3 <- rmse_ucsv / rmse_rw3

# --- Plot 1: Absolute RMSE ---
par(mfrow = c(1, 2))

plot(horizons, rmse_ucsv,
     type = "b", col = "blue", lwd = 2, pch = 16,
     ylim = range(c(rmse_rw2, rmse_rw3, rmse_ucsv)) * c(0.9, 1.1),
     xlab = "Forecast horizon (months)",
     ylab = "RMSE (m/m annualised, %)",
     main = "RMSE by horizon",
     xaxt = "n")
axis(1, at = horizons)
lines(horizons, rmse_rw2, type = "b", col = "red",        lwd = 2, pch = 16, lty = 2)
lines(horizons, rmse_rw3, type = "b", col = "darkorange", lwd = 2, pch = 16, lty = 3)
legend("topright",
       legend = c("UCSV", "RW (3m avg)", "RW (12m avg)"),
       col    = c("blue", "red", "darkorange"),
       lty    = c(1, 2, 3), lwd = 2, pch = 16, bty = "n")

# --- Plot 2: Relative RMSE ---
plot(horizons, ratio_rw2,
     type = "b", col = "red", lwd = 2, pch = 16, lty = 2,
     ylim = range(c(ratio_rw2, ratio_rw3)) * c(0.93, 1.05),
     xlab = "Forecast horizon (months)",
     ylab = "RMSE ratio (UCSV / benchmark)",
     main = "Relative RMSE: UCSV vs benchmarks",
     xaxt = "n")
axis(1, at = horizons)
lines(horizons, ratio_rw3, type = "b", col = "darkorange", lwd = 2, pch = 16, lty = 3)
abline(h = 1, col = "grey40", lty = 2, lwd = 1.5)
text(x = horizons, y = ratio_rw2,
     labels = round(ratio_rw2, 2), pos = 3, cex = 0.7, col = "red")
text(x = horizons, y = ratio_rw3,
     labels = round(ratio_rw3, 2), pos = 1, cex = 0.7, col = "darkorange")
legend("bottomleft",
       legend = c("UCSV / RW (3m avg)", "UCSV / RW (12m avg)"),
       col    = c("red", "darkorange"),
       lty    = c(2, 3), lwd = 2, pch = 16, bty = "n")

par(mfrow = c(1, 1))

# --- Plot 3: h=1 forecasts vs realised ---
fc_dates    <- dates[eval_start:(eval_start + n_origins - 1)]
realised_h1 <- pi_core_mom[(eval_start + 1):(eval_start + n_origins)]

plot(fc_dates, fc_ucsv[, "h1"],
     type = "l", col = "blue", lwd = 1.5,
     ylim = range(c(fc_ucsv[, "h1"], fc_rw2[, "h1"],
                    fc_rw3[, "h1"], realised_h1), na.rm = TRUE),
     main = "h=1 forecasts vs realised (m/m annualised)",
     ylab = "% (annualised)", xlab = "")
lines(fc_dates, fc_rw2[, "h1"], col = "red",        lwd = 1.2, lty = 2)
lines(fc_dates, fc_rw3[, "h1"], col = "darkorange", lwd = 1.2, lty = 3)
lines(fc_dates, realised_h1,    col = "grey40",     lwd = 1.5)
legend("topleft",
       legend = c("UCSV", "RW (3m avg)", "RW (12m avg)", "Realised"),
       col    = c("blue", "red", "darkorange", "grey40"),
       lty    = c(1, 2, 3, 1), lwd = 2, bty = "n")

# --- Plot 4: h=12 forecasts vs realised ---
realised_h12 <- numeric(n_origins)
for (i in 1:n_origins) {
  t_now           <- eval_start - 1 + i
  realised_h12[i] <- mean(pi_core_mom[(t_now + 1):(t_now + 12)], na.rm = TRUE)
}

plot(fc_dates, fc_ucsv[, "h12"],
     type = "l", col = "blue", lwd = 1.5,
     ylim = range(c(fc_ucsv[, "h12"], fc_rw2[, "h12"],
                    fc_rw3[, "h12"], realised_h12), na.rm = TRUE),
     main = "h=12 forecasts vs realised (m/m annualised avg)",
     ylab = "% (annualised)", xlab = "")
lines(fc_dates, fc_rw2[, "h12"], col = "red",        lwd = 1.2, lty = 2)
lines(fc_dates, fc_rw3[, "h12"], col = "darkorange", lwd = 1.2, lty = 3)
lines(fc_dates, realised_h12,    col = "grey40",     lwd = 1.5)
legend("topleft",
       legend = c("UCSV", "RW (3m avg)", "RW (12m avg)", "Realised"),
       col    = c("blue", "red", "darkorange", "grey40"),
       lty    = c(1, 2, 3, 1), lwd = 2, bty = "n")
 


