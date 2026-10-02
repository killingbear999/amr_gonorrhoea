# One-dimensional doxy-PEP efficacy sensitivity analysis.
# Based on the user's current working run_sensitivity_efficacy.R.
# Run from the folder containing the calibrated fit, or edit fit_file below.
# Dependencies: install.packages(c("deSolve", "ggplot2", "rstan", "patchwork"))
library(deSolve)
library(ggplot2)
library(patchwork)

# USER SETTINGS (efficacy and projection years confirmed by the user)
fit_file <- "fit_results_fixedinitialstate_UK_6years_allcases_covid_123.rds"
output_dir <- "sensitivity_doxypep_four_strategies"
# Vary doxy-PEP efficacy over the same 0--100% range and 10-point spacing.
efficacy_doxypep <- seq(0, 1, by = 0.1)
# Confirmed original 4CMenB vaccine efficacy: 40%.
efficacy_vaccine_fixed <- 0.40
p_d_combined <- 0.66
p_v_combined <- 0.66
n_iter <- 1000L
seed <- 42L
# Confirmed projection period: 2027--2041 (15 annual intervals).
# Intervention begins at the start of the projection, as in the original run_amr.
# Model times 9:24 and initial state from y[,8,] are unchanged.
first_intervention_year <- 2027L
target_year <- 2041L
n_years <- 15L
years <- first_intervention_year + seq_len(n_years) - 1L
idx_2041 <- match(target_year, years)
stopifnot(!is.na(idx_2041), n_iter > 0, n_iter == as.integer(n_iter),
          length(efficacy_vaccine_fixed) == 1L, is.finite(efficacy_vaccine_fixed),
          efficacy_vaccine_fixed >= 0, efficacy_vaccine_fixed <= 1,
          length(efficacy_doxypep) >= 2, all(is.finite(efficacy_doxypep)),
          all(efficacy_doxypep >= 0 & efficacy_doxypep <= 1))
efficacy_doxypep <- sort(unique(efficacy_doxypep))
ensure_output_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Cannot create output directory: ", path,
                              ". Choose a writable output_dir.")
  if (file.access(path, 2) != 0) stop("Output directory is not writable: ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}
# Resolve once so a later change in working directory cannot redirect output.
output_dir <- ensure_output_dir(output_dir)

# IMPORTANT: equations, isFixed = TRUE, initial-state seeding, and the original
# annual-incidence approximation are retained. Despite its trapezoidal comment,
# that approximation holds lambda at t+8 and averages endpoint U values;
# it is not exact integration of the infection flow over each year.
# These figures use that same cases_all definition for comparability.

# Helper functions
get_C <- function(E_N, A_N, S_N, E_D, A_D, S_D, E_V, A_V, S_V, E_M, A_M, S_M) {
  return(E_N + A_N + S_N + E_D + A_D + S_D + E_V + A_V + S_V + E_M + A_M + S_M)
}

get_N <- function(U_N, E_N_0, A_N_0, S_N_0, T_N_0, E_N_c, A_N_c, S_N_c, T_N_c, E_N_t, A_N_t, S_N_t, T_N_t, E_N_d2, A_N_d2, S_N_d2, T_N_d2, U_D, E_D_0, A_D_0, S_D_0, T_D_0, E_D_c, A_D_c, S_D_c, T_D_c, E_D_t, A_D_t, S_D_t, T_D_t, E_D_d2, A_D_d2, S_D_d2, T_D_d2, U_V, E_V_0, A_V_0, S_V_0, T_V_0, E_V_c, A_V_c, S_V_c, T_V_c, E_V_t, A_V_t, S_V_t, T_V_t, E_V_d2, A_V_d2, S_V_d2, T_V_d2, U_M, E_M_0, A_M_0, S_M_0, T_M_0, E_M_c, A_M_c, S_M_c, T_M_c, E_M_t, A_M_t, S_M_t, T_M_t, E_M_d2, A_M_d2, S_M_d2, T_M_d2) {
  return(U_N + E_N_0 + A_N_0 + S_N_0 + T_N_0 + E_N_c + A_N_c + S_N_c + T_N_c + E_N_t + A_N_t + S_N_t + T_N_t + E_N_d2 + A_N_d2 + S_N_d2 + T_N_d2 + U_D + E_D_0 + A_D_0 + S_D_0 + T_D_0 + E_D_c + A_D_c + S_D_c + T_D_c + E_D_t + A_D_t + S_D_t + T_D_t + E_D_d2 + A_D_d2 + S_D_d2 + T_D_d2 + U_V + E_V_0 + A_V_0 + S_V_0 + T_V_0 + E_V_c + A_V_c + S_V_c + T_V_c + E_V_t + A_V_t + S_V_t + T_V_t + E_V_d2 + A_V_d2 + S_V_d2 + T_V_d2 + U_M + E_M_0 + A_M_0 + S_M_0 + T_M_0 + E_M_c + A_M_c + S_M_c + T_M_c + E_M_t + A_M_t + S_M_t + T_M_t + E_M_d2 + A_M_d2 + S_M_d2 + T_M_d2)
}

get_pi <- function(c_target, N_target, c_remain, N_remain) {
  return((c_target * N_target) / (c_target * N_target + c_remain * N_remain))
}

get_lambda <- function(t, t_0, c, beta, phi_beta, epsilon, C_target, N_target, pi_target, C_remain, N_remain, pi_remain, isFixed) {
  if (isFixed && t > 8) {
    t <- 8
  }
  return(c * beta * (1 + phi_beta * (t - t_0)) * (epsilon * C_target / N_target + (1 - epsilon) * 
                                                    (pi_target * C_target / N_target + pi_remain * C_remain / N_remain)))
}

get_eta <- function(t, t_0, eta_H_init, phi_eta, isFixed) {
  if (isFixed && t > 8) {
    t <- 8
  }
  return(eta_H_init * (1 + phi_eta * (t - t_0)))
}

# ODE system
amr_model <- function(t, y, parameters) {
  with(as.list(c(y, parameters)), {
    # two cases: 1. the inferred trends in the time-varying behavioural parameters stabilise
    #            2. the trends continue until the end of the modelled period
    isFixed <- TRUE
    
    C_H_0 = get_C(E_N_H_0, A_N_H_0, S_N_H_0, E_D_H_0, A_D_H_0, S_D_H_0, E_V_H_0, A_V_H_0, S_V_H_0, E_M_H_0, A_M_H_0, S_M_H_0)
    C_L_0 = get_C(E_N_L_0, A_N_L_0, S_N_L_0, E_D_L_0, A_D_L_0, S_D_L_0, E_V_L_0, A_V_L_0, S_V_L_0, E_M_L_0, A_M_L_0, S_M_L_0)
    C_H_c = get_C(E_N_H_c, A_N_H_c, S_N_H_c, E_D_H_c, A_D_H_c, S_D_H_c, E_V_H_c, A_V_H_c, S_V_H_c, E_M_H_c, A_M_H_c, S_M_H_c)
    C_L_c = get_C(E_N_L_c, A_N_L_c, S_N_L_c, E_D_L_c, A_D_L_c, S_D_L_c, E_V_L_c, A_V_L_c, S_V_L_c, E_M_L_c, A_M_L_c, S_M_L_c)
    C_H_t = get_C(E_N_H_t, A_N_H_t, S_N_H_t, E_D_H_t, A_D_H_t, S_D_H_t, E_V_H_t, A_V_H_t, S_V_H_t, E_M_H_t, A_M_H_t, S_M_H_t)
    C_L_t = get_C(E_N_L_t, A_N_L_t, S_N_L_t, E_D_L_t, A_D_L_t, S_D_L_t, E_V_L_t, A_V_L_t, S_V_L_t, E_M_L_t, A_M_L_t, S_M_L_t)
    C_H_d2 = get_C(E_N_H_d2, A_N_H_d2, S_N_H_d2, E_D_H_d2, A_D_H_d2, S_D_H_d2, E_V_H_d2, A_V_H_d2, S_V_H_d2, E_M_H_d2, A_M_H_d2, S_M_H_d2)
    C_L_d2 = get_C(E_N_L_d2, A_N_L_d2, S_N_L_d2, E_D_L_d2, A_D_L_d2, S_D_L_d2, E_V_L_d2, A_V_L_d2, S_V_L_d2, E_M_L_d2, A_M_L_d2, S_M_L_d2)
    N_H <- get_N(U_N_H, E_N_H_0, A_N_H_0, S_N_H_0, T_N_H_0, E_N_H_c, A_N_H_c, S_N_H_c, T_N_H_c, E_N_H_t, A_N_H_t, S_N_H_t, T_N_H_t, E_N_H_d2, A_N_H_d2, S_N_H_d2, T_N_H_d2, 
                 U_D_H, E_D_H_0, A_D_H_0, S_D_H_0, T_D_H_0, E_D_H_c, A_D_H_c, S_D_H_c, T_D_H_c, E_D_H_t, A_D_H_t, S_D_H_t, T_D_H_t, E_D_H_d2, A_D_H_d2, S_D_H_d2, T_D_H_d2, 
                 U_V_H, E_V_H_0, A_V_H_0, S_V_H_0, T_V_H_0, E_V_H_c, A_V_H_c, S_V_H_c, T_V_H_c, E_V_H_t, A_V_H_t, S_V_H_t, T_V_H_t, E_V_H_d2, A_V_H_d2, S_V_H_d2, T_V_H_d2, 
                 U_M_H, E_M_H_0, A_M_H_0, S_M_H_0, T_M_H_0, E_M_H_c, A_M_H_c, S_M_H_c, T_M_H_c, E_M_H_t, A_M_H_t, S_M_H_t, T_M_H_t, E_M_H_d2, A_M_H_d2, S_M_H_d2, T_M_H_d2)
    N_L <- get_N(U_N_L, E_N_L_0, A_N_L_0, S_N_L_0, T_N_L_0, E_N_L_c, A_N_L_c, S_N_L_c, T_N_L_c, E_N_L_t, A_N_L_t, S_N_L_t, T_N_L_t, E_N_L_d2, A_N_L_d2, S_N_L_d2, T_N_L_d2, 
                 U_D_L, E_D_L_0, A_D_L_0, S_D_L_0, T_D_L_0, E_D_L_c, A_D_L_c, S_D_L_c, T_D_L_c, E_D_L_t, A_D_L_t, S_D_L_t, T_D_L_t, E_D_L_d2, A_D_L_d2, S_D_L_d2, T_D_L_d2, 
                 U_V_L, E_V_L_0, A_V_L_0, S_V_L_0, T_V_L_0, E_V_L_c, A_V_L_c, S_V_L_c, T_V_L_c, E_V_L_t, A_V_L_t, S_V_L_t, T_V_L_t, E_V_L_d2, A_V_L_d2, S_V_L_d2, T_V_L_d2, 
                 U_M_L, E_M_L_0, A_M_L_0, S_M_L_0, T_M_L_0, E_M_L_c, A_M_L_c, S_M_L_c, T_M_L_c, E_M_L_t, A_M_L_t, S_M_L_t, T_M_L_t, E_M_L_d2, A_M_L_d2, S_M_L_d2, T_M_L_d2)
    
    pi_H <- get_pi(c_H, N_H, c_L, N_L)
    pi_L <- get_pi(c_L, N_L, c_H, N_H)
    
    eta_H <- get_eta(t, t_0, eta_H_init, phi_eta, isFixed)
    eta_L <- omega * eta_H
    
    lambda_H_0 = get_lambda(t, t_0, c_H, beta, phi_beta, epsilon, C_H_0, N_H, pi_H, C_L_0, N_L, pi_L, isFixed)
    lambda_L_0 = get_lambda(t, t_0, c_L, beta, phi_beta, epsilon, C_L_0, N_L, pi_L, C_H_0, N_H, pi_H, isFixed)
    lambda_H_c = get_lambda(t, t_0, c_H, beta, phi_beta, epsilon, C_H_c, N_H, pi_H, C_L_c, N_L, pi_L, isFixed)
    lambda_L_c = get_lambda(t, t_0, c_L, beta, phi_beta, epsilon, C_L_c, N_L, pi_L, C_H_c, N_H, pi_H, isFixed)
    lambda_H_t = get_lambda(t, t_0, c_H, beta, phi_beta, epsilon, C_H_t, N_H, pi_H, C_L_t, N_L, pi_L, isFixed)
    lambda_L_t = get_lambda(t, t_0, c_L, beta, phi_beta, epsilon, C_L_t, N_L, pi_L, C_H_t, N_H, pi_H, isFixed)
    lambda_H_d2 = get_lambda(t, t_0, c_H, beta, phi_beta, epsilon, C_H_d2, N_H, pi_H, C_L_d2, N_L, pi_L, isFixed)
    lambda_L_d2 = get_lambda(t, t_0, c_L, beta, phi_beta, epsilon, C_L_d2, N_L, pi_L, C_H_d2, N_H, pi_H, isFixed)
    
    # ODEs
    # high-risk group
    # no intervention (N)
    dU_N_H = q_H * alpha + rho * (1 - w_c) * (T_N_H_0 + T_N_H_t) + rho * (1 - phi) * (T_N_H_c + T_N_H_d2) + nu * (A_N_H_0 + A_N_H_c + A_N_H_t + A_N_H_d2) - (lambda_H_0 + f_c * lambda_H_c + f_t * lambda_H_t + f_d2 * lambda_H_d2 + (p_d + p_v) * eta_H + 1/gamma) * U_N_H + xi_n * U_D_H + xi_d * U_M_H
    dE_N_H_0 = lambda_H_0 * U_N_H - (sigma + 1/gamma) * E_N_H_0 + xi_n * E_D_H_0 + xi_d * E_M_H_0
    dA_N_H_0 = sigma * (1 - psi) * E_N_H_0 - (nu + eta_H + 1/gamma) * A_N_H_0 + xi_n * A_D_H_0 + xi_d * A_M_H_0
    dS_N_H_0 = sigma * psi * E_N_H_0 - (mu + 1/gamma) * S_N_H_0 + xi_n * S_D_H_0 + xi_d * S_M_H_0
    dT_N_H_0 = eta_H * A_N_H_0 + mu * S_N_H_0 - (rho + 1/gamma) * T_N_H_0 + xi_n * T_D_H_0 + xi_d * T_M_H_0
    dE_N_H_c = f_c * lambda_H_c * U_N_H - (sigma + 1/gamma) * E_N_H_c + xi_n * E_D_H_c + xi_d * E_M_H_c
    dA_N_H_c = sigma * (1 - psi) * E_N_H_c - (nu + eta_H + 1/gamma) * A_N_H_c + phi * rho * T_N_H_c + xi_n * A_D_H_c + xi_d * A_M_H_c
    dS_N_H_c = sigma * psi * E_N_H_c - (mu + 1/gamma) * S_N_H_c + xi_n * S_D_H_c + xi_d * S_M_H_c
    dT_N_H_c = eta_H * A_N_H_c + mu * S_N_H_c - (rho + 1/gamma) * T_N_H_c + w_c * rho * T_N_H_0 + xi_n * T_D_H_c + xi_d * T_M_H_c
    dE_N_H_t = f_t * lambda_H_t * U_N_H - (sigma + 1/gamma) * E_N_H_t + xi_n * E_D_H_t + xi_d * E_M_H_t
    dA_N_H_t = sigma * (1 - psi) * E_N_H_t - (nu + eta_H + 1/gamma) * A_N_H_t + xi_n * A_D_H_t + xi_d * A_M_H_t
    dS_N_H_t = sigma * psi * E_N_H_t - (mu + 1/gamma) * S_N_H_t + xi_n * S_D_H_t + xi_d * S_M_H_t
    dT_N_H_t = eta_H * A_N_H_t + mu * S_N_H_t - (rho + 1/gamma) * T_N_H_t + xi_n * T_D_H_t + xi_d * T_M_H_t
    dE_N_H_d2 = f_d2 * lambda_H_d2 * U_N_H - (sigma + 1/gamma) * E_N_H_d2 + xi_n * E_D_H_d2 + xi_d * E_M_H_d2
    dA_N_H_d2 = sigma * (1 - psi) * E_N_H_d2 - (nu + eta_H + 1/gamma) * A_N_H_d2 + phi * rho * T_N_H_d2 + xi_n * A_D_H_d2 + xi_d * A_M_H_d2
    dS_N_H_d2 = sigma * psi * E_N_H_d2 - (mu + 1/gamma) * S_N_H_d2 + xi_n * S_D_H_d2 + xi_d * S_M_H_d2
    dT_N_H_d2 = eta_H * A_N_H_d2 + mu * S_N_H_d2 - (rho + 1/gamma) * T_N_H_d2 + w_c * rho * T_N_H_t + xi_n * T_D_H_d2 + xi_d * T_M_H_d2
    
    # doxy-PEP (D)
    dU_D_H = eta_H * p_d * U_N_H + rho * (1 - w_c) * (T_D_H_0 + T_D_H_t) + rho * (1 - phi) * (T_D_H_c + T_D_H_d2) + nu * (A_D_H_0 + A_D_H_c + A_D_H_t + A_D_H_d2) - (e_d * lambda_H_0 + e_d * f_c * lambda_H_c + f_t * lambda_H_t + f_d2 * lambda_H_d2 + xi_n + p_v * eta_H + 1/gamma) * U_D_H + xi_d * U_V_H
    dE_D_H_0 = e_d * lambda_H_0 * U_D_H - (sigma + xi_n + 1/gamma) * E_D_H_0 + xi_d * E_V_H_0
    dA_D_H_0 = (1 - w_t) * sigma * (1 - psi) * E_D_H_0 - (nu + eta_H + xi_n + 1/gamma) * A_D_H_0 + xi_d * A_V_H_0
    dS_D_H_0 = (1 - w_t) * sigma * psi * E_D_H_0 - (mu + xi_n + 1/gamma) * S_D_H_0 + xi_d * S_V_H_0
    dT_D_H_0 = eta_H * A_D_H_0 + mu * S_D_H_0 - (rho + xi_n + 1/gamma) * T_D_H_0 + xi_d * T_V_H_0
    dE_D_H_c = e_d * f_c * lambda_H_c * U_D_H - (sigma + xi_n + 1/gamma) * E_D_H_c + xi_d * E_V_H_c
    dA_D_H_c = (1 - w_t) * sigma * (1 - psi) * E_D_H_c - (nu + eta_H + xi_n + 1/gamma) * A_D_H_c + phi * rho * T_D_H_c + xi_d * A_V_H_c
    dS_D_H_c = (1 - w_t) * sigma * psi * E_D_H_c - (mu + xi_n + 1/gamma) * S_D_H_c + xi_d * S_V_H_c
    dT_D_H_c = eta_H * A_D_H_c + mu * S_D_H_c - (rho + xi_n + 1/gamma) * T_D_H_c + w_c * rho * T_D_H_0 + xi_d * T_V_H_c
    dE_D_H_t = f_t * lambda_H_t * U_D_H - (sigma + xi_n + 1/gamma) * E_D_H_t + xi_d * E_V_H_t
    dA_D_H_t = w_t * sigma * (1 - psi) * E_D_H_0 + sigma * (1 - psi) * E_D_H_t - (nu + eta_H + xi_n + 1/gamma) * A_D_H_t + xi_d * A_V_H_t
    dS_D_H_t = w_t * sigma * psi * E_D_H_0 + sigma * psi * E_D_H_t - (mu + xi_n + 1/gamma) * S_D_H_t + xi_d * S_V_H_t
    dT_D_H_t = eta_H * A_D_H_t + mu * S_D_H_t - (rho + xi_n + 1/gamma) * T_D_H_t + xi_d * T_V_H_t
    dE_D_H_d2 = f_d2 * lambda_H_d2 * U_D_H - (sigma + xi_n + 1/gamma) * E_D_H_d2 + xi_d * E_V_H_d2
    dA_D_H_d2 = w_t * sigma * (1 - psi) * E_D_H_c  + sigma * (1 - psi) * E_D_H_d2 - (nu + eta_H + xi_n + 1/gamma) * A_D_H_d2 + phi * rho * T_D_H_d2 + xi_d * A_V_H_d2
    dS_D_H_d2 = w_t * sigma * psi * E_D_H_c + sigma * psi * E_D_H_d2 - (mu + xi_n + 1/gamma) * S_D_H_d2 + xi_d * S_V_H_d2
    dT_D_H_d2 = eta_H * A_D_H_d2 + mu * S_D_H_d2 - (rho + xi_n + 1/gamma) * T_D_H_d2 + w_c * rho * T_D_H_t + xi_d * T_V_H_d2
    
    # doxy-PEP + 4CMenB (V)
    dU_V_H = eta_H * p_v * U_D_H + eta_H * p_d * U_M_H + rho * (1 - w_c) * (T_V_H_0 + T_V_H_t) + rho * (1 - phi) * (T_V_H_c + T_V_H_d2) + nu * (A_V_H_0 + A_V_H_c + A_V_H_t + A_V_H_d2) - (e_vd * lambda_H_0 + e_vd * f_c * lambda_H_c + e_v * f_t * lambda_H_t + e_v * f_d2 * lambda_H_d2 + xi_n + xi_d + 1/gamma) * U_V_H
    dE_V_H_0 = e_vd * lambda_H_0 * U_V_H - (sigma + xi_n + xi_d + 1/gamma) * E_V_H_0
    dA_V_H_0 = (1 - w_t) * sigma * (1 - psi) * E_V_H_0 - (nu + eta_H + xi_n + xi_d + 1/gamma) * A_V_H_0
    dS_V_H_0 = (1 - w_t) * sigma * psi * E_V_H_0 - (mu + xi_n + xi_d + 1/gamma) * S_V_H_0
    dT_V_H_0 = eta_H * A_V_H_0 + mu * S_V_H_0 - (rho + xi_n + xi_d + 1/gamma) * T_V_H_0
    dE_V_H_c = e_vd * f_c * lambda_H_c * U_V_H - (sigma + xi_n + xi_d + 1/gamma) * E_V_H_c
    dA_V_H_c = (1 - w_t) * sigma * (1 - psi) * E_V_H_c - (nu + eta_H + xi_n + xi_d + 1/gamma) * A_V_H_c + phi * rho * T_V_H_c
    dS_V_H_c = (1 - w_t) * sigma * psi * E_V_H_c - (mu + xi_n + xi_d + 1/gamma) * S_V_H_c
    dT_V_H_c = eta_H * A_V_H_c + mu * S_V_H_c - (rho + xi_n + xi_d + 1/gamma) * T_V_H_c + w_c * rho * T_V_H_0
    dE_V_H_t = e_v * f_t * lambda_H_t * U_V_H - (sigma + xi_n + xi_d + 1/gamma) * E_V_H_t
    dA_V_H_t = w_t * sigma * (1 - psi) * E_V_H_0 + sigma * (1 - psi) * E_V_H_t - (nu + eta_H + xi_n + xi_d + 1/gamma) * A_V_H_t
    dS_V_H_t = w_t * sigma * psi * E_V_H_0 + sigma * psi * E_V_H_t - (mu + xi_n + xi_d + 1/gamma) * S_V_H_t
    dT_V_H_t = eta_H * A_V_H_t + mu * S_V_H_t - (rho + xi_n + xi_d + 1/gamma) * T_V_H_t
    dE_V_H_d2 = e_v * f_d2 * lambda_H_d2 * U_V_H - (sigma + xi_n + xi_d + 1/gamma) * E_V_H_d2
    dA_V_H_d2 = w_t * sigma * (1 - psi) * E_V_H_c  + sigma * (1 - psi) * E_V_H_d2 - (nu + eta_H + xi_n + xi_d + 1/gamma) * A_V_H_d2 + phi * rho * T_V_H_d2
    dS_V_H_d2 = w_t * sigma * psi * E_V_H_c + sigma * psi * E_V_H_d2 - (mu + xi_n + xi_d + 1/gamma) * S_V_H_d2
    dT_V_H_d2 = eta_H * A_V_H_d2 + mu * S_V_H_d2 - (rho + xi_n + xi_d + 1/gamma) * T_V_H_d2 + w_c * rho * T_V_H_t
    
    # 4CMenB (M)
    dU_M_H = eta_H * p_v * U_N_H + xi_n * U_V_H + rho * (1 - w_c) * (T_M_H_0 + T_M_H_t) + rho * (1 - phi) * (T_M_H_c + T_M_H_d2) + nu * (A_M_H_0 + A_M_H_c + A_M_H_t + A_M_H_d2) - (e_v * lambda_H_0 + e_v * f_c * lambda_H_c + e_v * f_t * lambda_H_t + e_v * f_d2 * lambda_H_d2 + p_d * eta_H + xi_d + 1/gamma) * U_M_H
    dE_M_H_0 = e_v * lambda_H_0 * U_M_H - (sigma + 1/gamma) * E_M_H_0 + xi_n * E_V_H_0 - xi_d * E_M_H_0
    dA_M_H_0 = sigma * (1 - psi) * E_M_H_0 - (nu + eta_H + 1/gamma) * A_M_H_0 + xi_n * A_V_H_0 - xi_d * A_M_H_0
    dS_M_H_0 = sigma * psi * E_M_H_0 - (mu + 1/gamma) * S_M_H_0 + xi_n * S_V_H_0 - xi_d * S_M_H_0
    dT_M_H_0 = eta_H * A_M_H_0 + mu * S_M_H_0 - (rho + 1/gamma) * T_M_H_0 + xi_n * T_V_H_0 - xi_d * T_M_H_0
    dE_M_H_c = e_v * f_c * lambda_H_c * U_M_H - (sigma + 1/gamma) * E_M_H_c + xi_n * E_V_H_c - xi_d * E_M_H_c
    dA_M_H_c = sigma * (1 - psi) * E_M_H_c - (nu + eta_H + 1/gamma) * A_M_H_c + phi * rho * T_M_H_c + xi_n * A_V_H_c - xi_d * A_M_H_c
    dS_M_H_c = sigma * psi * E_M_H_c - (mu + 1/gamma) * S_M_H_c + xi_n * S_V_H_c - xi_d * S_M_H_c
    dT_M_H_c = eta_H * A_M_H_c + mu * S_M_H_c - (rho + 1/gamma) * T_M_H_c + w_c * rho * T_M_H_0 + xi_n * T_V_H_c - xi_d * T_M_H_c
    dE_M_H_t = e_v * f_t * lambda_H_t * U_M_H - (sigma + 1/gamma) * E_M_H_t + xi_n * E_V_H_t - xi_d * E_M_H_t
    dA_M_H_t = sigma * (1 - psi) * E_M_H_t - (nu + eta_H + 1/gamma) * A_M_H_t + xi_n * A_V_H_t - xi_d * A_M_H_t
    dS_M_H_t = sigma * psi * E_M_H_t - (mu + 1/gamma) * S_M_H_t + xi_n * S_V_H_t - xi_d * S_M_H_t
    dT_M_H_t = eta_H * A_M_H_t + mu * S_M_H_t - (rho + 1/gamma) * T_M_H_t + xi_n * T_V_H_t - xi_d * T_M_H_t
    dE_M_H_d2 = e_v * f_d2 * lambda_H_d2 * U_M_H - (sigma + 1/gamma) * E_M_H_d2 + xi_n * E_V_H_d2 - xi_d * E_M_H_d2
    dA_M_H_d2 = sigma * (1 - psi) * E_M_H_d2 - (nu + eta_H + 1/gamma) * A_M_H_d2 + phi * rho * T_M_H_d2 + xi_n * A_V_H_d2 - xi_d * A_M_H_d2
    dS_M_H_d2 = sigma * psi * E_M_H_d2 - (mu + 1/gamma) * S_M_H_d2 + xi_n * S_V_H_d2 - xi_d * S_M_H_d2
    dT_M_H_d2 = eta_H * A_M_H_d2 + mu * S_M_H_d2 - (rho + 1/gamma) * T_M_H_d2 + w_c * rho * T_M_H_t + xi_n * T_V_H_d2 - xi_d * T_M_H_d2
    
    # low-risk group
    # no intervention (N)
    dU_N_L = q_L * alpha + rho * (1 - w_c) * (T_N_L_0 + T_N_L_t) + rho * (1 - phi) * (T_N_L_c + T_N_L_d2) + nu * (A_N_L_0 + A_N_L_c + A_N_L_t + A_N_L_d2) - (lambda_L_0 + f_c * lambda_L_c + f_t * lambda_L_t + f_d2 * lambda_L_d2 + (p_d + p_v) * eta_L + 1/gamma) * U_N_L + xi_n * U_D_L + xi_d * U_M_L
    dE_N_L_0 = lambda_L_0 * U_N_L - (sigma + 1/gamma) * E_N_L_0 + xi_n * E_D_L_0 + xi_d * E_M_L_0
    dA_N_L_0 = sigma * (1 - psi) * E_N_L_0 - (nu + eta_L + 1/gamma) * A_N_L_0 + xi_n * A_D_L_0 + xi_d * A_M_L_0
    dS_N_L_0 = sigma * psi * E_N_L_0 - (mu + 1/gamma) * S_N_L_0 + xi_n * S_D_L_0 + xi_d * S_M_L_0
    dT_N_L_0 = eta_L * A_N_L_0 + mu * S_N_L_0 - (rho + 1/gamma) * T_N_L_0 + xi_n * T_D_L_0 + xi_d * T_M_L_0
    dE_N_L_c = f_c * lambda_L_c * U_N_L - (sigma + 1/gamma) * E_N_L_c + xi_n * E_D_L_c + xi_d * E_M_L_c
    dA_N_L_c = sigma * (1 - psi) * E_N_L_c - (nu + eta_L + 1/gamma) * A_N_L_c + phi * rho * T_N_L_c + xi_n * A_D_L_c + xi_d * A_M_L_c
    dS_N_L_c = sigma * psi * E_N_L_c - (mu + 1/gamma) * S_N_L_c + xi_n * S_D_L_c + xi_d * S_M_L_c
    dT_N_L_c = eta_L * A_N_L_c + mu * S_N_L_c - (rho + 1/gamma) * T_N_L_c + w_c * rho * T_N_L_0 + xi_n * T_D_L_c + xi_d * T_M_L_c
    dE_N_L_t = f_t * lambda_L_t * U_N_L - (sigma + 1/gamma) * E_N_L_t + xi_n * E_D_L_t + xi_d * E_M_L_t
    dA_N_L_t = sigma * (1 - psi) * E_N_L_t - (nu + eta_L + 1/gamma) * A_N_L_t + xi_n * A_D_L_t + xi_d * A_M_L_t
    dS_N_L_t = sigma * psi * E_N_L_t - (mu + 1/gamma) * S_N_L_t + xi_n * S_D_L_t + xi_d * S_M_L_t
    dT_N_L_t = eta_L * A_N_L_t + mu * S_N_L_t - (rho + 1/gamma) * T_N_L_t + xi_n * T_D_L_t + xi_d * T_M_L_t
    dE_N_L_d2 = f_d2 * lambda_L_d2 * U_N_L - (sigma + 1/gamma) * E_N_L_d2 + xi_n * E_D_L_d2 + xi_d * E_M_L_d2
    dA_N_L_d2 = sigma * (1 - psi) * E_N_L_d2 - (nu + eta_L + 1/gamma) * A_N_L_d2 + phi * rho * T_N_L_d2 + xi_n * A_D_L_d2 + xi_d * A_M_L_d2
    dS_N_L_d2 = sigma * psi * E_N_L_d2 - (mu + 1/gamma) * S_N_L_d2 + xi_n * S_D_L_d2 + xi_d * S_M_L_d2
    dT_N_L_d2 = eta_L * A_N_L_d2 + mu * S_N_L_d2 - (rho + 1/gamma) * T_N_L_d2 + w_c * rho * T_N_L_t + xi_n * T_D_L_d2 + xi_d * T_M_L_d2
    
    # doxy-PEP (D)
    dU_D_L = eta_L * p_d * U_N_L + rho * (1 - w_c) * (T_D_L_0 + T_D_L_t) + rho * (1 - phi) * (T_D_L_c + T_D_L_d2) + nu * (A_D_L_0 + A_D_L_c + A_D_L_t + A_D_L_d2) - (e_d * lambda_L_0 + e_d * f_c * lambda_L_c + f_t * lambda_L_t + f_d2 * lambda_L_d2 + xi_n + p_v * eta_L + 1/gamma) * U_D_L + xi_d * U_V_L
    dE_D_L_0 = e_d * lambda_L_0 * U_D_L - (sigma + xi_n + 1/gamma) * E_D_L_0 + xi_d * E_V_L_0
    dA_D_L_0 = (1 - w_t) * sigma * (1 - psi) * E_D_L_0 - (nu + eta_L + xi_n + 1/gamma) * A_D_L_0 + xi_d * A_V_L_0
    dS_D_L_0 = (1 - w_t) * sigma * psi * E_D_L_0 - (mu + xi_n + 1/gamma) * S_D_L_0 + xi_d * S_V_L_0
    dT_D_L_0 = eta_L * A_D_L_0 + mu * S_D_L_0 - (rho + xi_n + 1/gamma) * T_D_L_0 + xi_d * T_V_L_0
    dE_D_L_c = e_d * f_c * lambda_L_c * U_D_L - (sigma + xi_n + 1/gamma) * E_D_L_c + xi_d * E_V_L_c
    dA_D_L_c = (1 - w_t) * sigma * (1 - psi) * E_D_L_c - (nu + eta_L + xi_n + 1/gamma) * A_D_L_c + phi * rho * T_D_L_c + xi_d * A_V_L_c
    dS_D_L_c = (1 - w_t) * sigma * psi * E_D_L_c - (mu + xi_n + 1/gamma) * S_D_L_c + xi_d * S_V_L_c
    dT_D_L_c = eta_L * A_D_L_c + mu * S_D_L_c - (rho + xi_n + 1/gamma) * T_D_L_c + w_c * rho * T_D_L_0 + xi_d * T_V_L_c
    dE_D_L_t = f_t * lambda_L_t * U_D_L - (sigma + xi_n + 1/gamma) * E_D_L_t + xi_d * E_V_L_t
    dA_D_L_t = w_t * sigma * (1 - psi) * E_D_L_0 + sigma * (1 - psi) * E_D_L_t - (nu + eta_L + xi_n + 1/gamma) * A_D_L_t + xi_d * A_V_L_t
    dS_D_L_t = w_t * sigma * psi * E_D_L_0 + sigma * psi * E_D_L_t - (mu + xi_n + 1/gamma) * S_D_L_t + xi_d * S_V_L_t
    dT_D_L_t = eta_L * A_D_L_t + mu * S_D_L_t - (rho + xi_n + 1/gamma) * T_D_L_t + xi_d * T_V_L_t
    dE_D_L_d2 = f_d2 * lambda_L_d2 * U_D_L - (sigma + xi_n + 1/gamma) * E_D_L_d2 + xi_d * E_V_L_d2
    dA_D_L_d2 = w_t * sigma * (1 - psi) * E_D_L_c  + sigma * (1 - psi) * E_D_L_d2 - (nu + eta_L + xi_n + 1/gamma) * A_D_L_d2 + phi * rho * T_D_L_d2 + xi_d * A_V_L_d2
    dS_D_L_d2 = w_t * sigma * psi * E_D_L_c + sigma * psi * E_D_L_d2 - (mu + xi_n + 1/gamma) * S_D_L_d2 + xi_d * S_V_L_d2
    dT_D_L_d2 = eta_L * A_D_L_d2 + mu * S_D_L_d2 - (rho + xi_n + 1/gamma) * T_D_L_d2 + w_c * rho * T_D_L_t + xi_d * T_V_L_d2
    
    # doxy-PEP + 4CMenB (V)
    dU_V_L = eta_L * p_v * U_D_L + eta_L * p_d * U_M_L + rho * (1 - w_c) * (T_V_L_0 + T_V_L_t) + rho * (1 - phi) * (T_V_L_c + T_V_L_d2) + nu * (A_V_L_0 + A_V_L_c + A_V_L_t + A_V_L_d2) - (e_vd * lambda_L_0 + e_vd * f_c * lambda_L_c + e_v * f_t * lambda_L_t + e_v * f_d2 * lambda_L_d2 + xi_n + xi_d + 1/gamma) * U_V_L
    dE_V_L_0 = e_vd * lambda_L_0 * U_V_L - (sigma + xi_n + xi_d + 1/gamma) * E_V_L_0
    dA_V_L_0 = (1 - w_t) * sigma * (1 - psi) * E_V_L_0 - (nu + eta_L + xi_n + xi_d + 1/gamma) * A_V_L_0
    dS_V_L_0 = (1 - w_t) * sigma * psi * E_V_L_0 - (mu + xi_n + xi_d + 1/gamma) * S_V_L_0
    dT_V_L_0 = eta_L * A_V_L_0 + mu * S_V_L_0 - (rho + xi_n + xi_d + 1/gamma) * T_V_L_0
    dE_V_L_c = e_vd * f_c * lambda_L_c * U_V_L - (sigma + xi_n + xi_d + 1/gamma) * E_V_L_c
    dA_V_L_c = (1 - w_t) * sigma * (1 - psi) * E_V_L_c - (nu + eta_L + xi_n + xi_d + 1/gamma) * A_V_L_c + phi * rho * T_V_L_c
    dS_V_L_c = (1 - w_t) * sigma * psi * E_V_L_c - (mu + xi_n + xi_d + 1/gamma) * S_V_L_c
    dT_V_L_c = eta_L * A_V_L_c + mu * S_V_L_c - (rho + xi_n + xi_d + 1/gamma) * T_V_L_c + w_c * rho * T_V_L_0
    dE_V_L_t = e_v * f_t * lambda_L_t * U_V_L - (sigma + xi_n + xi_d + 1/gamma) * E_V_L_t
    dA_V_L_t = w_t * sigma * (1 - psi) * E_V_L_0 + sigma * (1 - psi) * E_V_L_t - (nu + eta_L + xi_n + xi_d + 1/gamma) * A_V_L_t
    dS_V_L_t = w_t * sigma * psi * E_V_L_0 + sigma * psi * E_V_L_t - (mu + xi_n + xi_d + 1/gamma) * S_V_L_t
    dT_V_L_t = eta_L * A_V_L_t + mu * S_V_L_t - (rho + xi_n + xi_d + 1/gamma) * T_V_L_t
    dE_V_L_d2 = e_v * f_d2 * lambda_L_d2 * U_V_L - (sigma + xi_n + xi_d + 1/gamma) * E_V_L_d2
    dA_V_L_d2 = w_t * sigma * (1 - psi) * E_V_L_c  + sigma * (1 - psi) * E_V_L_d2 - (nu + eta_L + xi_n + xi_d + 1/gamma) * A_V_L_d2 + phi * rho * T_V_L_d2
    dS_V_L_d2 = w_t * sigma * psi * E_V_L_c + sigma * psi * E_V_L_d2 - (mu + xi_n + xi_d + 1/gamma) * S_V_L_d2
    dT_V_L_d2 = eta_L * A_V_L_d2 + mu * S_V_L_d2 - (rho + xi_n + xi_d + 1/gamma) * T_V_L_d2 + w_c * rho * T_V_L_t
    
    # 4CMenB (M)
    dU_M_L = eta_L * p_v * U_N_L + xi_n * U_V_L + rho * (1 - w_c) * (T_M_L_0 + T_M_L_t) + rho * (1 - phi) * (T_M_L_c + T_M_L_d2) + nu * (A_M_L_0 + A_M_L_c + A_M_L_t + A_M_L_d2) - (e_v * lambda_L_0 + e_v * f_c * lambda_L_c + e_v * f_t * lambda_L_t + e_v * f_d2 * lambda_L_d2 + p_d * eta_L + xi_d + 1/gamma) * U_M_L
    dE_M_L_0 = e_v * lambda_L_0 * U_M_L - (sigma + 1/gamma) * E_M_L_0 + xi_n * E_V_L_0 - xi_d * E_M_L_0
    dA_M_L_0 = sigma * (1 - psi) * E_M_L_0 - (nu + eta_L + 1/gamma) * A_M_L_0 + xi_n * A_V_L_0 - xi_d * A_M_L_0
    dS_M_L_0 = sigma * psi * E_M_L_0 - (mu + 1/gamma) * S_M_L_0 + xi_n * S_V_L_0 - xi_d * S_M_L_0
    dT_M_L_0 = eta_L * A_M_L_0 + mu * S_M_L_0 - (rho + 1/gamma) * T_M_L_0 + xi_n * T_V_L_0 - xi_d * T_M_L_0
    dE_M_L_c = e_v * f_c * lambda_L_c * U_M_L - (sigma + 1/gamma) * E_M_L_c + xi_n * E_V_L_c - xi_d * E_M_L_c
    dA_M_L_c = sigma * (1 - psi) * E_M_L_c - (nu + eta_L + 1/gamma) * A_M_L_c + phi * rho * T_M_L_c + xi_n * A_V_L_c - xi_d * A_M_L_c
    dS_M_L_c = sigma * psi * E_M_L_c - (mu + 1/gamma) * S_M_L_c + xi_n * S_V_L_c - xi_d * S_M_L_c
    dT_M_L_c = eta_L * A_M_L_c + mu * S_M_L_c - (rho + 1/gamma) * T_M_L_c + w_c * rho * T_M_L_0 + xi_n * T_V_L_c - xi_d * T_M_L_c
    dE_M_L_t = e_v * f_t * lambda_L_t * U_M_L - (sigma + 1/gamma) * E_M_L_t + xi_n * E_V_L_t - xi_d * E_M_L_t
    dA_M_L_t = sigma * (1 - psi) * E_M_L_t - (nu + eta_L + 1/gamma) * A_M_L_t + xi_n * A_V_L_t - xi_d * A_M_L_t
    dS_M_L_t = sigma * psi * E_M_L_t - (mu + 1/gamma) * S_M_L_t + xi_n * S_V_L_t - xi_d * S_M_L_t
    dT_M_L_t = eta_L * A_M_L_t + mu * S_M_L_t - (rho + 1/gamma) * T_M_L_t + xi_n * T_V_L_t - xi_d * T_M_L_t
    dE_M_L_d2 = e_v * f_d2 * lambda_L_d2 * U_M_L - (sigma + 1/gamma) * E_M_L_d2 + xi_n * E_V_L_d2 - xi_d * E_M_L_d2
    dA_M_L_d2 = sigma * (1 - psi) * E_M_L_d2 - (nu + eta_L + 1/gamma) * A_M_L_d2 + phi * rho * T_M_L_d2 + xi_n * A_V_L_d2 - xi_d * A_M_L_d2
    dS_M_L_d2 = sigma * psi * E_M_L_d2 - (mu + 1/gamma) * S_M_L_d2 + xi_n * S_V_L_d2 - xi_d * S_M_L_d2
    dT_M_L_d2 = eta_L * A_M_L_d2 + mu * S_M_L_d2 - (rho + 1/gamma) * T_M_L_d2 + w_c * rho * T_M_L_t + xi_n * T_V_L_d2 - xi_d * T_M_L_d2
    
    list(c(dU_N_H, dE_N_H_0, dA_N_H_0, dS_N_H_0, dT_N_H_0, dE_N_H_c, dA_N_H_c, dS_N_H_c, dT_N_H_c, dE_N_H_t, dA_N_H_t, dS_N_H_t, dT_N_H_t, dE_N_H_d2, dA_N_H_d2, dS_N_H_d2, dT_N_H_d2,
           dU_D_H, dE_D_H_0, dA_D_H_0, dS_D_H_0, dT_D_H_0, dE_D_H_c, dA_D_H_c, dS_D_H_c, dT_D_H_c, dE_D_H_t, dA_D_H_t, dS_D_H_t, dT_D_H_t, dE_D_H_d2, dA_D_H_d2, dS_D_H_d2, dT_D_H_d2,
           dU_V_H, dE_V_H_0, dA_V_H_0, dS_V_H_0, dT_V_H_0, dE_V_H_c, dA_V_H_c, dS_V_H_c, dT_V_H_c, dE_V_H_t, dA_V_H_t, dS_V_H_t, dT_V_H_t, dE_V_H_d2, dA_V_H_d2, dS_V_H_d2, dT_V_H_d2,
           dU_M_H, dE_M_H_0, dA_M_H_0, dS_M_H_0, dT_M_H_0, dE_M_H_c, dA_M_H_c, dS_M_H_c, dT_M_H_c, dE_M_H_t, dA_M_H_t, dS_M_H_t, dT_M_H_t, dE_M_H_d2, dA_M_H_d2, dS_M_H_d2, dT_M_H_d2,
           dU_N_L, dE_N_L_0, dA_N_L_0, dS_N_L_0, dT_N_L_0, dE_N_L_c, dA_N_L_c, dS_N_L_c, dT_N_L_c, dE_N_L_t, dA_N_L_t, dS_N_L_t, dT_N_L_t, dE_N_L_d2, dA_N_L_d2, dS_N_L_d2, dT_N_L_d2,
           dU_D_L, dE_D_L_0, dA_D_L_0, dS_D_L_0, dT_D_L_0, dE_D_L_c, dA_D_L_c, dS_D_L_c, dT_D_L_c, dE_D_L_t, dA_D_L_t, dS_D_L_t, dT_D_L_t, dE_D_L_d2, dA_D_L_d2, dS_D_L_d2, dT_D_L_d2,
           dU_V_L, dE_V_L_0, dA_V_L_0, dS_V_L_0, dT_V_L_0, dE_V_L_c, dA_V_L_c, dS_V_L_c, dT_V_L_c, dE_V_L_t, dA_V_L_t, dS_V_L_t, dT_V_L_t, dE_V_L_d2, dA_V_L_d2, dS_V_L_d2, dT_V_L_d2,
           dU_M_L, dE_M_L_0, dA_M_L_0, dS_M_L_0, dT_M_L_0, dE_M_L_c, dA_M_L_c, dS_M_L_c, dT_M_L_c, dE_M_L_t, dA_M_L_t, dS_M_L_t, dT_M_L_t, dE_M_L_d2, dA_M_L_d2, dS_M_L_d2, dT_M_L_d2
    ))
  })
}

# Extract once; use identical posterior rows in every scenario.
if (!file.exists(fit_file)) stop("Missing calibrated fit: ", fit_file)
fit_amr_negbin <- readRDS(fit_file)
pars <- c("beta", "phi_beta", "epsilon", "sigma", "psi", "mu",
          "eta_H_init", "omega", "phi_eta", "rho", "nu", "phi",
          "f_c", "f_t", "f_d2", "w_c", "w_t", "kappa_T", "kappa_S")
posterior_df <- as.data.frame(rstan::extract(fit_amr_negbin, pars = pars,
                                           permuted = TRUE))
samples_y <- rstan::extract(fit_amr_negbin, pars = "y", permuted = TRUE)
stopifnot(length(dim(samples_y$y)) == 3, dim(samples_y$y)[2] >= 8,
          dim(samples_y$y)[3] >= 34, nrow(posterior_df) >= n_iter)
initial_medians <- vapply(seq_len(34), function(j) median(samples_y$y[, 8, j]),
                          numeric(1))
stopifnot(all(is.finite(initial_medians)), all(is.finite(as.matrix(posterior_df))))
set.seed(seed)
# Identical to sample(1:6000, ...) when the fit has the original 6000 draws.
random_integers <- sample(seq_len(nrow(posterior_df)), n_iter, replace = FALSE)
rm(samples_y, fit_amr_negbin)

# The original common median initial state is computed only once, not per draw.
    # initial conditions

    # high-risk group
    # no intervention (N)
    U_N_H = initial_medians[1]
    E_N_H_0 = initial_medians[2] - 5
    A_N_H_0 = initial_medians[3]
    S_N_H_0 = initial_medians[4]
    T_N_H_0 = initial_medians[5]
    E_N_H_c = initial_medians[6] + 3
    A_N_H_c = initial_medians[7]
    S_N_H_c = initial_medians[8]
    T_N_H_c = initial_medians[9]
    E_N_H_t = initial_medians[10]
    A_N_H_t = initial_medians[11]
    S_N_H_t = initial_medians[12]
    T_N_H_t = initial_medians[13]
    E_N_H_d2 = initial_medians[14] + 2
    A_N_H_d2 = initial_medians[15]
    S_N_H_d2 = initial_medians[16]
    T_N_H_d2 = initial_medians[17]
    # doxy-PEP (D)
    U_D_H = 0 # 18
    E_D_H_0 = 0 # 19
    A_D_H_0 = 0 # 20
    S_D_H_0 = 0 # 21
    T_D_H_0 = 0 # 22
    E_D_H_c = 0 # 23
    A_D_H_c = 0 # 24
    S_D_H_c = 0 # 25
    T_D_H_c = 0 # 26
    E_D_H_t = 0 # 27
    A_D_H_t = 0 # 28
    S_D_H_t = 0 # 29
    T_D_H_t = 0 # 30
    E_D_H_d2 = 0 # 31
    A_D_H_d2 = 0 # 32
    S_D_H_d2 = 0 # 33
    T_D_H_d2 = 0 # 34
    # doxy-PEP + 4CMenB (V)
    U_V_H = 0 # 35
    E_V_H_0 = 0 # 36
    A_V_H_0 = 0 # 37
    S_V_H_0 = 0 # 38
    T_V_H_0 = 0 # 39
    E_V_H_c = 0 # 40
    A_V_H_c = 0 # 41
    S_V_H_c = 0 # 42
    T_V_H_c = 0 # 43
    E_V_H_t = 0 # 44
    A_V_H_t = 0 # 45
    S_V_H_t = 0 # 46
    T_V_H_t = 0 # 47
    E_V_H_d2 = 0 # 48
    A_V_H_d2 = 0 # 49
    S_V_H_d2 = 0 # 50
    T_V_H_d2 = 0 # 51
    # 4CMenB (M)
    U_M_H = 0 # 52
    E_M_H_0 = 0 # 53
    A_M_H_0 = 0 # 54
    S_M_H_0 = 0 # 55
    T_M_H_0 = 0 # 56
    E_M_H_c = 0 # 57
    A_M_H_c = 0 # 58
    S_M_H_c = 0 # 59
    T_M_H_c = 0 # 60
    E_M_H_t = 0 # 61
    A_M_H_t = 0 # 62
    S_M_H_t = 0 # 63
    T_M_H_t = 0 # 64
    E_M_H_d2 = 0 # 65
    A_M_H_d2 = 0 # 66
    S_M_H_d2 = 0 # 67
    T_M_H_d2 = 0 # 68
    # low-risk group
    # no intervention (N)
    U_N_L = initial_medians[18]  # 69
    E_N_L_0 = initial_medians[19] - 2
    A_N_L_0 = initial_medians[20]
    S_N_L_0 = initial_medians[21]
    T_N_L_0 = initial_medians[22]
    E_N_L_c = initial_medians[23] + 1
    A_N_L_c = initial_medians[24]
    S_N_L_c = initial_medians[25]
    T_N_L_c = initial_medians[26]
    E_N_L_t = initial_medians[27]
    A_N_L_t = initial_medians[28]
    S_N_L_t = initial_medians[29]
    T_N_L_t = initial_medians[30]
    E_N_L_d2 = initial_medians[31] + 1
    A_N_L_d2 = initial_medians[32]
    S_N_L_d2 = initial_medians[33]
    T_N_L_d2 = initial_medians[34]
    # doxy-PEP (D)
    U_D_L = 0
    E_D_L_0 = 0
    A_D_L_0 = 0
    S_D_L_0 = 0
    T_D_L_0 = 0
    E_D_L_c = 0
    A_D_L_c = 0
    S_D_L_c = 0
    T_D_L_c = 0
    E_D_L_t = 0
    A_D_L_t = 0
    S_D_L_t = 0
    T_D_L_t = 0
    E_D_L_d2 = 0
    A_D_L_d2 = 0
    S_D_L_d2 = 0
    T_D_L_d2 = 0
    # doxy-PEP + 4CMenB (V)
    U_V_L = 0
    E_V_L_0 = 0
    A_V_L_0 = 0
    S_V_L_0 = 0
    T_V_L_0 = 0
    E_V_L_c = 0
    A_V_L_c = 0
    S_V_L_c = 0
    T_V_L_c = 0
    E_V_L_t = 0
    A_V_L_t = 0
    S_V_L_t = 0
    T_V_L_t = 0
    E_V_L_d2 = 0
    A_V_L_d2 = 0
    S_V_L_d2 = 0
    T_V_L_d2 = 0
    # 4CMenB (M)
    U_M_L = 0
    E_M_L_0 = 0
    A_M_L_0 = 0
    S_M_L_0 = 0
    T_M_L_0 = 0
    E_M_L_c = 0
    A_M_L_c = 0
    S_M_L_c = 0
    T_M_L_c = 0
    E_M_L_t = 0
    A_M_L_t = 0
    S_M_L_t = 0
    T_M_L_t = 0
    E_M_L_d2 = 0
    A_M_L_d2 = 0
    S_M_L_d2 = 0
    T_M_L_d2 = 0
    
    y0 = c(U_N_H=U_N_H, E_N_H_0=E_N_H_0, A_N_H_0=A_N_H_0, S_N_H_0=S_N_H_0, T_N_H_0=T_N_H_0, E_N_H_c=E_N_H_c, A_N_H_c=A_N_H_c, S_N_H_c=S_N_H_c, T_N_H_c=T_N_H_c, E_N_H_t=E_N_H_t, A_N_H_t=A_N_H_t, S_N_H_t=S_N_H_t, T_N_H_t=T_N_H_t, E_N_H_d2=E_N_H_d2, A_N_H_d2=A_N_H_d2, S_N_H_d2=S_N_H_d2, T_N_H_d2=T_N_H_d2,
           U_D_H=U_D_H, E_D_H_0=E_D_H_0, A_D_H_0=A_D_H_0, S_D_H_0=S_D_H_0, T_D_H_0=T_D_H_0, E_D_H_c=E_D_H_c, A_D_H_c=A_D_H_c, S_D_H_c=S_D_H_c, T_D_H_c=T_D_H_c, E_D_H_t=E_D_H_t, A_D_H_t=A_D_H_t, S_D_H_t=S_D_H_t, T_D_H_t=T_D_H_t, E_D_H_d2=E_D_H_d2, A_D_H_d2=A_D_H_d2, S_D_H_d2=S_D_H_d2, T_D_H_d2=T_D_H_d2,
           U_V_H=U_V_H, E_V_H_0=E_V_H_0, A_V_H_0=A_V_H_0, S_V_H_0=S_V_H_0, T_V_H_0=T_V_H_0, E_V_H_c=E_V_H_c, A_V_H_c=A_V_H_c, S_V_H_c=S_V_H_c, T_V_H_c=T_V_H_c, E_V_H_t=E_V_H_t, A_V_H_t=A_V_H_t, S_V_H_t=S_V_H_t, T_V_H_t=T_V_H_t, E_V_H_d2=E_V_H_d2, A_V_H_d2=A_V_H_d2, S_V_H_d2=S_V_H_d2, T_V_H_d2=T_V_H_d2,
           U_M_H=U_M_H, E_M_H_0=E_M_H_0, A_M_H_0=A_M_H_0, S_M_H_0=S_M_H_0, T_M_H_0=T_M_H_0, E_M_H_c=E_M_H_c, A_M_H_c=A_M_H_c, S_M_H_c=S_M_H_c, T_M_H_c=T_M_H_c, E_M_H_t=E_M_H_t, A_M_H_t=A_M_H_t, S_M_H_t=S_M_H_t, T_M_H_t=T_M_H_t, E_M_H_d2=E_M_H_d2, A_M_H_d2=A_M_H_d2, S_M_H_d2=S_M_H_d2, T_M_H_d2=T_M_H_d2,
           U_N_L=U_N_L, E_N_L_0=E_N_L_0, A_N_L_0=A_N_L_0, S_N_L_0=S_N_L_0, T_N_L_0=T_N_L_0, E_N_L_c=E_N_L_c, A_N_L_c=A_N_L_c, S_N_L_c=S_N_L_c, T_N_L_c=T_N_L_c, E_N_L_t=E_N_L_t, A_N_L_t=A_N_L_t, S_N_L_t=S_N_L_t, T_N_L_t=T_N_L_t, E_N_L_d2=E_N_L_d2, A_N_L_d2=A_N_L_d2, S_N_L_d2=S_N_L_d2, T_N_L_d2=T_N_L_d2,
           U_D_L=U_D_L, E_D_L_0=E_D_L_0, A_D_L_0=A_D_L_0, S_D_L_0=S_D_L_0, T_D_L_0=T_D_L_0, E_D_L_c=E_D_L_c, A_D_L_c=A_D_L_c, S_D_L_c=S_D_L_c, T_D_L_c=T_D_L_c, E_D_L_t=E_D_L_t, A_D_L_t=A_D_L_t, S_D_L_t=S_D_L_t, T_D_L_t=T_D_L_t, E_D_L_d2=E_D_L_d2, A_D_L_d2=A_D_L_d2, S_D_L_d2=S_D_L_d2, T_D_L_d2=T_D_L_d2,
           U_V_L=U_V_L, E_V_L_0=E_V_L_0, A_V_L_0=A_V_L_0, S_V_L_0=S_V_L_0, T_V_L_0=T_V_L_0, E_V_L_c=E_V_L_c, A_V_L_c=A_V_L_c, S_V_L_c=S_V_L_c, T_V_L_c=T_V_L_c, E_V_L_t=E_V_L_t, A_V_L_t=A_V_L_t, S_V_L_t=S_V_L_t, T_V_L_t=T_V_L_t, E_V_L_d2=E_V_L_d2, A_V_L_d2=A_V_L_d2, S_V_L_d2=S_V_L_d2, T_V_L_d2=T_V_L_d2,
           U_M_L=U_M_L, E_M_L_0=E_M_L_0, A_M_L_0=A_M_L_0, S_M_L_0=S_M_L_0, T_M_L_0=T_M_L_0, E_M_L_c=E_M_L_c, A_M_L_c=A_M_L_c, S_M_L_c=S_M_L_c, T_M_L_c=T_M_L_c, E_M_L_t=E_M_L_t, A_M_L_t=A_M_L_t, S_M_L_t=S_M_L_t, T_M_L_t=T_M_L_t, E_M_L_d2=E_M_L_d2, A_M_L_d2=A_M_L_d2, S_M_L_d2=S_M_L_d2, T_M_L_d2=T_M_L_d2)
    

run_amr <- function(p_d, p_v, efficacy_doxypep, efficacy_vaccine) {
  efficacy_combined = 1 - (1- efficacy_doxypep) * (1 - efficacy_vaccine)
  
  # Annual MSM population entrants (at age 15)
  alpha <- 12000
  
  # Proportion of the MSM population in group j
  q_H <- 0.15
  q_L <- 0.85
  
  # Annual rate of partner change in group j
  c_H <- 15.6
  c_L <- 0.6
  
  # Years spent in the sexually-active population
  gamma <- 50
  
  # Efficacy of doxycycline against infections
  e_d <- 1 - efficacy_doxypep
  
  # Efficacy of 4cmenb against infections
  e_v <- 1 - efficacy_vaccine
  
  # Combined efficacy of 4cmenb and doxypep against infections
  e_vd <- 1 - efficacy_combined
  
  # Discontinuation rate of doxy-PEP for syphilis
  xi_n <- 0.362
  
  # Duration of protection of 4cmenb
  xi_d <- 0.20
  
  cases_all <- matrix(NA, nrow = n_iter, ncol = n_years)
  cases_0 <- matrix(NA, nrow = n_iter, ncol = n_years)
  cases_c <- matrix(NA, nrow = n_iter, ncol = n_years)
  cases_t <- matrix(NA, nrow = n_iter, ncol = n_years)
  cases_d2 <- matrix(NA, nrow = n_iter, ncol = n_years)
  for (i in 1:n_iter) {
    # times
    t <- seq(8, 8+n_years+1, by = 1)
    t_0 = 0 
    t <- t[-1]
    
    # get index for posterior sample 
    idx <- random_integers[i]
    
    params <- list(
      q_H = q_H, c_H = c_H, c_L = c_L, q_L = q_L,
      t_0 = t_0, alpha = alpha, gamma = gamma,
      beta = posterior_df$beta[idx], phi_beta = posterior_df$phi_beta[idx], epsilon=posterior_df$epsilon[idx], sigma = posterior_df$sigma[idx], psi = posterior_df$psi[idx],
      mu=posterior_df$mu[idx], eta_H_init=posterior_df$eta_H_init[idx], omega=posterior_df$omega[idx], phi_eta=posterior_df$phi_eta[idx], rho=posterior_df$rho[idx],
      nu=posterior_df$nu[idx], phi=posterior_df$phi[idx], f_c=posterior_df$f_c[idx], f_t=posterior_df$f_t[idx], f_d2=posterior_df$f_d2[idx], w_c=posterior_df$w_c[idx], w_t=posterior_df$w_t[idx],
      e_d = e_d, e_v = e_v, e_vd = e_vd, xi_n = xi_n, xi_d = xi_d, p_d = p_d, p_v = p_v
    )
    
    # solve the system
    out <- ode(y = y0, times = t, func = amr_model, parms = params)
    if (nrow(out) != length(t) || any(!is.finite(out))) {
      stop("Incomplete or non-finite ODE output for draw ", idx)
    }
    out <- as.data.frame(out)
    
    # compute incidences and prescriptions
    incidence_all <- numeric(n_years)
    incidence_0 <- numeric(n_years)
    incidence_c <- numeric(n_years)
    incidence_t <- numeric(n_years)
    incidence_d2 <- numeric(n_years)
    
    for (t in 1:(n_years)) {
      isFixed = TRUE
      
      C_H_0 = get_C(out[t,3],  out[t,4],  out[t,5],
                    out[t,20], out[t,21], out[t,22],
                    out[t,37], out[t,38], out[t,39],
                    out[t,54], out[t,55], out[t,56])
      
      C_L_0 = get_C(out[t,71], out[t,72], out[t,73],
                    out[t,88], out[t,89], out[t,90],
                    out[t,105],out[t,106],out[t,107],
                    out[t,122],out[t,123],out[t,124])
      
      C_H_c = get_C(out[t,7],  out[t,8],  out[t,9],
                    out[t,24], out[t,25], out[t,26],
                    out[t,41], out[t,42], out[t,43],
                    out[t,58], out[t,59], out[t,60])
      
      C_L_c = get_C(out[t,75], out[t,76], out[t,77],
                    out[t,92], out[t,93], out[t,94],
                    out[t,109],out[t,110],out[t,111],
                    out[t,126],out[t,127],out[t,128])
      
      C_H_t = get_C(out[t,11], out[t,12], out[t,13],
                    out[t,28], out[t,29], out[t,30],
                    out[t,45], out[t,46], out[t,47],
                    out[t,62], out[t,63], out[t,64])
      
      C_L_t = get_C(out[t,79], out[t,80], out[t,81],
                    out[t,96], out[t,97], out[t,98],
                    out[t,113],out[t,114],out[t,115],
                    out[t,130],out[t,131],out[t,132])
      
      C_H_d2 = get_C(out[t,15], out[t,16], out[t,17],
                     out[t,32], out[t,33], out[t,34],
                     out[t,49], out[t,50], out[t,51],
                     out[t,66], out[t,67], out[t,68])
      
      C_L_d2 = get_C(out[t,83], out[t,84], out[t,85],
                     out[t,100],out[t,101],out[t,102],
                     out[t,117],out[t,118],out[t,119],
                     out[t,134],out[t,135],out[t,136])
      
      N_H <- get_N(
        out[t,2], out[t,3], out[t,4], out[t,5], out[t,6],
        out[t,7], out[t,8], out[t,9], out[t,10],
        out[t,11], out[t,12], out[t,13], out[t,14],
        out[t,15], out[t,16], out[t,17], out[t,18],
        
        out[t,19], out[t,20], out[t,21], out[t,22], out[t,23],
        out[t,24], out[t,25], out[t,26], out[t,27],
        out[t,28], out[t,29], out[t,30], out[t,31],
        out[t,32], out[t,33], out[t,34], out[t,35],
        
        out[t,36], out[t,37], out[t,38], out[t,39], out[t,40],
        out[t,41], out[t,42], out[t,43], out[t,44],
        out[t,45], out[t,46], out[t,47], out[t,48],
        out[t,49], out[t,50], out[t,51], out[t,52],
        
        out[t,53], out[t,54], out[t,55], out[t,56], out[t,57],
        out[t,58], out[t,59], out[t,60], out[t,61],
        out[t,62], out[t,63], out[t,64], out[t,65],
        out[t,66], out[t,67], out[t,68], out[t,69]
      )
      
      N_L <- get_N(
        out[t,70], out[t,71], out[t,72], out[t,73], out[t,74],
        out[t,75], out[t,76], out[t,77], out[t,78],
        out[t,79], out[t,80], out[t,81], out[t,82],
        out[t,83], out[t,84], out[t,85], out[t,86],
        
        out[t,87], out[t,88], out[t,89], out[t,90], out[t,91],
        out[t,92], out[t,93], out[t,94], out[t,95],
        out[t,96], out[t,97], out[t,98], out[t,99],
        out[t,100], out[t,101], out[t,102], out[t,103],
        
        out[t,104], out[t,105], out[t,106], out[t,107], out[t,108],
        out[t,109], out[t,110], out[t,111], out[t,112],
        out[t,113], out[t,114], out[t,115], out[t,116],
        out[t,117], out[t,118], out[t,119], out[t,120],
        
        out[t,121], out[t,122], out[t,123], out[t,124], out[t,125],
        out[t,126], out[t,127], out[t,128], out[t,129],
        out[t,130], out[t,131], out[t,132], out[t,133],
        out[t,134], out[t,135], out[t,136], out[t,137]
      ) 
      
      pi_H <- get_pi(c_H, N_H, c_L, N_L)
      pi_L <- get_pi(c_L, N_L, c_H, N_H)
      
      lambda_H_0 = get_lambda(t+8, t_0, c_H, params$beta, params$phi_beta, params$epsilon, C_H_0, N_H, pi_H, C_L_0, N_L, pi_L, isFixed)
      lambda_L_0 = get_lambda(t+8, t_0, c_L, params$beta, params$phi_beta, params$epsilon, C_L_0, N_L, pi_L, C_H_0, N_H, pi_H, isFixed)
      lambda_H_c = get_lambda(t+8, t_0, c_H, params$beta, params$phi_beta, params$epsilon, C_H_c, N_H, pi_H, C_L_c, N_L, pi_L, isFixed)
      lambda_L_c = get_lambda(t+8, t_0, c_L, params$beta, params$phi_beta, params$epsilon, C_L_c, N_L, pi_L, C_H_c, N_H, pi_H, isFixed)
      lambda_H_t = get_lambda(t+8, t_0, c_H, params$beta, params$phi_beta, params$epsilon, C_H_t, N_H, pi_H, C_L_t, N_L, pi_L, isFixed)
      lambda_L_t = get_lambda(t+8, t_0, c_L, params$beta, params$phi_beta, params$epsilon, C_L_t, N_L, pi_L, C_H_t, N_H, pi_H, isFixed)
      lambda_H_d2 = get_lambda(t+8, t_0, c_H, params$beta, params$phi_beta, params$epsilon, C_H_d2, N_H, pi_H, C_L_d2, N_L, pi_L, isFixed)
      lambda_L_d2 = get_lambda(t+8, t_0, c_L, params$beta, params$phi_beta, params$epsilon, C_L_d2, N_L, pi_L, C_H_d2, N_H, pi_H, isFixed)
      
      # Trapezoidal rule: (f(a) + f(b)) / 2 * (b - a)
      E_N_0 = 0.5 * lambda_H_0 * (out[t, 2] + out[t + 1, 2]) + 0.5 * lambda_L_0 * (out[t, 2+68] + out[t + 1, 2+68])
      E_D_0 = 0.5 * e_d * lambda_H_0 * (out[t, 19] + out[t + 1, 19]) + 0.5 * e_d * lambda_L_0 * (out[t, 19+68] + out[t + 1, 19+68])
      E_V_0 = 0.5 * e_vd * lambda_H_0 * (out[t, 36] + out[t + 1, 36]) +  0.5 * e_vd * lambda_L_0 * (out[t, 36+68] + out[t + 1, 36+68])
      E_M_0 = 0.5 * e_v * lambda_H_0 * (out[t, 53] + out[t + 1, 53]) + 0.5 * e_v * lambda_L_0 * (out[t, 53+68] + out[t + 1, 53+68])
      incidence_0[t] = E_N_0 + E_D_0 + E_V_0 + E_M_0
      
      E_N_c = 0.5 * params$f_c * lambda_H_c * (out[t, 2] + out[t + 1, 2]) + 0.5 * params$f_c * lambda_L_c * (out[t, 2+68] + out[t + 1, 2+68])
      E_D_c = 0.5 * e_d * params$f_c * lambda_H_c * (out[t, 19] + out[t + 1, 19]) + 0.5 * e_d * params$f_c * lambda_L_c * (out[t, 19+68] + out[t + 1, 19+68])
      E_V_c = 0.5 * e_vd * params$f_c * lambda_H_c * (out[t, 36] + out[t + 1, 36]) +  0.5 * e_vd * params$f_c * lambda_L_c * (out[t, 36+68] + out[t + 1, 36+68])
      E_M_c = 0.5 * e_v * params$f_c * lambda_H_c * (out[t, 53] + out[t + 1, 53]) + 0.5 * e_v * params$f_c * lambda_L_c * (out[t, 53+68] + out[t + 1, 53+68])
      incidence_c[t] = E_N_c + E_D_c + E_V_c + E_M_c
      
      E_N_t = 0.5 * params$f_t * lambda_H_t * (out[t, 2] + out[t + 1, 2]) + 0.5 * params$f_t * lambda_L_t * (out[t, 2+68] + out[t + 1, 2+68])
      E_D_t = 0.5 * params$f_t * lambda_H_t * (out[t, 19] + out[t + 1, 19]) + 0.5 * params$f_t * lambda_L_t * (out[t, 19+68] + out[t + 1, 19+68])
      E_V_t = 0.5 * e_v * params$f_t * lambda_H_t * (out[t, 36] + out[t + 1, 36]) +  0.5 * e_v * params$f_t * lambda_L_t * (out[t, 36+68] + out[t + 1, 36+68])
      E_M_t = 0.5 * e_v * params$f_t * lambda_H_t * (out[t, 53] + out[t + 1, 53]) + 0.5 * e_v * params$f_t * lambda_L_t * (out[t, 53+68] + out[t + 1, 53+68])
      incidence_t[t] = E_N_t + E_D_t + E_V_t + E_M_t
      
      E_N_d2 = 0.5 * params$f_d2 * lambda_H_d2 * (out[t, 2] + out[t + 1, 2]) + 0.5 * params$f_d2 * lambda_L_d2 * (out[t, 2+68] + out[t + 1, 2+68])
      E_D_d2 = 0.5 * params$f_d2 * lambda_H_d2 * (out[t, 19] + out[t + 1, 19]) + 0.5 * params$f_d2 * lambda_L_d2 * (out[t, 19+68] + out[t + 1, 19+68])
      E_V_d2 = 0.5 * e_v * params$f_d2 * lambda_H_d2 * (out[t, 36] + out[t + 1, 36]) +  0.5 * e_v * params$f_d2 * lambda_L_d2 * (out[t, 36+68] + out[t + 1, 36+68])
      E_M_d2 = 0.5 * e_v * params$f_d2 * lambda_H_d2 * (out[t, 53] + out[t + 1, 53]) + 0.5 * e_v * params$f_d2 * lambda_L_d2 * (out[t, 53+68] + out[t + 1, 53+68])
      incidence_d2[t] = E_N_d2 + E_D_d2 + E_V_d2 + E_M_d2
      
      incidence_all[t] = incidence_0[t] + incidence_c[t] + incidence_t[t] + incidence_d2[t]
    }
    
    cases_all[i,] <- incidence_all
    cases_0[i,] <- incidence_0
    cases_c[i,] <- incidence_c
    cases_t[i,] <- incidence_t
    cases_d2[i,] <- incidence_d2
  }
  
  return(list(
    cases_all = cases_all,
    cases_0 = cases_0,
    cases_c = cases_c,
    cases_t = cases_t,
    cases_d2 = cases_d2
  ))
}


validate_cases <- function(x, negative_tolerance = 1e-5) {
  if (!is.matrix(x) || !identical(dim(x), c(as.integer(n_iter), as.integer(n_years)))) {
    stop("Unexpected annual incidence dimensions: ", paste(dim(x), collapse = " x "))
  }
  if (any(!is.finite(x))) stop("Non-finite annual incidence values")
  if (any(x < -negative_tolerance)) {
    stop(sprintf("Annual incidence minimum %.9g is below tolerance -%.9g; inspect res before continuing.",
                 min(x), negative_tolerance))
  }
  # An absolute tolerance in infections, not a fraction of the largest count.
  # Handles negligible negative numerical error; larger negatives still stop.
  if (any(x < 0)) {
    message(sprintf("Set %d near-zero negative incidence values to zero (minimum %.9g).",
                    sum(x < 0), min(x)))
    x[x < 0] <- 0
  }
  x
}
summary3 <- function(x) {
  c(median = median(x), lower = unname(quantile(x, 0.025)),
    upper = unname(quantile(x, 0.975)))
}
# Every panel reports ANNUAL infections in 2041, not cumulative infections.
# Tet-R excludes Dual-R; each is a distinct compartment in the original model.
outcome_columns <- c(Total = "cases_all", Susceptible = "cases_0",
                     `Tet-R` = "cases_t", `Cef-R` = "cases_c", `Dual-R` = "cases_d2")
summarise_2041 <- function(res) {
  do.call(rbind, lapply(names(outcome_columns), function(outcome) {
    annual <- validate_cases(res[[outcome_columns[[outcome]]]])
    q <- summary3(annual[, idx_2041])
    data.frame(outcome = outcome, year = target_year,
               median = q[[1]], lower = q[[2]], upper = q[[3]])
  }))
}

# The same posterior rows and initial conditions are used in all four strategies.
strategy_levels <- c("No Intervention", "Doxy-PEP", "4CMenB", "Doxy-PEP + 4CMenB")
resume_saved_scenarios <- FALSE # Full runs use the current fitted model.
# Keep the incidence tolerance at 1e-5. Retry only posterior rows exceeding it.
# Ordinary near-zero cleanup remains in validate_cases.
retry_incidence_rows <- function(res, p_d, p_v, de, ve) {
  keys <- c("cases_all", "cases_0", "cases_t", "cases_c", "cases_d2")
  bad <- integer(0)
  for (key in keys) {
    x <- res[[key]]
    if (!is.matrix(x) || !identical(dim(x), c(as.integer(n_iter), as.integer(n_years))))
      stop("Unexpected dimensions in ", key, "; cannot retry automatically.")
    bad <- union(bad, which(rowSums(!is.finite(x) | x < -1e-5) > 0))
  }
  if (!length(bad)) return(res)
  message("Retrying ", length(bad), " affected posterior draw(s) with tighter ODE tolerances.")
  # Clone the existing runner; preserve all equations, initial conditions and
  # incidence calculations. Only numerical solver tolerances change.
  tightened <- run_amr
  replacements <- 0L
  replace_ode <- function(expr) {
    if (!is.call(expr)) return(expr)
    if (identical(expr[[1]], as.name("ode"))) {
      expr$rtol <- 1e-9
      expr$atol <- 1e-10
      expr$maxsteps <- 50000L
      replacements <<- replacements + 1L
      return(expr)
    }
    for (k in seq_along(expr)[-1L]) expr[[k]] <- replace_ode(expr[[k]])
    expr
  }
  body(tightened) <- replace_ode(body(tightened))
  if (replacements != 1L) stop("Could not identify the single ODE call in run_amr.")
  retry_environment <- new.env(parent = environment(run_amr))
  retry_environment$n_iter <- length(bad)
  retry_environment$random_integers <- random_integers[bad]
  environment(tightened) <- retry_environment
  rerun <- tightened(p_d = p_d, p_v = p_v,
                     efficacy_doxypep = de, efficacy_vaccine = ve)
  for (key in keys) {
    if (!is.matrix(rerun[[key]]) ||
        !identical(dim(rerun[[key]]), c(as.integer(length(bad)), as.integer(n_years))))
      stop("Unexpected retry dimensions in ", key)
    if (any(!is.finite(rerun[[key]])) || any(rerun[[key]] < -1e-5)) {
      stop("Tighter solve still fails in ", key,
           ". Do not widen the tolerance; inspect posterior rows: ",
           paste(random_integers[bad], collapse = ", "))
    }
  }
  for (key in keys) res[[key]][bad, ] <- rerun[[key]]
  attr(res, "retried_posterior_rows") <- random_integers[bad]
  res
}
run_and_summarise <- function(strategy, p_d, p_v, de, filename) {
  message(sprintf("%s: doxy-PEP efficacy %.0f%%; 4CMenB efficacy %.0f%%",
                  strategy, 100 * de, 100 * efficacy_vaccine_fixed))
  checkpoint <- file.path(output_dir, filename)
  if (resume_saved_scenarios && file.exists(checkpoint)) {
    saved <- readRDS(checkpoint)
    same <- function(a, b) isTRUE(all.equal(a, b, check.attributes = FALSE))
    if (!same(saved$years, years) || !same(saved$posterior_rows, random_integers) ||
        !same(saved$p_d, p_d) || !same(saved$p_v, p_v) ||
        !same(saved$efficacy_doxypep, de) ||
        !same(saved$efficacy_vaccine, efficacy_vaccine_fixed))
      stop("Saved scenario settings differ: ", checkpoint)
    res <- saved[c("cases_all", "cases_0", "cases_c", "cases_t", "cases_d2")]
    message("Reusing saved scenario: ", basename(checkpoint))
  } else {
    res <- run_amr(p_d = p_d, p_v = p_v,
                   efficacy_doxypep = de, efficacy_vaccine = efficacy_vaccine_fixed)
  }
  # Raw total and all strain trajectories are checkpointed before validation.
  # Inspect this file if validation stops; the numerical guard is unchanged.
  saveRDS(c(list(Strategy = strategy, p_d = p_d, p_v = p_v,
                 efficacy_doxypep = de, efficacy_vaccine = efficacy_vaccine_fixed,
                 years = years, posterior_rows = random_integers), res),
          file.path(output_dir, filename))
  corrected <- retry_incidence_rows(res, p_d, p_v, de, efficacy_vaccine_fixed)
  if (length(attr(corrected, "retried_posterior_rows"))) {
    backup <- paste0(checkpoint, ".before_retry.rds")
    if (!file.exists(backup) && !file.copy(checkpoint, backup))
      stop("Could not preserve original scenario before retry")
    res <- corrected
    saveRDS(c(list(Strategy = strategy, p_d = p_d, p_v = p_v,
                   efficacy_doxypep = de, efficacy_vaccine = efficacy_vaccine_fixed,
                   years = years, posterior_rows = random_integers,
                   retried_posterior_rows = attr(corrected, "retried_posterior_rows")), res),
            checkpoint)
  }
  ans <- summarise_2041(res)
  ans$Strategy <- strategy
  ans$efficacy_doxypep <- de
  ans$efficacy_vaccine <- efficacy_vaccine_fixed
  ans
}

# With no doxy-PEP uptake, its efficacy has no effect: run each fixed strategy once.
non_summary <- run_and_summarise("No Intervention", 0, 0, 0,
                                 "no_intervention_annual.rds")
vaccine_summary <- run_and_summarise("4CMenB", 0, p_v_combined, 0,
                                     "vaccine_only_annual.rds")
# Repeat summaries over x coordinates only, not simulations. This makes horizontal
# median curves and uncertainty bands for the two fixed strategies.
expand_fixed_strategy <- function(summary_df) {
  do.call(rbind, lapply(efficacy_doxypep, function(de) {
    ans <- summary_df
    ans$efficacy_doxypep <- de
    ans
  }))
}
results <- list(expand_fixed_strategy(non_summary),
                expand_fixed_strategy(vaccine_summary))
for (j in seq_along(efficacy_doxypep)) {
  de <- efficacy_doxypep[j]
  message(sprintf("Doxy-PEP efficacy value %d/%d", j, length(efficacy_doxypep)))
  results[[length(results) + 1L]] <- run_and_summarise(
    "Doxy-PEP", p_d_combined, 0, de, sprintf("doxypep_only_DE_%0.8f.rds", de))
  results[[length(results) + 1L]] <- run_and_summarise(
    "Doxy-PEP + 4CMenB", p_d_combined, p_v_combined, de,
    sprintf("combined_DE_%0.8f.rds", de))
}
final_df <- do.call(rbind, results)
final_df$Strategy <- factor(final_df$Strategy, levels = strategy_levels)
final_df <- final_df[order(match(final_df$outcome, names(outcome_columns)),
                           final_df$Strategy, final_df$efficacy_doxypep), ]
rownames(final_df) <- NULL
write.csv(final_df, file.path(output_dir, "sensitivity_summary.csv"), row.names = FALSE)
write.csv(non_summary, file.path(output_dir, "no_intervention_summary.csv"), row.names = FALSE)
write.csv(vaccine_summary, file.path(output_dir, "vaccine_only_summary.csv"), row.names = FALSE)

# Preserve theme, sizes, ribbons and stacked panel layout from the working file.
# Colours now identify strategies consistently across all five panels.
# A shared legend is necessary to identify the four strategy curves.
# load(file.path(
#   output_dir,
#   "workspace_fixed_sensitivity_doxypep_four_strategies_final.RData"
# ))

strategy_colours <- c("No Intervention" = "grey35", "Doxy-PEP" = "#0072B2",
                      "Vaccination" = "#009E73", "Doxy-PEP + Vaccination" = "#D55E00")
strategy_linetypes <- c("No Intervention" = "dashed", "Doxy-PEP" = "solid",
                        "Vaccination" = "dotdash", "Doxy-PEP + Vaccination" = "longdash")
make_incidence_plot <- function(outcome_name, y_label) {
  dat <- final_df[final_df$outcome == outcome_name, ]
  ggplot(dat, aes(x = 100 * efficacy_doxypep, y = median,
                  colour = Strategy, fill = Strategy, linetype = Strategy,
                  group = Strategy)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.18, colour = NA, show.legend = FALSE) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    scale_colour_manual(name = NULL, values = strategy_colours, drop = FALSE) +
    scale_fill_manual(name = NULL, values = strategy_colours, drop = FALSE) +
    scale_linetype_manual(name = NULL, values = strategy_linetypes, drop = FALSE) +
    scale_y_continuous(labels = scales::label_number(accuracy = 1, big.mark = ",")) +
    labs(x = "Doxy-PEP Efficacy (%)", y = y_label) +
    theme_classic(base_size = 7) + theme(legend.position = "bottom")
}
plot_total <- make_incidence_plot("Total", "Total Infections\nin 2041")
plot_susceptible <- make_incidence_plot("Susceptible", "Susceptible Infections\nin 2041")
plot_tetr <- make_incidence_plot("Tet-R", "Tet-R Infections\nin 2041")
plot_cefr <- make_incidence_plot("Cef-R", "Cef-R Infections\nin 2041")
plot_dualr <- make_incidence_plot("Dual-R", "Dual-R Infections\nin 2041")
plot_five_panels <- (
  plot_total / plot_susceptible / plot_tetr / plot_cefr / plot_dualr
) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 9), plot.tag.position = c(0, 1),
        legend.position = "right")
print(plot_five_panels)

# As in the supplied working file, figure exports are optional (commented out).
# Uncomment to export the five-panel figure.
# for (extension in c("pdf", "png")) {
#   ggsave(file.path(output_dir, paste0("figure_doxypep_four_strategies_five_panels.", extension)),
#          plot_five_panels, width = 7, height = 15, dpi = 300)
# }
# save.image(file.path(output_dir, "workspace_fixed_sensitivity_doxypep_four_strategies_final.RData"))
message("Saved summaries, annual posterior draws and workspace to: ", output_dir)
# Lines: posterior medians; ribbons: 95% posterior intervals.
# All draws share the original median initial state. Total medians are calculated
# from cases_all, not by adding the separate strain medians.
# Default cost: 1 No Intervention + 1 vaccine-only + 11 doxy-PEP-only + 11 combined.

# plot in 4 columns
strategy_levels <- c(
  "No Intervention",
  "Doxy-PEP",
  "Vaccination",
  "Doxy-PEP + Vaccination"
)

final_df$Strategy <- factor(
  final_df$Strategy,
  levels = strategy_levels
)

strategy_colours <- c(
  "No Intervention" = "grey35",
  "Doxy-PEP" = "#0072B2",
  "Vaccination" = "#009E73",
  "Doxy-PEP + Vaccination" = "#D55E00"
)

strategy_linetypes <- c(
  "No Intervention" = "dashed",
  "Doxy-PEP" = "solid",
  "Vaccination" = "dotdash",
  "Doxy-PEP + Vaccination" = "longdash"
)

make_incidence_plot <- function(outcome_name, y_label,
                                show_headers = FALSE) {
  
  dat <- final_df[final_df$outcome == outcome_name, ]
  
  p <- ggplot(
    dat,
    aes(
      x = 100 * efficacy_doxypep,
      y = median,
      colour = Strategy,
      fill = Strategy,
      linetype = Strategy,
      group = Strategy
    )
  ) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      alpha = 0.18,
      colour = NA,
      show.legend = FALSE
    ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    facet_grid(
      cols = vars(Strategy),
      scales = "fixed",
      drop = FALSE,
      labeller = labeller(
        Strategy = c(
          "No Intervention" = "No Intervention",
          "Doxy-PEP" = "Doxy-PEP",
          "Vaccination" = "Vaccination",
          "Doxy-PEP + Vaccination" = "Doxy-PEP +\nVaccination"
        )
      )
    ) +
    scale_colour_manual(values = strategy_colours, drop = FALSE) +
    scale_fill_manual(values = strategy_colours, drop = FALSE) +
    scale_linetype_manual(values = strategy_linetypes, drop = FALSE) +
    scale_x_continuous(
      breaks = c(0, 25, 50, 75, 100)
    ) +
    scale_y_continuous(
      labels = scales::label_number(
        accuracy = 1,
        big.mark = ","
      )
    ) +
    labs(
      x = "Doxy-PEP Efficacy (%)",
      y = y_label
    ) +
    theme_classic(base_size = 7) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text.x = element_text(size = 8, face = "bold"),
      panel.spacing.x = grid::unit(5, "mm")
    )
  
  # Show intervention headings only above the first row.
  if (!show_headers) {
    p <- p + theme(strip.text.x = element_blank())
  }
  
  p
}

plot_total <- make_incidence_plot(
  "Total", "Total Infections\nin 2041",
  show_headers = TRUE
)

plot_susceptible <- make_incidence_plot(
  "Susceptible", "Susceptible Infections\nin 2041"
)

plot_tetr <- make_incidence_plot(
  "Tet-R", "Tet-R Infections\nin 2041"
)

plot_cefr <- make_incidence_plot(
  "Cef-R", "Cef-R Infections\nin 2041"
)

plot_dualr <- make_incidence_plot(
  "Dual-R", "Dual-R Infections\nin 2041"
)

plot_five_panels <- (
  plot_total /
    plot_susceptible /
    plot_tetr /
    plot_cefr /
    plot_dualr
) +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(size = 9),
    plot.tag.position = c(0, 1),
    legend.position = "none"
  )

print(plot_five_panels)
