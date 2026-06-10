library(deSolve)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(writexl)

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
    isFixed <- FALSE
    
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

run_amr <- function(p_d, p_v) {
  print(p_d)
  print(p_v)
  
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
  e_d <- 1 - 0.55
  
  # Efficacy of 4cmenb against infections
  e_v <- 1 - 0.40
  
  # Combined efficacy of 4cmenb and doxypep against infections
  e_vd <- 1 - 0.73
  
  # Discontinuation rate of doxy-PEP for syphilis
  xi_n <- 0.362
  
  # Duration of protection of 4cmenb
  xi_d <- 0.20
  
  # load calibrated parameters
  fit_amr_negbin <- readRDS("fit_results_fixedinitialstate_UK_6years_allcases_covid_123.rds")
  
  # extract all posterior samples for these parameters as a list (permuted = TRUE merges chains)
  pars=c('beta', 'phi_beta', 'epsilon', 'sigma', 'psi', 'mu', 'eta_H_init', 'omega', 'phi_eta', 'rho', 'nu', 'phi', 'f_c', 'f_t', 'f_d2', 'w_c', 'w_t', 'kappa_T', 'kappa_S')
  posterior_samples <- rstan::extract(fit_amr_negbin, pars = pars, permuted = TRUE)
  
  # convert the list of arrays to a data frame where each row is a posterior draw
  posterior_df <- as.data.frame(posterior_samples)
  
  # run forward simulations
  set.seed(42) # for reproducibility
  n_iter <- 1000
  random_integers <- sample(1:6000, size = n_iter, replace = FALSE) # draw random integers without replacement
  n_years <- 15
  
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
    
    # initial conditions
    samples_y <- rstan::extract(fit_amr_negbin, pars = "y", permuted = TRUE)
    # high-risk group
    # no intervention (N)
    U_N_H = median(samples_y$y[, 8, 1])
    E_N_H_0 = median(samples_y$y[, 8, 2]) - 5
    A_N_H_0 = median(samples_y$y[, 8, 3])
    S_N_H_0 = median(samples_y$y[, 8, 4])
    T_N_H_0 = median(samples_y$y[, 8, 5])
    E_N_H_c = median(samples_y$y[, 8, 6]) + 3
    A_N_H_c = median(samples_y$y[, 8, 7])
    S_N_H_c = median(samples_y$y[, 8, 8])
    T_N_H_c = median(samples_y$y[, 8, 9])
    E_N_H_t = median(samples_y$y[, 8, 10])
    A_N_H_t = median(samples_y$y[, 8, 11])
    S_N_H_t = median(samples_y$y[, 8, 12])
    T_N_H_t = median(samples_y$y[, 8, 13])
    E_N_H_d2 = median(samples_y$y[, 8, 14]) + 2
    A_N_H_d2 = median(samples_y$y[, 8, 15])
    S_N_H_d2 = median(samples_y$y[, 8, 16])
    T_N_H_d2 = median(samples_y$y[, 8, 17])
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
    U_N_L = median(samples_y$y[, 8, 18])  # 69
    E_N_L_0 = median(samples_y$y[, 8, 19]) - 2
    A_N_L_0 = median(samples_y$y[, 8, 20])
    S_N_L_0 = median(samples_y$y[, 8, 21])
    T_N_L_0 = median(samples_y$y[, 8, 22])
    E_N_L_c = median(samples_y$y[, 8, 23]) + 1
    A_N_L_c = median(samples_y$y[, 8, 24])
    S_N_L_c = median(samples_y$y[, 8, 25])
    T_N_L_c = median(samples_y$y[, 8, 26])
    E_N_L_t = median(samples_y$y[, 8, 27])
    A_N_L_t = median(samples_y$y[, 8, 28])
    S_N_L_t = median(samples_y$y[, 8, 29])
    T_N_L_t = median(samples_y$y[, 8, 30])
    E_N_L_d2 = median(samples_y$y[, 8, 31]) + 1
    A_N_L_d2 = median(samples_y$y[, 8, 32])
    S_N_L_d2 = median(samples_y$y[, 8, 33])
    T_N_L_d2 = median(samples_y$y[, 8, 34])
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
    out <- as.data.frame(out)
    
    # compute incidences and prescriptions
    incidence_all <- numeric(n_years)
    incidence_0 <- numeric(n_years)
    incidence_c <- numeric(n_years)
    incidence_t <- numeric(n_years)
    incidence_d2 <- numeric(n_years)
    
    for (t in 1:(n_years)) {
      isFixed = FALSE
      
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

############################################################## get results ##############################################################
# Define Parameter Grid
pd_vals <- seq(0, 1, by = 0.1)
pv_vals <- seq(0, 1, by = 0.1)
param_grid <- expand.grid(p_d = pd_vals, p_v = pv_vals)

compute_result <- function(res, strat_name, p_d, p_v, idx_2040) {
  total <- median(res$cases_all[, idx_2040])
  
  data.frame(
    Strategy = strat_name,
    p_d = p_d,
    p_v = p_v,
    Susceptible = (median(res$cases_0[, idx_2040]) / total) * 100,
    Cef_R       = (median(res$cases_c[, idx_2040]) / total) * 100,
    Tet_R       = (median(res$cases_t[, idx_2040]) / total) * 100,
    Dual_R      = (median(res$cases_d2[, idx_2040]) / total) * 100
  )
}

idx_2040 <- 15

all_results <- list()

# -----------------------------
# 1. No Intervention (run once)
# -----------------------------
print("No Intervention")
res <- run_amr(p_d = 0, p_v = 0)

no_int <- expand.grid(p_d = pd_vals, p_v = pv_vals) %>%
  rowwise() %>%
  do(compute_result(res, "No Intervention", .$p_d, .$p_v, idx_2040))

all_results[["No Intervention"]] <- no_int

# -----------------------------
# 2. Doxy-PEP (vary p_d only)
# -----------------------------
print("Doxy-PEP")
doxy_list <- lapply(pd_vals, function(p_d) {
  res <- run_amr(p_d = p_d, p_v = 0)
  
  expand.grid(p_v = pv_vals) %>%
    rowwise() %>%
    do(compute_result(res, "Doxy-PEP", p_d, .$p_v, idx_2040))
})

all_results[["Doxy-PEP"]] <- bind_rows(doxy_list)

# -----------------------------
# 3. 4CMenB (vary p_v only)
# -----------------------------
print("4CMenB")
vac_list <- lapply(pv_vals, function(p_v) {
  res <- run_amr(p_d = 0, p_v = p_v)
  
  expand.grid(p_d = pd_vals) %>%
    rowwise() %>%
    do(compute_result(res, "4CMenB", .$p_d, p_v, idx_2040))
})

all_results[["4CMenB"]] <- bind_rows(vac_list)

# -----------------------------
# 4. Combined (full grid)
# -----------------------------
print("Doxy-PEP + 4CMenB")
both_list <- mapply(function(p_d, p_v) {
  res <- run_amr(p_d = p_d, p_v = p_v)
  compute_result(res, "Doxy-PEP + 4CMenB", p_d, p_v, idx_2040)
}, param_grid$p_d, param_grid$p_v, SIMPLIFY = FALSE)

all_results[["Doxy-PEP + 4CMenB"]] <- bind_rows(both_list)

# -----------------------------
# Final dataset
# -----------------------------
final_df <- bind_rows(all_results)
# save.image("workspace_fixed_heatmap_uptakes_final.RData")

load("workspace_fixed_heatmap_uptakes_final.RData")
final_df <- final_df %>%
  mutate(Strategy = case_match(Strategy,
                               "4CMenB"            ~ "Vaccination",
                               "Doxy-PEP + 4CMenB" ~ "Doxy-PEP + Vaccination",
                               .default            = Strategy
  ))

# 1. Reshape and Re-label Strains in the exact requested order
plot_data <- final_df %>%
  pivot_longer(
    cols = c(Susceptible, Tet_R, Cef_R, Dual_R), 
    names_to = "Strain",
    values_to = "Percentage"
  ) %>%
  mutate(
    Strain = factor(Strain, 
                    levels = c("Susceptible", "Tet_R", "Cef_R", "Dual_R"),
                    labels = c("Susceptible", "Tet-R", "Cef-R", "Dual-R")),
    Strategy = factor(Strategy, levels = c("No Intervention", "Doxy-PEP", "Vaccination", "Doxy-PEP + Vaccination"))
  )

# 2. Generate coordinates and sequential letters (a) to (p) for all 16 facets
facet_labels <- expand.grid(
  Strategy = factor(c("No Intervention", "Doxy-PEP", "Vaccination", "Doxy-PEP + Vaccination"),
                    levels = c("No Intervention", "Doxy-PEP", "Vaccination", "Doxy-PEP + Vaccination")),
  Strain = factor(c("Susceptible", "Tet-R", "Cef-R", "Dual-R"), 
                  levels = c("Susceptible", "Tet-R", "Cef-R", "Dual-R"))
) %>%
  mutate(
    p_d = -0.05,     # Maps to the absolute left edge of your x-axis (p_d)
    p_v = 1.08,  # Floats slightly above the absolute top edge of your y-axis (p_v)
    Label = paste0("(", letters[row_number()], ")") 
  )

# 3. Setup Plotting Parameters
size <- 7

# 4. Generate the Figure
plot_16_panels_vax <- ggplot(plot_data, aes(x = p_d, y = p_v)) +
  geom_tile(aes(fill = Percentage, color = Percentage)) + 
  
  # Apply pseudo-log scale to the border colors so they match the tile interiors
  scale_color_viridis_c(option = "magma", direction = 1, trans = "pseudo_log", guide = "none") + 
  
  # Overlay the 16 letter panel labels (a) to (p) outside the boxes
  geom_text(data = facet_labels, aes(label = Label), 
            color = "black", fontface = "bold", size = (size / 7) * 2.8, 
            hjust = 0, vjust = 0) +                             
  
  facet_grid(Strain ~ Strategy) + 
  
  # Applied pseudo-log and explicit pretty breaks to expose hidden low-value data
  scale_fill_viridis_c(
    option = "magma", 
    direction = 1,
    name = "Strain Prevalence\n(%)", 
    trans = "pseudo_log",
    breaks = c(0, 1, 5, 10, 25, 50, 100),
    limits = c(0, 100),
    na.value = "grey90"
  ) +
  labs(
    x = "Doxy-PEP Uptake Rate",
    y = "Vaccination Uptake Rate"
  ) +
  # Prevents ggplot from clipping the outside (a)-(p) labels
  coord_cartesian(clip = "off") + 
  theme_minimal(base_size = size, base_family = "Helvetica") +
  theme(
    text = element_text(family = "Helvetica", size = size),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = size),
    panel.border = element_blank(),
    
    axis.ticks = element_line(color = "black", linewidth = 0.25),
    axis.ticks.length = unit(1, "mm"),
    axis.line.x = element_line(color = "black", linewidth = 0.25),
    axis.line.y = element_line(color = "black", linewidth = 0.25),
    
    # Breathing room so labels don't collide with upper text panels
    panel.spacing.y = unit(9, "pt"), 
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )

plot_16_panels_vax