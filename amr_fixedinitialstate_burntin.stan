functions {
  real get_C(real E_N_0, real A_N_0, real S_N_0) {
    # compute strain-specific infectious population size of MSM in group j
    return(E_N_0 + A_N_0 + S_N_0);
  }
  real get_N(real U_N, real E_N_0, real A_N_0, real S_N_0, real T_N_0, real E_N_c, real A_N_c, real S_N_c, real T_N_c, real E_N_t, real A_N_t, real S_N_t, real T_N_t, real E_N_d2, real A_N_d2, real S_N_d2, real T_N_d2) {
    # compute population size of MSM in group j
    return(U_N + E_N_0 + A_N_0 + S_N_0 + T_N_0 + E_N_c + A_N_c + S_N_c + T_N_c + E_N_t + A_N_t + S_N_t + T_N_t + E_N_d2 + A_N_d2 + S_N_d2 + T_N_d2);
  }
  real get_pi(real c_target, real N_target, real c_remain, real N_remain) {
    # compute proportion of all partnerships in the population that involve a member of group j
    return((c_target * N_target) / (c_target * N_target + c_remain * N_remain));
  }
  real get_lambda(real t, real t_0, real c, real beta, real phi_beta, real epsilon, real C_target, real N_target, real pi_target, real C_remain, real N_remain, real pi_remain) {
    # compute force of infection in group j
    return(c * beta * (1 + phi_beta * (t - t_0)) * (epsilon * C_target / N_target + (1 - epsilon) * (pi_target * C_target / N_target + pi_remain * C_remain / N_remain)));
  }
  real get_eta(real t, real t_0, real eta_H_init, real phi_eta) {
    # compute rate of screening in the absence of symptoms in group H
    return(eta_H_init * (1 + phi_eta * (t - t_0)));
  }
  real[] syphilis_model(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
      # initial conditions
      real U_N_H = y[1];
      real E_N_H_0 = y[2];
      real A_N_H_0 = y[3];
      real S_N_H_0 = y[4];
      real T_N_H_0 = y[5];
      real E_N_H_c = y[6];
      real A_N_H_c = y[7];
      real S_N_H_c = y[8];
      real T_N_H_c = y[9];
      real E_N_H_t = y[10];
      real A_N_H_t = y[11];
      real S_N_H_t = y[12];
      real T_N_H_t = y[13];
      real E_N_H_d2 = y[14];
      real A_N_H_d2 = y[15];
      real S_N_H_d2 = y[16];
      real T_N_H_d2 = y[17];
      real U_N_L = y[18];
      real E_N_L_0 = y[19];
      real A_N_L_0 = y[20];
      real S_N_L_0 = y[21];
      real T_N_L_0 = y[22];
      real E_N_L_c = y[23];
      real A_N_L_c = y[24];
      real S_N_L_c = y[25];
      real T_N_L_c = y[26];
      real E_N_L_t = y[27];
      real A_N_L_t = y[28];
      real S_N_L_t = y[29];
      real T_N_L_t = y[30];
      real E_N_L_d2 = y[31];
      real A_N_L_d2 = y[32];
      real S_N_L_d2 = y[33];
      real T_N_L_d2 = y[34];
      
      # fixed real parameters
      real q_H = x_r[1]; # Proportion of the MSM population in group H
      real c_H = x_r[2]; # Annual rate of partner change in group H
      real c_L = x_r[3]; # Annual rate of partner change in group L
      real t_0 = x_r[4]; # Initial time
      real q_L = x_r[5]; # Proportion of the MSM population in group L
      
      # fixed integer parameters
      int alpha = x_i[1]; # Annual MSM population entrants (at age 15)
      int N_t0 = x_i[2]; # Initial population size of MSM
      int gamma = x_i[3]; # Years spent in the sexually-active population
      
      # model parameters
      real beta = theta[1]; # Probability of transmission per partnership
      real phi_beta = theta[2]; # Annual increase in transmission risk behaviour
      real epsilon = theta[3]; # Level of assortativity in sexual mixing
      real sigma = theta[4]; # Rate of leaving incubation period
      real psi = theta[5]; # Probability that incident infection is symptomatic
      real mu = theta[6]; # Rate of seeking treatment due to symptom
      real eta_H_init = theta[7]; # Initial rate of asymptomatic screening in group H
      real omega = theta[8]; # Ratio of screening rate in group L vs group H
      real phi_eta = theta[9]; # Annual increase in screening rate
      real rho = theta[10]; # Rate of recovery after treatment
      real nu = theta[11]; # Rate of naural recovery
      real phi = theta[12]; # Treatment failure rate for ceftriaxone-resistant infections treated with ceftriaxone
      real f_c = theta[13]; # Relative fitness of ceftriaxone-resistant bacteria, compared with susceptible
      real f_t = theta[14]; # Relative fitness of tetracycline-resistant bacteria, compared with susceptible
      real f_d2 = theta[15]; # Relative fitness of dual resistant bacteria, compared with susceptible
      real w_c = theta[16]; # Probability of emergence of ceftriaxone resistance upon treatment with ceftriaxone
      real w_t = theta[17]; # Probability of emergence of tetracycline resistance upon taking doxy-PEP for syphilis after gonorrhoea infection
      real kappa_T = theta[18]; # Shape parameter of GUMCAD dataset (annual number of gonorrhoea tests)
      real kappa_0 = theta[19]; # Shape parameter of GUMCAD dataset (annual number of gonorrhoea diagnoses, susceptible)
      real kappa_c = theta[20]; # Shape parameter of GUMCAD dataset (annual number of gonorrhoea diagnoses, ceftriaxone-resistant)
      real kappa_t = theta[21]; # Shape parameter of GUMCAD dataset (annual number of gonorrhoea diagnoses, tetracycline-resistant)
      real kappa_d2 = theta[22]; # Shape parameter of GUMCAD dataset (annual number of gonorrhoea diagnoses, dual-resistant)
      real kappa_S = theta[23]; # Shape parameter of GRASP dataset (annual number of symptomatic diagnoses)
      
      # time-dependent variables
      real C_H_0 = get_C(E_N_H_0, A_N_H_0, S_N_H_0); # Number of susceptible infectious individuals in group H
      real C_L_0 = get_C(E_N_L_0, A_N_L_0, S_N_L_0); # Number of susceptible infectious individuals in group L
      real C_H_c = get_C(E_N_H_c, A_N_H_c, S_N_H_c); # Number of ceftriaxone-resistant infectious individuals in group H
      real C_L_c = get_C(E_N_L_c, A_N_L_c, S_N_L_c); # Number of ceftriaxone-resistant infectious individuals in group L
      real C_H_t = get_C(E_N_H_t, A_N_H_t, S_N_H_t); # Number of tetracycline-resistant infectious individuals in group H
      real C_L_t = get_C(E_N_L_t, A_N_L_t, S_N_L_t); # Number of tetracycline-resistant infectious individuals in group L
      real C_H_d2 = get_C(E_N_H_d2, A_N_H_d2, S_N_H_d2); # Number of dual-resistant infectious individuals in group H
      real C_L_d2 = get_C(E_N_L_d2, A_N_L_d2, S_N_L_d2); # Number of dual-resistant infectious individuals in group L
      
      real N_H = get_N(U_N_H, I_N_H, P_N_H, S_N_H, E_N_H, L_N_H, T_N_H, R_N_H); # Total number of MSM population in group H
      real N_L = get_N(U_N_L, I_N_L, P_N_L, S_N_L, E_N_L, L_N_L, T_N_L, R_N_L); # Total number of MSM population in group L
      
      real pi_H = get_pi(c_H, N_H, c_L, N_L); # Proportion of all partnerships in the population that involve a member of group H
      real pi_L = get_pi(c_L, N_L, c_H, N_H); # Proportion of all partnerships in the population that involve a member of group L
      
      real eta_H = get_eta(t, t_0, eta_H_init, phi_eta); # Rate of screening in the absence of symptoms in group H
      real eta_L = omega * eta_H; # Rate of screening in the absence of symptoms in group L
      
      real lambda_H_0 = get_lambda(t, t_0, c_H, beta, phi_beta, epsilon, C_H_0, N_H, pi_H, C_L_0, N_L, pi_L); # Force of infection in group H for susceptible strain
      real lambda_L_0 = get_lambda(t, t_0, c_L, beta, phi_beta, epsilon, C_L_0, N_L, pi_L, C_H_0, N_H, pi_H); # Force of infection in group L for susceptible strain
      real lambda_H_c = get_lambda(t, t_0, c_H, beta, phi_beta, epsilon, C_H_c, N_H, pi_H, C_L_c, N_L, pi_L); # Force of infection in group H for ceftriaxone-resistant strain
      real lambda_L_c = get_lambda(t, t_0, c_L, beta, phi_beta, epsilon, C_L_c, N_L, pi_L, C_H_c, N_H, pi_H); # Force of infection in group L for ceftriaxone-resistant strain
      real lambda_H_t = get_lambda(t, t_0, c_H, beta, phi_beta, epsilon, C_H_t, N_H, pi_H, C_L_t, N_L, pi_L); # Force of infection in group H for tetracycline-resistant strain
      real lambda_L_t = get_lambda(t, t_0, c_L, beta, phi_beta, epsilon, C_L_t, N_L, pi_L, C_H_t, N_H, pi_H); # Force of infection in group L for tetracycline-resistant strain
      real lambda_H_d2 = get_lambda(t, t_0, c_H, beta, phi_beta, epsilon, C_H_d2, N_H, pi_H, C_L_d2, N_L, pi_L); # Force of infection in group H for dual-resistant strain
      real lambda_L_d2 = get_lambda(t, t_0, c_L, beta, phi_beta, epsilon, C_L_d2, N_L, pi_L, C_H_d2, N_H, pi_H); # Force of infection in group L for dual-resistant strain
      
      # ODEs
      # high-risk group
      # no intervention (N)
      real dU_N_H = q_H * alpha + rho * (1 - w_c) * (T_N_H_0 + T_N_H_t) + rho * (1 - phi) * (T_N_H_c + T_N_H_d2) + nu * (A_N_H_0 + A_N_H_c + A_N_H_t + A_N_H_d2) - (lambda_H_0 + f_c * lambda_H_c + f_t * lambda_H_t + f_d2 * lambda_H_d2 + 1/gamma) * U_N_H;
      real dE_N_H_0 = lambda_H_0 * U_N_H - (sigma + 1/gamma) * E_N_H_0;
      real dA_N_H_0 = sigma * (1 - psi) * E_N_H_0 - (nu + eta_H + 1/gamma) * A_N_H_0;
      real dS_N_H_0 = sigma * psi * E_N_H_0 - (mu + 1/gamma) * S_N_H_0;
      real dT_N_H_0 = eta_H * A_N_H_0 + mu * S_N_H_0 - (rho + 1/gamma) * T_N_H_0;
      real dE_N_H_c = f_c * lambda_H_c * U_N_H - (sigma + 1/gamma) * E_N_H_c;
      real dA_N_H_c = sigma * (1 - psi) * E_N_H_c - (nu + eta_H + 1/gamma) * A_N_H_c + phi * rho * T_N_H_c;
      real dS_N_H_c = sigma * psi * E_N_H_c - (mu + 1/gamma) * S_N_H_c;
      real dT_N_H_c = eta * A_N_H_c + mu * S_N_H_c - (rho + 1/gamma) * T_N_H_c + w_c * rho * T_N_H_0;
      real dE_N_H_t = f_c * lambda_H_t * U_N_H - (sigma + 1/gamma) * E_N_H_t;
      real dA_N_H_t = sigma * (1 - psi) * E_N_H_t - (nu + eta_H + 1/gamma) * A_N_H_t;
      real dS_N_H_t = sigma * psi * E_N_H_t - (mu + 1/gamma) * S_N_H_t;
      real dT_N_H_t = eta_H * A_N_H_t + mu * S_N_H_t - (rho + 1/gamma) * T_N_H_t;
      real dE_N_H_d2 = f_d2 * lambda_H_d2 * U_N_H - (sigma + 1/gamma) * E_N_H_d2;
      real dA_N_H_d2 = sigma * (1 - psi) * E_N_H_d2 - (nu + eta_H + 1/gamma) * A_N_H_d2 + phi * rho * T_N_H_d2;
      real dS_N_H_d2 = sigma * psi * E_N_H_d2 - (mu + 1/gamma) * S_N_H_d2;
      real dT_N_H_d2 = eta_H * A_N_H_d2 + mu * S_N_H_d2 - (rho + 1/gamma) * T_N_H_d2 + w_c * rho * T_N_H_t;

      # low-risk group
      # no intervention (N)
      real dU_N_L = q_L * alpha + rho * (1 - w_c) * (T_N_L_0 + T_N_L_t) + rho * (1 - phi) * (T_N_L_c + T_N_L_d2) + nu * (A_N_L_0 + A_N_L_c + A_N_L_t + A_N_L_d2) - (lambda_L_0 + f_c * lambda_L_c + f_t * lambda_L_t + f_d2 * lambda_L_d2 + 1/gamma) * U_N_L;
      real dE_N_L_0 = lambda_L_0 * U_N_L - (sigma + 1/gamma) * E_N_L_0;
      real dA_N_L_0 = sigma * (1 - psi) * E_N_L_0 - (nu + eta_L + 1/gamma) * A_N_L_0;
      real dS_N_L_0 = sigma * psi * E_N_L_0 - (mu + 1/gamma) * S_N_L_0;
      real dT_N_L_0 = eta_L * A_N_L_0 + mu * S_N_L_0 - (rho + 1/gamma) * T_N_L_0;
      real dE_N_L_c = f_c * lambda_L_c * U_N_L - (sigma + 1/gamma) * E_N_L_c;
      real dA_N_L_c = sigma * (1 - psi) * E_N_L_c - (nu + eta_L + 1/gamma) * A_N_L_c + phi * rho * T_N_L_c;
      real dS_N_L_c = sigma * psi * E_N_L_c - (mu + 1/gamma) * S_N_L_c;
      real dT_N_L_c = eta * A_N_L_c + mu * S_N_L_c - (rho + 1/gamma) * T_N_L_c + w_c * rho * T_N_L_0;
      real dE_N_L_t = f_c * lambda_L_t * U_N_L - (sigma + 1/gamma) * E_N_L_t;
      real dA_N_L_t = sigma * (1 - psi) * E_N_L_t - (nu + eta_L + 1/gamma) * A_N_L_t;
      real dS_N_L_t = sigma * psi * E_N_L_t - (mu + 1/gamma) * S_N_L_t;
      real dT_N_L_t = eta_L * A_N_L_t + mu * S_N_L_t - (rho + 1/gamma) * T_N_L_t;
      real dE_N_L_d2 = f_d2 * lambda_L_d2 * U_N_L - (sigma + 1/gamma) * E_N_L_d2;
      real dA_N_L_d2 = sigma * (1 - psi) * E_N_L_d2 - (nu + eta_L + 1/gamma) * A_N_L_d2 + phi * rho * T_N_L_d2;
      real dS_N_L_d2 = sigma * psi * E_N_L_d2 - (mu + 1/gamma) * S_N_L_d2;
      real dT_N_L_d2 = eta_L * A_N_L_d2 + mu * S_N_L_d2 - (rho + 1/gamma) * T_N_L_d2 + w_c * rho * T_N_L_t;

      return {dU_N_H, dE_N_H_0, dA_N_H_0, dS_N_H_0, dT_N_H_0, dE_N_H_c, dA_N_H_c, dS_N_H_c, dT_N_H_c, dE_N_H_t, dA_N_H_t, dS_N_H_t, dT_N_H_t, dE_N_H_d2, dA_N_H_d2, dS_N_H_d2, dT_N_H_d2, dU_N_L, dE_N_L_0, dA_N_L_0, dS_N_L_0, dT_N_L_0, dE_N_L_c, dA_N_L_c, dS_N_L_c, dT_N_L_c, dE_N_L_t, dA_N_L_t, dS_N_L_t, dT_N_L_t, dE_N_L_d2, dA_N_L_d2, dS_N_L_d2, dT_N_L_d2};
  }
}
data {
  int<lower=1> n_years; # Number of years modelling
  real<lower=0> y0[34]; # Initial conditions for all compartments
  real ts[n_years + 6]; # Sequences of time steps
  real t_0;
  real q_H;
  real c_H;
  real c_L;
  real q_L;
  int tests[n_years]; # Annual gonorrhoea tests
  int cases_0[n_years]; # Annual susceptible gonorrhoea cases
  int cases_c[n_years]; # Annual ceftriaxone-resistant gonorrhoea cases
  int cases_t[n_years]; # Annual tetracycline-resistant gonorrhoea cases
  int cases_d2[n_years]; # Annual dual-resistant gonorrhoea cases
  int cases_symptomatic[n_years]; # Annual symptomatic gonorrhoea cases
  int alpha;
  int N_t0;
  int gamma;
}
transformed data {
  real x_r[5];
  int x_i[3];
  
  # assign values to x_r
  x_r[1] = q_H;
  x_r[2] = c_H;
  x_r[3] = c_L;
  x_r[4] = t_0;
  x_r[5] = q_L;
  
  # assign values to x_i
  x_i[1] = alpha;
  x_i[2] = N_t0;
  x_i[3] = gamma;
}
parameters {
  real<lower=0, upper=1> beta;
  real<lower=0, upper=0.5> phi_beta;
  real<lower=0, upper=1> epsilon;
  real<lower=5, upper=200> sigma;
  real<lower=0, upper=1> psi;
  real<lower=5, upper=300> mu;
  real<lower=0, upper=2> eta_H_init;
  real<lower=0.1, upper=1> omega;
  real<lower=0, upper=0.3> phi_eta;
  real<lower=50, upper=500> rho;
  real<lower=0.5, upper=20> nu;
  real<lower=0, upper=0.2> phi;
  real<lower=0.7, upper=1.1> f_c;
  real<lower=0.7, upper=1.1> f_t;
  real<lower=0.5, upper=1.0> f_d2;
  real<lower=1e-10, upper=1e-5> w_c;
  real<lower=1e-6, upper=1e-2> w_t;
  real<lower=0.05, upper=1> kappa_T;
  real<lower=0.05, upper=1> kappa_0;
  real<lower=0.05, upper=1> kappa_c;
  real<lower=0.05, upper=1> kappa_t;
  real<lower=0.05, upper=1> kappa_d2;
  real<lower=0.05, upper=1> kappa_S;
}
transformed parameters{
  real y[n_years+6, 16];
  real incidence[n_years+5];
  real theta[13];
  theta[1] = beta;
  theta[2] = phi_beta;
  theta[3] = epsilon;
  theta[4] = rho;
  theta[5] = eta_H_init;
  theta[6] = phi_eta;
  theta[7] = omega;
  theta[8] = sigma;
  theta[9] = mu;
  theta[10] = psi_S;
  theta[11] = psi_E;
  theta[12] = psi_L;
  theta[13] = kappa_D;

  y = integrate_ode_bdf(syphilis_model, y0, t_0, ts, theta, x_r, x_i);
  for (t in 1:(n_years+5)) {
    # Trapezoidal rule: (f(a) + f(b)) / 2 * (b - a)
    incidence[t] = 0.5 * rho * (y[t, 8] + y[t + 1, 8] + y[t, 16] + y[t + 1, 16]);
  }
}
model {
  # priors
  beta ~ uniform(0,1);
  phi_beta ~ uniform(0,1);
  epsilon ~ uniform(0,1);
  rho ~ lognormal(log(20), 0.3); # 20.94 (11.1, 35.9)
  eta_H_init ~ uniform(0,4);
  phi_eta ~ uniform(0,1);
  omega ~ lognormal(-0.87,0.39); # 0.451 (0.22, 0.80)
  sigma ~ lognormal(log(15), 0.5); # 17 (6.6, 34.1)
  mu ~ lognormal(log(200), 0.6); # 239.5 (74.6, 535.9)
  psi_S ~ lognormal(log(10), 0.4); # 10.83 (5.18, 19.19)
  psi_E ~ lognormal(log(2), 0.4); # 2.17 (1.04, 3.86)
  psi_L ~ lognormal(log(0.5), 0.4); # 0.542 (0.26, 0.97)
  kappa_D ~ uniform(0,1);
  
  # sampling distribution
  for (i in 1:(n_years)) {
    cases[i] ~ neg_binomial_2(incidence[i + 5], kappa_D);
  }
}
generated quantities {
  real pred_cases[n_years];
  for (t in 1:n_years) {
    pred_cases[t] = neg_binomial_2_rng(incidence[t+5], kappa_D);
  }
}
