library(tidyverse)
library(rstan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(patchwork)
rstan_options (auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE)) # use all available cores by default when sampling

# time series of MSM gonorrhoea data
tests = c(86472, 99558, 114282, 146991, 162071, 170936, 190941, 228103, 268605, 323012, 364015, 372371) 
cases_0 = as.integer(c(22042, 17297, 21209, 26748, 33634, 36933) - c(41.6, 45.1, 54.4, 63.6, 75.0, 68.5) / 100 * c(22042, 17297, 21209, 26748, 33634, 36933))
cases_c = c(0, 0, 0, 0, 0, 0)
cases_t = as.integer(c(41.6, 45.1, 54.4, 63.6, 75.0, 68.5) / 100 * c(22042, 17297, 21209, 26748, 33634, 36933))
cases_d2 = c(0, 0, 0, 0, 0, 0) 
cases_symptomatic = c(402, 573, 454, 366, 355, 361, 408, 367, 406, 437) 
samples = c(1090, 1073, 1226, 821, 839, 799, 971, 852, 1052, 971) 

# Initial population size of MSM
N_t0 <- 600000

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

# times
# burnt in 10 years
n_years <- length(tests)
n_years_cases <- length(cases_0)
n_years_symptomatic <- length(cases_symptomatic)
n_burntin <- 0
t <- seq(0, n_years+n_burntin+1, by = 1)
t_0 = 0 
t <- t[-1]

# initial conditions
U_N_H = 86400
E_N_H_0 = 126
A_N_H_0 = 1386
S_N_H_0 = 756
T_N_H_0 = 252
E_N_H_c = 1e-8
A_N_H_c = 1e-8
S_N_H_c = 1e-8
T_N_H_c = 1e-8
E_N_H_t = 54
A_N_H_t = 594
S_N_H_t = 324
T_N_H_t = 108
E_N_H_d2 = 1e-8
A_N_H_d2 = 1e-8
S_N_H_d2 = 1e-8
T_N_H_d2 = 1e-8

U_N_L = 504900
E_N_L_0 = 179
A_N_L_0 = 1964
S_N_L_0 = 1071
T_N_L_0 = 357
E_N_L_c = 1e-8
A_N_L_c = 1e-8
S_N_L_c = 1e-8
T_N_L_c = 1e-8
E_N_L_t = 77
A_N_L_t = 842
S_N_L_t = 459
T_N_L_t = 153
E_N_L_d2 = 1e-8
A_N_L_d2 = 1e-8
S_N_L_d2 = 1e-8
T_N_L_d2 = 1e-8

y0 = c(U_N_H=U_N_H, E_N_H_0=E_N_H_0, A_N_H_0=A_N_H_0, S_N_H_0=S_N_H_0, T_N_H_0=T_N_H_0, E_N_H_c=E_N_H_c, A_N_H_c=A_N_H_c, S_N_H_c=S_N_H_c, T_N_H_c=T_N_H_c, E_N_H_t=E_N_H_t, A_N_H_t=A_N_H_t, S_N_H_t=S_N_H_t, T_N_H_t=T_N_H_t, E_N_H_d2=E_N_H_d2, A_N_H_d2=A_N_H_d2, S_N_H_d2=S_N_H_d2, T_N_H_d2=T_N_H_d2, 
       U_N_L=U_N_L, E_N_L_0=E_N_L_0, A_N_L_0=A_N_L_0, S_N_L_0=S_N_L_0, T_N_L_0=T_N_L_0, E_N_L_c=E_N_L_c, A_N_L_c=A_N_L_c, S_N_L_c=S_N_L_c, T_N_L_c=T_N_L_c, E_N_L_t=E_N_L_t, A_N_L_t=A_N_L_t, S_N_L_t=S_N_L_t, T_N_L_t=T_N_L_t, E_N_L_d2=E_N_L_d2, A_N_L_d2=A_N_L_d2, S_N_L_d2=S_N_L_d2, T_N_L_d2=T_N_L_d2)

# data for Stan
data_amr <- list(n_years = n_years, n_years_cases = n_years_cases, n_years_symptomatic = n_years_symptomatic, n_burntin = n_burntin, y0 = y0, ts = t, t_0 = t_0, q_H = q_H, c_H = c_H, c_L = c_L, q_L = q_L, tests = tests, cases_0 = cases_0, cases_c = cases_c, cases_t = cases_t, cases_d2 = cases_d2, cases_symptomatic = cases_symptomatic, samples = samples, alpha = alpha, N_t0 = N_t0, gamma = gamma)

# run MCMC
model <- stan_model("amr_fixedinitialstate_burntin.stan")
fit_amr_negbin <- sampling(model,
                                data = data_amr,
                                iter = 100,
                                warmup = 50,
                                control = list(adapt_delta = 0.995, max_treedepth = 15),
                                chains = 4,
                                cores = 4,
                                seed = 43, # seed = 43 for UK
                                verbose=TRUE)

saveRDS(fit_amr_negbin, file = "fit_results_fixedinitialstate_UK.rds")
fit_amr_negbin <- readRDS("fit_results_fixedinitialstate_UK.rds")

# print the data
options(max.print = 1000000)
sink("output_fixedinitialstate_uk.txt")
print(fit_amr_negbin)
sink()

# get the parameter median and 95% CI
samples <- rstan::extract(fit_amr_negbin, pars = "kappa_T")$kappa_T
posterior_median <- median(samples)
ci_lower <- quantile(samples, 0.025)
ci_upper <- quantile(samples, 0.975)
cat(sprintf("Posterior median: %.5f\n95%% credible interval: [%.5f, %.5f]\n",
            posterior_median, ci_lower, ci_upper))

# print the mcmc results
pars=c('beta', 'phi_beta', 'epsilon', 'sigma', 'psi', 'mu', 'eta_H_init', 'omega', 'phi_eta', 'rho', 'nu', 'phi', 'f_c', 'f_t', 'f_d2', 'w_c', 'w_t', 'kappa_T', 'kappa_0', 'kappa_c', 'kappa_t', 'kappa_d2', 'kappa_S')
print(fit_amr_negbin, pars = pars)

# trace plots to assess mixing of a chain
traceplot(fit_amr_negbin, pars = pars)

# marginal posterior densities
stan_dens(fit_amr_negbin, pars = pars, separate_chains = TRUE)

########################################################## plot England results ##########################################################
years_pred <- c(2011:2019, 2022:2024)
years_obs  <- c(2015:2019, 2022)

# prepare plotting data for susceptible cases
smr_pred_0 <- as.data.frame(summary(
  fit_amr_negbin,
  pars = "pred_cases_0",
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
)$summary)
smr_pred_0$year <- years_pred

df_obs <- data.frame(
  year = years_obs,
  observation = cases_0
)

df_0 <- merge(
  smr_pred_0,
  df_obs,
  by = "year",
  all.x = TRUE
)

df_0 <- transform(
  df_0,
  group  = factor(year),
  ymin   = X2.5.,
  lower  = X25.,
  middle = X50.,
  upper  = X75.,
  ymax   = X97.5.
)

# prepare plotting data for ceftriaxone-resistant cases
smr_pred_c <- as.data.frame(summary(
  fit_amr_negbin,
  pars = "pred_cases_c",
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
)$summary)
smr_pred_c$year <- years_pred

df_obs <- data.frame(
  year = years_obs,
  observation = cases_c
)

df_c <- merge(
  smr_pred_c,
  df_obs,
  by = "year",
  all.x = TRUE
)

df_c <- transform(
  df_c,
  group  = factor(year),
  ymin   = X2.5.,
  lower  = X25.,
  middle = X50.,
  upper  = X75.,
  ymax   = X97.5.
)

# prepare plotting data for tetracycline-resistant cases
smr_pred_c <- as.data.frame(summary(
  fit_amr_negbin,
  pars = "pred_cases_t",
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
)$summary)
smr_pred_t$year <- years_pred

df_obs <- data.frame(
  year = years_obs,
  observation = cases_t
)

df_t <- merge(
  smr_pred_t,
  df_obs,
  by = "year",
  all.x = TRUE
)

df_t <- transform(
  df_t,
  group  = factor(year),
  ymin   = X2.5.,
  lower  = X25.,
  middle = X50.,
  upper  = X75.,
  ymax   = X97.5.
)

# prepare plotting data for dual-resistant cases
smr_pred_d2 <- as.data.frame(summary(
  fit_amr_negbin,
  pars = "pred_cases_d2",
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
)$summary)
smr_pred_d2$year <- years_pred

df_obs <- data.frame(
  year = years_obs,
  observation = cases_d2
)

df_d2 <- merge(
  smr_pred_d2,
  df_obs,
  by = "year",
  all.x = TRUE
)

df_d2 <- transform(
  df_d2,
  group  = factor(year),
  ymin   = X2.5.,
  lower  = X25.,
  middle = X50.,
  upper  = X75.,
  ymax   = X97.5.
)

# posterior predictive check on number of cases
smr_pred_0 <- cbind(as.data.frame(summary(
  fit_amr_negbin, pars = "pred_cases_0", probs = c(0.025, 0.25, 0.5, 0.75, 0.975))$summary), t[1:(n_years)], cases_0[1:(n_years_cases)])
colnames(smr_pred_0) <- make.names(colnames(smr_pred_0)) # to remove % in the col names

# posterior predictive check on number of cases
smr_pred_c <- cbind(as.data.frame(summary(
  fit_amr_negbin, pars = "pred_cases_c", probs = c(0.025, 0.25, 0.5, 0.75, 0.975))$summary), t[1:(n_years)], cases_c[1:(n_years_cases)])
colnames(smr_pred_c) <- make.names(colnames(smr_pred_c)) # to remove % in the col names

# posterior predictive check on number of cases
smr_pred_t <- cbind(as.data.frame(summary(
  fit_amr_negbin, pars = "pred_cases_t", probs = c(0.025, 0.25, 0.5, 0.75, 0.975))$summary), t[1:(n_years)], cases_t[1:(n_years_cases)])
colnames(smr_pred_t) <- make.names(colnames(smr_pred_t)) # to remove % in the col names

# posterior predictive check on number of cases
smr_pred_d2 <- cbind(as.data.frame(summary(
  fit_amr_negbin, pars = "pred_cases_d2", probs = c(0.025, 0.25, 0.5, 0.75, 0.975))$summary), t[1:(n_years)], cases_d2[1:(n_years_cases)])
colnames(smr_pred_d2) <- make.names(colnames(smr_pred_d2)) # to remove % in the col names

selected_years <- c(2011:2019, 2022:2024)
all_years <- 2011:2024
year_indices <- which(all_years %in% selected_years)

# box plot, output size (8, 6)
df_0 <- data.frame(
  group = factor(all_years[year_indices]),
  ymin = smr_pred_0$X2.5.[1:(n_years)],  # minimum
  lower = smr_pred_0$X25.[1:(n_years)],  # Q1
  middle = smr_pred_0$X50.[1:(n_years)], # median
  upper = smr_pred_0$X75.[1:(n_years)],  # Q3
  ymax = smr_pred_0$X97.5.[1:(n_years)], # maximum
  observation = cases_0[1:(n_years)],    # observed
  year = all_years[year_indices]          # year
)

size = 35

p1 <- ggplot(
  df_0,
  aes(
    x = group,
    ymin   = ymin,
    lower  = lower,
    middle = middle,
    upper  = upper,
    ymax   = ymax,
    color  = "Predicted"
  )
) +
  geom_boxplot(
    stat = "identity",
    fill = "salmon",
    width = 0.6
  ) +
  geom_point(
    data = subset(df_0, !is.na(observation)),
    aes(x = group, y = observation, color = "Observed"),
    shape = 18,
    size = 4
  ) +
  labs(
    title = "England",
    x = "Year",
    y = "Annual Number of Susceptible Diagnosed Cases"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Predicted" = "darkred",
      "Observed"  = "deepskyblue4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, NA)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    legend.position = "inside",
    legend.position.inside = c(0.20, 0.85),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
    legend.box.background = element_rect(color = "black"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(),
    axis.text.x  = element_text(),
    axis.text.y  = element_text(),
    
    axis.ticks.y = element_blank(),
    axis.line.x  = element_line(color = "black", linewidth = 0.5),
    axis.line.y  = element_line(color = "black", linewidth = 0.5)
  )

p2 <- ggplot(
  df_c,
  aes(
    x = group,
    ymin   = ymin,
    lower  = lower,
    middle = middle,
    upper  = upper,
    ymax   = ymax,
    color  = "Predicted"
  )
) +
  geom_boxplot(
    stat = "identity",
    fill = "salmon",
    width = 0.6
  ) +
  geom_point(
    data = subset(df_c, !is.na(observation)),
    aes(x = group, y = observation, color = "Observed"),
    shape = 18,
    size = 4
  ) +
  labs(
    title = "England",
    x = "Year",
    y = "Annual Number of Ceftriaxone-Resistant Diagnosed Cases"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Predicted" = "darkred",
      "Observed"  = "deepskyblue4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, NA)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    legend.position = "inside",
    legend.position.inside = c(0.20, 0.85),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
    legend.box.background = element_rect(color = "black"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(),
    axis.text.x  = element_text(),
    axis.text.y  = element_text(),
    
    axis.ticks.y = element_blank(),
    axis.line.x  = element_line(color = "black", linewidth = 0.5),
    axis.line.y  = element_line(color = "black", linewidth = 0.5)
  )

p3 <- ggplot(
  df_t,
  aes(
    x = group,
    ymin   = ymin,
    lower  = lower,
    middle = middle,
    upper  = upper,
    ymax   = ymax,
    color  = "Predicted"
  )
) +
  geom_boxplot(
    stat = "identity",
    fill = "salmon",
    width = 0.6
  ) +
  geom_point(
    data = subset(df_t, !is.na(observation)),
    aes(x = group, y = observation, color = "Observed"),
    shape = 18,
    size = 4
  ) +
  labs(
    title = "England",
    x = "Year",
    y = "Annual Number of Tetracycline-Resistant Diagnosed Cases"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Predicted" = "darkred",
      "Observed"  = "deepskyblue4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, NA)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    legend.position = "inside",
    legend.position.inside = c(0.20, 0.85),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
    legend.box.background = element_rect(color = "black"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(),
    axis.text.x  = element_text(),
    axis.text.y  = element_text(),
    
    axis.ticks.y = element_blank(),
    axis.line.x  = element_line(color = "black", linewidth = 0.5),
    axis.line.y  = element_line(color = "black", linewidth = 0.5)
  )

p4 <- ggplot(
  df_d2,
  aes(
    x = group,
    ymin   = ymin,
    lower  = lower,
    middle = middle,
    upper  = upper,
    ymax   = ymax,
    color  = "Predicted"
  )
) +
  geom_boxplot(
    stat = "identity",
    fill = "salmon",
    width = 0.6
  ) +
  geom_point(
    data = subset(df_d2, !is.na(observation)),
    aes(x = group, y = observation, color = "Observed"),
    shape = 18,
    size = 4
  ) +
  labs(
    title = "England",
    x = "Year",
    y = "Annual Number of Dual-Resistant Diagnosed Cases"
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Predicted" = "darkred",
      "Observed"  = "deepskyblue4"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, NA)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    legend.position = "inside",
    legend.position.inside = c(0.20, 0.85),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
    legend.box.background = element_rect(color = "black"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(),
    axis.text.x  = element_text(),
    axis.text.y  = element_text(),
    
    axis.ticks.y = element_blank(),
    axis.line.x  = element_line(color = "black", linewidth = 0.5),
    axis.line.y  = element_line(color = "black", linewidth = 0.5)
  )
