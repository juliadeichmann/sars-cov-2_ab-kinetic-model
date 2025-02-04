library(cmdstanr)
library(bayesplot)
library(dplyr)
library(tidyr)
library(readr)


### create stan data input ###

data <- readRDS(file.path(here::here(), "data_synthetic", "synthetic_data.RDS"))

N = nrow(data)                                           # total number of observations
np = length(unique(data$ID))                             # number of individuals
nv = 2*np                                                # number of vaccination events

nd = (data %>% group_by(ID) %>% summarize(n=n()))$n      # number of observations per person
nb = rep(2, np)                                          # number of doses per person

tmp <- data %>% 
  group_by(ID) %>% 
  slice(1)

nc = 10                                                  # number of covariates
X = tmp[, 6:15]                                          # covariate matrix
X$age = log(X$age/43)
X$BMI = log(X$BMI/25)

t = data$time                                            # observation time points
obs = data$obs                                           # observed antibody levels

tvacc = (tmp[,c("V1","V2")] %>% pivot_longer(cols=everything()))$value     # vaccination time points

sigma_prior = 0.2                                        # sigma prior for random effects

data_stan <- list(N=N,
                  np=np,
                  nv=nv,
                  nd=nd,
                  nb=nb,
                  nc=nc,
                  X=X,
                  t=t,
                  obs=obs,
                  tvacc=tvacc,
                  sigma_prior=sigma_prior)

saveRDS(data_stan, file.path(here::here(), "data_synthetic", "data_stan.RDS"))


### sample ###

warmups <- 1000      
total_iterations <- 3000      
n_chains <- 4 

mod <- cmdstan_model(file.path(here::here(), "stan", "model.stan"))

fit <- mod$sample(
  data = data_stan,
  chains = n_chains,
  parallel_chains = n_chains,
  iter_warmup = warmups,
  iter_sampling = total_iterations-warmups,
  seed = 55,     
  refresh = 50,
  adapt_delta = 0.95,
  max_treedepth = 14
)

fit$save_object(file = file.path(here::here(), "results", "fit.RDS"))


### fit diagnostics ###

fit$diagnostic_summary()
post <- fit$draws(format="df")

pop_params <- c("log_A0_pop", "log_gB_pop[1]", "log_gB_pop[2]",
                "rho_pop", "log_r", "log_c_pop[1]", "log_c_pop[2]", "sigma")
sigma_params <- c("A0_sigma", "gB1_sigma", "gB2_sigma", "rho_sigma",
                  "c2_sigma", "c1_sigma")
beta_params <- names(post %>% select(matches("beta")))

fit$summary(variables=pop_params)
fit$summary(variables=beta_params)
fit$summary(variables=sigma_params)

mcmc_trace(post, pars=pop_params)
mcmc_trace(post, pars=beta_params)
mcmc_trace(post, pars=sigma_params)


### parameter posteriors ###

post_pop_pars <- post[pop_params] %>%
  tidyr::pivot_longer(cols=everything(), names_to="parameter")
post_beta_pars <- post[beta_params] %>%
  tidyr::pivot_longer(cols=everything(), names_to="parameter")
post_sigma_pars <- post[sigma_params] %>%
  tidyr::pivot_longer(cols=everything(), names_to="parameter")

saveRDS(post_pop_pars, file.path(here::here(), "results", "posterior_pop-params.RDS"))
saveRDS(post_beta_pars, file.path(here::here(), "results", "posterior_beta-params.RDS"))
saveRDS(post_sigma_pars, file.path(here::here(), "results", "posterior_sigma-params.RDS"))


### posterior prediction ###

pp_summary <- post %>%                     
  select(matches("ypred")) %>%
  pivot_longer(cols=everything()) %>%
  mutate(i=parse_number(name)) %>%
  group_by(i) %>%
  summarize(
    mean = mean(value),
    lower = quantile(value, .025),
    upper = quantile(value, .975),
    GMT_pred = exp(mean(log(value)))
  ) %>%
  mutate(ID=data$ID) %>%
  mutate(time=data_stan$t) 

saveRDS(pp_summary[,2:7], file.path(here::here(), "results", "posterior-summary_y.RDS"))
