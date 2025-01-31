library(cmdstanr)
library(dplyr)
library(tidyr)


### create stan data input ###

data <- readRDS("data_synthetic/synthetic_data.RDS")

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

saveRDS(data_stan, "data_synthetic/data_stan.RDS")
# data <- readRDS("data_synthetic/data_stan.RDS")


### sample ###

warmups <- 1000      
total_iterations <- 3000      
n_chains <- 4 

file <- "model.stan"
mod <- cmdstan_model(file)

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

fit$save_object(file = "results/fit.RDS")   

