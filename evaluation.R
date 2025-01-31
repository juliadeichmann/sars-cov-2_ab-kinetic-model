library(bayesplot)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(gridExtra)


obs <- readRDS("data_synthetic/synthetic_data.RDS")
data_stan <- readRDS("data_synthetic/data_stan.RDS")

fit <- readRDS("results/fit.RDS")  
post <- fit$draws(format="df")


### fit diagnostics ###

fit$diagnostic_summary()


pop_params <- c("log_A0_pop", "log_gB_pop[1]", "log_gB_pop[2]",
                "rho_pop", "log_r", "log_c_pop[1]", "log_c_pop[2]", "sigma")
sigma_params <- c("A0_sigma", "gB1_sigma", "gB2_sigma", "rho_sigma",
                  "c2_sigma", "c1_sigma")

beta_gB1_params <- names(post %>% select(matches("beta_gB1")))
beta_gB2_params <- names(post %>% select(matches("beta_gB2")))
beta_c1_params <- names(post %>% select(matches("beta_c1")))
beta_c2_params <- names(post %>% select(matches("beta_c2")))
beta_params <- c(beta_gB1_params, beta_gB2_params, beta_c1_params, beta_c2_params)

fit$summary(variables=pop_params)
fit$summary(variables=beta_params)
fit$summary(variables=sigma_params)

mcmc_trace(post, pars=pop_params)
mcmc_trace(post, pars=beta_params)
mcmc_trace(post, pars=sigma_params)


### parameters ###

tmp_gB1 <- fit$summary(variables = beta_gB1_params, ~quantile(.x, probs = c(0.025, 0.5, 0.975)))
tmp_gB1$name <- "gB1"
tmp_gB2 <- fit$summary(variables = beta_gB2_params, ~quantile(.x, probs = c(0.025, 0.5, 0.975)))
tmp_gB2$name <- "gB2"
tmp_c2 <- fit$summary(variables = beta_c2_params, ~quantile(.x, probs = c(0.025, 0.5, 0.975)))
tmp_c2$name <- "cs"
tmp_c1 <- fit$summary(variables = beta_c1_params, ~quantile(.x, probs = c(0.025, 0.5, 0.975)))
tmp_c1$name <- "cl"

beta_summary <- rbind(tmp_gB1, tmp_gB2, tmp_c2, tmp_c1)

p_beta <- ggplot(beta_summary, aes(x = exp(`50%`), y = rep(10:1, 4))) + 
  geom_vline(aes(xintercept = 1), size = .25, color="grey80") +
  geom_errorbarh(aes(xmax = exp(`97.5%`), xmin = exp(`2.5%`)), linewidth = .25, 
                 height = 0) +
  geom_point(size = 1) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y=element_blank(),
        axis.title=element_text(size=9, color="grey20"),
        text = element_text(size=9),
        strip.background=element_rect(fill="white")) +
  facet_wrap(~factor(name, c("gB1", "gB2", "cs", "cl")), ncol=4) +
  scale_y_continuous(breaks=10:1, labels=colnames(data_stan$X)) +
  xlim(exp(-1.5), exp(0.85)) +
  ylab("") +
  xlab("exp(beta)")
p_beta

ggsave(
  filename = file.path(here::here(), "figures", "Fig1c_posterior_beta-params.png"),
  p_beta,
  width = 5.5, height = 2.5, units = "in", dpi = 450)


fit_pars <- post[pop_params] %>%
  tidyr::pivot_longer(cols=everything())

post_pop <- ggplot(data=fit_pars, aes(value)) + 
  geom_histogram(color="black", fill="#680369FF", size=.25) +
  facet_wrap(~name, scales = "free_x", nrow=2) +   
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=9, color="grey20"),
        text = element_text(size=9),
        strip.background=element_rect(fill="white")) +
  xlab("parameter value")
post_pop

ggsave(
  filename = file.path(here::here(), "figures", "FigS1a_posterior_pop-params.png"),
  post_pop,
  width = 5.5, height = 3, units = "in", dpi = 450)


tmp_sigma <- fit$summary(variables = sigma_params, ~quantile(.x, probs = c(0.025, 0.5, 0.975)))
p_sigma <- ggplot(tmp_sigma, aes(x = `50%`, y = 6:1)) +
  geom_vline(aes(xintercept = 0), size = .25, color="grey80") +
  geom_errorbarh(aes(xmax = `97.5%`, xmin = `2.5%`), linewidth = .25, 
                 height = 0) +
  geom_point(size = 1) +
  scale_y_continuous(breaks=6:1, labels=tmp_sigma$variable) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size=9, color="grey20"),
        text = element_text(size=9),
        axis.title.y = element_blank(),
        legend.position="") +
  xlab("parameter value")
p_sigma

ggsave(
  filename = file.path(here::here(), "figures", "FigS1c_posterior_sigma-params.png"),
  p_sigma,
  width = 2, height = 3, units = "in", dpi = 450)


### residuals ###

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
  mutate(ID=obs$ID) %>%
  mutate(time=data_stan$t) 
  
pp_summary <- merge(pp_summary[,c("ID", "time", "GMT_pred")], 
                    obs[,c("ID", "time", "obs")], by=c("ID", "time")) %>%
  mutate(res_log=log(obs)-log(GMT_pred))

p_res <- ggplot(pp_summary, aes(x=time, y=res_log)) +
  geom_point(size=.25) +
  geom_hline(yintercept=0, color="darkorange3", lwd=.5, linetype="dashed") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title=element_text(size=9, color="grey20"),
        text = element_text(size=9),
        plot.margin=unit(c(0.2,-0.2,0.2,0.2), "cm")) +
  ylim(-1.6, 1.6) +
  xlab("time (days)") + ylab("residual")     

p_hist <- ggplot(data=pp_summary, aes(y=res_log)) +
  geom_histogram(fill="white", color="black",
                 binwidth=0.1, size=.25) +
  geom_hline(yintercept=0, color="darkorange3", lwd=.5, linetype="dashed") +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        axis.title=element_text(size=9, color="grey20"),
        text = element_text(size=9),
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_line(),
        plot.margin=unit(c(0.2,0.2,0.2,-0.2), "cm")) +
  xlab("count") + ylab("") +
  scale_x_continuous(breaks = c(0, 1500, 3000)) +
  ylim(-1.6, 1.6) 

p1 <- grid.arrange(p_res, p_hist, widths = c(3.1,0.9))
p1

ggsave(
  filename = file.path(here::here(), "figures", "FigS1b_residuals.png"),
  p1,
  width = 5.1, height = 3, units = "in", dpi = 450)


p2 <- ggplot(pp_summary, aes(x=log(GMT_pred), y=log(obs))) +
  geom_abline(intercept = 0, slope = 1, color="darkorange3", linewidth=.5) +
  geom_point(size=.25) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title=element_text(size=9, color="grey20"),
        text = element_text(size=9)) +
  xlim(-6, 7) + ylim(-6, 7) +
  xlab("estimated log IgG") +     
  ylab("observed log IgG")    
p2

ggsave(
  filename = file.path(here::here(), "figures", "FigS1b_obs-vs-pred.png"),
  p2,
  width = 3, height = 3, units = "in", dpi = 450)


### ab trajectories ###

pp_summary <- pp_summary %>%
  mutate(week = (time+3)%/%7)

month_start <- c(39, 64, 92, 120, 148, 176, 204)
pp_summary$month=0
for (i in 1:(length(month_start)-1)){
  pp_summary[pp_summary$time>=month_start[i] & 
             pp_summary$time<month_start[i+1],]$month = i
}

p_traj <- ggplot() +
  geom_hline(yintercept=0.62, linetype="dashed", lwd=.25, color="grey80") +
  geom_jitter(data=subset(pp_summary, week<=5), aes(x=week*7-1.5, y=log(obs)), color="slategray3",
              alpha=.5, size=.2, width=1) +
  geom_boxplot(data=subset(pp_summary, week<=5), aes(x=week*7-1.5, y=log(obs), group=factor(week)),
               lwd=.5, outlier.shape = NA, colour="black", width=1.4) +   
  geom_jitter(data=subset(pp_summary, week<=5), aes(x=week*7+1.5, y=log(GMT_pred)), color="darkorange3",
              alpha=.5, size=.2, width=1) +
  geom_boxplot(data=subset(pp_summary, week<=5), aes(x=week*7+1.5, y=log(GMT_pred), group=factor(week)),
               lwd=.5, outlier.shape = NA, colour="black", width=1.4) +  
  
  geom_jitter(data=subset(pp_summary, month>0), aes(x=month*28+21-3, y=log(obs)), color="slategray3",
              alpha=.5, size=.2, width=2) +
  geom_boxplot(data=subset(pp_summary, month>0), aes(x=month*28+21-3, y=log(obs), group=factor(month)),
               lwd=.5, outlier.shape = NA, colour="black", width=2.5) +
  geom_jitter(data=subset(pp_summary, month>0), aes(x=month*28+21+3, y=log(GMT_pred)), color="darkorange3",
              alpha=.5, size=.2, width=2) +
  geom_boxplot(data=subset(pp_summary, month>0), aes(x=month*28+21+3, y=log(GMT_pred), group=factor(month)),
               lwd=.5, outlier.shape = NA, colour="black", width=2.5) +
  
  scale_x_continuous(breaks=c(0, 7, 14, 21, 28, 35, 49, 77, 105, 133, 161, 189),
                     labels=c("w0","w1","w2","w3","w4","w5", "M1", "M2", "M3", "M4", "M5", "M6"),
                     name="") +
  scale_y_continuous(name="log IgG") +
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.title=element_text(size=9, color="grey20"))
p_traj

ggsave(
  filename = file.path(here::here(), "figures", "Fig1a_fit_ab-trajectories.png"),
  p_traj,
  width = 5.5, height = 3, units = "in", dpi = 450)

