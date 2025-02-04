# Code that accompanies the following manuscript:
#
# Title:   Predicting antibody kinetics and duration of protection against 
#          SARS-CoV-2 following vaccination from sparse serological data
# Authors: Julia Deichmann, Noam Barda, Michal Canetti, Yovel Peretz, Yael Ottolenghi,
#          Yaniv Lustig, Gili Regev-Yochay, Marc Lipsitch
#
# Date:    February 4, 2025
# Author:  Julia Deichmann <jdeichmann@hsph.harvard.edu>

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)


obs <- readRDS(file.path(here::here(), "data_synthetic", "synthetic_data.RDS"))
data_stan <- readRDS(file.path(here::here(), "data_synthetic", "data_stan.RDS"))

post_pop_pars <- readRDS(file.path(here::here(), "results", "posterior_pop-params.RDS"))
post_beta_pars <- readRDS(file.path(here::here(), "results", "posterior_beta-params.RDS"))
post_sigma_pars <- readRDS(file.path(here::here(), "results", "posterior_sigma-params.RDS"))

y_pp <- readRDS(file.path(here::here(), "results", "posterior-summary_y.RDS"))


### parameters ###

beta_summary <- post_beta_pars %>%
  mutate(parameter=factor(parameter, levels=unique(parameter))) %>%
  group_by(parameter) %>%
  summarize(q025=quantile(value, 0.025),
            median=quantile(value, 0.5),
            q975=quantile(value, 0.975))

beta_summary$name = NA
beta_summary[grepl("gB1", beta_summary$parameter),]$name = "gB1"
beta_summary[grepl("gB2", beta_summary$parameter),]$name = "gB2"
beta_summary[grepl("c2", beta_summary$parameter),]$name = "cs"
beta_summary[grepl("c1", beta_summary$parameter),]$name = "cl"


p_beta <- ggplot(beta_summary, aes(x = exp(median), y = rep(10:1, 4))) + 
  geom_vline(aes(xintercept = 1), linewidth = .25, color="grey80") +
  geom_errorbarh(aes(xmax = exp(q975), xmin = exp(q025)), linewidth = .25, 
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


post_pop <- ggplot(data=post_pop_pars, aes(value)) + 
  geom_histogram(color="black", fill="#680369FF", size=.25) +
  facet_wrap(~parameter, scales = "free_x", nrow=2) +   
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


sigma_summary <- post_sigma_pars %>%
  mutate(parameter=factor(parameter, levels=unique(parameter))) %>%
  group_by(parameter) %>%
  summarize(q025=quantile(value, 0.025),
            median=quantile(value, 0.5),
            q975=quantile(value, 0.975))

p_sigma <- ggplot(sigma_summary, aes(x = median, y = 6:1)) +
  geom_vline(aes(xintercept = 0), linewidth = .25, color="grey80") +
  geom_errorbarh(aes(xmax = q975, xmin = q025), linewidth = .25, 
                 height = 0) +
  geom_point(size = 1) +
  scale_y_continuous(breaks=6:1, labels=sigma_summary$parameter) +
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
  
pp_summary <- merge(y_pp[,c("ID", "time", "GMT_pred")], 
                    obs[,c("ID", "time", "obs")], by=c("ID", "time")) %>%
  mutate(res_log=log(obs)-log(GMT_pred))

p_res <- ggplot(pp_summary, aes(x=time, y=res_log)) +
  geom_point(size=.25) +
  geom_hline(yintercept=0, color="darkorange3", linewidth=.5, linetype="dashed") +
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
  geom_hline(yintercept=0.62, linetype="dashed", linewidth=.25, color="grey80") +
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

