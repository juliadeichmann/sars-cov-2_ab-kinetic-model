library(cmdstanr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)


obs <- readRDS(file.path(here::here(), "data_synthetic", "synthetic_data.RDS")) 
data_stan <- readRDS(file.path(here::here(), "data_synthetic", "data_stan.RDS"))

# fit <- readRDS(file.path(here::here(), "results", "fit.RDS"))

mod <- cmdstan_model(file.path(here::here(), "stan", "model_gq.stan"))


### run posterior predictive simulation ###

data_stan$tpred <- seq(0, 200, 5)
data_stan$nt <- length(data_stan$tpred)

# gq <- mod$generate_quantities(fit, data=data_stan, parallel_chains=4)
# 
# pred_1 <- gq$summary()
# pred_2 <- gq$summary(quantiles = ~ quantile(., probs = c(0.025, 0.975)),
#                      geom_mean = ~ exp(mean(log(.)))) %>%
#   rename(q025 = `2.5%`, q975 = `97.5%`)
# pred <- cbind(pred_1, pred_2[, 2:4]) %>%
#   as.data.frame()

# saveRDS(pred, file.path(here::here(), "results", "pp_simulation.RDS"))
pred <- readRDS(file.path(here::here(), "results", "pp_simulation.RDS"))


### plot individual Ab trajectories ###

ids <- rep(unique(obs$ID), each=data_stan$nt)
ypred <- pred %>% 
  filter(str_detect(variable,'pp_yhat')) %>%
  mutate(ID=ids) %>%
  mutate(time = rep(data_stan$tpred, data_stan$np)) %>%
  mutate(yhat = geom_mean) %>%    
  select(ID, time, yhat, q5, q95)


id <- c(5, 233, 488, 503, 745, 914)  
p <- ggplot() +
  geom_hline(yintercept=log(0.62), 
             linetype="dashed", linewidth=.25, color='grey80') +
  geom_ribbon(
    data = subset(ypred, ID %in% id),
    mapping = aes(time, log(yhat), ymin = log(q5), ymax = log(q95)),
    fill = "darkorange3", alpha=.3
  ) +
  geom_line(
    data = subset(ypred, ID %in% id),
    mapping = aes(time, log(yhat)),
    color = "darkorange3"
  ) +
  geom_point(
    data = subset(obs, ID %in% id), 
    mapping = aes(time, log(obs)),    
    size = 1, color="black"   
  ) +
  facet_wrap(~ factor(ID, levels = id), nrow = 6) + 
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    strip.text = element_blank(), 
    strip.background = element_blank(),
    axis.title=element_text(size=9, color="grey20")
  ) +
  scale_y_continuous(breaks = c(-2, 0, 2, 4)) +
  labs(
    x = "time (days)",
    y = "log IgG"
  ) 
p

ggsave(
  filename = file.path(here::here(), "figures", "Fig1b_fit_individuals.png"),
  p,
  width = 1.6, height = 4.7, units = "in", dpi = 450)


### analyze maximum Ab level ###

data_stan$X$age <- 43 * exp(data_stan$X$age)
data_stan$X$BMI <- 25 * exp(data_stan$X$BMI)

Amax <- pred %>%
  filter(str_detect(variable,'pp_ymax')) %>%
  cbind(data_stan$X[,]) %>%
  mutate(age1 = ifelse(age<30, 1, 0),
         age2 = ifelse(age>=30 & age<45, 1, 0),
         age3 = ifelse(age>=45 & age<60, 1, 0),
         age4 = ifelse(age>=60, 1, 0)) %>%
  mutate(bmi1 = ifelse(BMI<=25, 1, 0),
         bmi2 = ifelse(BMI>25 & BMI<30, 1, 0),
         bmi3 = ifelse(BMI>=30, 1, 0))

lm.Amax <- lm(log(geom_mean) ~ age2 + age3 + age4 + sex +  
                bmi2 + bmi3 + hypertension + dyslipidemia + 
                autoimmune_disease + diabetes + lung_disease +            
                heart_disease + immunosuppression,
              data = Amax)
summary(lm.Amax)

Amax_rm <- data.frame(mean=exp(lm.Amax$coefficients), exp(confint(lm.Amax)), 
                      type="maximum IgG")


### analyze duration of protection ###

tprot <- pred %>%
  filter(str_detect(variable,'t_prot')) %>%
  cbind(data_stan$X[,]) %>%
  mutate(age1 = ifelse(age<30, 1, 0),
         age2 = ifelse(age>=30 & age<45, 1, 0),
         age3 = ifelse(age>=45 & age<60, 1, 0),
         age4 = ifelse(age>=60, 1, 0)) %>%
  mutate(bmi1 = ifelse(BMI<=25, 1, 0),
         bmi2 = ifelse(BMI>25 & BMI<30, 1, 0),
         bmi3 = ifelse(BMI>=30, 1, 0))

tprot$unprotected <- 1
tprot[tprot$median>0,]$unprotected <- 0

# linear regression on duration of protection
lm.tprot <- lm(log(median) ~ age2 + age3 + age4 + sex +
                bmi2 + bmi3 + hypertension + dyslipidemia + 
                autoimmune_disease + diabetes + lung_disease +            
                heart_disease + immunosuppression,
              data = subset(tprot, median>0))
summary(lm.tprot)

t_prot_rm <- data.frame(mean=exp(lm.tprot$coefficients), exp(confint(lm.tprot)), 
                        type="duration of protection")

# logistic regression and odds ratio of being unprotected
glm.prot <- glm(unprotected ~ age2 + age3 + age4 + sex + bmi2 + bmi3 + 
                  hypertension + dyslipidemia + autoimmune_disease + diabetes + 
                  lung_disease + heart_disease + immunosuppression,
                data = tprot, family = "binomial")
odds <- data.frame(OR = coef(glm.prot), confint(glm.prot), type="unprotected")


### plot analysis ###

labels_tmp <- c("", "30 to <45", "45 to <60", ">=60", "male", "25 to 30", ">=30",
                "Hypertension", "dyslipidemia", "Autoimmune dis.", "diabetes", 
                "Lung disease", "Heart disease", "Immunosuppr.")

Amax_rm$label <- labels_tmp
Amax_rm$ref <- 1   
t_prot_rm$label <- labels_tmp
t_prot_rm$ref <- 1
odds$mean <- odds$OR
odds$label <- labels_tmp
odds$ref <- 0     

df_lm <- rbind(Amax_rm[2:14,], t_prot_rm[2:14,], odds[2:14, c("mean", "X2.5..", "X97.5..", "type", "label", "ref")])
df_lm$type <- factor(df_lm$type, 
                     levels=c("maximum IgG", 
                              "duration of protection",
                              "unprotected"),
                     labels=c("maximum IgG", 
                              "duration of protection",
                              "unprotected"))
df_lm$label <- factor(df_lm$label,
                      levels=labels_tmp[14:1],
                      labels=labels_tmp[14:1])

p <- ggplot(df_lm, aes(x = mean, y = label)) + 
  geom_vline(aes(xintercept = ref), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = X97.5.., xmin = X2.5..), linewidth = .375, 
                 height = 0, color = "#680369FF") +
  geom_point(size = 1.25, color = "#680369FF") +
  facet_wrap(~type, scales = "free_x") +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title=element_text(size=9, color="grey20"),
        strip.background=element_rect(fill="white")) + 
  ylab("") +
  xlab("ratio of means / log odds (95% CI)") +
  scale_y_discrete(labels=rev(c("30 to <45", "45 to <60", expression("">=60), 
                                "male", 
                                "25 to 30", expression("">=30),
                                "Hypertension", "dyslipidemia", "Autoimmune dis.", "diabetes", 
                                "Lung disease", "Heart disease", "Immunosuppr.")))
p

ggsave(
  filename = file.path(here::here(), "figures", "FigS3_lm_maxAb-tprot.png"),
  p,
  width = 5.5, height = 3.5, units = "in", dpi = 450)


##################### plot per group ###########################

ypred <- merge(ypred, cbind(ID=unique(obs$ID), data_stan$X), by="ID") %>%      
  mutate(
    age_group = dplyr::case_when(
      age < 30             ~ "<30",
      age >= 30 & age < 45 ~ "30 to <45",
      age >= 45 & age < 60 ~ "45 to <60",
      age >= 60 ~ ">=60"
    ),
    age_group = factor(
      age_group,
      level = c("<30", "30 to <45", "45 to <60", ">=60")
    )
  ) %>%
  mutate(
    BMI_group = dplyr::case_when(
      BMI <= 25 ~ "<=25",
      BMI > 25 & BMI < 30 ~ "25 to 30",
      BMI >= 30 ~ ">=30"
    ),
    BMI_group = factor(
      BMI_group,
      level = c("<=25", "25 to 30", ">=30")
    )
  )


covariates_cont <- c("age", "BMI")
p_cont <- list()
for (i in 1:2) {
  tmp_factor <- covariates_cont[i]
  tmp_y <- ypred %>%
    group_by(time, .data[[paste0(tmp_factor, "_group")]]) %>%
    summarize(mean=mean(log(yhat)), se=sd(log(yhat))/sqrt(n()))
  
  lvls <- levels(tmp_y[[paste0(tmp_factor, "_group")]])
  n_lvls <- length(lvls)
  
  p1 <- ggplot(data=tmp_y) +
    geom_ribbon(aes(x=time, y=mean, ymin=mean-se, ymax=mean+se, 
                    fill=.data[[paste0(tmp_factor, "_group")]]), alpha=.2) +
    geom_line(aes(x=time, y=mean, color=.data[[paste0(tmp_factor, "_group")]]), 
              linewidth=.75) +
    scale_color_manual(values = paletteer_d("MoMAColors::Flash")[seq(1,7, 6/(n_lvls-1))]) +
    scale_fill_manual(values = paletteer_d("MoMAColors::Flash")[seq(1,7, 6/(n_lvls-1))]) +
    geom_hline(yintercept=log(0.62), color='grey80', linetype="dashed") +
    xlab("time (days)") +
    ylab("log IgG") +
    ylim(-3.35, 4.7) +
    theme_bw() + theme(panel.grid = element_blank(),
                       legend.position = c(0.57, 0.2),
                       legend.title = element_blank(),
                       legend.text=element_text(size=8, color="grey20"),
                       legend.key.size = unit(0.5, "cm"),
                       legend.spacing.y = unit(-0.2, 'cm'),
                       legend.background = element_rect(fill='transparent'),
                       axis.title=element_text(size=9, color="grey20")) +   
    guides(fill=FALSE,
           color=guide_legend(ncol=2))
  
  p2 <- ggplot(data=Amax, aes(x=.data[[tmp_factor]], y=log(median))) +
    geom_point(aes(color=.data[[tmp_factor]]), size=.75) +
    scale_color_gradientn(colors=paletteer_d("MoMAColors::Flash")) +  
    geom_smooth(color='black', fill='grey60') +
    xlab("") + ylab("log IgG") +
    theme_bw() + theme(panel.grid = element_blank(), legend.position="",
                       axis.title=element_text(size=9, color="grey20"))

  p3 <- ggplot(data=tprot, aes(x=.data[[tmp_factor]], y=median)) +
    geom_point(aes(color=.data[[tmp_factor]]), size=.75) +  
    scale_color_gradientn(colors=paletteer_d("MoMAColors::Flash")) +   
    geom_smooth(data=subset(tprot, median>0), color='black', fill='grey60') +
    xlab("") + ylab("time (days)") +
    theme_bw() + theme(panel.grid = element_blank(), legend.position="",
                       axis.title=element_text(size=9, color="grey20"))
  
  p_cont[[i]] <- grid.arrange(p1, p2, p3, nrow=3,
                              top = textGrob(tmp_factor,
                                             gp = gpar(fontsize = 9, fontface = "bold"),
                                             x = 0.55, y = 0.5))
}

p_cont <- grid.arrange(grobs=p_cont, nrow=2)


covariates_disc <- c("sex", "hypertension", "dyslipidemia", "autoimmune_disease", 
                     "diabetes", "lung_disease", "heart_disease", "immunosuppression")
p_disc <- list()
for (i in 1:8) {
  tmp_factor <- covariates_disc[i]
  tmp_y <- ypred %>%
    group_by(time, .data[[tmp_factor]]) %>%
    summarize(mean=mean(log(yhat)), se=sd(log(yhat))/sqrt(n()))
  tmp_y[[tmp_factor]] <- as.factor(tmp_y[[tmp_factor]])
  
  if (i==1) {
    lvls <- c("female", "male")
  }
  else {
    lvls <- c("no", "yes")
  }
  
  p1 <- ggplot(data=tmp_y) +
    geom_ribbon(aes(x=time, y=mean, ymin=mean-se, ymax=mean+se, fill=.data[[tmp_factor]]), alpha=.2) +
    geom_line(aes(x=time, y=mean, group=.data[[tmp_factor]], color=.data[[tmp_factor]]), linewidth=.75) +
    scale_color_manual(values=c("grey75", "#680369FF")) +
    scale_fill_manual(values=c("grey75", "#680369FF")) +
    geom_hline(yintercept=log(0.62), color='grey80', linetype="dashed") +
    xlab("time (days)") + ylab("log IgG") +
    ylim(-3.35, 4.7) +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "",
                       axis.title=element_text(size=9, color="grey20"))
  
  p2 <- ggplot(data=Amax, aes(y=log(median), group=.data[[tmp_factor]], color=as.factor(.data[[tmp_factor]]))) +
    geom_boxplot() +
    scale_color_manual(values=c("grey75", "#680369FF")) +  
    scale_x_continuous(breaks=c(-0.2, 0.2),
                       labels=lvls,
                       name="") + 
    ylab("log IgG") +
    theme_bw() + theme(panel.grid = element_blank(), legend.position="",
                       axis.title=element_text(size=9, color="grey20")) 
  
  p3 <- ggplot(data=tprot, aes(y=median, group=.data[[tmp_factor]], color=as.factor(.data[[tmp_factor]]))) +
    geom_boxplot() +
    scale_color_manual(values=c("grey75", "#680369FF")) +  
    scale_x_continuous(breaks=c(-0.2, 0.2),
                       labels=lvls,
                       name="") +    # tmp_factor 
    ylab("time (days)") +
    theme_bw() + theme(panel.grid = element_blank(), legend.position="",
                       axis.title=element_text(size=9, color="grey20"),
                       axis.title.x = element_text(face = "bold"))
  
  p_disc[[i]] <- grid.arrange(p1, p2, p3, nrow=3, 
                              top = textGrob(tmp_factor,
                                             gp = gpar(fontsize = 9, fontface = "bold"),
                                             x = 0.55, y = 0.5))
}

p_disc <- grid.arrange(grobs=p_disc, nrow=2)
p_full <- grid.arrange(p_cont, p_disc, nrow=1, widths=c(1,4))

ggsave(
  filename = file.path(here::here(), "figures", "FigS2_evaluation_Ab-response.png"),
  p_full,
  width = 14, height = 12, units = "in", dpi = 450)
