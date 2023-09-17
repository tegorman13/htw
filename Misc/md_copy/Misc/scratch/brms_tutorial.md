---
title: "Some Brms Testing based on Tutorials"
subtitle: "Fitting mixed effects models"
date: last-modified
categories: [Analysis, R, Bayesian]
code-fold: true
---


```{r}
pacman::p_load(tidyverse,tidybayes,brms,bayesplot,broom,broom.mixed,lme4,emmeans,here,knitr,kableExtra,gt,gghalves,patchwork,ggforce,ggdist)
e1 <- readRDS(here("data/e1_08-04-23.rds"))
source(here("Functions/Display_Functions.R"))
test <- e1 |> filter(expMode2 == "Test")  
train <- e1 |> filter(expMode2 == "Train")  


dolphin <- aida::data_MT
project_colors = cspplot::list_colors() |> pull(hex)

options(brms.backend="cmdstanr",mc.cores=4)





```


```{r}

dolphin <- dolphin %>% 
  filter(correct == 1) 

# plotting the data
ggplot(data = dolphin, 
       aes(x = MAD, 
           color = condition, fill = condition)) + 
  geom_density(alpha = 0.3, linewidth = 0.4, trim = F) +
  facet_grid(~group) + xlab("MAD")


model1_noInnteraction_F <- brm(MAD ~ group + condition,data=dolphin)

model2_interaction_FE <- brm(MAD ~ group*condition,data=dolphin)

model3_interaction_RandSlopes <- brm(MAD ~ group*condition + (1|subject_id),data=dolphin)

model4_interaction_MaxRE = brm(
  MAD ~ condition * group + (1 + group | exemplar) + (1 | subject_id), 
  data = dolphin
) 

loo_comp1 <- loo_compare(loo(model1_noInnteraction_F), loo(model2_interaction_FE))
loo_comp1
1 - pnorm(-loo_comp1[2,1], loo_comp1[2,2])

loo_comp2 <- loo_compare(loo(model3_interaction_RandSlopes), loo(model4_interaction_MaxRE))
loo_comp2
1 - pnorm(-loo_comp2[2,1], loo_comp2[2,2])


lc3 <- loo_compare(
  loo(model1_noInnteraction_F),
  loo(model2_interaction_FE),
  loo(model3_interaction_RandSlopes),
  loo(model4_interaction_MaxRE)
)

```









```{r}
data_forget <- tibble(
  y = c(.94, .77, .40, .26, .24, .16),
  t = c(  1,   3,   6,   9,  12,  18),
  N = 100,
  k = y * N
)

data_forget |>
  ggplot(aes(x = t, y = y)) +
  geom_line(color = project_colors[6]) +
  geom_point(color = project_colors[2], size = 2.5) +
  ylim(0,1) +
  xlab("time after memorization") +
  ylab("recall rate") 

fit_exponential <- brms::brm(
    formula = brms::bf(k | trials(N) ~ a * exp(-b * t), 
                       a + b ~ 1, 
                       nl=TRUE),
    data    = data_forget,
    prior   = prior(lognormal(0,0.5), nlpar = "a", lb = 0) + 
              prior(lognormal(0,0.5), nlpar = "b", lb = 0),
    family  = binomial(link = "identity"),
    control = list(adapt_delta = 0.99)
  )

fit_exponential

```



Linear learning model
```{r}

mmLearnLin <- brm(dist ~ 1 + (gt.train | id),data=train)

mmLearnLin2 <- brm(dist ~ 1 + gt.train+ (gt.train | id),data=train)
conditional_effects(mmLearnLin2)

mmLearnLin3 <- brm(dist ~ 1 + gt.train + bandInt + (1 + gt.train | id),data=train)



mmLearnLin3 <- brm(dist ~ 1 + condit+ gt.train + bandInt + (1 + gt.train | id),data=train)
conditional_effects(mmLearnLin3)


pp_check(mmLearnLin3)
pp_check(mmLearnLin3,type="stat_grouped",ndraws=500, group="bandInt",stat="mean")

bayesplot::ppc_dens_overlay_grouped(train$dist,yrep=posterior_predict(mmLearnLin3,ndraws=100),group=train$bandInt)
bayesplot::ppc_ribbon_grouped(train$dist,yrep=posterior_predict(mmLearnLin3,ndraws=100),x=train$gt.train,group=train$bandInt)

bayesplot::ppc_error_hist_grouped(train$dist,yrep=posterior_predict(mmLearnLin3,ndraws=100),x=train$gt.train,group=train$bandInt)
bayesplot::ppc_intervals_grouped(train$dist,yrep=posterior_predict(mmLearnLin3,ndraws=100),x=train$gt.train,group=train$bandInt)

bayesplot::ppc_intervals_grouped(train$dist,yrep=posterior_predict(mmLearnLin3,ndraws=100),x=train$gt.train,group=train$bandInt)

trainS <- train %>% filter(id %in% 1:10)



```


                
https://osf.io/rb2vh

s3_training<-brm(accuracy~0+Condition*Set+Block+Condition:Set:Block+(1+Set:Block|gr(subject,by=Condition))+(1|stimulus), 
                 data=P2_3_training_mixed, family = bernoulli(),
                 save_all_pars = TRUE, save_model = TRUE, prior=prior, inits=0,
                 control = list(adapt_delta =.99, max_treedepth = 15), 
                 sample_prior=TRUE,  cores=5, chains=5, iter=25000)
mep_correct<-marginal_effects(s2_training, effects="Block:Condition",
                              condition=list(Set=c(1,2)), plot=F)
```{r}

tl1 <- brm(dist ~ 0 + condit*bandInt + (1+bandInt|gr(id,by=condit)),data=train,
           save_model=TRUE,chains=2)



tl2 <- brm(dist ~ 0 + condit*bandInt + (1+bandInt|gr(id,by=condit)),data=train,
           save_model=TRUE,chains=2,iter=1000)

samples1 <- posterior_samples(tl2)
summary(tl2)


tl3 <- brm(dist ~ 0 + condit*bandInt+gt.train + (1+bandInt:gt.train|gr(id,by=condit)),data=train,
           save_model=TRUE,chains=2,iter=1000)


Set1_base_vs_R<-(samples5$`b_ConditionBase:Set1:Block`-samples5$`b_ConditionCon1:Set1:Block`)
hist(Set1_base_vs_R) 


conditional_effects(tl3)


tl4 <- brm(dist ~ 0 + condit+bandInt+gt.train + (1+gt.train|id),data=train,
           save_model=TRUE,chains=2,iter=1000)


samples2 <- posterior_samples(tl4)
s2 <- as_draws_df(tl4)

ps2 <- as_tibble(posterior_samples(tl4))

tl4 |> epred_draws(newdata= crossing(id = unique(train$id), gt.train = unique(train$gt.train),condit=unique(train$condit),bandInt=unique(train$bandInt)) |> filter((condit=="Constant" & bandInt==800) | condit=="Varied"),ndraws=5)

train_C_vs_V <- samples2$b_conditConstant - samples2$b_conditVaried
hist(train_C_vs_V)


  
conditional_effects(tl4,groups="condit")


```









set.seed(123)
data <- tibble(
  subject = rep(1:50, each = 5),
  condition = rep(c("base","conA","conB"), each = 50*5),
  stimulus = rep(sample(1:12, size = 50*5, replace = TRUE), 3),
  block1 = rep(rep(1:5, each = 10), 3),
  log_rt = rnorm(50*5*3, mean = 7, sd = 0.4)
)

# Filtering and preprocessing
data <- data %>%
  mutate(bl = round(block1/2 + 0.1), 
         condition = factor(condition, levels = c("base","conA","conB")),
         set = ifelse(stimulus %in% c(5, 8), "Set1", "Set2"))

# Using `group_by` and `summarise` to get means and confidence intervals
dts <- data %>%
  group_by(condition, bl, set) %>%
  summarise(mean_log_rt = mean(log_rt, na.rm = TRUE), 
            CI_low = mean_log_rt - qt(0.975, df = n() - 1)*sd(log_rt)/sqrt(n()), 
            CI_high = mean_log_rt + qt(0.975, df = n() - 1)*sd(log_rt)/sqrt(n())) 

# Plotting using ggplot
ggplot(dts, aes(x = bl, y = mean_log_rt, color = condition, group = condition)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  facet_wrap(~set) +
  scale_color_manual(values = c("black", "red", "blue")) +
  labs(x = "Trainingblock", y = "Mean Log-RT") +
  theme_minimal()



```{r}

texp1 <- fit_exponential <- brms::brm(
    formula = brms::bf(dist~ a * exp(-b * gt.train), 
                       a + b ~ 1, 
                       nl=TRUE),
    data    = train,
    family  = gaussian(),
     prior   = prior(normal(0,10), nlpar = "a", lb = 0) + 
              prior(normal(5,100), nlpar = "b", lb = 0),
  )

texp1
conditional_effects(texp1)


texp2 <- fit_exponential <- brms::brm(
    formula = brms::bf(dist~ a * exp(-b * gt.train)+b2*bandInt, 
                       a + b +b2 ~ 1, 
                       nl=TRUE),
    data    = train,
    family  = gaussian(),
     prior   = prior(normal(0,10), nlpar = "a", lb = 0) + 
              prior(normal(5,100), nlpar = "b", lb = 0)+
              prior(normal(0,5), nlpar = "b2"),
  )

texp2
conditional_effects(texp2)



texp3 <- fit_exponential <- brms::brm(
    formula = brms::bf(dist~ a * exp(-b * gt.train)+b2*bandInt, 
                       a + b +b2 ~ 1 + (1 | id), 
                       nl=TRUE),
    data    = train,
    family  = gaussian(),
    prior   = prior(normal(0,100), nlpar = "a", lb = 0) + 
              prior(normal(5,100), nlpar = "b", lb = 0)+
              prior(normal(0,5), nlpar = "b2"),
  )
texp3


```


