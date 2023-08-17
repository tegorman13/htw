---
title: "E1 Testing - Bayesian Mixed Models"
subtitle: "Fitting mixed effects models"
date: last-modified
categories: [Analysis, R, Bayesian]
code-fold: true
execute: 
  warning: false
  eval: false
---



```{r}
pacman::p_load(tidyverse,tidybayes,brms,bayesplot,bayestestR,
               broom.mixed,lme4,emmeans,here,knitr,kableExtra,gt)
walk(c(here("Functions/Display_Functions.R"), here("Functions/org_functions.R")), source)

test <- readRDS(here("data/e1_08-04-23.rds")) |> 
  filter(expMode2 == "Test") |> 
  mutate(distS2 = custom_scale(dist))

options(brms.backend="cmdstanr",mc.cores=4)


```



## Random Effects - Interaction condit x vb
```{r}

modelName <- "e1_testDistBT_RF"
e1_testDistBT_RF <- brm(dist ~ (condit * vb) + bandType + (1 + vb +bandType|id),
                      data=test,file=paste0(here::here("data/model_cache",modelName)),
                      iter=2000,chains=4,silent=0,
                      control=list(adapt_delta=0.92, max_treedepth=13))
```

## Random Effects - Interaction condit x vb
```{r}

modelName <- "e1_testDistRF2"
e1_testDistRF2 <- brm(dist ~ condit * vb + (1 + vb|id),
                      data=test,file=paste0(here::here("data/model_cache",modelName)),
                      iter=5000,chains=4)

GetModelStats(e1_testDistRF2)

brms_eq_tidy <-tidyMCMC(e1_testDistRF2, conf.int = TRUE, conf.level = 0.95,
           estimate.method = "median", conf.method = "HPDinterval")

bayes_R2(e1_testDistRF2)
sjPlot::tab_model(e1_testDistRF2)


fixef(e1_testDistRF2)
newdat <-data.frame(crossing(condit=c("Constant","Varied"), vb = unique(test$vb)))
preds <- fitted(e1_testDistRF2, re_formula = NA, newdata = newdat, probs = c(0.025, 0.975))


#shinystan::launch_shinystan(e1_testDistRF2)

```




## 0 intercept;  Random Effects - Interaction condit x vb
```{r}

modelName <- "e1_testDistRF2_0"
e1_testDistRF2_0 <- brm(dist ~ 0 + (condit * vb) + (0 + vb|id),
                        data=test,file=paste0(here::here("data/model_cache",modelName)), 
                        iter=5000,chains=4)
GetModelStats(e1_testDistRF2_0)

kable(fixef(e1_testDistRF2_0))
newdat <-data.frame(crossing(condit=c("Constant","Varied"), vb = unique(test$vb)))
preds <- fitted(e1_testDistRF2_0, re_formula = NA, newdata = newdat, probs = c(0.025, 0.975))


kable(tidy(e1_testDistRF2_0, effects="fixed",ess=TRUE))
conditional_effects(e1_testDistRF2_0,"condit:vb",method="pp_expect",points=TRUE)


draws_fit <- as_draws_df(e1_testDistRF2_0, variable = "^b_", regex = TRUE)


mcmc_plot(e1_testDistRF2_0, type="trace",variable="^b_",regex=TRUE)
mcmc_plot(e1_testDistRF2_0, type="intervals",variable="^b_",regex=TRUE)
mcmc_plot(e1_testDistRF2_0, type="areas",variable="^b_",regex=TRUE)

mcmc_hist(e1_testDistRF2_0,pars=c("b_conditConstant","b_conditVaried"))
plot(e1_testDistRF2_0,variable=c("b_conditConstant","b_conditVaried"))


bayesplot::ppc_dens_overlay_grouped(test$dist,yrep=posterior_predict(e1_testDistRF2_0,ndraws=200),group=test$vb)


mcmc_hist(e1_testDistRF2_0,prob=.5,regex_pars=c("^r_id\\[1,.*\\]"))


model_parameters(e1_testDistRF2_0,effects="random",keep="^r_id\\[3,.*\\]")
plot_subject_fits(e1_testDistRF2_0,3)


indvFit <- GetIndvFits(e1_testDistRF2_0)

indvFit |> ggplot(aes(x=condit,y=Median,fill=condit))+stat_halfeye()+facet_wrap(~vb)

```






## Random Effects
```{r}
#| results: asis
modelName <- "e1_testDistRF"
e1_testDistRF <- brm(dist ~ condit + vb + (1 + vb|id),
                     data=test,file=paste0(here::here("data/model_cache",modelName)))


pt <- posterior_table(e1_testDistRF) |> select(-CI)
kable(pt)
brms_posterior_checks(e1_testDistRF,dist,vb)

map_estimate(e1_testDistRF)


hdi(e1_testDistRF$fit, ci = c(0.5, 0.75, 0.89, 0.95))


plot(pt)

intPlot <- plot(conditional_effects(e1_testDistRF,effects="vb:condit"))


md <- tidy_draws(e1_testDistRF) |> select(b_Intercept:`b_vb1200M1400`)
plot(bayestestR::hdi(md, ci = c(.89, .95)))

plot(bayestestR::bayesfactor_parameters(e1_testDistRF, null = c(-.5, .5)))

p1 <- GetModelStats(e1_testDistRF)

kable(p1) |> column_spec(1:9,width="5em")

#GetBrmsModelStats(e1_testDistRF)



e1_testDistRF %>%
  spread_draws(b_Intercept, r_condition[condition,])


draws1 <- e1_testDistRF |> spread_draws()

```

## Interaction - Grouped Random Effects
```{r}

modelName <- "e1_testConditVb_Dist_Gr"
e1_testConditVb_Dist_Gr <- brm(dist ~ condit * vb + (1 + vb|gr(id,by=condit)),
                      data=test,file=paste0(here::here("data/model_cache",modelName)),
                      iter=5000,chains=4)
fixef(e1_testConditVb_Dist_Gr)
head(coef(e1_testConditVb_Dist_Gr)$id)
bayes_R2(e1_testConditVb_Dist_Gr)

GetModelStats(e1_testConditVb_Dist_Gr)

# coef(e1_testConditVb_Dist_Gr)$id %>% as_tibble(rownames="id") %>% select(id, starts_with("Est")) |> print(n=10)

individual_coefs <- coef(e1_testConditVb_Dist_Gr)$id %>%
    as_tibble(rownames = "id") %>%
    select(id, starts_with("Estim")) %>%
    pivot_longer(cols = -id, names_to = "variable", values_to = "value") %>%
    separate(variable, c("type", "effect"), sep = "\\.")

merged_data <- test |> group_by(id, condit) |> summarise(n=n()) %>%
               left_join(individual_coefs, by = "id") |>  filter(effect == "Intercept" | grepl("^vb", effect))


head(merged_data)



ggplot(merged_data, aes(x = effect, y = value, color = condit)) +
  geom_point(position = position_jitterdodge()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Individual Differences in Velocity Production Task (Intercept & vb Effects)",
    x = "Velocity Band Effect",
    y = "Estimated Value",
    color = "Condition"
  )



 epred_draws(e1_testConditVb_Dist_Gr, newdata = test,ndraws=50) |>
ggplot(aes(x = vb, y = .epred, color = condit)) +
    gghalves::geom_half_violin(position=position_dodge(),alpha=.7) +
    gghalves::geom_half_boxplot(position=position_dodge(),alpha=.5) +
    #geom_jitter(width = 0.2, height = 0) +
    labs(
        title = "Posterior Predictive Distribution",
        x = "Effect Category",
        y = "Expected Predicted Value",
        color = "Condition"
    )


(r_fit <- e1_testConditVb_Dist_Gr %>% 
  tidy() %>% filter(effect=="fixed") |> select(-effect,-component, -group) |> 
  mutate(term = janitor::make_clean_names(term)) |>
    mutate(across(where(is.numeric), \(x) round(x, 0))) |> kbl() |>
    column_spec(1:5, width = "5em") ) 
 
cat(r_fit$term)
paste(r_fit$term)
```




## 0 intercept;  Random Effects - Interaction condit x vb
```{r}

modelName <- "e1_testDistRF2_02"
e1_testDistRF2_02 <- brm(dist ~ 0 + condit:vb + (0 + vb|id),data=test,file=paste0(here::here("data/model_cache",modelName)))
brms_eq_tidy <-tidyMCMC(e1_testDistRF2_02, conf.int = TRUE, conf.level = 0.95,
           estimate.method = "median", conf.method = "HPDinterval")

GetModelStats(e1_testDistRF2_02)

bayes_R2(e1_testDistRF2_02)
sjPlot::tab_model(e1_testDistRF2_02)


fixef(e1_testDistRF2_02)
newdat <-data.frame(crossing(condit=c("Constant","Varied"), vb = unique(test$vb)))
preds <- fitted(e1_testDistRF2_02, re_formula = NA, newdata = newdat, probs = c(0.025, 0.975))


#shinystan::launch_shinystan(e1_testDistRF2_02)

```





```{r}

modelName <- "e1_test_Dist_Int"
e1_testDist_Int <- brm(dist ~ 0 + Intercept + bandInt + condit + (1|id), 
                       data=test,
                       iter=1000, chains=4, silent=0,
                       file=paste0(here::here("data/model_cache",modelName)))

e1_testDist_Int

pp_check(e1_testDist_Int,group="bandInt")
bayesplot::ppc_dens_overlay_grouped(test$dist,yrep=posterior_predict(e1_testDist_Int,ndraws=100),group=test$bandInt)
pp_check(e1_testDist_Int,type="stat_grouped",ndraws=500, group="bandInt",stat="mean")
coef(e1_testDist_Int)$id

```


```{r}

modelName <- "e1Test_conditBand_RS1"
e1Test_conditBand_RS1 <- brm(dist ~ 1 + bandInt + condit + (1+ bandInt|id), 
                       data=test,
                       iter=1000, chains=4, silent=0,
                       file=paste0(here::here("data/model_cache",modelName)))

e1Test_conditBand_RS1
fixef(e1Test_conditBand_RS1)

pp_check(e1Test_conditBand_RS1,group="bandInt")
bayesplot::ppc_dens_overlay_grouped(test$dist,yrep=posterior_predict(e1Test_conditBand_RS1,ndraws=100),group=test$bandInt)
pp_check(e1Test_conditBand_RS1,type="stat_grouped",ndraws=500, group="bandInt",stat="mean")
coef(e1Test_conditBand_RS1)$id


e1Test_conditBand_RS1 |> 
  tidy_draws() |> 
  select(starts_with("b_"),.chain,.iteration,.draw) 
  
  
e1Test_conditBand_RS1 |> 
  spread_draws(b_Intercept,b_conditVaried) 

```


```{r}

modelName <- "e1Test_conditVb_RS1"
e1Test_conditVb_RS1 <- brm(dist ~ 1 + vb + condit + (1+ vb|id), 
                       data=test,
                       iter=1000, chains=4, silent=0,
                       file=paste0(here::here("data/model_cache",modelName)))

e1Test_conditVb_RS1

pp_check(e1Test_conditVb_RS1,group="bandInt")
bayesplot::ppc_dens_overlay_grouped(test$dist,yrep=posterior_predict(e1Test_conditVb_RS1,ndraws=100),group=test$bandInt)
pp_check(e1Test_conditVb_RS1,type="stat_grouped",ndraws=500, group="vb",stat="mean")
coef(e1Test_conditVb_RS1)$id

```



```{r}




```





## Fixed Effects Only
```{r}

modelName <- "e1_testDist"
e1_testDist <- brm(dist ~ condit,data=test,file=paste0(here::here("data/model_cache",modelName)))
brms_eq_tidy <-tidyMCMC(e1_testDist, conf.int = TRUE, conf.level = 0.95,
           estimate.method = "median", conf.method = "HPDinterval")


modelName <- "e1_testDist0"
e1_testDist0 <- brm(dist ~ 0+condit,data=test,file=paste0(here::here("data/model_cache",modelName)))
brms_eq_tidy <-tidyMCMC(e1_testDist0, conf.int = TRUE, conf.level = 0.95,
           estimate.method = "median", conf.method = "HPDinterval")


modelName <- "e1_testDist_uneq"

e1_testDist_uneq <- brm(bf(dist ~ condit,sigma~condit),data=test,file=paste0(here::here("data/model_cache",modelName)))

brms_eq_tidy_uneq <-tidyMCMC(e1_testDist_uneq, conf.int = TRUE, conf.level = 0.95,
           estimate.method = "median", conf.method = "HPDinterval")

brms_eq_tidy_uneq |> mutate_at(vars(estimate, std.error, conf.low, conf.high),
            funs(ifelse(str_detect(term, "sigma"), exp(.), .)))

modelName <- "e1_testDist_uneq_robust"
brms_uneq_robust <- brm(
  bf(dist ~ condit, sigma ~ condit), data=test,
  family = student,file=paste0(here::here("data/model_cache",modelName)))

brms_uneq_robust_tidy <- tidyMCMC(brms_uneq_robust, conf.int = TRUE, conf.level = 0.95,estimate.method = "median", conf.method = "HPDinterval") %>% mutate_at(vars(estimate, std.error, conf.low, conf.high),
            funs(ifelse(str_detect(term, "sigma"), exp(.), .)))
brms_uneq_robust_tidy



#brm(dist ~ (0+vb)+(1+condit),data=test)


bayes_R2(e1_testDist)
sjPlot::tab_model(e1_testDist)

```





## Vx



```{r}
modelName <- "e1_testVxRF"
e1_testDistRF <- brm(vx ~ condit + bandInt+ (1 + vb|id),
                     data=test,file=paste0(here::here("data/model_cache",modelName)))

modelName <- "e1_testVxRF2"
e1_testVxRF2 <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                      data=test,file=paste0(here::here("data/model_cache",modelName)),
                      iter=5000,chains=4)


modelName <- "e1_test_vx_Int"
e1_testDist_Int <- brm(vx ~ 0 + Intercept + bandInt + condit + (1|id), 
                       data=test,
                       iter=2000, chains=4, silent=0,
                       file=paste0(here::here("data/model_cache",modelName)))


modelName <- "e1_testVxRF2_02"
e1_testDistRF2_02 <- brm(vx ~ 0 + condit:bandInt + (0 + bandInt|id),data=test,file=paste0(here::here("data/model_cache",modelName)))
brms_eq_tidy <-tidyMCMC(e1_testDistRF2_02, conf.int = TRUE, conf.level = 0.95,
           estimate.method = "median", conf.method = "HPDinterval")



modelName <- "e1_testConditBand_vx_Gr"
e1_testConditBand_vx_Gr <- brm(vx ~ condit * bandInt + (1 + bandInt|gr(id,by=condit)),
                      data=test,file=paste0(here::here("data/model_cache",modelName)),
                      iter=5000,chains=4)

modelName <- "e1_testConditBand0_vx_Gr"
e1_testConditBand_vx_Gr <- brm(vx ~ 0 + condit * bandInt + (0 + bandInt|gr(id,by=condit)),
                      data=test,file=paste0(here::here("data/model_cache",modelName)),
                      iter=5000,chains=4)

modelName <- "e1_testConditBand02_vx_Gr"
e1_testConditBand_vx_Gr <- brm(vx ~ 1 + condit * bandInt + (0 + bandInt|gr(id,by=condit)),
                      data=test,file=paste0(here::here("data/model_cache",modelName)),
                      iter=5000,chains=4,silent=0)


```



















```{r}

bform1 <- bf(mvbind(vx, dist) ~ condit * vb + (1|p|id)) + set_rescor(TRUE)
mv1 <- brm(bform1, data = test, chains = 2, cores = 2)

conditional_effects(mv1)
pp_check(mv1,resp="dist")
pp_check(mv1,resp="vx")


ppc_dens_overlay_grouped(test$vx,posterior_predict(mv1,ndraws=200),test$vb)
bayes_R2(mv1)


bf_dist <- bf(dist|trunc(lb=-70) ~ condit * vb + (0 + vb|id) ) + gaussian()
bf_vx <- bf(vx|trunc(lb=0) ~ condit * vb + (0 + vb|id)) + gaussian()

mv3 <- brm(bf_dist + bf_vx + set_rescor(FALSE),
            data = test, chains = 2, cores = 2,iter=1300, silent=0,
           file=here::here("data/model_cache/mv_trunc2"))
bayes_R2(mv3)
pp_check(mv3,resp="dist")

ppc_dens_overlay_grouped(test$vx,posterior_predict(mv3,ndraws=200),test$vb)



bform2 <- bf(mvbind(vx, vy) ~ condit * vb + (0+vb|id)) + set_rescor(FALSE)
mv3 <- brm(bform2, data = test, chains = 2, cores = 2,silent=0)
bayes_R2(mv3)
pp_check(mv3,resp="vx")
pp_check(mv3,resp="vy")



fitted(mv3) %>%
  as_tibble() %>%
  bind_cols(test %>% select(id,condit,vx,vy,vb,bandInt)) %>%
  ggplot(aes(x = vx, y = Estimate.vx)) +
  geom_abline(linetype = 2, color = "grey50", linewidth = .5) +  
  geom_point(size = 1.5, aes(color=), alpha = 3/4) +
  geom_linerange(aes(ymin = Q2.5.vx, ymax = Q97.5.vx),
                 size = 1/4, color = "firebrick4") +
  facet_wrap(~vb)



```



```{r}

conditRF <- brm(vx ~ vb + condit + (1+vb|condit),data=test,iter=1000,chains=2)

conditRFX <- brm(vx ~ vb * condit + (1+vb|condit),data=test,iter=1000,chains=2,silent=0,
                 file=here("data/model_cache/e1_conditRFX"))


conditRFX %>%
  spread_draws(b_Intercept, r_condit[condit,]) %>%
  mutate(condition_mean = b_Intercept + r_condit) %>%
  ggplot(aes(y = condit, x = condition_mean)) +
  stat_halfeye() + facet_wrap(~vb)

```


```{r}

vx_mm4 <- brm(vx ~ (1+vb|id), data=test,chains=2,silent=0)

mcmc_areas(vx_mm4,prob=.5,regex_pars=c("^r_id\\[1,.*\\]"),regex=T)


vx_mm5<- brm(vx ~ 0 + (0+vb|id), data=test,chains=2,silent=0,iter=1200)


plot_subject_fits <- function(model, subject_code) {
  pattern <- glue("^r_id\\[{subject_code},.*\\]")
  plot <- mcmc_areas(model, prob = .5, regex_pars = c(pattern)) +
            ggtitle(glue("fit for subject #{subject_code}:"))
  return(plot)
}

plot_subject_fits(vx_mm5,86)

```


individual_coefs <- coef(gt_vx)$id %>%
    as_tibble(rownames = "id") %>%
    select(id, starts_with("Estim")) %>%
    pivot_longer(cols = -id, names_to = "variable", values_to = "value") %>%
    separate(variable, c("type", "effect"), sep = "\\.")


dist_gauss_trunc0 <- brm(dist|trunc(lb=0) ~ 1+vb*condit + (0+vb|id),data=testExtrap,family=gaussian(), iter=2000, control=list(adapt_delta=0.92, max_treedepth=13), chains=3, file=paste0(here::here("data/model_cache","gauss_t_vbCondit_dist_extrap_ml6")), silent=0); 


gauss_trunc0 <- brm(vx|trunc(lb=0) ~ 1+vb*condit + (0+vb|id),data=testExtrap,family=gaussian(), iter=2000, control=list(adapt_delta=0.92, max_treedepth=13), chains=3, file=paste0(here::here("data/model_cache","gauss_t_vbCondit_vx_extrap_ml6")), silent=0); 


intercept condit_varied vb350m550 vb600m800 vb800m1000 vb1000m1200 vb1200m1400 condit_varied_vb350m550 condit_varied_vb600m800 condit_varied_vb800m1000 condit_varied_vb1000m1200 condit_varied_vb1200m1400

