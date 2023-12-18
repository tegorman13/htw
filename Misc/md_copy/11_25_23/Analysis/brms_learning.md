---
title: "Bayesian Mixed Effects Models - Training Data"
subtitle: "Fitting mixed effects models"
date: last-modified
categories: [Analysis, R, Bayesian]
code-fold: true
---






```{r}
pacman::p_load(tidyverse,tidybayes,brms,bayesplot,broom,broom.mixed,lme4,emmeans,here,knitr,kableExtra,gt,gghalves,patchwork,ggforce,ggdist)
e1 <- readRDS(here("data/e1_08-21-23.rds"))
source(here("Functions/Packages.R"))
train <- e1 |> filter(expMode2 == "Train")  
train$bandIntS <- scale(train$bandInt)
train$distS <- scale(train$dist)

options(brms.backend="cmdstanr",mc.cores=4)


```




```{r}

mod_train_expo3_dist <-bf(dist ~ betaMu + (alphaMu - betaMu) * exp(-exp(gammaMu) * gt.train),
                 betaMu ~ 1 + (1|id), 
                 alphaMu ~ 1 +   (1|id), 
                 gammaMu ~ 1 +  (1|id), 
                 nl = TRUE)

e1_train_expo3_dist <- brm(mod_train_expo3_dist,
                            chains=4, iter=2000, silent=0,
                            file=here::here("data/model_cache/e1_train_expo3_dist"))


mod_train_expo3_condit_dist <- bf(dist ~ betaMu + (alphaMu - betaMu) * exp(-exp(gammaMu) * gt.train),
                 betaMu ~ 1 + condit + (1|id), 
                 alphaMu ~ 1 + condit +   (1|id), 
                 gammaMu ~ 1 + condit + (1|id), 
                 nl = TRUE)

e1_train_expo3_condit_dist <- brm(mod_train_expo3_condit_dist,
                            chains=4, iter=2000, silent=0,
                            file=here::here("data/model_cache/e1_train_expo3_condit_dist"))



mod_train_expo3_conditBandit_dist <- bf(dist ~ betaMu + (alphaMu - betaMu) * exp(-exp(gammaMu) * gt.train),
                 betaMu ~ 0 + condit + (1|id) + (1|id:bandInt), 
                 alphaMu ~ 0 + condit + (1|id) + (1|id:bandInt), 
                 gammaMu ~ 0 + condit  +  (1|id) + (1|id:bandInt), 
                 nl = TRUE)



e1_train_expo3_conditBanditS_dist <- brm(mod_train_expo3_conditBandit_dist,
                            chains=4, iter=2000, silent=0,
                            data=train,
                            file=here::here("data/model_cache/e1_train_expo3_conditBanditSdist"))




mod_train_expo3_band_distS <- bf(distS ~ betaMu + (alphaMu - betaMu) * exp(-exp(gammaMu) * gt.train),
                 betaMu ~ 1 + bandIntS + (1|id), 
                 alphaMu ~ 1 + bandIntS +   (1|id), 
                 gammaMu ~ 1 + bandIntS + (1|id), 
                 nl = TRUE)

e1_train_expo3_band_distS <- brm(mod_train_expo3_band_distS,
                            chains=4, iter=2000, silent=0,
                            data=train,
                            file=here::here("data/model_cache/e1_train_expo3_band_distS"))
coef(e1_train_expo3_band_distS)$id |> as_tibble(rownames="id") |> select(starts_with("Esti")) |> print(n=15)




 mod_train_expo3_bandCondit_distS <- bf(distS ~ betaMu + (alphaMu - betaMu) * exp(-exp(gammaMu) * gt.train),
                 betaMu ~ 1 + bandIntS + condit + (1|id), 
                 alphaMu ~ 1 + bandIntS + condit +  (1|id), 
                 gammaMu ~ 1 + bandIntS + condit + (1|id), 
                 nl = TRUE)

e1_train_expo3_bandCondit_distS <- brm(mod_train_expo3_bandCondit_distS,
                            chains=4, iter=2000, silent=0,
                            data=train,
                            file=here::here("data/model_cache/e1_train_expo3_bandCondit_distS"))                           
coef(e1_train_expo3_bandCondit_distS)$id |> as_tibble(rownames="id") |> select(starts_with("Esti")) |> print(n=15)


```





b_mod3 <- bf(dist ~ betaMu + (alphaMu - betaMu) * exp(-exp(gammaMu) * gt.train),
                 betaMu ~ 1 + condit + bandInt + (1|id), 
                 alphaMu ~ 1 + condit + bandInt +   (1|id), 
                 gammaMu ~ 1 + condit + bandInt + (1|id), 
                 nl = TRUE)
