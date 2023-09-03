---
title: "Distributional_Explorations"
subtitle: "Fitting mixed effects models"
date: last-modified
categories: [Analysis, R, Bayesian]
code-fold: true
execute: 
  warning: false
  eval: false
---

```{r}
pacman::p_load(tidyverse,tidybayes,brms,bayesplot,broom,broom.mixed,lme4,emmeans,here,knitr,kableExtra,gt,gghalves,patchwork,ggforce,ggdist,moments)
e1 <- readRDS(here("data/e1_08-21-23.rds"))
source(here("Functions/Display_Functions.R"))
source(here("Functions/org_functions.R"))
test <- e1 |> 
  filter(expMode2 == "Test") |> 
  group_by(vb) |> 
  mutate(distS = as.numeric(scale(dist, scale = FALSE)), 
         distS2 = custom_scale(dist))

options(brms.backend="cmdstanr",mc.cores=4)

test %>% group_by(condit) |> summarise(mean=mean(dist),sd=sd(dist),sk=moments::skewness(dist)) 
test %>% group_by(vb) |> summarise(mean=mean(dist),sd=sd(dist),sk=moments::skewness(dist)) 


```


```{r fig.width=12}

# test |> ggplot(aes(x=dist))+geom_histogram() + facet_wrap(~vb) + ggtitle("empirical_dist")  +
# test |> ggplot(aes(x=distS))+geom_histogram() + facet_wrap(~vb) + ggtitle("centered dist")  +
# test |> ggplot(aes(x=distS2))+geom_histogram() + facet_wrap(~vb) + ggtitle("scaled dist")

test |> filter(id %in% 1:5) |> ggplot(aes(x=dist))+geom_histogram() + facet_wrap(id~vb) + ggtitle("empirical_dist")  
plot(density(test$vx))
plot(density(test$dist))


test |> ggplot(aes(x=dist))+geom_density() + facet_wrap(~vb) + ggtitle("empirical_dist") 
test |> ggplot(aes(x=vx))+geom_density() + facet_wrap(~vb) + ggtitle("empirical_dist") 




```





```{r}

sk1 <- brm(dist ~ vb,family=skew_normal(),data=test,iter=800,chains=2)
g1 <- brm(dist ~ vb,family=gaussian(),data=test,iter=800,chains=2)
ga2 <- brm(dist+.01 ~ vb,family=Gamma(),data=test,iter=800,chains=2)
ln1 <- brm(dist+.01 ~ vb,family=lognormal(),data=test,iter=800,chains=2)
ln2 <- brm(dist+.0001 ~ vb,family=lognormal(link="inverse"),data=test,iter=800,chains=2)

g2 <- brm(bf(dist+.001|trunc(lb=0) ~ vb),data=testS,iter=800,chains=2,family=gaussian())
bayesplot::ppc_dens_overlay_grouped(testS$dist,yrep=posterior_predict(g2,ndraws=200),group=testS$vb)
pp_check(g2,type="stat_grouped",ndraws=200, group="vb",stat="mean")




bayesplot::ppc_dens_overlay_grouped(test$dist,yrep=posterior_predict(sk1,ndraws=200),group=test$vb)
bayesplot::ppc_dens_overlay_grouped(test$dist,yrep=posterior_predict(g1,ndraws=200),group=test$vb)
bayesplot::ppc_dens_overlay_grouped(test$dist,yrep=posterior_predict(ga2,ndraws=200),group=test$vb)
bayesplot::ppc_dens_overlay_grouped(test$dist,yrep=posterior_predict(ln1,ndraws=200),group=test$vb)



bayes_R2(g1)
bayes_R2(ga1)
bayes_R2(sk1)
bayes_R2(sk1)
bayes_R2(g1S)


#testS <- test %>% group_by(id,vb) |> filter(id %in% c("1","2","139")) |> select(id,vb,gt.stage,trial,condit,vx,dist,distS,distS2)
testS <- test %>% group_by(id,vb) |> filter(id %in% 1:15) |> select(id,vb,gt.stage,trial,condit,vx,dist,distS,distS2)

testS |> ggplot(aes(x=trial,y=dist,col=vb))+geom_line()+facet_wrap(~id)
testS |> ggplot(aes(x=trial,y=distS2,col=vb))+geom_line()+facet_wrap(~id)



g1S <- brm(distS ~ 0+ vb + condit,family=gaussian(),data=test,iter=1000,chains=4)
bayesplot::ppc_dens_overlay_grouped(test$distS,yrep=posterior_predict(g1S,ndraws=200),group=test$vb)

g1S2 <- brm(distS2 ~ vb,family=gaussian(),data=test,iter=1000,chains=4)
bayesplot::ppc_dens_overlay_grouped(test$distS2,yrep=posterior_predict(g1S2,ndraws=200),group=test$vb)


g1S_F <- brm(distS ~ 0+ vb + condit + (0+vb|id),family=gaussian(),data=test,iter=1000,chains=4,
             file=paste0(here::here("data/model_cache","e1_test_centeredDistS")))

bayesplot::ppc_dens_overlay_grouped(test$distS,yrep=posterior_predict(g1S_F,ndraws=200),group=test$vb)
bayesplot::ppc_dens_overlay_grouped(test$distS,yrep=posterior_predict(g1S_F,ndraws=200),group=test$vb)



g1S_gamma <- brm(dist+.001 ~ 1+ vb + condit + (0+vb|id),family=Gamma(),data=test,iter=1000,chains=4,
             file=paste0(here::here("data/model_cache","e1_test_Gamma")))
bayesplot::ppc_dens_overlay_grouped(test$dist,yrep=posterior_predict(g1S_gamma,ndraws=200),group=test$vb)

pp_check(g1S_weibull)+ xlim(c(-300,300))


g1S_exG <- brm(dist ~ 0 + vb + condit + (0+vb|id),family=exgaussian(),data=test,iter=1000,chains=4,
             file=paste0(here::here("data/model_cache","e1_test_exgauss")))
bayesplot::ppc_dens_overlay_grouped(test$dist,yrep=posterior_predict(g1S_exG,ndraws=200),group=test$vb)


```


```{r}



```


