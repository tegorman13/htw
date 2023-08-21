---
title: Individual Differenes
date: last-modified
categories: [Analysis, R]
page-layout: full
toc: false
execute: 
  warning: false
  eval: false
---



```{r setup, include=FALSE}
source(here::here("Functions", "packages.R"))
test <- readRDS(here("data/e1_08-04-23.rds")) |> filter(expMode2 == "Test") 
options(brms.backend="cmdstanr",mc.cores=4)
e1Sbjs <- test |> group_by(id,condit) |> summarise(n=n())

testAvg <- test %>% group_by(id, condit, vb, bandInt,bandType,tOrder) %>%
  summarise(nHits=sum(dist==0),vx=mean(vx),dist=mean(dist),sdist=mean(sdist),n=n(),Percent_Hit=nHits/n)

e1_vxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                        data=test,file=paste0(here::here("data/model_cache", "e1_testVxBand_RF_5k_Ml1")),
                        iter=5000,chains=4,silent=0,
                        control=list(adapt_delta=0.94, max_treedepth=13))
mt2 <-GetModelStats(e1_vxBMM ) |> kable(escape=F,booktabs=T)


new_data_grid=map_dfr(1, ~data.frame(unique(test[,c("id","condit","bandInt")]))) |> 
  dplyr::arrange(id,bandInt) |> 
  mutate(condit_dummy = ifelse(condit == "Varied", 1, 0)) 


nTop <- 5

indv_coefs <- coef(e1_vxBMM)$id |> 
  as_tibble(rownames="id") |> 
  select(id, starts_with("Est")) |>
  left_join(e1Sbjs, by=join_by(id) ) |> 
  group_by(condit) |> 
  mutate(rank = rank(desc(Estimate.bandInt)),
         intErrorRank=rank((Est.Error.Intercept)),
         bandErrorRank=rank((Est.Error.bandInt)),
         nCond = n()) |> arrange(intErrorRank) |> ungroup()

indv_coefs %>% group_by(condit) |> top_n(n=5, wt=Est.Error.Intercept)
indv_coefs |> slice_max(n=5, Est.Error.Intercept,by=condit)


top_indv_coefs <- indv_coefs %>%
  filter(rank <=nTop) |>
  select(id,condit,rank,Est.BandInt=Estimate.bandInt, Est.Intercept=Estimate.Intercept) 

small_slope <- indv_coefs %>%
  filter(rank >=nCond - nTop) |>
  select(id,condit,rank,Est.BandInt=Estimate.bandInt, Est.Intercept=Estimate.Intercept) 



well_fit <- indv_coefs %>%
  filter(intErrorRank<=nTop)

```


```{r}

indv_coefs |> ggplot(aes(y=id, x=Estimate.bandInt)) + geom_pointrange()

fixed_effects <- e1_vxBMM |> 
  spread_draws(`^b_.*`,regex=TRUE) |> arrange(.chain,.draw,.iteration)


random_effects <- e1_vxBMM |> 
  gather_draws(`^r_id.*$`, regex = TRUE, ndraws = 2000) |> 
  separate(.variable, into = c("effect", "id", "term"), sep = "\\[|,|\\]") |> 
  mutate(id = factor(id,levels=levels(test$id))) |> 
  pivot_wider(names_from = term, values_from = .value) |> arrange(id,.chain,.draw,.iteration)

cd <- left_join(random_effects, fixed_effects, by = join_by(".chain", ".iteration", ".draw")) |> 
  rename(bandInt_RF = bandInt) |>
  mutate(Slope=bandInt_RF+b_bandInt) |> group_by(id) 

cdMed <- cd |> group_by(id) |> median_qi(Slope) |> mutate(rankSlope=rank(Slope)) |> arrange(rankSlope) |> left_join(e1Sbjs, by=join_by(id))

cdMed %>% ggplot(aes(y=rankSlope, x=Slope,fill=condit,color=condit)) + 
  geom_pointrange(aes(xmin=.lower , xmax=.upper)) + 
  labs(x="Estimated Slope", y="Participant") # + facet_wrap(~condit) + 

cd |>  ggplot(aes(y=id, x=Slope)) + geom_pointrange()
```


```{r}




coef(e1_vxBMM)$id %>% as_tibble(rownames="id") %>% select(id, starts_with("Est")) |> print(n=15)


n_lines <- 20
f <- fitted(e1_vxBMM, newdata = new_data_grid, re_formula = NA, summary=F, ndraws=n_lines) |> 
  as_tibble() %>% 
  mutate(draw = 1:n()) %>% 
  pivot_longer(-draw)




combined_df <- left_join(random_effects, fixed_effects, by = join_by(".chain", ".iteration", ".draw")) |> 
  rename(bandInt_RF = bandInt) |>
  right_join(new_data_grid, by = join_by("id")) |> 
  mutate(
    fixed_effects = b_Intercept +
      (bandInt * b_bandInt) +
      condit_dummy * b_conditVaried + 
      (bandInt * condit_dummy) * `b_conditVaried:bandInt` ,  
    random_effects = Intercept + bandInt_RF,  # Assuming random effects for intercept and bandInt
    estimate = fixed_effects + random_effects
  )


combined_df |> 
  filter(id %in% top_indv_coefs$id) |> 
  group_by(id, bandInt) |>
  sample_n(50) |>
  ggplot(aes(x=bandInt,y=estimate)) + 
  geom_abline(aes(intercept=Intercept+b_Intercept, slope=bandInt_RF+b_bandInt), color="grey50") +
  geom_abline(data=indv_coefs |> filter(id %in% top_indv_coefs$id),aes(intercept=Estimate.Intercept,slope=Estimate.bandInt,color="red")) +
  stat_halfeye() +
  stat_halfeye(data=testAvg |> filter(id %in% top_indv_coefs$id), aes(x=bandInt,y=vx),color="blue") +
  facet_wrap(~id)


combined_df |> 
  filter(id %in% (indv_coefs |> slice_min(Est.Error.Intercept,n=5,by=condit) |> pull(id))) |> 
  group_by(id, bandInt) |>
  sample_n(50) |>
  ggplot(aes(x=bandInt,y=estimate)) + 
  geom_abline(aes(intercept=Intercept+b_Intercept, slope=bandInt_RF+b_bandInt), color="grey50") +
  geom_abline(data=indv_coefs |> slice_min(Est.Error.Intercept,n=5,by=condit),aes(intercept=Estimate.Intercept,slope=Estimate.bandInt,color="red")) +
  stat_halfeye() +
  stat_halfeye(data=testAvg |> filter(id %in% well_fit$id), aes(x=bandInt,y=vx),color="blue") +
  ggh4x::facet_nested_wrap(vars(condit,id))



combined_df |> 
  filter(id %in% (indv_coefs |> filter(rank<=2) |> pull(id)) )  |> 
  group_by(id, bandInt) |>
  sample_n(50) |>
  ggplot(aes(x=bandInt,y=estimate)) + 
  geom_abline(aes(intercept=Intercept+b_Intercept, slope=bandInt_RF+b_bandInt), color="grey50") +
  geom_abline(data=indv_coefs |> filter(id %in% (indv_coefs |> filter(rank<=2) |> pull(id))),aes(intercept=Estimate.Intercept,slope=Estimate.bandInt,color="red")) +
  stat_halfeye() +
  stat_halfeye(data=testAvg |> filter(id %in% (indv_coefs |> filter(rank<=2) |> pull(id))), aes(x=bandInt,y=vx),color="blue") +
  ggh4x::facet_nested_wrap(vars(condit,id))

```

```{r}


indv_coefs |> 
  filter(id %in% top_indv_coefs$id) |> ggplot(aes(x=bandInt,y=estimate)) + 
  geom_abline(aes(intercept=Intercept+b_Intercept, slope=bandInt_RF+b_bandInt), color="grey50") +
  stat_halfeye()



all_effects <- e1_vxBMM |> 
  gather_draws(b_Intercept, b_conditVaried, b_bandInt, `b_conditVaried:bandInt`, `^r_id.*$`, regex = TRUE, ndraws = 1)

head(all_effects)


 nd <- e1_vxBMM |>  
   spread_draws(b_Intercept, b_bandInt, `b_conditVaried`,`b_conditVaried:bandInt`,  r_id[id,term], ndraws=10)  
 
 

```



```{r}

df_pred <- 
  posterior_predict(e1_distBMM, ndraws = 500) |>
  array_branch(margin=1) |> 
   map_dfr( 
    function(yrep_iter) {
      test %>%
        mutate(dist_pred = yrep_iter)
    },
    .id = 'iter'
  ) |>
  mutate(iter = as.numeric(iter))

df_pred %>% filter(id %in% 1:162 ) |>
  filter(iter < 100) %>%
  ggplot(aes(dist_pred, group = iter)) +
  geom_line(alpha = .05, stat = 'density', color = 'blue') +
  geom_density(data = test |> filter(id %in% 1:162),
               aes(dist,col=vb),
               inherit.aes = FALSE,
               size = 0.8) + # 1
  facet_grid(condit ~ vb) +
  xlab('vx')

df_pred %>% filter(id %in% 1:10 ) |>
  filter(iter < 100) %>%
  ggplot(aes(dist_pred, group = iter)) +
  geom_line(alpha = .05, stat = 'density', color = 'blue') +
  geom_density(data = test |> filter(id %in% 1:10),
               aes(dist),
               inherit.aes = FALSE,
               size = 0.5) + # 1
  facet_wrap(id ~ .) +
  xlab('vx')
```

```{r}

vx_pred <- 
  posterior_predict(e1_vxBMM, ndraws = 500) |> 
  array_branch(margin=1) |> 
   map_dfr( 
    function(yrep_iter) {
      test %>%
        mutate(vx_pred = yrep_iter)
    },
    .id = 'iter'
  ) %>%
  mutate(iter = as.numeric(iter))

vx_pred %>% filter(id %in% 1:2 ) |>
  filter(iter < 100) %>%
  ggplot(aes(vx_pred, group = iter)) +
  geom_line(alpha = .05, stat = 'density', color = 'blue') +
  geom_density(data = test |> filter(id %in% 1:2),
               aes(vx,col=vb),
               inherit.aes = FALSE,
               size = 0.5) + # 1
  facet_grid(id ~ vb) +
  xlab('vx')

vx_pred  |>
  filter(iter < 100) %>%
  ggplot(aes(vx_pred, group = iter)) +
  geom_line(alpha = .05, stat = 'density', color = 'blue') +
  geom_density(data = test,
               aes(vx,col=vb),
               inherit.aes = FALSE,
               size = 0.5) + # 1
  facet_grid(condit ~ vb) +
  xlab('vx')

vx_pred %>% filter(id %in% 1:9 ) |>
  filter(iter < 100) %>%
  ggplot(aes(vx_pred, group = iter)) +
  geom_line(alpha = .05, stat = 'density', color = 'blue') +
  geom_density(data = test |> filter(id %in% 1:9),
               aes(vx),
               inherit.aes = FALSE,
               size = 0.5) + # 1
  facet_wrap(id ~ .,ncol=3) +
  xlab('vx')

```
```


```{r}


indv_coefs |> 
  filter(id %in% top_indv_coefs$id) |> ggplot(aes(x=bandInt,y=estimate)) + 
  geom_abline(aes(intercept=Intercept+b_Intercept, slope=bandInt_RF+b_bandInt), color="grey50") +
  stat_halfeye()



all_effects <- e1_vxBMM |> 
  gather_draws(b_Intercept, b_conditVaried, b_bandInt, `b_conditVaried:bandInt`, `^r_id.*$`, regex = TRUE, ndraws = 1)

head(all_effects)


 nd <- e1_vxBMM |>  
   spread_draws(b_Intercept, b_bandInt, `b_conditVaried`,`b_conditVaried:bandInt`,  r_id[id,term], ndraws=10)  
 
 

```



```{r}

df_posterior <- e1_vxBMM |> as_tibble()

ggplot(df_posterior) + 
  aes(x = b_Intercept, y = b_bandInt) + 
  stat_density_2d( geom = "polygon")


df_posterior |> sample_n(100) |> ggplot() + 
  aes(x = b_Intercept, y = b_bandInt) + 
  geom_raster(interpolate = T)

df_posterior |> sample_n(100) |> ggplot(aes(x = b_Intercept, y = b_bandInt)) + 
  stat_ellipse(geom = "polygon", level = 0.1, alpha = 1/2)

df_posterior |> 
  sample_n(1000) |>
  ggplot() + 
  aes(x = b_Intercept, y = b_bandInt) + 
  stat_dist_dotsinterval()

```
