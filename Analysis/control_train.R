pacman::p_load(dplyr,purrr,tidyr,tibble,ggplot2,
  brms,tidybayes, rstanarm,emmeans,broom,bayestestR,
  stringr, here,conflicted, patchwork, knitr,kableExtra)

walk(c("brms","dplyr","bayestestR"), conflict_prefer_all, quiet = TRUE)
walk(c("Display_Functions","org_functions"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
options(scipen = 999)


e1 <- readRDS(here("data/e1_08-21-23.rds")) 
e1Sbjs <- e1 |> group_by(id,condit) |> summarise(n=n())
testE1 <- e1 |> filter(expMode2 == "Test")
nbins=5
trainE1 <-  e1 |> filter(expMode2=="Train") |> group_by(id,condit, vb) |> 
    mutate(Trial_Bin = cut( gt.train, breaks = seq(1, max(gt.train),length.out=nbins+1),include.lowest = TRUE, labels=FALSE)) 
trainE1_max <- trainE1 |> filter(Trial_Bin == nbins, bandInt==800)
trainE1_avg <- trainE1_max |> group_by(id,condit) |> summarise(avg = mean(dist))


nbins=5
trainE1 <-  e1 |> filter(expMode2=="Train") |> group_by(id,condit, vb) |> 
    mutate(Trial_Bin = cut( gt.bandStage, breaks = seq(1, max(gt.bandStage),length.out=nbins+1),include.lowest = TRUE, labels=FALSE)) 
trainE1_max <- trainE1 |> filter(Trial_Bin == nbins, bandInt==800)
trainE1_avg <- trainE1_max |> group_by(id,condit) |> summarise(train_end = mean(dist))


trainE1 |> select(id,condit,Trial_Bin,trial,vb,bandInt,dist,vx,gt.bandStage) |> 
  group_by(id,condit,vb,Trial_Bin) |> 
  summarise(mean_dist=mean(dist),mean_vx=mean(vx),n=n()) 


testE1 <- testE1 |> left_join(trainE1_avg, by=c("id","condit")) |>
  select(id,condit,bandType,bandInt,vb,vx,dist,train_end)

bmtd2 <- brm(dist ~ condit * train_end + (1|bandInt) + (1|id), 
    data=testE1,
    file=paste0(here::here("data/model_cache","e1_trainEnd_RF2")), 
    iter=1000,chains=2, control = list(adapt_delta = .92, max_treedepth = 11))
summary(bmtd2)
bayestestR::describe_posterior(bmtd2)


# summary(bmtd2)
#  Family: gaussian 
#   Links: mu = identity; sigma = identity 
# Formula: dist ~ condit * train_end + (1 | bandInt) + (1 | id) 
#    Data: testE1 (Number of observations: 9491) 
#   Draws: 2 chains, each with iter = 1000; warmup = 500; thin = 1;
#          total post-warmup draws = 1000

# Group-Level Effects: 
# ~bandInt (Number of levels: 6) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)    60.40     26.44    30.21   127.38 1.01      269      496

# ~id (Number of levels: 156) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)   140.16      8.69   124.84   159.52 1.00      125      292

# Population-Level Effects: 
#                        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                113.62     35.51    53.37   190.69 1.01      118      217
# conditVaried              -5.26     41.53   -86.65    81.12 1.04       62       92
# train_end                  0.82      0.21     0.41     1.20 1.01       90      235
# conditVaried:train_end     0.09      0.26    -0.38     0.56 1.02       80      263

# Family Specific Parameters: 
#       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma   243.38      1.61   240.38   246.59 1.00     1831      639

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).
# Warning message:
# There were 1 divergent transitions after warmup. Increasing adapt_delta above 0.92 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 

# r$> bayestestR::describe_posterior(bmtd2)
# Summary of Posterior Distribution 

# Parameter              | Median |           95% CI |     pd |            ROPE | % in ROPE |  Rhat |    ESS
# ----------------------------------------------------------------------------------------------------------
# (Intercept)            | 111.72 | [ 53.37, 190.69] | 99.70% | [-29.55, 29.55] |        0% | 1.012 | 116.00
# conditVaried           |  -7.29 | [-86.65,  81.12] | 56.00% | [-29.55, 29.55] |    55.47% | 1.038 |  64.00
# train_end              |   0.82 | [  0.41,   1.20] |   100% | [-29.55, 29.55] |      100% | 1.011 |  86.00
# conditVaried:train_end |   0.10 | [ -0.38,   0.56] | 62.30% | [-29.55, 29.55] |      100% | 1.017 |  77.00






bmtd3 <- brm(dist ~ condit * bandType * train_end + (1|bandInt) + (1|id), 
    data=testE1, 
    file=paste0(here::here("data/model_cache","e1_trainEnd_BT_RF2")),
    iter=1000,chains=2, control = list(adapt_delta = .92, max_treedepth = 11))
summary(bmtd3)
bayestestR::describe_posterior(bmtd3)


# summary(bmtd3)
#  Family: gaussian 
#   Links: mu = identity; sigma = identity 
# Formula: dist ~ condit * bandType * train_end + (1 | bandInt) + (1 | id) 
#    Data: testE1 (Number of observations: 9491) 
#   Draws: 2 chains, each with iter = 1000; warmup = 500; thin = 1;
#          total post-warmup draws = 1000

# Group-Level Effects: 
# ~bandInt (Number of levels: 6) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)    73.85     31.99    35.91   157.68 1.00      215      429

# ~id (Number of levels: 156) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)   140.36      8.62   125.43   159.07 1.03      111      142

# Population-Level Effects: 
#                                              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                                       43.41     49.99   -62.85   136.23 1.01      140      255
# conditVaried                                    67.01     49.92   -32.25   154.69 1.05       71      147
# bandTypeExtrapolation                           97.92     27.14    45.56   151.52 1.00      290      382
# train_end                                        1.02      0.28     0.45     1.51 1.01      115      221
# conditVaried:bandTypeExtrapolation             -76.93     27.48  -132.46   -22.04 1.00      273      487
# conditVaried:train_end                          -0.56      0.31    -1.11     0.11 1.03       87      205
# bandTypeExtrapolation:train_end                 -0.26      0.17    -0.58     0.09 1.01      205      310
# conditVaried:bandTypeExtrapolation:train_end     0.89      0.18     0.52     1.24 1.01      207      432

# Family Specific Parameters: 
#       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma   240.85      1.82   237.19   244.52 1.00     1029      512

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).
# Summary of Posterior Distribution 

# Parameter                                    | Median |            95% CI |     pd |            ROPE | % in ROPE |  Rhat |    ESS
# ---------------------------------------------------------------------------------------------------------------------------------
# (Intercept)                                  |  41.58 | [ -62.85, 136.23] | 82.20% | [-29.55, 29.55] |    33.16% | 1.014 | 142.00
# conditVaried                                 |  69.52 | [ -32.25, 154.69] | 90.30% | [-29.55, 29.55] |    18.84% | 1.046 |  70.00
# bandTypeExtrapolation                        |  97.03 | [  45.56, 151.52] |   100% | [-29.55, 29.55] |        0% | 1.000 | 280.00
# train_end                                    |   1.04 | [   0.45,   1.51] | 99.90% | [-29.55, 29.55] |      100% | 1.010 | 108.00
# conditVaried:bandTypeExtrapolation           | -76.11 | [-132.46, -22.04] | 99.90% | [-29.55, 29.55] |     2.00% | 1.000 | 268.00
# conditVaried:train_end                       |  -0.57 | [  -1.11,   0.11] | 95.10% | [-29.55, 29.55] |      100% | 1.029 |  84.00
# bandTypeExtrapolation:train_end              |  -0.27 | [  -0.58,   0.09] | 92.00% | [-29.55, 29.55] |      100% | 1.006 | 204.00
# conditVaried:bandTypeExtrapolation:train_end |   0.89 | [   0.52,   1.24] |   100% | [-29.55, 29.55] |      100% | 1.005 | 206.00



condEffects <- function(m,xvar){
  m |> ggplot(aes(x = {{xvar}}, y = .value, color = condit, fill = condit)) + 
  stat_dist_pointinterval() + 
  stat_halfeye(alpha=.1, height=.5) +
  theme(legend.title=element_blank(),axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) 
  
}

 bmtd3 |> emmeans( ~condit * bandType * train_end) |>
  gather_emmeans_draws() |>
 condEffects(bandType) + labs(y="Absolute Deviation From Band", x="Band Type")

plot(conditional_effects(bmtd3),points=FALSE)
ce_bmtd3 <- plot(conditional_effects(bmtd3),points=FALSE)
wrap_plots(ce_bmtd3)



bmt6d <- brm(dist ~ condit * bandInt * train_end + (1|id), 
    data=testE1, 
    file=paste0(here::here("data/model_cache","e1_trainEnd_B6_RF2")),
    iter=2000,chains=2, control = list(adapt_delta = .92, max_treedepth = 11))

summary(bmt6d)
bayestestR::describe_posterior(bmt6d)

# summary(bmt6d)
#     bayestestR::describe_posterior(bmt6d)
#  Family: gaussian 
#   Links: mu = identity; sigma = identity 
# Formula: dist ~ condit * bandInt * train_end + (1 | id) 
#    Data: testE1 (Number of observations: 9491) 
#   Draws: 2 chains, each with iter = 2000; warmup = 1000; thin = 1;
#          total post-warmup draws = 2000

# Group-Level Effects: 
# ~id (Number of levels: 156) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)   139.82      8.68   123.12   156.77 1.02      362      557

# Population-Level Effects: 
#                                Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                        109.11     33.97    40.25   173.24 1.06       34      120
# conditVaried                      -8.52     52.43  -113.21    88.68 1.11       18       61
# bandInt                            0.04      0.02     0.00     0.07 1.01      502     1132
# train_end                          0.91      0.25     0.43     1.42 1.04       45      163
# conditVaried:bandInt              -0.02      0.03    -0.07     0.03 1.01      258      819
# conditVaried:train_end             0.51      0.31    -0.11     1.11 1.09       28       82
# bandInt:train_end                 -0.00      0.00    -0.00    -0.00 1.01      549     1338
# conditVaried:bandInt:train_end    -0.00      0.00    -0.00    -0.00 1.01      343     1213

# Family Specific Parameters: 
#       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma   243.66      1.73   240.38   247.07 1.00     2046     1378

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).
# Warning message:
# Parts of the model have not converged (some Rhats are > 1.05). Be careful when analysing the results! We recommend running more iterations and/or setting stronger priors. 
# Summary of Posterior Distribution 

# Parameter                      |    Median |            95% CI |     pd |            ROPE | % in ROPE |  Rhat |    ESS
# ----------------------------------------------------------------------------------------------------------------------
# (Intercept)                    |    110.96 | [  40.25, 173.24] |   100% | [-29.55, 29.55] |        0% | 1.063 |  32.00
# conditVaried                   |     -9.93 | [-113.21,  88.68] | 56.35% | [-29.55, 29.55] |    38.58% | 1.112 |  18.00
# bandInt                        |      0.04 | [   0.00,   0.07] | 98.05% | [-29.55, 29.55] |      100% | 1.007 | 417.00
# train_end                      |      0.90 | [   0.43,   1.42] |   100% | [-29.55, 29.55] |      100% | 1.041 |  42.00
# conditVaried:bandInt           |     -0.02 | [  -0.07,   0.03] | 77.45% | [-29.55, 29.55] |      100% | 1.013 | 210.00
# conditVaried:train_end         |      0.53 | [  -0.11,   1.11] | 94.50% | [-29.55, 29.55] |      100% | 1.086 |  23.00
# bandInt:train_end              | -3.19e-04 | [  -0.00,  -0.00] | 98.10% | [-29.55, 29.55] |      100% | 1.005 | 471.00
# conditVaried:bandInt:train_end | -6.13e-04 | [  -0.00,  -0.00] |   100% | [-29.55, 29.55] |      100% | 1.009 | 292.00


 bmt6d |> emmeans( ~condit * bandInt*train_end, at=list(bandInt=c(100,350,600,800,1000,1200),train_end=c(249,144,70, 30) )) |>
  gather_emmeans_draws() |>
   ggplot(aes(x=bandInt,y=.value,color=condit,fill=condit)) +
  stat_bar + facet_wrap(~train_end) + labs(y="Absolute Deviation From Band", x="Band Intensity")
plot(conditional_effects(bmt6d),points=FALSE)






bmt6d2 <- brm(dist ~ condit * bandInt * train_end + (1 + bandInt|id), 
    data=testE1, 
    file=paste0(here::here("data/model_cache","e1_trainEnd_B6_RFint")),
    iter=2000,chains=4, control = list(adapt_delta = .92, max_treedepth = 12))

plot(conditional_effects(bmt6d2),points=FALSE)

# summary(bmt6d2)
#  Family: gaussian 
#   Links: mu = identity; sigma = identity 
# Formula: dist ~ condit * bandInt * train_end + (1 + bandInt | id) 
#    Data: testE1 (Number of observations: 9491) 
#   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#          total post-warmup draws = 4000

# Group-Level Effects: 
# ~id (Number of levels: 156) 
#                        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)            282.43     16.13   251.45   315.45 1.01      519      779
# sd(bandInt)                0.30      0.02     0.27     0.34 1.01      599     1199
# cor(Intercept,bandInt)    -0.95      0.01    -0.97    -0.93 1.00     1019     1838

# Population-Level Effects: 
#                                Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                        104.91     60.76   -19.48   221.88 1.02      107      186
# conditVaried                      -6.65     92.55  -183.06   199.80 1.05       72       47
# bandInt                            0.04      0.07    -0.08     0.17 1.02      125      242
# train_end                          0.93      0.46     0.08     1.89 1.02      161      228
# conditVaried:bandInt              -0.02      0.10    -0.23     0.18 1.04       78       64
# conditVaried:train_end             0.50      0.57    -0.73     1.52 1.02       92       76
# bandInt:train_end                 -0.00      0.00    -0.00     0.00 1.02      173      303
# conditVaried:bandInt:train_end    -0.00      0.00    -0.00     0.00 1.03      102      111

# Family Specific Parameters: 
#       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma   220.06      1.62   216.89   223.23 1.00     5469     2623

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).
# Warning message:
# Parts of the model have not converged (some Rhats are > 1.05). Be careful when analysing the results! We recommend running more iterations and/or setting stronger priors. 



bmtd3_mv <- brm(mvbind(dist,train_end) ~ condit * bandType + (1|bandInt) + (1|id), 
    data=testE1, 
    file=paste0(here::here("data/model_cache","e1_trainEnd_BT_mv_RF2")),
    iter=1000,chains=2, control = list(adapt_delta = .92, max_treedepth = 11))



# summary(bmtd3_mv)
#  Family: MV(gaussian, gaussian) 
#   Links: mu = identity; sigma = identity
#          mu = identity; sigma = identity 
# Formula: dist ~ condit * bandType + (1 | bandInt) + (1 | id) 
#          train_end ~ condit * bandType + (1 | bandInt) + (1 | id) 
#    Data: testE1 (Number of observations: 9491) 
#   Draws: 2 chains, each with iter = 1000; warmup = 500; thin = 1;
#          total post-warmup draws = 1000

# Group-Level Effects: 
# ~bandInt (Number of levels: 6) 
#                        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(dist_Intercept)        92.83     49.67    43.18   142.48 2.48        2       11
# sd(trainend_Intercept)    36.52     36.54     0.00    73.05 2.36        3       11

# ~id (Number of levels: 156) 
#                        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(dist_Intercept)       108.02     48.62    59.42   156.63 2.32        3       21
# sd(trainend_Intercept)    91.70     28.80    62.91   120.48 2.23        3       NA

# Population-Level Effects: 
#                                             Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# dist_Intercept                                  7.84      5.28     2.56    13.12 2.50        2       NA
# trainend_Intercept                             21.03      6.86    14.17    27.88 2.62        2        3
# dist_conditVaried                              -0.19      0.10    -0.29    -0.10 2.24        3       15
# dist_bandTypeExtrapolation                     -0.41      1.11    -1.52     0.71 2.48        2       21
# dist_conditVaried:bandTypeExtrapolation         6.08      6.27    -0.19    12.35 2.37        3       NA
# trainend_conditVaried                          27.67     22.96     4.72    50.62 2.46        2       13
# trainend_bandTypeExtrapolation                  0.00      0.00    -0.00     0.00 1.83        3      123
# trainend_conditVaried:bandTypeExtrapolation     0.00      0.00    -0.00     0.00 1.84        3       65

# Family Specific Parameters: 
#                Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma_dist       264.51     22.04   242.47   286.59 2.28        3       13
# sigma_trainend     0.00      0.00     0.00     0.00 3.05        2       11

# Residual Correlations: 
#                       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# rescor(dist,trainend)    -0.08      0.11    -0.19     0.02 2.89        2       11

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).
# Warning message:
# Parts of the model have not converged (some Rhats are > 1.05). Be careful when analysing the results! We recommend running more iterations and/or setting stronger priors. 













e1 <- readRDS(here("data/e1_08-21-23.rds")) 
e2 <- readRDS(here("data/e2_08-04-23.rds")) 
e3 <- readRDS(here("data/e3_08-04-23.rds")) 

# combine all 3 experiments
d <- rbind(e1,e2,e3)

d <- d |> 
    mutate(trainCon=case_when(
    bandOrder=="Original" ~ "800",
    bandOrder=="Reverse" ~ "600",
    TRUE ~ NA_character_
    ), trainCon=as.numeric(trainCon)) 


nbins=5
train <-  d |> filter(expMode2=="Train") |> group_by(id,condit,fb,bandOrder, vb) |> 
    mutate(Trial_Bin = cut( gt.train, breaks = seq(1, max(gt.train),length.out=nbins+1),include.lowest = TRUE, labels=FALSE)) 
train_max <- train |> filter(Trial_Bin == nbins, bandInt==trainCon)
train_avg <- train_max |> group_by(id,condit,fb,bandOrder) |> summarise(train_end = mean(dist))
train_avg2 <- train_max |> select(id,condit,fb,bandOrder,expMode2,vb,bandInt,dist,vx)

test2 <- d |> filter(expMode2=="Test") |> 
  select(id,condit,fb,bandOrder,expMode2,vb,bandInt,dist,vx) |> 
  rbind(train_avg2) |>
  left_join(train_avg, by=c("id","condit", "fb", "bandOrder")) 

test <- d |> filter(expMode2=="Test") |> left_join(train_avg, by=c("id","condit", "fb", "bandOrder")) |>
  select(id,condit,bandType,bandInt,vb,vx,dist,train_end,fb,bandOrder)




  ltest|> group_by(id,condit) |> pivot_longer(c("dist","train_end"),names_to="var",values_to="value") |> 
  ggplot(aes(x=var,y=value, fill=condit)) + stat_bar + facet_wrap(~var)


test2 |> ggplot(aes(x=bandInt,y=dist,fill=expMode2)) + stat_bar + facet_wrap(~condit+fb+bandOrder)

test2 |> mutate(train_end_q = ntile(train_end,4)) |>  
  ggplot(aes(x=bandInt,y=dist,fill=expMode2)) + stat_bar + facet_wrap(~condit+train_end_q+bandOrder)

test2 |> mutate(train_end_q = ntile(train_end,4)) |>  
  ggplot(aes(x=condit,y=dist,fill=vb)) + stat_bar + 
  facet_wrap(~expMode2+train_end_q+bandOrder,scales="free")

test2 |> mutate(train_end_q = ntile(train_end,4)) |>  
  ggplot(aes(x=vb,y=dist,fill=condit)) + stat_bar + 
  facet_wrap(~expMode2+train_end_q+bandOrder,scales="free")

test2 |> mutate(train_end_q = ntile(train_end,4)) |>  
  ggplot(aes(x=train_end_q,y=dist,fill=condit)) + stat_bar + 
  facet_wrap(~expMode2+vb+bandOrder,scales="free")


test |> ggplot(aes(x=train_end,y=dist,fill=condit)) + 
  stat_summary(geom = "line", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(.9), fun.data = mean_se, width = .4, alpha = .7) + 
    facet_wrap(~vb) +
  labs(x="Band", y="Deviation From Target")

test |> filter(bandType=="Extrapolation") |>
  ggplot(aes(x=train_end,y=dist,fill=condit,col=condit)) + 
  #geom_point() +
  geom_smooth(method="loess") +
    facet_nested_wrap(~bandOrder+vb+fb) +
  labs(x="Band", y="Deviation From Target")

test |> filter(bandType=="Extrapolation") |>
  ggplot(aes(x=train_end,y=dist,fill=condit,col=condit)) + 
  #geom_point() +
  geom_smooth() +
    facet_nested_wrap(~bandOrder+vb+fb) +
  labs(x="Band", y="Deviation From Target")

test |> filter(bandType=="Extrapolation") |>
  ggplot(aes(x=train_end,y=dist,fill=condit,col=condit)) + 
  #geom_point() +
  geom_smooth(method="lm") +
    facet_nested_wrap(~bandOrder+vb) +
  labs(x="Band", y="Deviation From Target")


# create quartiles for train_end
test |> group_by(condit,vb,fb,bandOrder) |>  filter(bandType=="Extrapolation") |>
  mutate(train_end_q = ntile(train_end,4)) |> 
  ggplot(aes(x=train_end_q,y=dist,fill=condit)) + 
 stat_bar + 
    ggh4x::facet_wrap2(bandOrder~vb~fb) +
  labs(x="Band", y="Deviation From Target")


test |> group_by(condit,vb,fb,bandOrder) |>  #filter(bandType=="Extrapolation") |>
  mutate(train_end_q = ntile(train_end,4)) |> 
  ggplot(aes(x=condit,y=vx,fill=vb)) +
  stat_bar + 
    ggh4x::facet_wrap2(~train_end_q+bandOrder+fb) 