[[ReadItLater]] [[Article]]

# [HTW Project - HTW E2 Testing](http://localhost:4200/Analysis/e2_test.html)

[Figure 1](http://localhost:4200/Analysis/e2_test.html#fig-design-e2) illustrates the design of Experiment 2. The stages of the experiment (i.e. training, testing no-feedback, test with feedback), are identical to that of Experiment 1. The only change is that Experiment 2 participants train, and then test, on bands in the reverse order of Experiment 1 (i.e. training on the softer bands; and testing on the harder bands).

Code

```
digraph {
  graph [layout = dot, rankdir = LR]

  // define the global styles of the nodes
  node [shape = rectangle, style = filled]

  data1 [label = " Varied Training \n100-300\n350-550\n600-800", fillcolor = "#FF0000"]
  data2 [label = " Constant Training \n600-800", fillcolor = "#00A08A"]
  Test3 [label = "    Final Test \n  Novel With Feedback  \n800-1000\n1000-1200\n1200-1400", fillcolor = "#ECCBAE"]

  // edge definitions with the node IDs
  data1 -> Test1
  data2 -> Test1
  subgraph cluster {
    label = "Test Phase \n(Counterbalanced Order)"
    Test1 [label = "Test  \nNovel Bands  \n800-1000\n1000-1200\n1200-1400", fillcolor = "#ECCBAE"]
    Test2 [label = "  Test \n  Varied Training Bands  \n100-300\n350-550\n600-800", fillcolor = "#ECCBAE"]
    Test1 -> Test2
  }

  Test2 -> Test3
}
```

cluster Test Phase (Counterbalanced Order)data1 Varied Training 100-300350-550600-800Test1 Test  Novel Bands  800-10001000-12001200-1400data1->Test1 data2 Constant Training 600-800data2->Test1 Test3    Final Test  Novel With Feedback  800-10001000-12001200-1400Test2  Test  Varied Training Bands  100-300350-550600-800Test1->Test2 Test2->Test3

Figure 1: Experiment 2 Design. Constant and Varied participants complete different training conditions. The training and testing bands are the reverse of Experiment 1.

## Results

### Testing Phase - No feedback.

In the first part of the testing phase, participants are tested from each of the velocity bands, and receive no feedback after each throw.

#### Deviation From Target Band

Descriptive summaries testing deviation data are provided in [Table 1](http://localhost:4200/Analysis/e2_test.html#tbl-e2-test-nf-deviation) and [Figure 2](http://localhost:4200/Analysis/e2_test.html#fig-e2-test-dev). To model differences in accuracy between groups, we used Bayesian mixed effects regression models to the trial level data from the testing phase. The primary model predicted the absolute deviation from the target velocity band (dist) as a function of training condition (condit), target velocity band (band), and their interaction, with random intercepts and slopes for each participant (id).

\\\[\\begin{equation} dist\_{ij} = \\beta\_0 + \\beta\_1 \\cdot condit\_{ij} + \\beta\_2 \\cdot band\_{ij} + \\beta\_3 \\cdot condit\_{ij} \\cdot band\_{ij} + b\_{0i} + b\_{1i} \\cdot band\_{ij} + \\epsilon\_{ij} \\end{equation}\\\]

Code

````
```{r}
#| label: tbl-e2-test-nf-deviation
#| tbl-cap: "Testing Deviation - Empirical Summary"
#| tbl-subcap: ["Constant Testing - Deviation", "Varied Testing - Deviation"]

result <- test_summary_table(testE2, "dist","Deviation", mfun = list(mean = mean, median = median, sd = sd))
result$constant 
result$varied 
# make kable table with smaller font size
# result$constant |> kbl(caption="Constant Testing - Deviation",booktabs=T,escape=F) |> kable_styling(font_size = 7)
```
````

Table 1: Testing Deviation - Empirical Summary

Summary of Deviation - Constant
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Extrapolation | 206 | 48 | 317 |
| 350-550 | Extrapolation | 194 | 86 | 268 |
| 600-800 | Trained | 182 | 112 | 240 |
| 800-1000 | Extrapolation | 200 | 129 | 233 |
| 1000-1200 | Extrapolation | 238 | 190 | 234 |
| 1200-1400 | Extrapolation | 311 | 254 | 288 |

Summary of Deviation - Varied
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Trained | 153 | 25 | 266 |
| 350-550 | Trained | 138 | 53 | 233 |
| 600-800 | Trained | 160 | 120 | 183 |
| 800-1000 | Extrapolation | 261 | 207 | 257 |
| 1000-1200 | Extrapolation | 305 | 258 | 273 |
| 1200-1400 | Extrapolation | 363 | 314 | 297 |

Code

```
testE2 |>  ggplot(aes(x = vb, y = dist,fill=condit)) +
    stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(.9), fun.data = mean_se, width = .4, alpha = .7) + 
  labs(x="Band", y="Deviation From Target")
```

![](md_extract/assets/fig-e2-test-dev-1.png)

Figure 2: E2. Deviations from target band during testing without feedback stage.

Code

```
#contrasts(test$condit) 
# contrasts(testE2$vb)

modelName <- "e2_testDistBand_RF_5K"
e2_distBMM <- brm(dist ~ condit * bandInt + (1 + bandInt|id),
                      data=testE2,file=paste0(here::here("data/model_cache",modelName)),
                      iter=5000,chains=4)
mp2 <- GetModelStats(e2_distBMM) |> kable(escape=F,booktabs=T)
mp2
```

Table 2: Experiment 2. Bayesian Mixed Model predicting absolute deviation as a function of condition (Constant vs. Varied) and Velocity Band

| Term | Estimate | 95% CrI Lower | 95% CrI Upper | pd |
| --- | --- | --- | --- | --- |
| Intercept | 151.71 | 90.51 | 215.86 | 1.00 |
| conditVaried | \-70.33 | \-156.87 | 16.66 | 0.94 |
| Band | 0.10 | 0.02 | 0.18 | 1.00 |
| condit\*Band | 0.12 | 0.02 | 0.23 | 0.99 |

Code

```
e2_distBMM |> 
  emmeans("condit",by="bandInt",at=list(bandInt=c(100,350,600,800,1000,1200)),
          epred = TRUE, re_formula = NA) |> 
  pairs() |> gather_emmeans_draws()  |> 
   summarize(median_qi(.value),pd=sum(.value>0)/n()) |>
   select(contrast,Band=bandInt,value=y,lower=ymin,upper=ymax,pd) |> 
   mutate(across(where(is.numeric), \(x) round(x, 2)),
          pd=ifelse(value<0,1-pd,pd)) |>
   kbl(caption="Contrasts")
```

Contrasts
| contrast | Band | value | lower | upper | pd |
| --- | --- | --- | --- | --- | --- |
| Constant - Varied | 100 | 57.57 | \-20.48 | 135.32 | 0.93 |
| Constant - Varied | 350 | 26.60 | \-30.93 | 83.84 | 0.83 |
| Constant - Varied | 600 | \-4.30 | \-46.73 | 38.52 | 0.58 |
| Constant - Varied | 800 | \-29.30 | \-69.38 | 11.29 | 0.92 |
| Constant - Varied | 1000 | \-54.62 | \-101.06 | \-5.32 | 0.98 |
| Constant - Varied | 1200 | \-79.63 | \-139.47 | \-15.45 | 0.99 |

Code

```
coef_details <- get_coef_details(e2_distBMM, "conditVaried")
```

The model predicting absolute deviation showed a modest tendency for the varied training group to have lower deviation compared to the constant training group (β = -70.33, 95% CI \[-156.87, 16.66\]),with 94% of the posterior distribution being less than 0. This suggests a potential benefit of training with variation, though the evidence is not definitive.

(SHOULD PROBABLY DO ALTERNATE ANALYSIS THAT ONLY CONSIDERS THE NOVEL EXTRAPOLATION BANDS)

Code

```
e2_distBMM |> emmeans( ~condit + bandInt, 
                       at = list(bandInt = c(100, 350, 600, 800, 1000, 1200))) |>
  gather_emmeans_draws() |>
 condEffects()+
  scale_x_continuous(breaks = c(100, 350, 600, 800, 1000, 1200), 
                     labels = levels(testE2$vb), 
                     limits = c(0, 1400)) 
``` 

![](md_extract/assets/fig-e2-bmm-dist-1.png)

Figure 3: E2. Conditioinal Effect of Training Condition and Band. Ribbon indicated 95% Credible Intervals.

#### Discrimination between Velocity Bands

In addition to accuracy/deviation. We also assessed the ability of participants to reliably discriminate between the velocity bands (i.e. responding differently when prompted for band 600-800 than when prompted for band 150-350). [Table 3](http://localhost:4200/Analysis/e2_test.html#tbl-e2-test-nf-vx) shows descriptive statistics of this measure, and Figure 1 visualizes the full distributions of throws for each combination of condition and velocity band. To quantify discrimination, we again fit Bayesian Mixed Models as above, but this time the dependent variable was the raw x velocity generated by participants.

\\\[\\begin{equation} vx\_{ij} = \\beta\_0 + \\beta\_1 \\cdot condit\_{ij} + \\beta\_2 \\cdot bandInt\_{ij} + \\beta\_3 \\cdot condit\_{ij} \\cdot bandInt\_{ij} + b\_{0i} + b\_{1i} \\cdot bandInt\_{ij} + \\epsilon\_{ij} \\end{equation}\\\]

Code

```
testE2 %>% group_by(id,vb,condit) |> plot_distByCondit()
```

![](md_extract/assets/fig-e2-test-vx-1.png)

Figure 4: E2 testing x velocities. Translucent bands with dash lines indicate the correct range for each velocity band.

Code

````
```{r}
#| label: tbl-e2-test-nf-vx
#| tbl-cap: "Testing vx - Empirical Summary"
#| tbl-subcap: ["Constant Testing - vx", "Varied Testing - vx"]
#| layout-ncol: 2

result <- test_summary_table(testE2, "vx","X Velocity" ,mfun = list(mean = mean, median = median, sd = sd))
result$constant 
result$varied 
```
````

Table 3: Testing vx - Empirical Summary

Summary of X Velocity - Constant
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Extrapolation | 457 | 346 | 354 |
| 350-550 | Extrapolation | 597 | 485 | 368 |
| 600-800 | Trained | 728 | 673 | 367 |
| 800-1000 | Extrapolation | 953 | 913 | 375 |
| 1000-1200 | Extrapolation | 1064 | 1012 | 408 |
| 1200-1400 | Extrapolation | 1213 | 1139 | 493 |

Summary of X Velocity - Varied
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Trained | 410 | 323 | 297 |
| 350-550 | Trained | 582 | 530 | 303 |
| 600-800 | Trained | 696 | 641 | 316 |
| 800-1000 | Extrapolation | 910 | 848 | 443 |
| 1000-1200 | Extrapolation | 1028 | 962 | 482 |
| 1200-1400 | Extrapolation | 1095 | 1051 | 510 |

Code

```
e2_vxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                        data=testE2,file=paste0(here::here("data/model_cache", "e2_testVxBand_RF_5k")),
                        iter=5000,chains=4,silent=0,
                        control=list(adapt_delta=0.94, max_treedepth=13))
mt3 <-GetModelStats(e2_vxBMM ) |> kable(escape=F,booktabs=T)
mt3
```

Table 4: Experiment 2. Bayesian Mixed Model Predicting Vx as a function of condition (Constant vs. Varied) and Velocity Band

| Term | Estimate | 95% CrI Lower | 95% CrI Upper | pd |
| --- | --- | --- | --- | --- |
| Intercept | 362.64 | 274.85 | 450.02 | 1.00 |
| conditVaried | \-8.56 | \-133.97 | 113.98 | 0.55 |
| Band | 0.71 | 0.58 | 0.84 | 1.00 |
| condit\*Band | \-0.06 | \-0.24 | 0.13 | 0.73 |

Code

```
cd1 <- get_coef_details(e2_vxBMM, "conditVaried")
sc1 <- get_coef_details(e2_vxBMM, "bandInt")
intCoef1 <- get_coef_details(e2_vxBMM, "conditVaried:bandInt")
```

See [Table 4](http://localhost:4200/Analysis/e2_test.html#tbl-e2-bmm-vx) for the full model results.

When examining discrimination ability using the model predicting raw x-velocity, the results were less clear than those of the absolute deviation analysis. The slope on Velocity Band (β = 0.71, 95% CrI \[0.58, 0.84\]) indicates that participants showed good discrimination between bands overall. However, the interaction term suggested this effect was not modulated by training condition (β = -0.06, 95% CrI \[-0.24, 0.13\]) Thus, while varied training may provide some advantage for accuracy, it does not appear to impair the ability to discriminate between velocity bands relative to constant training.

Code

```
e2_vxBMM |> emmeans( ~condit + bandInt, 
                       at = list(bandInt = c(100, 350, 600, 800, 1000, 1200))) |>
  gather_emmeans_draws() |>
  ggplot(aes(x = bandInt, y = .value, color = condit, fill = condit)) + 
  stat_dist_pointinterval() +
  stat_lineribbon(alpha = .25, size = 1, .width = c(.95)) +
  ylab("Predicted X Velocity") + xlab("Band")+
  scale_x_continuous(breaks = c(100, 350, 600, 800, 1000, 1200), 
                     labels = levels(testE2$vb), 
                     limits = c(0, 1400)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) 
``` 

![](md_extract/assets/fig-e2-bmm-vx-1.png)

Figure 5: Conditional effect of training condition and Band. Ribbons indicate 95% HDI.

Code

```
new_data_grid=map_dfr(1, ~data.frame(unique(testE2[,c("id","condit","bandInt")]))) |> 
  dplyr::arrange(id,bandInt) |> 
  mutate(condit_dummy = ifelse(condit == "Varied", 1, 0)) 

indv_coefs <- coef(e2_vxBMM)$id |> 
  as_tibble(rownames="id") |> 
  select(id, starts_with("Est")) |>
  left_join(e2Sbjs, by=join_by(id) ) |> 
  group_by(condit) |> 
  mutate(rank = rank(desc(Estimate.bandInt)),
         intErrorRank=rank((Est.Error.Intercept)),
         bandErrorRank=rank((Est.Error.bandInt)),
         nCond = n()) |> arrange(intErrorRank)

fixed_effects <- e2_vxBMM |> 
  spread_draws(`^b_.*`,regex=TRUE) |> arrange(.chain,.draw,.iteration)


random_effects <- e2_vxBMM |> 
  gather_draws(`^r_id.*[[ReadItLater]] [[Article]]

# [HTW Project - HTW E2 Testing](http://localhost:4200/Analysis/e2_test.html)

, regex = TRUE, ndraws = 2000) |> 
  separate(.variable, into = c("effect", "id", "term"), sep = "\\[|,|\\]") |> 
  mutate(id = factor(id,levels=levels(testE2$id))) |> 
  pivot_wider(names_from = term, values_from = .value) |> arrange(id,.chain,.draw,.iteration)

cd <- left_join(random_effects, fixed_effects, by = join_by(".chain", ".iteration", ".draw")) |> 
  rename(bandInt_RF = bandInt) |>
  mutate(Slope=bandInt_RF+b_bandInt) |> group_by(id) 

cdMed <- cd |> group_by(id) |> median_qi(Slope)  |> 
  left_join(e2Sbjs, by=join_by(id)) |> group_by(condit) |>
  mutate(rankSlope=rank(Slope)) |> arrange(rankSlope)

cdMed %>% ggplot(aes(y=rankSlope, x=Slope,fill=condit,color=condit)) + 
  geom_pointrange(aes(xmin=.lower , xmax=.upper)) + 
  labs(x="Estimated Slope", y="Participant")  + facet_wrap(~condit)  

# cdMed |>  ggplot(aes(x = condit, y = Slope,fill=condit)) +
#     stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
#     stat_summary(geom = "errorbar", position=position_dodge(.9), fun.data = mean_se, width = .4, alpha = .7) + 
#   geom_jitter()
#   labs(x="Band", y="Deviation From Target")
```

Code

```
new_data_grid=map_dfr(1, ~data.frame(unique(testE2[,c("id","condit","bandInt")]))) |> 
  dplyr::arrange(id,bandInt) |> 
  mutate(condit_dummy = ifelse(condit == "Varied", 1, 0)) 

indv_coefs <- as_tibble(coef(e2_vxBMM)$id, rownames="id")|> 
  select(id, starts_with("Est")) |>
  left_join(e2Sbjs, by=join_by(id) ) 


fixed_effects <- e2_vxBMM |> 
  spread_draws(`^b_.*`,regex=TRUE) |> arrange(.chain,.draw,.iteration)


random_effects <- e2_vxBMM |> 
  gather_draws(`^r_id.*[[ReadItLater]] [[Article]]

# [HTW Project - HTW E2 Testing](http://localhost:4200/Analysis/e2_test.html)

, regex = TRUE, ndraws = 1500) |> 
  separate(.variable, into = c("effect", "id", "term"), sep = "\\[|,|\\]") |> 
  mutate(id = factor(id,levels=levels(testE2$id))) |> 
  pivot_wider(names_from = term, values_from = .value) |> arrange(id,.chain,.draw,.iteration)


 indvDraws <- left_join(random_effects, fixed_effects, by = join_by(".chain", ".iteration", ".draw")) |> 
  rename(bandInt_RF = bandInt,RF_Intercept=Intercept) |>
  right_join(new_data_grid, by = join_by("id")) |> 
  mutate(
    Slope = bandInt_RF+b_bandInt,
    Intercept= RF_Intercept + b_Intercept,
    estimate = (b_Intercept + RF_Intercept) + (bandInt*(b_bandInt+bandInt_RF)) + (bandInt * condit_dummy) * `b_conditVaried:bandInt`,
    SlopeInt = Slope + (`b_conditVaried:bandInt`*condit_dummy)
  ) 

  indvSlopes <- indvDraws |> group_by(id) |> median_qi(Slope,SlopeInt, Intercept,b_Intercept,b_bandInt) |>
  left_join(e2Sbjs, by=join_by(id)) |> group_by(condit) |>
    select(id,condit,Intercept,b_Intercept,starts_with("Slope"),b_bandInt, n) |>
  mutate(rankSlope=rank(Slope)) |> arrange(rankSlope)   |> ungroup()
 
  
  indvSlopes |> mutate(Condition=condit) |>  group_by(Condition) |> 
    reframe(enframe(quantile(SlopeInt, c(0.0,0.25, 0.5, 0.75,1)), "quantile", "SlopeInt")) |> 
  pivot_wider(names_from=quantile,values_from=SlopeInt,names_prefix="Q_") |>
  group_by(Condition) |>
  summarise(across(starts_with("Q"), list(mean = mean))) |> kbl()
```

Table 5: Slope coefficients by quartile, per condition

| Condition | Q\_0%\_mean | Q\_25%\_mean | Q\_50%\_mean | Q\_75%\_mean | Q\_100%\_mean |
| --- | --- | --- | --- | --- | --- |
| Constant | \-0.2874640 | 0.3911363 | 0.6988131 | 1.075805 | 1.617691 |
| Varied | \-0.2427695 | 0.3050573 | 0.6919588 | 0.949496 | 1.811897 |

[Figure 6](http://localhost:4200/Analysis/e2_test.html#fig-e2-bmm-bx2) shows the distributions of estimated slopes relating velocity band to x velocity for each participant, ordered from lowest to highest within condition. Slope values are lower overall for varied training compared to constant training. Figure Xb plots the density of these slopes for each condition. The distribution for varied training has more mass at lower values than the constant training distribution. Both figures illustrate the model’s estimate that varied training resulted in less discrimination between velocity bands, evidenced by lower slopes on average.

Code

```
  indvSlopes |> ggplot(aes(y=rankSlope, x=SlopeInt,fill=condit,color=condit)) + 
  geom_pointrange(aes(xmin=SlopeInt.lower , xmax=SlopeInt.upper)) + 
  labs(x="Estimated Slope", y="Participant")  + facet_wrap(~condit) +
   ggplot(indvSlopes, aes(x = SlopeInt, color = condit)) + 
  geom_density() + labs(x="Slope Coefficient",y="Density")
```

![](md_extract/assets/fig-e2-bmm-bx2-1.png)

(a) Slope estimates by participant - ordered from lowest to highest within each condition.

Figure 6: Slope distributions between condition

Code

```
nSbj <- 3
indvDraws  |> indv_model_plot(indvSlopes, testE2Avg, SlopeInt,rank_variable=Slope,n_sbj=nSbj,"max")
indvDraws |> indv_model_plot(indvSlopes, testE2Avg,SlopeInt, rank_variable=Slope,n_sbj=nSbj,"min")
```

![](md_extract/assets/fig-e2-indv-slopes-1.png)

(a) subset with largest slopes

![](md_extract/assets/fig-e2-indv-slopes-2.png)

(b) subset with smallest slopes

Figure 7: Subset of Varied and Constant Participants with the smallest and largest estimated slope values. Red lines represent the best fitting line for each participant, gray lines are 200 random samples from the posterior distribution. Colored points and intervals at each band represent the empirical median and 95% HDI.