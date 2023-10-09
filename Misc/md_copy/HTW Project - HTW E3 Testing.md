[[ReadItLater]] [[Article]]

# [HTW Project - HTW E3 Testing](http://localhost:4200/Analysis/e3_test.html)

The major manipulation adjustment of experiment 3 is for participants to receive ordinal feedback during training, in contrast to the continuous feedback of the earlier experiments. Ordinal feedback informs participants whether a throw was too soft, too hard, or fell within the target velocity range. Experiment 3 participants were randomly assigned to both a training condition (Constant vs. Varied) and a Band Order condition (original order used in Experiment 1, or the Reverse order of Experiment 2).

## Results

### Testing Phase - No feedback.

In the first part of the testing phase, participants are tested from each of the velocity bands, and receive no feedback after each throw. Note that these no-feedback testing trials are identical to those of Experiment 1 and 2, as the ordinal feedback only occurs during the training phase, and final testing phase, of Experiment 3.

#### Deviation From Target Band

Descriptive summaries testing deviation data are provided in [Table 1](http://localhost:4200/Analysis/e3_test.html#tbl-e3-test-nf-deviation) and [Figure 1](http://localhost:4200/Analysis/e3_test.html#fig-e3-test-dev). To model differences in accuracy between groups, we fit Bayesian mixed effects regression models to the trial level data from the testing phase. The primary model predicted the absolute deviation from the target velocity band (dist) as a function of training condition (condit), target velocity band (band), and their interaction, with random intercepts and slopes for each participant (id).

Code

````
```{r}
#| label: tbl-e3-test-nf-deviation
#| tbl-cap: "Testing Deviation - Empirical Summary"
#| tbl-subcap: ["Constant Testing - Deviation", "Varied Testing - Deviation"]
#| layout-ncol: 2

resultOrig <- test_summary_table(testE3 |> filter(bandOrder=="Original"), "dist","Deviation", mfun = list(mean = mean, median = median, sd = sd))
resultOrig$constant 
resultOrig$varied 

resultRev <- test_summary_table(testE3 |> filter(bandOrder=="Reverse"), "dist","Deviation", mfun = list(mean = mean, median = median, sd = sd))
resultRev$constant 
resultRev$varied 
```
````

Table 1: Testing Deviation - Empirical Summary

Summary of Deviation - Constant
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Extrapolation | 396 | 325 | 350 |
| 350-550 | Extrapolation | 278 | 176 | 299 |
| 600-800 | Extrapolation | 173 | 102 | 215 |
| 800-1000 | Trained | 225 | 126 | 284 |
| 1000-1200 | Extrapolation | 253 | 192 | 271 |
| 1200-1400 | Extrapolation | 277 | 210 | 262 |

Summary of Deviation - Varied
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Extrapolation | 383 | 254 | 385 |
| 350-550 | Extrapolation | 287 | 154 | 318 |
| 600-800 | Extrapolation | 213 | 140 | 244 |
| 800-1000 | Trained | 199 | 142 | 209 |
| 1000-1200 | Trained | 222 | 163 | 221 |
| 1200-1400 | Trained | 281 | 227 | 246 |

Summary of Deviation - Constant
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Extrapolation | 403 | 334 | 383 |
| 350-550 | Extrapolation | 246 | 149 | 287 |
| 600-800 | Trained | 155 | 82 | 209 |
| 800-1000 | Extrapolation | 207 | 151 | 241 |
| 1000-1200 | Extrapolation | 248 | 220 | 222 |
| 1200-1400 | Extrapolation | 322 | 281 | 264 |

Summary of Deviation - Varied
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Trained | 153 | 0 | 307 |
| 350-550 | Trained | 147 | 55 | 258 |
| 600-800 | Trained | 159 | 107 | 192 |
| 800-1000 | Extrapolation | 221 | 160 | 235 |
| 1000-1200 | Extrapolation | 244 | 185 | 235 |
| 1200-1400 | Extrapolation | 324 | 264 | 291 |

Code

```
testE3 |>  ggplot(aes(x = vb, y = dist,fill=condit)) +
    stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(.9), fun.data = mean_se, width = .4, alpha = .7) + 
  labs(x="Band", y="Deviation From Target") + facet_wrap(~bandOrder)
```

![](md_extract/assets/fig-e3-test-dev-1.png)

Figure 1: e3. Deviations from target band during testing without feedback stage.

Code

```
#contrasts(test$condit) 

# contrasts(testE3$vb)

modelName <- "e3_testDistBand_RF_5K"
e3_distBMM <- brm(dist ~ condit * bandInt + (1 + bandInt|id),
                      data=testE3,file=paste0(here::here("data/model_cache",modelName)),
                      iter=5000,chains=4)
mp3 <- GetModelStats(e3_distBMM) |> kable(escape=F,booktabs=T)
mp3
```

Table 2: Experiment 3. Bayesian Mixed Model predicting absolute deviation as a function of condition (Constant vs. Varied) and Velocity Band

| Term | Estimate | 95% CrI Lower | 95% CrI Upper | pd |
| --- | --- | --- | --- | --- |
| Intercept | 306.47 | 243.89 | 368.75 | 1.00 |
| conditVaried | \-90.65 | \-182.79 | 3.75 | 0.97 |
| Band | \-0.07 | \-0.13 | 0.00 | 0.97 |
| condit\*Band | 0.09 | \-0.01 | 0.19 | 0.96 |

Code

```
cd1 <- get_coef_details(e3_distBMM, "conditVaried")
sc1 <- get_coef_details(e3_distBMM, "bandInt")
intCoef1 <- get_coef_details(e3_distBMM, "conditVaried:bandInt")
```

The effect of training condition in Experiment 3 showed a similar pattern to Experiment 2, with the varied group tending to have lower deviation than the constant group (β = -90.65, 95% CrI \[-182.79, 3.75\]), with 97% of the posterior distribution falling under 0.

(NEED TO CONTROL FOR BAND ORDER HERE)

Code

```
e3_distBMM |> emmeans( ~condit + bandInt, 
                       at = list(bandInt = c(100, 350, 600, 800, 1000, 1200))) |>
  gather_emmeans_draws() |>
  ggplot(aes(x = bandInt, y = .value, color = condit, fill = condit)) + 
  stat_dist_pointinterval() +
  stat_lineribbon(alpha = .25, size = 1, .width = c(.95)) +
    ylab("Predicted Deviation") + xlab("Velocity Band")+
  scale_x_continuous(breaks = c(100, 350, 600, 800, 1000, 1200), 
                     labels = levels(testE3$vb), 
                     limits = c(0, 1400)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) 
``` 

![](md_extract/assets/fig-e3-bmm-dist-1.png)

Figure 2: e3. Conditioinal Effect of Training Condition and Band. Ribbon indicated 95% Credible Intervals.

#### Discrimination between Velocity Bands

In addition to accuracy/deviation. We also assessed the ability of participants to reliably discriminate between the velocity bands (i.e. responding differently when prompted for band 600-800 than when prompted for band 150-350). [Table 3](http://localhost:4200/Analysis/e3_test.html#tbl-e3-test-nf-vx) shows descriptive statistics of this measure, and Figure 1 visualizes the full distributions of throws for each combination of condition and velocity band. To quantify discrimination, we again fit Bayesian Mixed Models as above, but this time the dependent variable was the raw x velocity generated by participants.

\\\[\\begin{equation} vx\_{ij} = \\beta\_0 + \\beta\_1 \\cdot condit\_{ij} + \\beta\_2 \\cdot bandInt\_{ij} + \\beta\_3 \\cdot condit\_{ij} \\cdot bandInt\_{ij} + b\_{0i} + b\_{1i} \\cdot bandInt\_{ij} + \\epsilon\_{ij} \\end{equation}\\\]

Code

```
# testE3 |> filter(bandOrder=="Original")|> group_by(id,vb,condit) |> plot_distByCondit()
# testE3 |> filter(bandOrder=="Reverse")|> group_by(id,vb,condit) |> plot_distByCondit() +ggtitle("test")

testE3 |> group_by(id,vb,condit,bandOrder) |> plot_distByCondit() + 
  facet_wrap(bandOrder~condit,scale="free_x") 
``` 

![](md_extract/assets/fig-e3-test-vx-1.png)

Figure 3: e3 testing x velocities. Translucent bands with dash lines indicate the correct range for each velocity band.

Code

````
```{r}
#| label: tbl-e3-test-nf-vx
#| tbl-cap: "Testing vx - Empirical Summary"
#| tbl-subcap: ["Constant Testing - vx", "Varied Testing - vx"]
#| layout-ncol: 2

resultOrig <- test_summary_table(testE3 |> filter(bandOrder=="Original"), "vx","X Velocity", mfun = list(mean = mean, median = median, sd = sd))
resultOrig$constant 
resultOrig$varied 

resultRev <- test_summary_table(testE3 |> filter(bandOrder=="Reverse"), "vx","X Velocity", mfun = list(mean = mean, median = median, sd = sd))
resultRev$constant 
resultRev$varied 
```
````

Table 3: Testing vx - Empirical Summary

Summary of X Velocity - Constant
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Extrapolation | 680 | 625 | 370 |
| 350-550 | Extrapolation | 771 | 716 | 357 |
| 600-800 | Extrapolation | 832 | 786 | 318 |
| 800-1000 | Trained | 1006 | 916 | 417 |
| 1000-1200 | Extrapolation | 1149 | 1105 | 441 |
| 1200-1400 | Extrapolation | 1180 | 1112 | 443 |

Summary of X Velocity - Varied
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Extrapolation | 667 | 554 | 403 |
| 350-550 | Extrapolation | 770 | 688 | 383 |
| 600-800 | Extrapolation | 869 | 814 | 358 |
| 800-1000 | Trained | 953 | 928 | 359 |
| 1000-1200 | Trained | 1072 | 1066 | 388 |
| 1200-1400 | Trained | 1144 | 1093 | 426 |

Summary of X Velocity - Constant
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Extrapolation | 684 | 634 | 406 |
| 350-550 | Extrapolation | 729 | 679 | 350 |
| 600-800 | Trained | 776 | 721 | 318 |
| 800-1000 | Extrapolation | 941 | 883 | 387 |
| 1000-1200 | Extrapolation | 1014 | 956 | 403 |
| 1200-1400 | Extrapolation | 1072 | 1014 | 442 |

Summary of X Velocity - Varied
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Trained | 392 | 270 | 343 |
| 350-550 | Trained | 540 | 442 | 343 |
| 600-800 | Trained | 642 | 588 | 315 |
| 800-1000 | Extrapolation | 943 | 899 | 394 |
| 1000-1200 | Extrapolation | 1081 | 1048 | 415 |
| 1200-1400 | Extrapolation | 1185 | 1129 | 500 |

Code

```
e3_vxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                        data=testE3,file=paste0(here::here("data/model_cache", "e3_testVxBand_RF_5k")),
                        iter=5000,chains=4,silent=0,
                        control=list(adapt_delta=0.94, max_treedepth=13))
mt4 <-GetModelStats(e3_vxBMM ) |> kable(escape=F,booktabs=T)
mt4
```

Table 4: Experiment 3. Bayesian Mixed Model Predicting Vx as a function of condition (Constant vs. Varied) and Velocity Band

| Term | Estimate | 95% CrI Lower | 95% CrI Upper | pd |
| --- | --- | --- | --- | --- |
| Intercept | 607.67 | 536.02 | 679.87 | 1 |
| conditVaried | \-167.76 | \-277.14 | \-64.08 | 1 |
| Band | 0.44 | 0.35 | 0.52 | 1 |
| condit\*Band | 0.18 | 0.06 | 0.31 | 1 |

Code

```
cd1 <- get_coef_details(e3_vxBMM, "conditVaried")
sc1 <- get_coef_details(e3_vxBMM, "bandInt")
intCoef1 <- get_coef_details(e3_vxBMM, "conditVaried:bandInt")
```

See [Table 4](http://localhost:4200/Analysis/e3_test.html#tbl-e3-bmm-vx) for the full model results.

Slope estimates for experiment 3 suggest that participants were capable of distinguishing between velocity bands even when provided only ordinal feedback during training (β = 0.44, 95% CrI \[0.35, 0.52\]). Unlike the previous two experiments, the posterior distribution for the interaction between condition and band was consistently positive, suggestive of superior discrimination for the varied participants β = 0.18, 95% CrI \[0.06, 0.31\].

Code

```
e3_vxBMM |> emmeans( ~condit + bandInt, 
                       at = list(bandInt = c(100, 350, 600, 800, 1000, 1200))) |>
  gather_emmeans_draws() |>
  ggplot(aes(x = bandInt, y = .value, color = condit, fill = condit)) + 
  stat_dist_pointinterval() +
  stat_lineribbon(alpha = .25, size = 1, .width = c(.95)) +
  ylab("Predicted X Velocity") + xlab("Band")+
  scale_x_continuous(breaks = c(100, 350, 600, 800, 1000, 1200), 
                     labels = levels(testE3$vb), 
                     limits = c(0, 1400)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) 
``` 

![](md_extract/assets/fig-e3-bmm-vx-1.png)

Figure 4: Conditional effect of training condition and Band. Ribbons indicate 95% HDI.

Code

```
new_data_grid=map_dfr(1, ~data.frame(unique(testE3[,c("id","condit","bandInt")]))) |> 
  dplyr::arrange(id,bandInt) |> 
  mutate(condit_dummy = ifelse(condit == "Varied", 1, 0)) 

indv_coefs <- as_tibble(coef(e3_vxBMM)$id, rownames="id")|> 
  select(id, starts_with("Est")) |>
  left_join(e3Sbjs, by=join_by(id) ) 


fixed_effects <- e3_vxBMM |> 
  spread_draws(`^b_.*`,regex=TRUE) |> arrange(.chain,.draw,.iteration)


random_effects <- e3_vxBMM |> 
  gather_draws(`^r_id.*[[ReadItLater]] [[Article]]

# [HTW Project - HTW E3 Testing](http://localhost:4200/Analysis/e3_test.html)

, regex = TRUE, ndraws = 1500) |> 
  separate(.variable, into = c("effect", "id", "term"), sep = "\\[|,|\\]") |> 
  mutate(id = factor(id,levels=levels(testE3$id))) |> 
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
  left_join(e3Sbjs, by=join_by(id)) |> group_by(condit) |>
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
| Constant | \-0.3546215 | 0.1158112 | 0.4250213 | 0.6400295 | 1.935485 |
| Varied | \-0.3138630 | 0.2473060 | 0.5963830 | 0.9144855 | 1.753327 |

[Figure 5](http://localhost:4200/Analysis/e3_test.html#fig-e3-bmm-bx2) shows the distributions of estimated slopes relating velocity band to x velocity for each participant, ordered from lowest to highest within condition. Slope values are lower overall for varied training compared to constant training. Figure Xb plots the density of these slopes for each condition. The distribution for varied training has more mass at lower values than the constant training distribution. Both figures illustrate the model’s estimate that varied training resulted in less discrimination between velocity bands, evidenced by lower slopes on average.

Code

```
indvSlopes |> ggplot(aes(y=rankSlope, x=SlopeInt,fill=condit,color=condit)) + 
  geom_pointrange(aes(xmin=SlopeInt.lower , xmax=SlopeInt.upper)) + 
  labs(x="Estimated Slope", y="Participant")  + facet_wrap(~condit) +
ggplot(indvSlopes, aes(x = SlopeInt, color = condit)) + 
geom_density() + labs(x="Slope Coefficient",y="Density") 
``` 

![](md_extract/assets/fig-e3-bmm-bx2-1.png)

(a) Slope estimates by participant - ordered from lowest to highest within each condition.

Figure 5: Slope distributions between condition

Code

```
nSbj <- 3
indvDraws  |> indv_model_plot(indvSlopes, testE3Avg, SlopeInt,rank_variable=Slope,n_sbj=nSbj,"max")
indvDraws |> indv_model_plot(indvSlopes, testE3Avg,SlopeInt, rank_variable=Slope,n_sbj=nSbj,"min")
```

![](md_extract/assets/fig-e3-indv-slopes-1.png)

(a) subset with largest slopes

![](md_extract/assets/fig-e3-indv-slopes-2.png)

(b) subset with smallest slopes

Figure 6: Subset of Varied and Constant Participants with the smallest and largest estimated slope values. Red lines represent the best fitting line for each participant, gray lines are 200 random samples from the posterior distribution. Colored points and intervals at each band represent the empirical median and 95% HDI.