[[ReadItLater]] [[Article]]

# [HTW Project - HTW E1 Testing](http://localhost:4200/Analysis/e1_test.html)

## Analyses Strategy

All data processing and statistical analyses were performed in R version 4.31 Team ([2020](http://localhost:4200/Analysis/e1_test.html#ref-rcoreteamLanguageEnvironmentStatistical2020)). To assess differences between groups, we used Bayesian Mixed Effects Regression. Model fitting was performed with the brms package in R Bürkner ([2017](http://localhost:4200/Analysis/e1_test.html#ref-burknerBrmsPackageBayesian2017)), and descriptive stats and tables were extracted with the BayestestR package Makowski et al. ([2019](http://localhost:4200/Analysis/e1_test.html#ref-makowskiBayestestRDescribingEffects2019a)). Mixed effects regression enables us to take advantage of partial pooling, simultaneously estimating parameters at the individual and group level. Our use of Bayesian, rather than frequentist methods allows us to directly quantify the uncertainty in our parameter estimates, as well as circumventing convergence issues common to the frequentist analogues of our mixed models. For each model, we report the median values of the posterior distribution, and 95% credible intervals.

Each model was set to run with 4 chains, 5000 iterations per chain, with the first 2500 of which were discarded as warmup chains. Rhat values were generally within an acceptable range, with values <=1.02 (see appendix for diagnostic plots). We used uninformative priors for the fixed effects of the model (condition and velocity band), and weakly informative Student T distributions for for the random effects.

We compared varied and constant performance across two measures, deviation and discrimination. Deviation was quantified as the absolute deviation from the nearest boundary of the velocity band, or set to 0 if the throw velocity fell anywhere inside the target band. Thus, when the target band was 600-800, throws of 400, 650, and 1100 would result in deviation values of 200, 0, and 300, respectively. Discrimination was measured by fitting a linear model to the testing throws of each subjects, with the lower end of the target velocity band as the predicted variable, and the x velocity produced by the participants as the predictor variable. Participants who reliably discriminated between velocity bands tended to have positive slopes with values ~1, while participants who made throws irrespective of the current target band would have slopes ~0.

Code

```
# Create the data frame for the table
table_data <- data.frame(
  Type = c(
    rep("Population-Level Effects", 4),
    rep("Group-Level Effects", 2),
    "Family Specific Parameters"
  ),
  Parameter = c(
    "\\(\\beta_0\\)", "\\(\\beta_1\\)", "\\(\\beta_2\\)", "\\(\\beta_3\\)",
    "\\(\\sigma_{\\text{Intercept}}\\)", "\\(\\sigma_{\\text{bandInt}}\\)", "\\(\\sigma_{\\text{Observation}}\\)"
  ),
  Term = c(
    "(Intercept)", "conditVaried", "bandInt", "conditVaried:bandInt",
    "sd__(Intercept)", "sd__bandInt", "sd__Observation"
  ),
  Description = c(
    "Intercept representing the baseline deviation", "Effect of condition (Varied vs. Constant) on deviation", "Effect of target velocity band (bandInt) on deviation", "Interaction effect between training condition and target velocity band on deviation",
    "Standard deviation for (Intercept)", "Standard deviation for bandInt", "Standard deviation for Gaussian Family"
  )
) |>   mutate(
    Term = glue::glue("<code>{Term}</code>")
  ) 

# Create the table
kable_out <- table_data %>%
  kbl(format = 'html', escape = FALSE, booktabs = TRUE, 
      #caption = '<span style = "color:black;"><center><strong>Table 1: General Model Structure Information</strong></center></span>',
      col.names = c("Type", "Parameter", "Term", "Description")) %>%
  kable_styling(position="left", bootstrap_options = c("hover"), full_width = FALSE) %>%
  column_spec(1, bold = FALSE, border_right = TRUE) %>%
  column_spec(2, width = '4cm') %>%
  column_spec(3, width = '4cm') %>%
  row_spec(c(4, 7), extra_css = "border-bottom: 2px solid black;") %>%
  pack_rows("", 1, 4, bold = FALSE, italic = TRUE) %>%
  pack_rows("", 5, 6, bold = FALSE, italic = TRUE) %>%
  pack_rows("", 7, 7, bold = FALSE, italic = TRUE)

kable_out
```

Table 1: Mixed model structure and coefficient descriptions

| Type | Parameter | Term | Description |
| --- | --- | --- | --- |
|  |  |  |  |
| Population-Level Effects | \\(\\beta\_0\\) | `(Intercept)` | Intercept representing the baseline deviation |
| Population-Level Effects | \\(\\beta\_1\\) | `conditVaried` | Effect of condition (Varied vs. Constant) on deviation |
| Population-Level Effects | \\(\\beta\_2\\) | `bandInt` | Effect of target velocity band (bandInt) on deviation |
| Population-Level Effects | \\(\\beta\_3\\) | `conditVaried:bandInt` | Interaction effect between training condition and target velocity band on deviation |
|  |  |  |  |
| Group-Level Effects | \\(\\sigma\_{\\text{Intercept}}\\) | `sd__(Intercept)` | Standard deviation for (Intercept) |
| Group-Level Effects | \\(\\sigma\_{\\text{bandInt}}\\) | `sd__bandInt` | Standard deviation for bandInt |
|  |  |  |  |
| Family Specific Parameters | \\(\\sigma\_{\\text{Observation}}\\) | `sd__Observation` | Standard deviation for Gaussian Family |

## Results

### Testing Phase - No feedback.

In the first part of the testing phase, participants are tested from each of the velocity bands, and receive no feedback after each throw.

#### Deviation From Target Band

Descriptive summaries testing deviation data are provided in [Table 2](http://localhost:4200/Analysis/e1_test.html#tbl-e1-test-nf-deviation) and [Figure 1](http://localhost:4200/Analysis/e1_test.html#fig-e1-test-dev). To model differences in accuracy between groups, we used Bayesian mixed effects regression models to the trial level data from the testing phase. The primary model predicted the absolute deviation from the target velocity band (dist) as a function of training condition (condit), target velocity band (band), and their interaction, with random intercepts and slopes for each participant (id).

\\\[\\begin{equation} dist\_{ij} = \\beta\_0 + \\beta\_1 \\cdot condit\_{ij} + \\beta\_2 \\cdot band\_{ij} + \\beta\_3 \\cdot condit\_{ij} \\cdot band\_{ij} + b\_{0i} + b\_{1i} \\cdot band\_{ij} + \\epsilon\_{ij} \\end{equation}\\\] Code

````
```{r}
#| label: tbl-e1-test-nf-deviation
#| tbl-cap: "Testing Deviation - Empirical Summary"
#| tbl-subcap: ["Constant Testing - Deviation", "Varied Testing - Deviation"]
#| layout-ncol: 2

result <- test_summary_table(test, "dist","Deviation", mfun = list(mean = mean, median = median, sd = sd))
result$constant
result$varied
```
````

Table 2: Testing Deviation - Empirical Summary

Summary of Deviation - Constant
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Extrapolation | 254 | 148 | 298 |
| 350-550 | Extrapolation | 191 | 110 | 229 |
| 600-800 | Extrapolation | 150 | 84 | 184 |
| 800-1000 | Trained | 184 | 106 | 242 |
| 1000-1200 | Extrapolation | 233 | 157 | 282 |
| 1200-1400 | Extrapolation | 287 | 214 | 290 |

Summary of Deviation - Varied
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Extrapolation | 386 | 233 | 426 |
| 350-550 | Extrapolation | 285 | 149 | 340 |
| 600-800 | Extrapolation | 234 | 144 | 270 |
| 800-1000 | Trained | 221 | 149 | 248 |
| 1000-1200 | Trained | 208 | 142 | 226 |
| 1200-1400 | Trained | 242 | 182 | 235 |

Code

```
test |>  ggplot(aes(x = vb, y = dist,fill=condit)) +
    stat_summary(geom = "bar", position=position_dodge(), fun = mean) +
    stat_summary(geom = "errorbar", position=position_dodge(.9), fun.data = mean_se, width = .4, alpha = .7) + 
  labs(x="Band", y="Deviation From Target")
```

![](md_extract/assets/fig-e1-test-dev-1.png)

Figure 1: E1. Deviations from target band during testing without feedback stage.

Code

```
#contrasts(test$condit) 
#contrasts(test$vb)
modelName <- "e1_testDistBand_RF_5K"
e1_distBMM <- brm(dist ~ condit * bandInt + (1 + bandInt|id),
                      data=test,file=paste0(here::here("data/model_cache",modelName)),
                      iter=5000,chains=4)
GetModelStats(e1_distBMM) |> kable(escape=F,booktabs=T,caption="Model Coefficients")
```

Table 3: Experiment 1. Bayesian Mixed Model predicting absolute deviation as a function of condition (Constant vs. Varied) and Velocity Band

Model Coefficients
| Term | Estimate | 95% CrI Lower | 95% CrI Upper | pd |
| --- | --- | --- | --- | --- |
| Intercept | 205.09 | 136.86 | 274.06 | 1.00 |
| conditVaried | 157.44 | 60.53 | 254.90 | 1.00 |
| Band | 0.01 | \-0.07 | 0.08 | 0.57 |
| condit\*Band | \-0.16 | \-0.26 | \-0.06 | 1.00 |

Code

```
e1_distBMM |> 
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
| Constant - Varied | 100 | \-141.49 | \-229.19 | \-53.83 | 1.00 |
| Constant - Varied | 350 | \-101.79 | \-165.62 | \-36.32 | 1.00 |
| Constant - Varied | 600 | \-62.02 | \-106.21 | \-14.77 | 1.00 |
| Constant - Varied | 800 | \-30.11 | \-65.08 | 6.98 | 0.94 |
| Constant - Varied | 1000 | 2.05 | \-33.46 | 38.41 | 0.54 |
| Constant - Varied | 1200 | 33.96 | \-11.94 | 81.01 | 0.92 |

Code

```
coef_details <- get_coef_details(e1_distBMM, "conditVaried")
```

The model predicting absolute deviation (dist) showed clear effects of both training condition and target velocity band (Table X). Overall, the varied training group showed a larger deviation relative to the constant training group (β = 157.44, 95% CI \[60.53, 254.9\]). Deviation also depended on target velocity band, with lower bands showing less deviation. See [Table 3](http://localhost:4200/Analysis/e1_test.html#tbl-e1-bmm-dist) for full model output.

Code

```
condEffects <- function(m){
  m |> ggplot(aes(x = bandInt, y = .value, color = condit, fill = condit)) + 
  stat_dist_pointinterval() + stat_halfeye(alpha=.2) +
  stat_lineribbon(alpha = .25, size = 1, .width = c(.95)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
  ylab("Predicted X Velocity") + xlab("Band")
}

e1_distBMM |> emmeans( ~condit + bandInt, 
                       at = list(bandInt = c(100, 350, 600, 800, 1000, 1200))) |>
  gather_emmeans_draws() |>
 condEffects()+
  scale_x_continuous(breaks = c(100, 350, 600, 800, 1000, 1200), 
                     labels = levels(test$vb), 
                     limits = c(0, 1400)) 
``` 

![](md_extract/assets/fig-e1-bmm-dist-1.png)

Figure 2: E1. Conditioinal Effect of Training Condition and Band. Ribbon indicated 95% Credible Intervals.

#### Discrimination between bands

In addition to accuracy/deviation, we also assessed the ability of participants to reliably discriminate between the velocity bands (i.e. responding differently when prompted for band 600-800 than when prompted for band 150-350). [Table 4](http://localhost:4200/Analysis/e1_test.html#tbl-e1-test-nf-vx) shows descriptive statistics of this measure, and Figure 1 visualizes the full distributions of throws for each combination of condition and velocity band. To quantify discrimination, we again fit Bayesian Mixed Models as above, but this time the dependent variable was the raw x velocity generated by participants on each testing trial.

\\\[\\begin{equation} vx\_{ij} = \\beta\_0 + \\beta\_1 \\cdot condit\_{ij} + \\beta\_2 \\cdot bandInt\_{ij} + \\beta\_3 \\cdot condit\_{ij} \\cdot bandInt\_{ij} + b\_{0i} + b\_{1i} \\cdot bandInt\_{ij} + \\epsilon\_{ij} \\end{equation}\\\]

Code

```
test %>% group_by(id,vb,condit) |> plot_distByCondit()
```

![](md_extract/assets/fig-e1-test-vx-1-3.png)

Figure 3: E1 testing x velocities. Translucent bands with dash lines indicate the correct range for each velocity band.

Code

````
```{r}
#| label: tbl-e1-test-nf-vx
#| tbl-cap: "Testing vx - Empirical Summary"
#| tbl-subcap: ["Constant Testing - vx", "Varied Testing - vx"]
#| layout-ncol: 2

result <- test_summary_table(test, "vx","X Velocity", mfun = list(mean = mean, median = median, sd = sd))
result$constant 
result$varied 
```
````

Table 4: Testing vx - Empirical Summary

Summary of X Velocity - Constant
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Extrapolation | 524 | 448 | 327 |
| 350-550 | Extrapolation | 659 | 624 | 303 |
| 600-800 | Extrapolation | 770 | 724 | 300 |
| 800-1000 | Trained | 1001 | 940 | 357 |
| 1000-1200 | Extrapolation | 1167 | 1104 | 430 |
| 1200-1400 | Extrapolation | 1283 | 1225 | 483 |

Summary of X Velocity - Varied
| Band | Band Type | Mean | Median | Sd |
| --- | --- | --- | --- | --- |
| 100-300 | Extrapolation | 664 | 533 | 448 |
| 350-550 | Extrapolation | 768 | 677 | 402 |
| 600-800 | Extrapolation | 876 | 813 | 390 |
| 800-1000 | Trained | 1064 | 1029 | 370 |
| 1000-1200 | Trained | 1180 | 1179 | 372 |
| 1200-1400 | Trained | 1265 | 1249 | 412 |

Code

````
```{r}
#| label: tbl-e1-bmm-vx
#| tbl-cap: "Experiment 1. Bayesian Mixed Model Predicting Vx as a function of condition (Constant vs. Varied) and Velocity Band"
#| tbl-subcap: ["Model fit to all 6 bands", "Model fit to 3 extrapolation bands"]
#| layout-ncol: 2
e1_vxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                        data=test,file=paste0(here::here("data/model_cache", "e1_testVxBand_RF_5k")),
                        iter=5000,chains=4,silent=0,
                        control=list(adapt_delta=0.94, max_treedepth=13))
GetModelStats(e1_vxBMM ) |> kable(escape=F,booktabs=T, caption="Fit to all 6 bands")

cd1 <- get_coef_details(e1_vxBMM, "conditVaried")
sc1 <- get_coef_details(e1_vxBMM, "bandInt")
intCoef1 <- get_coef_details(e1_vxBMM, "conditVaried:bandInt")


modelName <- "e1_extrap_testVxBand"
e1_extrap_VxBMM <- brm(vx ~ condit * bandInt + (1 + bandInt|id),
                  data=test |>
                    filter(expMode=="test-Nf"),file=paste0(here::here("data/model_cache",modelName)),
                  iter=5000,chains=4)
GetModelStats(e1_extrap_VxBMM ) |> kable(escape=F,booktabs=T, caption="Fit to 3 extrapolation bands")


sc2 <- get_coef_details(e1_extrap_VxBMM, "bandInt")
intCoef2 <- get_coef_details(e1_extrap_VxBMM, "conditVaried:bandInt")
```
````

Table 5: Experiment 1. Bayesian Mixed Model Predicting Vx as a function of condition (Constant vs. Varied) and Velocity Band

Fit to all 6 bands
| Term | Estimate | 95% CrI Lower | 95% CrI Upper | pd |
| --- | --- | --- | --- | --- |
| Intercept | 408.55 | 327.00 | 490.61 | 1.00 |
| conditVaried | 164.05 | 45.50 | 278.85 | 1.00 |
| Band | 0.71 | 0.62 | 0.80 | 1.00 |
| condit\*Band | \-0.14 | \-0.26 | \-0.01 | 0.98 |

Fit to 3 extrapolation bands
| Term | Estimate | 95% CrI Lower | 95% CrI Upper | pd |
| --- | --- | --- | --- | --- |
| Intercept | 478.47 | 404.00 | 551.45 | 1.00 |
| conditVaried | 142.04 | 37.17 | 247.59 | 1.00 |
| Band | 0.50 | 0.42 | 0.57 | 1.00 |
| condit\*Band | \-0.07 | \-0.17 | 0.04 | 0.89 |

See [Table 5](http://localhost:4200/Analysis/e1_test.html#tbl-e1-bmm-vx) for the full model results. The estimated coefficient for training condition (β = 164.05, 95% CrI \[45.5, 278.85\]) suggests that the varied group tends to produce harder throws than the constant group, but is not in and of itself useful for assessing discrimination. Most relevant to the issue of discrimination is the slope on Velocity Band (β = 0.71, 95% CrI \[0.62, 0.8\]). Although the median slope does fall underneath the ideal of value of 1, the fact that the 95% credible interval does not contain 0 provides strong evidence that participants exhibited some discrimination between bands. The estimate for the interaction between slope and condition (β = -0.14, 95% CrI \[-0.26, -0.01\]), suggests that the discrimination was somewhat modulated by training condition, with the varied participants showing less senitivity between vands than the constant condition. This difference is depicted visually in [Figure 4](http://localhost:4200/Analysis/e1_test.html#fig-e1-bmm-vx).@tbl-e1-slope-quartile shows the average slope coefficients for varied and constant participants separately for each quartile. The constant participant participants appear to have larger slopes across quartiles, but the difference between conditions may be less pronounced for the top quartiles of subjects who show the strongest discrimination. Figure [Figure 5](http://localhost:4200/Analysis/e1_test.html#fig-e1-bmm-bx2) shows the distributions of slope values for each participant, and the compares the probability density of slope coefficients between training conditions. [Figure 6](http://localhost:4200/Analysis/e1_test.html#fig-e1-indv-slopes)

The second model, which focused solely on extrapolation bands, revealed similar patterns. The Velocity Band term (β = 0.5, 95% CrI \[0.42, 0.57\]) still demonstrates a high degree of discrimination ability. However, the posterior distribution for interaction term (β = -0.07, 95% CrI \[-0.17, 0.04\] ) does across over 0, suggesting that the evidence for decreased discrimination ability for the varied participants is not as strong when considering only the three extrapolation bands.

Code

```
e1_vxBMM |> emmeans( ~condit + bandInt, 
                       at = list(bandInt = c(100, 350, 600, 800, 1000, 1200))) |>
  gather_emmeans_draws() |> 
  condEffects() +
  scale_x_continuous(breaks = c(100, 350, 600, 800, 1000, 1200), 
                     labels = levels(test$vb), 
                     limits = c(0, 1400))

e1_extrap_VxBMM |> emmeans( ~condit + bandInt, 
                       at = list(bandInt = c(100, 350, 600))) |>
  gather_emmeans_draws() |>
  condEffects() +
  scale_x_continuous(breaks = c(100, 350, 600), 
                     labels = levels(test$vb)[1:3], 
                     limits = c(0, 1000)) 
``` 

![](md_extract/assets/fig-e1-bmm-vx-1.png)

(a) Model fit to all 6 bands

![](md_extract/assets/fig-e1-bmm-vx-2.png)

(b) Model fit to only 3 extrapolation bands

Figure 4: Conditional effect of training condition and Band. Ribbons indicate 95% HDI. The steepness of the lines serves as an indicator of how well participants discriminated between velocity bands.

Code

```
new_data_grid=map_dfr(1, ~data.frame(unique(test[,c("id","condit","bandInt")]))) |> 
  dplyr::arrange(id,bandInt) |> 
  mutate(condit_dummy = ifelse(condit == "Varied", 1, 0)) 

indv_coefs <- as_tibble(coef(e1_vxBMM)$id, rownames="id")|> 
  select(id, starts_with("Est")) |>
  left_join(e1Sbjs, by=join_by(id) ) 


fixed_effects <- e1_vxBMM |> 
  spread_draws(`^b_.*`,regex=TRUE) |> arrange(.chain,.draw,.iteration)


random_effects <- e1_vxBMM |> 
  gather_draws(`^r_id.*[[ReadItLater]] [[Article]]

# [HTW Project - HTW E1 Testing](http://localhost:4200/Analysis/e1_test.html)

, regex = TRUE, ndraws = 1500) |> 
  separate(.variable, into = c("effect", "id", "term"), sep = "\\[|,|\\]") |> 
  mutate(id = factor(id,levels=levels(test$id))) |> 
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
  left_join(e1Sbjs, by=join_by(id)) |> group_by(condit) |>
    select(id,condit,Intercept,b_Intercept,starts_with("Slope"),b_bandInt, n) |>
  mutate(rankSlope=rank(Slope)) |> arrange(rankSlope)   |> ungroup()
 
  
  indvSlopes |> mutate(Condition=condit) |>  group_by(Condition) |> 
    reframe(enframe(quantile(SlopeInt, c(0.0,0.25, 0.5, 0.75,1)), "quantile", "SlopeInt")) |> 
  pivot_wider(names_from=quantile,values_from=SlopeInt,names_prefix="Q_") |>
  group_by(Condition) |>
  summarise(across(starts_with("Q"), list(mean = mean))) |> kbl()
```

Table 6: Slope coefficients by quartile, per condition

| Condition | Q\_0%\_mean | Q\_25%\_mean | Q\_50%\_mean | Q\_75%\_mean | Q\_100%\_mean |
| --- | --- | --- | --- | --- | --- |
| Constant | \-0.1025710 | 0.4796476 | 0.6891841 | 0.9284214 | 1.394802 |
| Varied | \-0.2015313 | 0.2673139 | 0.5888388 | 0.8980842 | 1.295650 |

[Figure 5](http://localhost:4200/Analysis/e1_test.html#fig-e1-bmm-bx2) shows the distributions of estimated slopes relating velocity band to x velocity for each participant, ordered from lowest to highest within condition. Slope values are lower overall for varied training compared to constant training. Figure Xb plots the density of these slopes for each condition. The distribution for varied training has more mass at lower values than the constant training distribution. Both figures illustrate the model’s estimate that varied training resulted in less discrimination between velocity bands, evidenced by lower slopes on average.

Code

```
  indvSlopes |> ggplot(aes(y=rankSlope, x=SlopeInt,fill=condit,color=condit)) + 
  geom_pointrange(aes(xmin=SlopeInt.lower , xmax=SlopeInt.upper)) + 
  labs(x="Estimated Slope", y="Participant")  + facet_wrap(~condit)

   ggplot(indvSlopes, aes(x = SlopeInt, color = condit)) + 
  geom_density() + labs(x="Slope Coefficient",y="Density")
```

![](md_extract/assets/fig-e1-bmm-bx2-1-3.png)

(a) Slope estimates by participant - ordered from lowest to highest within each condition.

![](md_extract/assets/fig-e1-bmm-bx2-2-2.png)

(b) Destiny of slope coefficients by training group

Figure 5: Slope distributions between condition

Code

```
nSbj <- 3
indvDraws  |> indv_model_plot(indvSlopes, testAvg, SlopeInt,rank_variable=Slope,n_sbj=nSbj,"max")
indvDraws |> indv_model_plot(indvSlopes, testAvg,SlopeInt, rank_variable=Slope,n_sbj=nSbj,"min")
```

![](md_extract/assets/fig-e1-indv-slopes-1-3.png)

(a) subset with largest slopes

![](md_extract/assets/fig-e1-indv-slopes-2-3.png)

(b) subset with smallest slopes

Figure 6: Subset of Varied and Constant Participants with the smallest and largest estimated slope values. Red lines represent the best fitting line for each participant, gray lines are 200 random samples from the posterior distribution. Colored points and intervals at each band represent the empirical median and 95% HDI.

#### E1 Results Discussion

NEEDS TO BE WRITTEN


