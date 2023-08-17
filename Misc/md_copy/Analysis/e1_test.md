---
title: "E1 Testing"
categories: [Analysis, R, Bayesian]
---

## Analyses Strategy

To assess differences between groups, we used Bayesian Mixed Effects Regression. Analyses were done using brms package in R @burknerBrmsPackageBayesian2017 as well as  @makowskiBayestestRDescribingEffects2019a. 

Unless otherwise stated, we ran each model with 4 chains, 5000 iterations per chain, the first 2500 of which were discarded as warmup chains. Rhat values were generally within an acceptable range, with values <=1.02. 


We made use of two separate performance metrics, deviation and discrimination. Deviation was quantified as the absolute deviation from the nearest boundary of the velocity band, or set to 0 if the throw velocity fell anywhere inside the target band. Thus, when the target band was 600-800, throws of 400, 650, and 1100 would result in deviation values of 200, 0, and 300, respectively. Discrimination was measured by fitting a linear model to the testing throws of each subjects, with the lower end of the target velocity band as the predicted variable, and the x velocity produced by the participants as the predictor variable. Participants who reliably discriminated between velocity bands tended to have positive slopes with values ~1, while participants who made throws irrespective of the current target band would have slopes ~0.  



## Results

### Testing Phase - No feedback. 

In the first part of the testing phase, participants are tested from each of the velocity bands, and receive no feedback after each throw. 

```{r setup, include=FALSE}
source(here::here("Functions", "packages.R"))
test <- readRDS(here("data/e1_08-04-23.rds")) |> 
  filter(expMode2 == "Test") 
options(brms.backend="cmdstanr",mc.cores=4)
```



#### Deviation From Target Band


To assess differences in accuracy between groups, we used Bayesian mixed effects regression models implemented via brms in R. The primary model predicted the absolute deviation from the target velocity band (dist) as a function of training condition (condit), target velocity band (vb), and their interaction, with random intercepts and slopes for each participant (id). 

```{r}
#contrasts(test$condit) 
#  Varied
# Constant      0
# Varied        1

# contrasts(test$vb)
#           350-550 600-800 800-1000 1000-1200 1200-1400
# 100-300         0       0        0         0         0
# 350-550         1       0        0         0         0
# 600-800         0       1        0         0         0
# 800-1000        0       0        1         0         0
# 1000-1200       0       0        0         1         0
# 1200-1400       0       0        0         0         1

modelName <- "e1_testDistRF2"
e1_testDistRF2 <- brm(dist ~ condit * vb + (1 + vb|id),
                      data=test,file=paste0(here::here("data/model_cache",modelName)),
                      iter=5000,chains=4)
mp1 <- model_parameters(e1_testDistRF2,effects="fixed",diagnostic=c("ESS"),verbose=FALSE)
```



The model predicting absolute deviation (dist) showed clear effects of both training condition and target velocity band (Table X). Overall, the varied training group showed a larger deviation relative to the constant training group (β = 119, 95% CI [37, 202]). Deviation also depended on target velocity band, with lower bands showing less deviation. 



#### Discrimination between Velocity Bands

In addition to accuracy/deviation. We also assessed the ability of participants to reliably discriminate between the velocity bands (i.e. responding differently when prompted for band 600-800 than when prompted for band 150-350). 

NEEDS TO BE WRITTEN


#### E1 Results Discussion



```{r}
#| label: Band Deviation
#| tbl-cap: "Deviations from target"
#| layout: [[70, 30],[-5]]



```






## Misc 

{{< include "../Misc/alm_table.qmd" >}}







Here is a footnote reference,[^1] and another.[^2]

[^1]: Here is the footnote.
[^2]: Here's one with multiple blocks. Subsequent paragraphs are
    indented to show that they belong to the previous footnote.
    ```         
    { some.code }
    ```
Here is an inline note.[^3]

[^3]: Inlines notes are easier to write, since you don't have to pick an
    identifier and move down to type the note.


## References

::: {#refs}
:::
