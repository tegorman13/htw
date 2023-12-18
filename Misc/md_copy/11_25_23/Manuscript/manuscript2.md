---
title: "HTW"
short-title: "Variability and Extrapolation"
date: 10-11-2023
author:
- name: Thomas Gorman
  affiliation: Indiana University
  url: www.tegorman13.github.io
  email: tegorman@iu.edu
  orcid: 0000-0001-5366-5442
- name: Rob Goldstone
  affiliation: Indiana University
  url: www.pc.cogs.indiana.edu 
  email: rgoldsto@indiana.edu
  orcid: 0000-0001-8357-8358
abstract: |
  In project 1, we applied model-based techniques to quantify and control for the similarity between training and testing experience, which in turn enabled us to account for the difference between varied and constant training via an extended version of a similarity based generalization model. In project 2, we will go a step further, implementing a full process model capable of both 1) producing novel responses and 2) modeling behavior in both the learning and testing stages of the experiment. Project 2 also places a greater emphasis on extrapolation performance following training - as varied training has often been purported to be particularly beneficial in such situations. 
keywords:
  - Learning Generalization
  - Function Learning
  - Visuomotor learning
  - Training Variability
code-repo: "Access the code, data, and analysis at <https://github.com/tegorman13/htw>"
bibliography: ../Assets/Bib/Dissertation.bib
link-citations: true
---


```{r}
#| label: tbl-example
#| tbl-cap: "Example"
#| tbl-subcap: 
#|   - "Cars1"
#|   - "Pressure1"
#| echo: fenced
library(knitr)
kable(head(cars))
kable(head(pressure))

pacman::p_load(tidyverse,tidybayes,brms, lme4, bayesplot,bayestestR,parameters,marginaleffects,
                emmeans, equatiomatic, here, pacman,  broom,
               broom.mixed,lme4,emmeans,here,knitr)





walk(c(here("Functions/Display_Functions.R"), here("Functions/org_functions.R"), 
       here("Functions/Table_Functions.R")), source)
```





# Methods

## Participants

Data was collected from 647 participants (after exclusions). The results shown below consider data from subjects in our initial experiment, which consisted of 196 participants (106 constant, 90 varied). The follow-up experiments entailed minor manipulations: 1) reversing the velocity bands that were trained on vs. novel during testing; 2) providing ordinal rather than numerical feedback during training (e.g. correct, too low, too high). The data from these subsequent experiments are largely consistently with our initial results shown below.


## Procedure


```{r}
#| label: tbl-examplelllll
#| tbl-cap: "Example"
#| tbl-subcap: 
#|   - "Cars1l"
#|   - "Pressure1l"
#| echo: fenced
library(knitr)
kable(head(cars))
kable(head(pressure))
```


## extra

```{r}
pacman::p_load(tidyverse,knitr,tidybayes,brms, lme4, bayesplot,bayestestR,parameters,marginaleffects,
                emmeans, equatiomatic, here, pacman,  broom,
               broom.mixed,lme4,emmeans,here,
                wesanderson,glue, ggdist,ggforce,patchwork,gghalves,ggh4x,
                install = TRUE,
                update = FALSE
               )
```


```{r}
#| label: tbl-example2vvv
#| tbl-cap: "Example"
#| tbl-subcap: 
#|   - "Cars2vv"
#|   - "Pressure2vv"
#| echo: fenced
library(knitr)
kable(head(cars))
kable(head(pressure))
```


## Analyses Strategy

All data processing and statistical analyses were performed in R version 4.31 @rcoreteamLanguageEnvironmentStatistical2020. To assess differences between groups, we used Bayesian Mixed Effects Regression. Model fitting was performed with the brms package in R @burknerBrmsPackageBayesian2017, and descriptive stats and tables were extracted with the BayestestR package @makowskiBayestestRDescribingEffects2019a. Mixed effects regression enables us to take advantage of partial pooling, simultaneously estimating parameters at the individual and group level. Our use of Bayesian, rather than frequentist methods allows us to directly quantify the uncertainty in our parameter estimates, as well as circumventing convergence issues common to the frequentist analogues of our mixed models. For each model, we report the median values of the posterior distribution, and 95% credible intervals.

## load data
```{r}
library(here)
library(dplyr)
here::set_here(path='..')
test <- readRDS(here("data/e1_08-21-23.rds")) |> filter(expMode2 == "Test") 
#options(brms.backend="cmdstanr",mc.cores=4)
e1Sbjs <- test |> group_by(id,condit) |> summarise(n=n())
testAvg <- test %>% group_by(id, condit, vb, bandInt,bandType,tOrder) %>%
  summarise(nHits=sum(dist==0),vx=mean(vx),dist=mean(dist),sdist=mean(sdist),n=n(),Percent_Hit=nHits/n)

```

```{r}
#| label: tbl-htw
#| tbl-cap: "Example"
#| tbl-subcap: 
#|   - "Cars2"
#|   - "Pressure2"
#| echo: fenced
library(knitr)
kable(head(e1Sbjs))
kable(head(testAvg))
```








## source functions



## Test tables

```{r}
#| label: tbl-e1-test-nf-deviation
#| tbl-cap: "Testing Deviation - Empirical Summary"
#| tbl-subcap:
#|   - "Full datasets"
#|   - "Intersection of samples with all labels available"
#| eval: false


result <- test_summary_table(test, "dist","Deviation", mfun = list(mean = mean, median = median, sd = sd))


result$constant |>kable(booktabs = TRUE,
                        linesep = "\\addlinespace[0.5em]")
                        #caption = paste("Summary of Deviation- Constant"))
# |>
#   kable_styling(font_size = ifelse(fmt_out == "latex", 8.5, NA))

result$varied |> kable(booktabs = TRUE,
                        linesep = "\\addlinespace[0.5em]")
                        #caption = paste("Summary of Deviation- Varied"))


```


```{r}
# library(gt)
```

```{r}
#| label: tbl-example2bbb
#| tbl-cap: "Example"
#| tbl-subcap: 
#|   - "Cars2bb"
#|   - "Pressure2bb"
#| echo: fenced
# library(knitr)
# kable(head(cars))
# kable(head(pressure))
```

