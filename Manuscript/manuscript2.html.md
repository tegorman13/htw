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


::: {#tbl-example .cell tbl-cap='Example' tbl-subcap='["Cars1","Pressure1"]'}

````{.cell-code}
```{{r}}
#| label: tbl-example
#| tbl-cap: "Example"
#| tbl-subcap:
#|   - "Cars1"
#|   - "Pressure1"
library(knitr)
kable(head(cars))
kable(head(pressure))

pacman::p_load(tidyverse,tidybayes,brms, lme4, bayesplot,bayestestR,parameters,marginaleffects,
                emmeans, equatiomatic, here, pacman,  broom,
               broom.mixed,lme4,emmeans,here,knitr)





walk(c(here("Functions/Display_Functions.R"), here("Functions/org_functions.R"), 
       here("Functions/Table_Functions.R")), source)
```
````

::: {.cell-output-display}


| speed| dist|
|-----:|----:|
|     4|    2|
|     4|   10|
|     7|    4|
|     7|   22|
|     8|   16|
|     9|   10|



| temperature| pressure|
|-----------:|--------:|
|           0|   0.0002|
|          20|   0.0012|
|          40|   0.0060|
|          60|   0.0300|
|          80|   0.0900|
|         100|   0.2700|


:::
:::








# Methods

## Participants

Data was collected from 647 participants (after exclusions). The results shown below consider data from subjects in our initial experiment, which consisted of 196 participants (106 constant, 90 varied). The follow-up experiments entailed minor manipulations: 1) reversing the velocity bands that were trained on vs. novel during testing; 2) providing ordinal rather than numerical feedback during training (e.g. correct, too low, too high). The data from these subsequent experiments are largely consistently with our initial results shown below.


## Procedure





::: {#tbl-examplelllll .cell tbl-cap='Example' tbl-subcap='["Cars1l","Pressure1l"]'}

````{.cell-code}
```{{r}}
#| label: tbl-examplelllll
#| tbl-cap: "Example"
#| tbl-subcap:
#|   - "Cars1l"
#|   - "Pressure1l"
library(knitr)
kable(head(cars))
kable(head(pressure))
```
````

::: {.cell-output-display}


| speed| dist|
|-----:|----:|
|     4|    2|
|     4|   10|
|     7|    4|
|     7|   22|
|     8|   16|
|     9|   10|



| temperature| pressure|
|-----------:|--------:|
|           0|   0.0002|
|          20|   0.0012|
|          40|   0.0060|
|          60|   0.0300|
|          80|   0.0900|
|         100|   0.2700|


:::
:::





## extra




::: {.cell}

```{.r .cell-code}
pacman::p_load(tidyverse,knitr,tidybayes,brms, lme4, bayesplot,bayestestR,parameters,marginaleffects,
                emmeans, equatiomatic, here, pacman,  broom,
               broom.mixed,lme4,emmeans,here,
                wesanderson,glue, ggdist,ggforce,patchwork,gghalves,ggh4x,
                install = TRUE,
                update = FALSE
               )
```
:::

::: {#tbl-example2vvv .cell tbl-cap='Example' tbl-subcap='["Cars2vv","Pressure2vv"]'}

````{.cell-code}
```{{r}}
#| label: tbl-example2vvv
#| tbl-cap: "Example"
#| tbl-subcap:
#|   - "Cars2vv"
#|   - "Pressure2vv"
library(knitr)
kable(head(cars))
kable(head(pressure))
```
````

::: {.cell-output-display}


| speed| dist|
|-----:|----:|
|     4|    2|
|     4|   10|
|     7|    4|
|     7|   22|
|     8|   16|
|     9|   10|



| temperature| pressure|
|-----------:|--------:|
|           0|   0.0002|
|          20|   0.0012|
|          40|   0.0060|
|          60|   0.0300|
|          80|   0.0900|
|         100|   0.2700|


:::
:::





## Analyses Strategy

All data processing and statistical analyses were performed in R version 4.31 @rcoreteamLanguageEnvironmentStatistical2020. To assess differences between groups, we used Bayesian Mixed Effects Regression. Model fitting was performed with the brms package in R @burknerBrmsPackageBayesian2017, and descriptive stats and tables were extracted with the BayestestR package @makowskiBayestestRDescribingEffects2019a. Mixed effects regression enables us to take advantage of partial pooling, simultaneously estimating parameters at the individual and group level. Our use of Bayesian, rather than frequentist methods allows us to directly quantify the uncertainty in our parameter estimates, as well as circumventing convergence issues common to the frequentist analogues of our mixed models. For each model, we report the median values of the posterior distribution, and 95% credible intervals.

## load data



::: {.cell}

```{.r .cell-code}
library(here)
library(dplyr)
here::set_here(path='..')
```

::: {.cell-output .cell-output-stderr}

```
File .here already exists in /Users/thomasgorman/Library/CloudStorage/GoogleDrive-tegorman13@gmail.com/My Drive/HTW/gl
```


:::

```{.r .cell-code}
test <- readRDS(here("data/e1_08-21-23.rds")) |> filter(expMode2 == "Test") 
#options(brms.backend="cmdstanr",mc.cores=4)
e1Sbjs <- test |> group_by(id,condit) |> summarise(n=n())
testAvg <- test %>% group_by(id, condit, vb, bandInt,bandType,tOrder) %>%
  summarise(nHits=sum(dist==0),vx=mean(vx),dist=mean(dist),sdist=mean(sdist),n=n(),Percent_Hit=nHits/n)
```
:::

::: {#tbl-htw .cell tbl-cap='Example' tbl-subcap='["Cars2","Pressure2"]'}

````{.cell-code}
```{{r}}
#| label: tbl-htw
#| tbl-cap: "Example"
#| tbl-subcap:
#|   - "Cars2"
#|   - "Pressure2"
library(knitr)
kable(head(e1Sbjs))
kable(head(testAvg))
```
````

::: {.cell-output-display}


|id |condit   |  n|
|:--|:--------|--:|
|1  |Varied   | 63|
|2  |Varied   | 55|
|3  |Constant | 58|
|4  |Varied   | 63|
|5  |Constant | 63|
|6  |Varied   | 63|



|id |condit |vb        | bandInt|bandType      |tOrder     | nHits|        vx|     dist|     sdist|  n| Percent_Hit|
|:--|:------|:---------|-------:|:-------------|:----------|-----:|---------:|--------:|---------:|--:|-----------:|
|1  |Varied |100-300   |     100|Extrapolation |trainFirst |     0|  564.8623| 264.8623| 264.86235| 15|   0.0000000|
|1  |Varied |350-550   |     350|Extrapolation |trainFirst |     1|  791.2064| 243.6982| 243.69821| 15|   0.0666667|
|1  |Varied |600-800   |     600|Extrapolation |trainFirst |     2|  984.5210| 233.5615| 209.09513| 15|   0.1333333|
|1  |Varied |800-1000  |     800|Trained       |trainFirst |     2| 1060.7720| 155.5877| 102.60470|  6|   0.3333333|
|1  |Varied |1000-1200 |    1000|Trained       |trainFirst |     1| 1187.2458| 296.8225| 118.81533|  6|   0.1666667|
|1  |Varied |1200-1400 |    1200|Trained       |trainFirst |     2| 1378.7755| 136.2681|  83.08371|  6|   0.3333333|


:::
:::











## source functions



## Test tables




::: {#tbl-e1-test-nf-deviation .cell tbl-cap='Testing Deviation - Empirical Summary' tbl-subcap='["Full datasets","Intersection of samples with all labels available"]'}

```{.r .cell-code}
result <- test_summary_table(test, "dist","Deviation", mfun = list(mean = mean, median = median, sd = sd))


result$constant |>kable(booktabs = TRUE,
                        linesep = "\\addlinespace[0.5em]")
```

::: {.cell-output-display}


|Band      |Band Type     | Mean| Median|  Sd|
|:---------|:-------------|----:|------:|---:|
|100-300   |Extrapolation |  254|    148| 298|
|350-550   |Extrapolation |  191|    110| 229|
|600-800   |Extrapolation |  150|     84| 184|
|800-1000  |Trained       |  184|    106| 242|
|1000-1200 |Extrapolation |  233|    157| 282|
|1200-1400 |Extrapolation |  287|    214| 290|


:::

```{.r .cell-code}
                        #caption = paste("Summary of Deviation- Constant"))
# |>
#   kable_styling(font_size = ifelse(fmt_out == "latex", 8.5, NA))

result$varied |> kable(booktabs = TRUE,
                        linesep = "\\addlinespace[0.5em]")
```

::: {.cell-output-display}


|Band      |Band Type     | Mean| Median|  Sd|
|:---------|:-------------|----:|------:|---:|
|100-300   |Extrapolation |  386|    233| 426|
|350-550   |Extrapolation |  285|    149| 340|
|600-800   |Extrapolation |  234|    144| 270|
|800-1000  |Trained       |  221|    149| 248|
|1000-1200 |Trained       |  208|    142| 226|
|1200-1400 |Trained       |  242|    182| 235|


:::

```{.r .cell-code}
                        #caption = paste("Summary of Deviation- Varied"))
```
:::

::: {.cell}

```{.r .cell-code}
library(gt)
```
:::

::: {#tbl-example2bbb .cell tbl-cap='Example' tbl-subcap='["Cars2bb","Pressure2bb"]'}

````{.cell-code}
```{{r}}
#| label: tbl-example2bbb
#| tbl-cap: "Example"
#| tbl-subcap:
#|   - "Cars2bb"
#|   - "Pressure2bb"
library(knitr)
kable(head(cars))
kable(head(pressure))
```
````

::: {.cell-output-display}


| speed| dist|
|-----:|----:|
|     4|    2|
|     4|   10|
|     7|    4|
|     7|   22|
|     8|   16|
|     9|   10|



| temperature| pressure|
|-----------:|--------:|
|           0|   0.0002|
|          20|   0.0012|
|          40|   0.0060|
|          60|   0.0300|
|          80|   0.0900|
|         100|   0.2700|


:::
:::
