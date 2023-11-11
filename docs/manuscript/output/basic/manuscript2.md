# HTW
Thomas GormanRob Goldstone
2023-10-11

```` markdown
```{r}
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

<div class="cell-output-display">

| speed | dist |
|------:|-----:|
|     4 |    2 |
|     4 |   10 |
|     7 |    4 |
|     7 |   22 |
|     8 |   16 |
|     9 |   10 |

| temperature | pressure |
|------------:|---------:|
|           0 |   0.0002 |
|          20 |   0.0012 |
|          40 |   0.0060 |
|          60 |   0.0300 |
|          80 |   0.0900 |
|         100 |   0.2700 |

</div>

# Methods

## Participants

Data was collected from 647 participants (after exclusions). The results
shown below consider data from subjects in our initial experiment, which
consisted of 196 participants (106 constant, 90 varied). The follow-up
experiments entailed minor manipulations: 1) reversing the velocity
bands that were trained on vs. novel during testing; 2) providing
ordinal rather than numerical feedback during training (e.g. correct,
too low, too high). The data from these subsequent experiments are
largely consistently with our initial results shown below.

## Procedure

```` markdown
```{r}
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

<div class="cell-output-display">

| speed | dist |
|------:|-----:|
|     4 |    2 |
|     4 |   10 |
|     7 |    4 |
|     7 |   22 |
|     8 |   16 |
|     9 |   10 |

| temperature | pressure |
|------------:|---------:|
|           0 |   0.0002 |
|          20 |   0.0012 |
|          40 |   0.0060 |
|          60 |   0.0300 |
|          80 |   0.0900 |
|         100 |   0.2700 |

</div>

## extra

``` r
pacman::p_load(tidyverse,knitr,tidybayes,brms, lme4, bayesplot,bayestestR,parameters,marginaleffects,
                emmeans, equatiomatic, here, pacman,  broom,
               broom.mixed,lme4,emmeans,here,
                wesanderson,glue, ggdist,ggforce,patchwork,gghalves,ggh4x,
                install = TRUE,
                update = FALSE
               )
```

```` markdown
```{r}
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

<div class="cell-output-display">

| speed | dist |
|------:|-----:|
|     4 |    2 |
|     4 |   10 |
|     7 |    4 |
|     7 |   22 |
|     8 |   16 |
|     9 |   10 |

| temperature | pressure |
|------------:|---------:|
|           0 |   0.0002 |
|          20 |   0.0012 |
|          40 |   0.0060 |
|          60 |   0.0300 |
|          80 |   0.0900 |
|         100 |   0.2700 |

</div>

## Analyses Strategy

All data processing and statistical analyses were performed in R version
4.31 Team ([2020](#ref-rcoreteamLanguageEnvironmentStatistical2020)). To
assess differences between groups, we used Bayesian Mixed Effects
Regression. Model fitting was performed with the brms package in R
Bürkner ([2017](#ref-burknerBrmsPackageBayesian2017)), and descriptive
stats and tables were extracted with the BayestestR package
([**makowskiBayestestRDescribingEffects2019a?**](#ref-makowskiBayestestRDescribingEffects2019a)).
Mixed effects regression enables us to take advantage of partial
pooling, simultaneously estimating parameters at the individual and
group level. Our use of Bayesian, rather than frequentist methods allows
us to directly quantify the uncertainty in our parameter estimates, as
well as circumventing convergence issues common to the frequentist
analogues of our mixed models. For each model, we report the median
values of the posterior distribution, and 95% credible intervals.

## load data

``` r
library(here)
library(dplyr)
here::set_here(path='..')
```

    File .here already exists in /Users/thomasgorman/Library/CloudStorage/GoogleDrive-tegorman13@gmail.com/My Drive/HTW/gl

``` r
test <- readRDS(here("data/e1_08-21-23.rds")) |> filter(expMode2 == "Test") 
#options(brms.backend="cmdstanr",mc.cores=4)
e1Sbjs <- test |> group_by(id,condit) |> summarise(n=n())
testAvg <- test %>% group_by(id, condit, vb, bandInt,bandType,tOrder) %>%
  summarise(nHits=sum(dist==0),vx=mean(vx),dist=mean(dist),sdist=mean(sdist),n=n(),Percent_Hit=nHits/n)
```

```` markdown
```{r}
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

<div class="cell-output-display">

| id  | condit   |   n |
|:----|:---------|----:|
| 1   | Varied   |  63 |
| 2   | Varied   |  55 |
| 3   | Constant |  58 |
| 4   | Varied   |  63 |
| 5   | Constant |  63 |
| 6   | Varied   |  63 |

| id  | condit | vb        | bandInt | bandType      | tOrder     | nHits |        vx |     dist |     sdist |   n | Percent_Hit |
|:----|:-------|:----------|--------:|:--------------|:-----------|------:|----------:|---------:|----------:|----:|------------:|
| 1   | Varied | 100-300   |     100 | Extrapolation | trainFirst |     0 |  564.8623 | 264.8623 | 264.86235 |  15 |   0.0000000 |
| 1   | Varied | 350-550   |     350 | Extrapolation | trainFirst |     1 |  791.2064 | 243.6982 | 243.69821 |  15 |   0.0666667 |
| 1   | Varied | 600-800   |     600 | Extrapolation | trainFirst |     2 |  984.5210 | 233.5615 | 209.09513 |  15 |   0.1333333 |
| 1   | Varied | 800-1000  |     800 | Trained       | trainFirst |     2 | 1060.7720 | 155.5877 | 102.60470 |   6 |   0.3333333 |
| 1   | Varied | 1000-1200 |    1000 | Trained       | trainFirst |     1 | 1187.2458 | 296.8225 | 118.81533 |   6 |   0.1666667 |
| 1   | Varied | 1200-1400 |    1200 | Trained       | trainFirst |     2 | 1378.7755 | 136.2681 |  83.08371 |   6 |   0.3333333 |

</div>

## source functions

## Test tables

``` r
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

``` r
# library(gt)
```

```` markdown
```{r}
#| label: tbl-example2bbb
#| tbl-cap: "Example"
#| tbl-subcap:
#|   - "Cars2bb"
#|   - "Pressure2bb"
# library(knitr)
# kable(head(cars))
# kable(head(pressure))
```
````

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0" line-spacing="2">

<div id="ref-burknerBrmsPackageBayesian2017" class="csl-entry">

Bürkner, P.-C. (2017). Brms: An R Package for Bayesian Multilevel Models
Using Stan. *Journal of Statistical Software*, *80*, 1–28.
<https://doi.org/10.18637/jss.v080.i01>

</div>

<div id="ref-rcoreteamLanguageEnvironmentStatistical2020"
class="csl-entry">

Team, R. C. (2020). *R: A Language and Environment for Statistical
Computing*. R: A Language and Environment for Statistical Computing.

</div>

</div>
