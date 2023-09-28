

#source(here::here("Functions", "packages.R"))
library(pacman)
pacman::p_load(tidyverse,tidybayes,brms, lme4, bayesplot,bayestestR,parameters,marginaleffects,
                emmeans, equatiomatic, here, pacman,  broom,
               broom.mixed,lme4,emmeans,here,knitr,kableExtra,gt,
                wesanderson,glue, ggdist,ggforce,patchwork,gghalves,ggh4x,
                pander,
                install = TRUE,
                update = FALSE
               )

test <- readRDS(here("data/e1_08-21-23.rds")) |> filter(expMode2 == "Test") 
options(brms.backend="cmdstanr",mc.cores=4)
e1Sbjs <- test |> group_by(id,condit) |> summarise(n=n())
testAvg <- test %>% group_by(id, condit, vb, bandInt,bandType,tOrder) %>%
  summarise(nHits=sum(dist==0),vx=mean(vx),dist=mean(dist),sdist=mean(sdist),n=n(),Percent_Hit=nHits/n)


result <- test_summary_table(test, "dist", mfun = list(mean = mean, median = median, sd = sd))
result$constant
result$varied