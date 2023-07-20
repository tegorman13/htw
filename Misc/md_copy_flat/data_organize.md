---
title: Organization
---


```{r}
pacman::p_load(tidyverse,lme4,emmeans,here,knitr,kableExtra)
options(dplyr.summarise.inform=FALSE)
d <- readRDS(here("data/dPrune-05-21-23.rds")) %>% ungroup()
d$stage <- factor(d$stage, ordered=TRUE, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23"))
d$tOrder <- factor(d$tOrder, levels=c("testFirst","trainFirst"),labels=c("Test First","Train First"))
d$bandOrder <- factor(d$catOrder,levels=c("orig","rev"),labels=c("Original","Reverse"))
d$fb <- factor(d$fb,levels=c("continuous","ordinal"),labels=c("Continuous","Ordinal"))

d <- d |> mutate(bandType=case_when(
  vb %in% c("1000-1200","1200-1400") & condit=="Constant" & bandOrder=="Original" ~ "Extrapolation",
  vb %in% c("100-300","350-550") & condit=="Constant" & bandOrder=="Reverse" ~ "Extrapolation",
  (vb %in% c("100-300","350-550","600-800") & bandOrder=="Original") |
  (vb %in% c("1200-1400","1000-1200","800-1000") & bandOrder=="Reverse") ~ "Extrapolation",
  
  (vb %in% c("800-1000","1000-1200","1200-1400") & bandOrder=="Original") |
  (vb %in% c("100-300","350-550","600-800") & bandOrder=="Reverse") ~ "Trained",
  TRUE ~ NA_character_
))
d$bandType <- factor(d$bandType,levels=c("Trained","Extrapolation"))

d<- d %>% mutate(id=as.factor(d$id),fullCond=as.factor(fullCond)) %>% 
  select(-catOrder,-mode) %>% relocate(bandOrder,bandType,.after=fb)

# test <- d |> filter(expMode %in% c("test-Nf","test-train-nf"))
# rtest <- test |> filter(bandOrder=="Reverse")
# rtest |> group_by(vb,bandType,condit,tOrder) |> summarise(n=n()) 

```

head(d)
colnames(d)
d %>% select_if(is.factor) %>% colnames()


```{r}
d %>% select_if(is.factor) %>% select(-sbjCode,-id) %>% map(levels)
```


```{r}
d %>% select_if(is.numeric) |> colnames()

```


```{r}
#| eval: false
e1 <- d |> filter(fb=="Continuous" & bandOrder=="Original")
e2 <- d |> filter(fb=="Continuous" & bandOrder=="Reverse")
e3 <- d |> filter(fb=="Ordinal")

saveRDS(d,here("data/dAll_07_04_23.rds"))
saveRDS(e1,here("data/e1_07_04_23.rds"))
saveRDS(e2,here("data/e2_07_04_23.rds"))
saveRDS(e3,here("data/e3_07_04_23.rds"))
# save csv versions
d %>% write_csv(here("data/dAll_07_04_23.csv"))
e1 %>% write_csv(here("data/e1_07_04_23.csv"))
e2 %>% write_csv(here("data/e2_07_04_23.csv"))
e3 %>% write_csv(here("data/e3_07_04_23.csv"))
```






