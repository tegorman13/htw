#
#
#
#
#

pacman::p_load(tidyverse,here,knitr,kableExtra,reactable)
select <- dplyr::select; mutate <- dplyr::mutate 
options(dplyr.summarise.inform=FALSE)
d <- readRDS(here("data/dPrune-07-27-23.rds")) %>% ungroup()


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
d<- d %>% relocate(bandOrder,bandType,.after=bandInt)


d <- d |> mutate(Exp=case_when(
  fb=="Continuous" & bandOrder=="Original" ~ "E1",
  fb=="Continuous" & bandOrder=="Reverse" ~ "E2",
  fb=="Ordinal" ~ "E3",
  TRUE ~ NA_character_
))

alt <- d |> filter(expMode2=="Test") |> group_by(id,condit,Exp) |> summarise(n=n_distinct(bandInt)) |> 
  filter(n<=4) |> pull(id) |> droplevels()
d <- d |> filter(!(id %in% alt))

d <- d |> mutate(vxC = ifelse(vx>1800,1800,vx))
d <- d %>% relocate(vxC,.after=vx)
d <- d |> select(-goodThrow, -trainVec, -fullCond)

# test <- d |> filter(expMode %in% c("test-Nf","test-train-nf"))
# rtest <- test |> filter(bandOrder=="Reverse")
# rtest |> group_by(vb,bandType,condit,tOrder) |> summarise(n=n()) 



#
#
#
#| eval: false
e1 <- d |> filter(fb=="Continuous" & bandOrder=="Original") |> mutate(id=factor(id,levels=unique(id)))
e2 <- d |> filter(fb=="Continuous" & bandOrder=="Reverse") |> mutate(id=factor(id,levels=unique(id)))
e3 <- d |> filter(fb=="Ordinal") |> mutate(id=factor(id,levels=unique(id)))

date.append="08-21-23"

# saveRDS(d,here(paste0("data/dAll_",date.append,".rds")))
# saveRDS(e1,here(paste0("data/e1_",date.append,".rds")))
# saveRDS(e2,here(paste0("data/e2_",date.append,".rds")))
# saveRDS(e3,here(paste0("data/e3_",date.append,".rds")))

# save csv versions
# d %>% write_csv(here(paste0("data/dAll_",date.append,".csv")))
# e1 %>% write_csv(here(paste0("data/e1_",date.append,".csv")))
# e2 %>% write_csv(here(paste0("data/e2_",date.append,".csv")))
# e3 %>% write_csv(here(paste0("data/e3_",date.append,".csv")))

#
#
#
#
#
#
#
#
#
#



d <- readRDS(here::here("data/e1_08-21-23.rds"))
levels(d$condit)
dtest <- d %>% filter(expMode %in% c("test-Nf","test-train-nf")) %>% group_by(id,lowBound) %>% 
  mutate(nBand=n(),band=bandInt,id=factor(id)) %>% group_by(id) %>% mutate(nd=n_distinct(lowBound))
dtest <- dtest %>% group_by(id,lowBound) %>% filter(nBand>=5 & nd==6)
dtest <- dtest %>% group_by(id) %>% filter(!id %in% unique(dtest$id[dtest$nBand<5]))

ds <- d %>% filter(expMode2 %in% c("Train","Test")) |> 
  filter(!id %in% unique(dtest$id[dtest$nBand<5])) |>   
  group_by(id,condit,expMode2) |> 
  mutate(input=bandInt,x=bandInt, y=vx, tr= row_number()) |>
  select(id,condit,expMode2,tr,x,y) 

#saveRDS(ds,here::here("data/e1_md_11-06-23.rds"))

#
#
#
#
#
#
#
#

d %>% select_if(is.factor) %>% select(-sbjCode,-id) %>% map(levels)
d %>% select_if(is.factor) %>% select(-sbjCode,-id) %>% droplevels %>% map(levels)

#
#
#
#
d %>% select_if(is.numeric) |> colnames()
#
#
#
#
#
d %>%
  distinct(id, condit, fb, bandOrder, tOrder) %>%
  group_by(condit, fb, bandOrder, tOrder) %>%
  summarise(n = n()) %>%
  kable()
#
#
#
#
#
#

# Average trials per subject by condition  
d %>%
  group_by(condit, fb, bandOrder, tOrder, id) %>%
  summarise(n = n()) %>%
  group_by(condit, fb, bandOrder, tOrder) %>%
  summarise(mean_trials = mean(n)) %>%
  kable()


d %>%
  group_by(condit, fb, bandOrder, tOrder, id, expMode) %>% 
  summarise(n = n()) %>%
  pivot_wider(names_from = expMode, values_from = n, names_prefix = "n_") %>%
  group_by(condit, fb, bandOrder, tOrder) %>%
  summarise(across(starts_with("n_"), ~mean(., na.rm = TRUE))) %>%
  kable()

#
#
#
#
#| column: page-right
#| 
# Define column defs function 
col_defs <- function(data) {
  
  cols <- colnames(data)
  
  defs <- lapply(cols, function(x) {
    
    if(is.factor(data[[x]])) {
      colDef(sortable = TRUE, filterable = TRUE,minWidth=108) 
    } else {
      colDef(sortable = TRUE, filterable = FALSE,minWidth=90)
    }
    
  })
  
  setNames(defs, cols)
}

d %>%
  group_by(id,condit, fb, bandOrder, tOrder, expMode) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = expMode, values_from = n, names_prefix = "n_") %>%
  
  # Pass original data too 
  reactable(columns = col_defs(.), 
            highlight = TRUE,
            defaultPageSize = 25)
#
#
#
#
#
#
d %>% filter(nGoodTrial==1,nTotal>100) %>% ggplot(aes(nTotal)) + geom_histogram() + facet_wrap(~condit)
d %>% filter(nGoodTrial==1,nTotal>100) %>% ggplot(aes(nTestNf)) + geom_histogram() + facet_wrap(~condit)
d %>% filter(nGoodTrial==1,nTotal>100) %>% ggplot(aes(nTrain)) + geom_histogram() + facet_wrap(~condit)

#
#
#
#
