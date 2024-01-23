---
title: "`r params$title`"
author: Thomas Gorman
date: "`r Sys.Date()`"
format: 
  html:
    page-layout: full
code-fold: true
code-tools: true
execute: 
  warning: false
  eval: true
params:
    file_name: file_name
    title: title
   #ind_fits: ind_fits
   
---




```{r}
#| results: asis

pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted, stringr, knitr, gt, glue)
conflict_prefer_all("dplyr", quiet = TRUE)

#walk(c("Display_Functions"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
ids <- c(1,2,4,5,6,7,8, 10,11,12,13)
ids2 <- c(1,66,36)
ids3 <- c(20,71,101,4,76,192)
idsBad <- c(76,192, 101)


#source(here::here("Functions/fun_indv_fit.R"))

source("../../../Functions/fun_indv_fit.R")
source("../../../Functions/Display_Functions.R")

raw_name <- tools::file_path_sans_ext(basename(params$file_name))

# print(list.files("../../data"))
# print(list.files("../../../data"))
# print(getwd())

# cat(glue::glue('## {raw_name}'))

# print(params$file_name)

```


## `{r} raw_name`

#| eval: false

```{r}
#| eval: true

ind_fits <- readRDS(params$file_name)


ind_fits_df <- imap_dfr(ind_fits |> list_assign(runInfo = zap()), ~imap_dfr(.x, ~{
  teter_data <- .x[["teter_results"]] %>% mutate(result_type = "teter_results")
  te_data <- .x[["te_results"]] %>% mutate(result_type = "te_results")
  tr_data <- .x[["tr_results"]] %>% mutate(result_type = "tr_results")
  combined_data <- bind_rows(teter_data, te_data, tr_data)
}, .id = "id"))

ids2 <- c(1,66,36)
ids3 <- c(20,71,101,4,76,192)
idsBad <- c(76,192, 101)

#head(ind_fits$runInfo$systemInfo)

repair_names <- function(names) {
  dupes <- names[duplicated(names)]
  if ("rank" %in% dupes) {
    first_rank <- which(names == "rank")[1]
    names <- make.unique(names, sep = "_")
    names[first_rank] <- "rank"
  }
  names
}

# Single solution for both cases
post_dat <- ind_fits_df %>% 
  select(sim_dat, rank) %>% 
  unnest(sim_dat, names_repair = repair_names) |> as.data.table()


#post_dat <- ind_fits_df |> select(sim_dat,rank) |> unnest(sim_dat,names_repair = "universal") |> as.data.table()


post_dat_l <- setorder(melt(setcolorder(post_dat[expMode2 == "Test"][, `:=`(avg_y = mean(y)), 
        by = .(id, condit, x)], c("id", "condit", "Model", "Fit_Method", 
    "expMode2", "tr", "rank", "c", "lr", "x", "y", "avg_y", "pred", 
    "resid")), measure.vars = c("pred", "y"), variable.name = "Resp", 
        value.name = "val", variable.factor = FALSE)[, `:=`(Resp = fcase(Resp == 
        "y", "Observed", Model == "ALM", "ALM", Model == "EXAM", 
        "EXAM")), by = .(id, condit, x)], tr, id, Resp, na.last = TRUE)



```


##  Posterior Average Table: `{r} raw_name`
```{r fig.width=12, fig.height=17}
#| eval: true
post_tabs <- abc_tables(post_dat_l) 

#post_tabs$et_sum |> select(condit,Fit_Method,Avg_ALM_error, Avg_EXAM_error) |> filter(Fit_Method=="Test")

post_tabs$et_sum |> gt::gt()


```


## Posterior Predictive: `{r} raw_name`
```{r, fig.width=12, fig.height=14}
#| eval: true
group_predictive_plots(post_dat_l)
group_best_plots(post_dat_l)
```



## Individual Plots: `{r} raw_name`
```{r, fig.width=12, fig.height=17}
#| eval: true
indv_best_plots(post_dat_l)
indv_predictive_plots(post_dat_l, ids2)
indv_predictive_plots(post_dat_l, idsBad)

```


## Distributions: `{r} raw_name`
```{r, fig.width=12, fig.height=12}
#| eval: true
plot_sampled_posterior(ind_fits)
plot_indv_posterior(ind_fits_df)

```




## Bigger Tables: `{r} raw_name`
```{r}
#| eval: true

post_tabs$et_sum_x |> kable()
post_tabs$et_sum_x_indv |> kable()
```


```{r}


extract_info <- function(raw_names) {
  full_name <- raw_names
  n_samp <- str_extract(raw_names, "(?<=_)\\d+(?=_)")
  ng_value <- str_extract(raw_names, "(?<=ng)\\d+")
  buf_value <- str_replace(str_extract(raw_names, "(?<=buf)[0-9p]+"), "p", ".")
  run_type <- ifelse(str_detect(raw_names, "ss"), "ss", "trial")

  data.frame(full_name, n_samp, ng_value, buf_value, run_type)
}

#extract_info(raw_names)
post_tabs <- map(post_tabs, function(df) {
  # Assuming raw_names is a column in df, change to actual column name if different
  df_info <- extract_info(raw_name)

  # Mutating the dataframe with the new columns
  mutate(df, 
         n_samp = as.integer(df_info$n_samp),
         ng_value = as.integer(df_info$ng_value),
         buf_value = as.numeric(df_info$buf_value),
         run_type = df_info$run_type, 
         full_name = df_info$full_name)
})


path1 = "../../../data/abc_tabs"
saveRDS(post_tabs, paste0(path1, "/", raw_name,"_post_tab", ".rds"))
saveRDS(post_tabs, paste0(raw_name,"_post_tab", ".rds"))


```

