---
title: Individual ABC Fits
author: Thomas Gorman
date: "`r Sys.Date()`"
code-fold: true
code-tools: true
execute: 
  warning: false
  eval: true
---

```{r}
pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted, 
               stringr,future,furrr, knitr, reactable, flextable,ggstance)
conflict_prefer_all("dplyr", quiet = TRUE)
options(scipen = 999)
walk(c("Display_Functions","fun_alm","fun_indv_fit","fun_model"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
input_layer <<- output_layer <<-  c(100,350,600,800,1000,1200)

ids <- c(1,2,4,5,6,7,8, 10,11,12,13)
ids2 <- c(1,66,36)
ids3 <- c(20,71,101,4,76,192)
idsBad <- c(76,192, 101)


#file_name <- "n_iter_300_ntry_3000_0800"
#file_name <- "n_iter_100_ntry_200_4509"
file_name <- "n_iter_100_ntry_400_3247"

# list.files(here('data/abc_reject'))
# (grep("Train",list.files(here(paste0('data/abc_reject/',file_name)),
#                                            pattern="EXAM_Test",full.names = TRUE),
#                                 invert=TRUE, value=TRUE))


ind_fits <- map(list.files(here(paste0('data/abc_reject/'),file_name),full.names=TRUE), readRDS)
# ind_fits <- map(list.files(here('data/abc_reject/n_iter_2000_ntry_10_2918'),full.names=TRUE), readRDS)
ind_fits_df <- ind_fits |> map(~list(dat=.x[[1]], Model = .x[["Model"]], Fit_Method=.x[["Fit_Method"]]))
ind_fits_df <- ind_fits_df |> map(~rbindlist(.x$dat) |> mutate(Model = .x$Model, Fit_Method = .x$Fit_Method)) |> rbindlist() 


#run_params <- tibble::lst(file_name,cMean=ind_fits[[1]]$cMean, cSig=ind_fits[[1]]$cSig, 
                          
run_params <- keep(ind_fits[[1]], ~any(class(.x) %in% c("character", "numeric")))
run_params <- tibble(!!!run_params) |> select(-Model,-Fit_Method)

extract_info <- function(raw_names, run_params) {
  fname <- tools::file_path_sans_ext(basename(raw_names))
  n_samp <- str_extract(raw_names, "(?<=_)\\d+(?=_)")
  ntry_ <- str_extract(raw_names, "(?<=ntry_)\\d+")
  type <- 'ss'
  data.frame(run_params, n_samp, ntry_, type, fname)
}

run_info <- extract_info(file_name,run_params)

```


```{r}

generate_data <- function(Model, post_samples, data, num_samples = 1, return_dat = "train_data, test_data") {
  # Filter data for the specific id without invalidating selfref
  sbj_data <- copy(data[id == post_samples$id[1]])
  simulation_function <- ifelse(Model == "EXAM", full_sim_exam, full_sim_alm)

  target_data <- switch(return_dat,
                        "test_data" = copy(sbj_data[expMode2 == "Test"]),
                        "train_data" = copy(sbj_data[expMode2 == "Train"]),
                        "train_data, test_data" = copy(sbj_data[expMode2 %in% c("Test", "Train")]))
  
  post_samples <- post_samples[order(mean_error)][1:num_samples, .(c, lr, mean_error, rank = .I)]

  simulated_data_list <- lapply(1:nrow(post_samples), function(i) {
    params <- post_samples[i]
    sim_data <- simulation_function(sbj_data, params$c, params$lr, input_layer = input_layer, 
                                    output_layer = output_layer, return_dat = return_dat)
    sim_data_dt <- data.table(id = sbj_data$id[1], condit = sbj_data$condit[1], 
                              expMode2 = target_data$expMode2, Model = Model,tr=target_data$tr,
                              y = target_data$y, x = target_data$x, c = params$c, 
                              lr = params$lr, mean_error = params$mean_error, rank = i,
                              pred = sim_data)
    return(sim_data_dt)
  })
  
  result_dt <- rbindlist(simulated_data_list)
  setcolorder(result_dt, c("id", "condit", "expMode2","tr", "c", "lr", "x", "y", "pred"))
  return(result_dt)
}

future::plan(multisession)

nestSbjModelFit <- ind_fits_df %>% nest(.by=c(id,Model,Fit_Method))

post_dat <- nestSbjModelFit |> mutate(pp=furrr::future_pmap(list(id,Model,Fit_Method,data), ~{
   generate_data(..2, ..4 |> mutate(id=..1), ds, num_samples = 50, return_dat="test_data")
   })) |> 
  select(Fit_Method,pp,-data) |>  
  unnest(pp) |> filter(expMode2=="Test") |> as.data.table()

post_dat_avg <- post_dat |> group_by(id, condit, Model, Fit_Method, x, c, lr, rank) |> 
  summarise(y = mean(y), pred = mean(pred), error = y - pred) |> as.data.table()

setorder(post_dat_avg, id, x, rank)
post_dat_l <- melt(post_dat_avg, id.vars = c("id", "condit", "Model", "Fit_Method", "x", "c", "lr", "rank","error"),
                   measure.vars = c("pred", "y"), variable.name = "Resp", value.name = "val")
post_dat_l[, Resp := fifelse(Resp == "y", "Observed",
                             fifelse(Model == "ALM", "ALM", "EXAM"))]
setorder(post_dat_l, id, Resp)
rm(post_dat_avg)
plan(sequential)
# AIC and BIC
# mc2 <- post_dat |> filter(rank==1) |>
#   group_by(id,condit,Model,Fit_Method) |>
#   mutate(e2=abs(y-pred),n=n(),k=2) |>
#   group_by(id,condit,Model,Fit_Method,x) |>
#   summarise(y=mean(y), pred=mean(pred), e2=mean(e2), mean_error=first(mean_error),n=first(n),k=first(k)) |>
#   group_by(id,condit,Model,Fit_Method) |>
#   summarise(e2=mean(e2), mean_error=first(mean_error),n=first(n),k=first(k)) |>
#   mutate(AIC=2*k+n*mean_error, BIC=log(n)*k+n*mean_error)


# Single run extraction:
# exam_test <- ind_fits_df |> filter(Model == "EXAM", Fit_Method == "Test")
# post_dat_trial <- exam_test %>% split(f =c(.$id), drop=TRUE) |> 
#   map(~generate_data(.x$Model, .x, ds, num_samples = 15, return_dat="train_data, test_data")) |> 
#   rbindlist() |> 
#   filter(expMode2=="Test")
# 
# post_dat_avg <- post_dat_trial |> group_by(id,condit,x,c,lr,rank) |> 
#   summarise(y = mean(y),pred = mean(pred)) |> arrange(id,x,rank)



#post_dat_l |> filter(Model == "EXAM", Fit_Method == "Test",id==1) |> arrange(rank)

# head(post_dat_l)

```





#### `{r} kable(run_info)`


##  Posterior Average Table: 
```{r fig.width=12, fig.height=17}
#| eval: true
#| tbl-cap: "Posterior Distribution"
#| tbl-subcap: 
#| - "Full Posterior"
#| - "Best Posterior"

# post_tabs <- abc_tables(post_dat_l)
# post_tabs$et_sum |> gt::gt()

post_tabs <- abc_tables(post_dat,post_dat_l)


post_tabs$agg_full |> flextable::tabulator(rows=c("Fit_Method","Model"), columns=c("condit"), 
                       `ME` = as_paragraph(mean_error)) |> as_flextable() |>
  set_caption("Mean from full_posterior") 


post_tabs$agg_best |> flextable::tabulator(rows=c("Fit_Method","Model"), columns=c("condit"), 
                       `ME` = as_paragraph(mean_error)) |> as_flextable() |>
  set_caption("Mean from best parameters") 

post_tabs$agg_x_full |> 
  left_join(post_tabs$agg_x_best |> rename("best_error"=mean_error), 
            by=c("condit","Model","Fit_Method","x")) |> 
  arrange(desc(condit),x,Fit_Method) |>
  reactable::reactable() 

```


## Posterior Predictive: 
```{r, fig.width=12, fig.height=15}
#| eval: true
group_predictive_plots(post_dat_l) 
```
### Alt grouping of Posterior Predictive
```{r, fig.width=12, fig.height=15}
 group_predictive_plots2(post_dat_l) 
```




```{r fig.width=11, fig.height=12}
plot_indv_posterior(ind_fits_df |> mutate(Group=condit))

plot_indv_posterior(post_dat |> filter(rank==1) |> mutate(Group=condit))

# plot join density of c and lr
# post_dat |> filter(rank==1, c<.001) |> 
#   ggplot(aes(x=c, y=lr, color=condit)) + geom_point() + facet_wrap(~Model+Fit_Method)

post_dat |>  filter(c<.001) |> ggplot(aes(x=c,y=condit,fill=condit)) + geom_boxploth(width=.5) + facet_wrap(~Model+Fit_Method,scales="free")

```

## Model Comparison:
```{r, fig.width=12, fig.height=13}
#| eval: true
indv_best_plots(post_dat_l)

post_dat_l |> group_by(condit,Model,Fit_Method,rank) |> 
   summarise(mean_error=mean(abs(error)), n=n()) |> 
   ggplot(aes(x=rank,y=mean_error,fill=Model))+geom_col()+facet_wrap(~Fit_Method,scales="free")

post_dat_l |> group_by(condit,Model,Fit_Method,rank,x) |> 
   summarise(mean_error=mean(abs(error)), n=n()) |> 
   ggplot(aes(x=rank,y=mean_error,fill=Model))+geom_col()+facet_wrap(Fit_Method~x,ncol=6,scales="free")
 

group_best_plots(post_dat_l)

```

## Individual Predictive Plots

```{r fig.width=11, fig.height=12}
indv_predictive_plots(post_dat_l, ids2)
indv_predictive_plots(post_dat_l, idsBad)
```


### Subject 1
```{r fig.width=11, fig.height=12}
 
indv_predictive_dist((post_dat_l |> filter(rank<=200)),ind_fits_df, sbj=list(1))

abc_tables(post_dat |> filter(id==1))$agg_full |> flextable() |> set_caption("Mean from full posterior")
abc_tables(post_dat |> filter(id==1))$agg_best |> flextable() |> set_caption("Mean from best parameters")

post_dat_l |> filter(id==1) |> group_by(condit,Model,Fit_Method,rank) |> 
   summarise(mean_error=mean(abs(error)), n=n()) |> 
   ggplot(aes(x=rank,y=mean_error,fill=Model))+geom_col()+facet_wrap(~Fit_Method)

# plot_indv_posterior(ind_fits_df |> filter(id==1))
# ind_fits_df |> filter(id==1, Fit_Method=="Test Only", Model=="EXAM") |> pull(c) |> unique()



```



### Subject 36
```{r fig.width=11, fig.height=12}
 
indv_predictive_dist(post_dat_l,ind_fits_df, sbj=list(36))
#abc_tables(post_dat_l |> filter(id==36))$et_sum |> gt::gt()

abc_tables(post_dat |> filter(id==36))$agg_full |> flextable() |> set_caption("Mean from full posterior")
abc_tables(post_dat |> filter(id==36))$agg_best |> flextable() |> set_caption("Mean from best parameters")
# plot_indv_posterior(ind_fits_df |> filter(id==1))
# ind_fits_df |> filter(id==1, Fit_Method=="Test Only", Model=="EXAM") |> pull(c) |> unique()
```





## Learning Curves
```{r fig.width=11, fig.height=16}
#| fig-cap: "Learning Curves"
#| fig-width: 11
#| fig-height: 13


pd_train <- nestSbjModelFit |> mutate(pp=furrr::future_pmap(list(id,Model,Fit_Method,data), ~{
   generate_data(..2, ..4 |> mutate(id=..1), ds, num_samples = 10, return_dat="train_data, test_data")
   })) |> 
  select(Fit_Method,pp,-data) |>  
  unnest(pp) |> as.data.table() |> filter(expMode2=="Train") 



nbins <- 3
pd_train <- pd_train |> group_by(id,condit,Model,Fit_Method) |>
  mutate(Block=cut(tr,breaks=seq(1,max(tr), length.out=nbins+1),include.lowest=TRUE,labels=FALSE))
setorder(pd_train, id, x,Block, rank)

pd_train_l <- melt(pd_train, id.vars = c("id", "condit", "Model","Block", "Fit_Method", "x", "c", "lr", "rank"),
                   measure.vars = c("pred", "y"), variable.name = "Resp", value.name = "val") |> as.data.table()
pd_train_l[, Resp := fifelse(Resp == "y", "Observed",
                             fifelse(Model == "ALM", "ALM", "EXAM"))] 
setorder(pd_train_l, id,Block, Resp) 

pd_train_l <- pd_train_l |> mutate(dist=abs(x-val))

pd_train_l |> ggplot(aes(x=Block,y=dist,col=Resp))+
  geom_point( stat = "summary", fun = mean) + 
  stat_summary( geom = "line", fun = mean) +
  stat_summary(geom="errorbar",fun.data=mean_se,width=.1,alpha=.4) +
  #facet_wrap(condit~x,scales="free")
  ggh4x::facet_nested_wrap(Fit_Method~condit~x,scales="free") +
  #coord_cartesian(ylim = c(600, 1300)) +
  labs(title="Full Posterior")


pd_train_l |>filter(rank==1) |> ggplot(aes(x=Block,y=val,col=Resp))+
  geom_point( stat = "summary", fun = mean) + 
  stat_summary( geom = "line", fun = mean) +
  stat_summary(geom="errorbar",fun.data=mean_se,width=.1,alpha=.4) +
  #facet_wrap(condit~x,scales="free")
  ggh4x::facet_nested_wrap(Fit_Method~condit~x,scales="free") +
  coord_cartesian(ylim = c(600, 1300)) +labs(title="Best Parameters")

pd_train_l |> filter(Fit_Method=="Train") |> 
  mutate(x=as.factor(x), Resp=as.factor(Resp), Block=as.factor(Block)) |>
  ggplot(aes(x=x,y=val,fill=Block))+
  stat_bar +
  #facet_wrap(condit~x,scales="free")
  facet_wrap(condit~Resp,scales="free",ncol=3) +
  coord_cartesian(ylim = c(600, 1300))


# pd_train_l |> filter(id %in% c(1,33,66, 36,76), Fit_Method=="Train") |> 
#   ggplot(aes(x=Block,y=val,col=Resp))+
#   geom_point( stat = "summary", fun = mean) + 
#   stat_summary( geom = "line", fun = mean) +
#   stat_summary(geom="errorbar",fun.data=mean_se,width=.1,alpha=.4) +
#   #facet_wrap(condit~x,scales="free")
#   facet_wrap(id~x,scales="free",ncol=3) +
#   coord_cartesian(ylim = c(600, 1300))


pd_train_l |> 
  #filter(id %in% c(1,33,66, 36,76), Fit_Method=="Train") |> 
  mutate(x=as.factor(x), Resp=as.factor(Resp), Block=as.factor(Block)) |>
  ggplot(aes(x=x,y=dist,fill=Block))+
  stat_bar +
  #facet_wrap(condit~x,scales="free")
   #coord_cartesian(ylim = c(600, 1300)) +
  facet_wrap(id~Resp,scales="free",ncol=3) 
 


```



```{r}
#| eval: false
post_dat |> ggplot(aes(x=Model,y=mean_error,fill=Model))+geom_col()+facet_wrap(~Fit_Method)
post_dat |> ggplot(aes(x=Model,y=mean_error,fill=Model))+geom_col()+facet_wrap(condit~Fit_Method)


```



```{r}

# post_tabs <- map(post_tabs, ~cbind(.x,run_info) |> mutate(fn=file_name))
# # path1 = "../../../data/abc_tabs"
# saveRDS(post_tabs, paste0(path1, "/", tools::file_path_sans_ext(basename(file_name)),"_post_tab", ".rds"))
# saveRDS(post_tabs, paste0(tools::file_path_sans_ext(basename(file_name)),"_post_tab", ".rds"))
```


