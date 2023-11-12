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
#
#
# load and view data
pacman::p_load(tidyverse,patchwork,here, pander, latex2exp)
purrr::walk(here::here(c("Functions/Display_Functions.R", "Functions/alm_core.R","Functions/misc_model_funs.R")),source)
select <- dplyr::select; mutate <- dplyr::mutate 

ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))
dsAvg <- ds |> group_by(condit,expMode2,tr, x) |> 
  summarise(y=mean(y),.groups="keep") 

vAvg <- dsAvg |> filter(condit=="Varied")
cAvg <- dsAvg |> filter(condit=="Constant")

i1 <- ds |> filter(id=="3")

input.layer <- c(100,350,600,800,1000,1200)
output.layer <- c(100,350,600,800,1000,1200)


purrr::walk(c("con_group_exam_fits", "var_group_exam_fits", "hybrid_group_exam_fits"), 
            ~ list2env(readRDS(here::here(paste0("data/model_cache/", .x, ".rds"))), 
            envir = .GlobalEnv))
#
#
#
#
#
#
#| label: fig-alm-diagram
#| fig.cap: The basic structure of the ALM model. 

alm_plot()

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
#
#
#

almParamsV <- cbind(Model="ALM Test Only",pluck(a_te_v, "Fit"), pluck(a_te_v, "test") %>% summarise(Test_RMSE=RMSE(y,pred)) ) |>
  rbind(cbind(Model="ALM Test & Train", pluck(a_tetr_v,"Fit"), pluck(a_tetr_v, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  rbind(cbind(Model="ALM Train Only", pluck(a_tr_v, "Fit"), pluck(a_tr_v, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

examParamsV <- cbind(Model="EXAM Test Only",pluck(ex_te_v, "Fit"), pluck(ex_te_v, "test") %>% summarise(Test_RMSE=RMSE(y,pred))) |>
  rbind(cbind(Model="EXAM Test & Train", pluck(ex_tetr_v,"Fit"), pluck(ex_tetr_v, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  rbind(cbind(Model="EXAM Train Only", pluck(ex_tr_v, "Fit"), pluck(ex_tr_v, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))


hybridParamsV <-cbind(Model="Hybrid Test Only",pluck(hybrid_te_v, "Fit"), pluck(hybrid_te_v, "test") %>% summarise(Test_RMSE=RMSE(y,pred))) |>
  rbind(cbind(Model="Hybrid Test & Train", pluck(hybrid_tetr_v,"Fit"), pluck(hybrid_tetr_v, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  rbind(cbind(Model="Hybrid Train Only", pluck(hybrid_tr_v, "Fit"), pluck(hybrid_tr_v, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

almParamsC <- cbind(Model="ALM Test Only",pluck(a_te_c, "Fit"), pluck(a_te_c, "test") %>% summarise(Test_RMSE=RMSE(y,pred)) ) |>
  rbind(cbind(Model="ALM Test & Train", pluck(a_tetr_c,"Fit"), pluck(a_tetr_c, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  rbind(cbind(Model="ALM Train Only", pluck(a_tr_c, "Fit"), pluck(a_tr_c, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

examParamsC <- cbind(Model="EXAM Test Only",pluck(ex0_te_c, "Fit"), pluck(ex0_te_c, "test") %>% summarise(Test_RMSE=RMSE(y,pred))) |>
  rbind(cbind(Model="EXAM Test & Train", pluck(ex0_tetr_c,"Fit"), pluck(ex0_tetr_c, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  rbind(cbind(Model="EXAM Train Only", pluck(ex0_tr_c, "Fit"), pluck(ex0_tr_c, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))


hybridParamsC <-cbind(Model="Hybrid Test Only",pluck(hybrid_te_c, "Fit"), pluck(hybrid_te_c, "test") %>% summarise(Test_RMSE=RMSE(y,pred))) |>
  rbind(cbind(Model="Hybrid Test & Train", pluck(hybrid_tetr_c,"Fit"), pluck(hybrid_tetr_c, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  rbind(cbind(Model="Hybrid Train Only", pluck(hybrid_tr_c, "Fit"), pluck(hybrid_tr_c, "test") %>% summarise(Test_RMSE=RMSE(y,pred)))) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

#
#
#
#
#
#
#| label: tbl-e1-model-fits2V
#| tbl-cap: Varied Group - Fit Parameters and Model RMSE

pander(almParamsV, caption="ALM"); pander(examParamsV, caption="EXAM"); pander(hybridParamsV,caption="Hybrid") 
#
#
#
#| label: tbl-e1-model-fitsC
#| tbl-cap: Constant Group - Fit Parameters and Model RMSE

pander(almParamsC, caption="ALM"); pander(examParamsC, caption="EXAM"); pander(hybridParamsC,caption="Hybrid")
#
#
#
#
#
#
#
#| label: fig-model-preds-varied
#| fig-cap: Varied Group - Mean Model predictions vs. observations
#| fig-height: 12
#| fig-width: 14
#| column: screen-inset-right

####

vte <-  pluck(a_te_v, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex_te_v, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_te_v, "test") |> pull(pred)) |>  
  pivot_longer(Observed:Hybrid, names_to="Model", values_to = "vx") |> 
  ggplot(aes(x,vx,fill=Model, group=Model)) +geom_bar(position="dodge",stat="identity") +
  scale_fill_manual(values=col_themes$wes2)+
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(ds$x)))+ylim(0,1500) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Fit to Test Only")

vtetr <-  pluck(a_tetr_v, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex_tetr_v, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_tetr_v, "test") |> pull(pred)) |>  
  pivot_longer(Observed:Hybrid, names_to="Model", values_to = "vx") |> 
  ggplot(aes(x,vx,fill=Model, group=Model)) +geom_bar(position="dodge",stat="identity") + 
  scale_fill_manual(values=col_themes$wes2)+
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(ds$x)))+ylim(0,1500) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Fit to Test and Train")

vtr <-  pluck(a_tr_v, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex_tr_v, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_tr_v, "test") |> pull(pred)) |>  
  pivot_longer(Observed:Hybrid, names_to="Model", values_to = "vx") |> 
  ggplot(aes(x,vx,fill=Model, group=Model)) +geom_bar(position="dodge",stat="identity") +
  scale_fill_manual(values=col_themes$wes2)+
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(ds$x)))+ylim(0,1500) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Fit to Train Only")
  

  


vte/vtetr/vtr

#
#
#
#
#
#
#
#| label: tbl-e1-predsV
#| tbl-cap: Varied group - mean model predictions vs. observations
#| 
tvte<- pluck(a_te_v, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex_te_v, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_te_v, "test") |> pull(pred))

tvtetr<-pluck(a_tetr_v, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex_tetr_v, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_tetr_v, "test") |> pull(pred))

tvtr<- pluck(a_tr_v, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex_tr_v, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_tr_v, "test") |> pull(pred))

pander(tvte, caption="Varied fit to test only")
pander(tvtetr,caption="Varied fit to train and test")
pander(tvtr,caption="Varied fit to train only")

#
#
#
#
#
#
#| label: fig-model-preds-constant
#| fig-cap: Constant Group - Mean Model predictions vs. observations
#| fig-height: 12
#| fig-width: 14
#| column: screen-inset-right

####

cte <-  pluck(a_te_c, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex0_te_c, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_te_c, "test") |> pull(pred)) |>  
  pivot_longer(Observed:Hybrid, names_to="Model", values_to = "vx") |> 
  ggplot(aes(x,vx,fill=Model, group=Model)) +geom_bar(position="dodge",stat="identity") +
  scale_fill_manual(values=col_themes$wes2)+
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(ds$x)))+ylim(0,1500) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Fit to Test Only")

ctetr <-  pluck(a_tetr_c, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex0_tetr_c, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_tetr_c, "test") |> pull(pred)) |>  
  pivot_longer(Observed:Hybrid, names_to="Model", values_to = "vx") |> 
  ggplot(aes(x,vx,fill=Model, group=Model)) +geom_bar(position="dodge",stat="identity") + 
  scale_fill_manual(values=col_themes$wes2)+
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(ds$x)))+ylim(0,1500) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Fit to Test and Train")

ctr <-  pluck(a_tr_c, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex0_tr_c, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_tr_c, "test") |> pull(pred)) |>  
  pivot_longer(Observed:Hybrid, names_to="Model", values_to = "vx") |> 
  ggplot(aes(x,vx,fill=Model, group=Model)) +geom_bar(position="dodge",stat="identity") +
  scale_fill_manual(values=col_themes$wes2)+
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(ds$x)))+ylim(0,1500) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Fit to Train Only")
  
cte/ctetr/ctr

#
#
#
#
#
#| label: tbl-e1-predsC
#| tbl-cap: Constant group - mean model predictions vs. observations
#| 
tcte<- pluck(a_te_c, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex0_te_c, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_te_c, "test") |> pull(pred))

tctetr<-pluck(a_tetr_c, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex0_tetr_c, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_tetr_c, "test") |> pull(pred))

tctr<- pluck(a_tr_c, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex0_tr_c, "test") |> pull(pred)) %>%
  cbind(., Hybrid=pluck(hybrid_tr_c, "test") |> pull(pred))

pander(tcte, caption="Constant fit to test only")
pander(tctetr,caption="Constant fit to train and test")
pander(tctr,caption="Constant fit to train only")

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
#
#
#| eval: false
pluck(a_te_v, "train") |> pivot_longer(y:almResp, names_to="Resp", values_to = "vx") |> 
  mutate(dev=x-vx,abs_dev=abs(x-vx)) |>
  learn_curve_plot(tr, vx, Resp,facet_var=x, groupVec=Resp,nbins=8)

pluck(a_te_v, "train") |> pivot_longer(y:almResp, names_to="Resp", values_to = "vx") |> 
  mutate(dev=x-vx,abs_dev=abs(x-vx)) |>
  ungroup() %>%
  gather(key = "variable", value = "y_value", dev, abs_dev, vx) %>%
  group_by(variable) %>%
  group_map(~ learn_curve_plot(.x, x_var = tr, y_var = y_value, color_var = Resp, facet_var = x, groupVec = Resp, nbins = 8, y_label = .y$variable), .keep = TRUE)

list(a_tr_v, a_te_v,a_tetr_v) |> map( ~{pluck(.x, "train") |> pivot_longer(y:almResp, names_to="Resp", values_to = "vx") |> mutate(dev=x-vx,abs_dev=abs(x-vx)) |>
  ungroup() %>%
  gather(key = "variable", value = "y_value", dev, abs_dev, vx) %>%
  group_by(variable) %>%
  group_map(~ learn_curve_plot(.x, x_var = tr, y_var = y_value, color_var = Resp, facet_var = x, groupVec = Resp, nbins = 8, y_label = .y$variable), .keep = TRUE) })

#
#
#
#
#
#
#| eval: false
pluck(ex_te_v, "train") |> pivot_longer(y:almResp, names_to="Resp", values_to = "vx") |> 
  mutate(dev=x-vx,abs_dev=abs(x-vx)) |>
  learn_curve_plot(tr, vx, Resp,facet_var=x, groupVec=Resp,nbins=8)

pluck(ex_te_v, "train") |> pivot_longer(y:almResp, names_to="Resp", values_to = "vx") |> 
  mutate(dev=x-vx,abs_dev=abs(x-vx)) |>
  ungroup() %>%
  gather(key = "variable", value = "y_value", dev, abs_dev, vx) %>%
  group_by(variable) %>%
  group_map(~ learn_curve_plot(.x, x_var = tr, y_var = y_value, color_var = Resp, facet_var = x, groupVec = Resp, nbins = 8, y_label = .y$variable), .keep = TRUE)

#
#
#
#
#
#| eval: false


optimize_params_weighted_individual <- function(ds, c_values, lr_values, weight_exam_values, input.layer, output.layer) {
    all_results <- list()
    
    # Loop through each unique id
    for (individual in unique(ds$id)) {
        indiv_data <- ds[ds$id == individual, ]
        
        # Run the optimization function for the individual's data
        result <- optimize_params_weighted(indiv_data, c_values, lr_values, weight_exam_values, input.layer, output.layer)
        
        all_results[[as.character(individual)]] <- result
    }
    
    all_results
}

dss <- ds |> filter(id %in% c(1,2))

all_results_weighted_hybrid <- readRDS(here::here('data/model_cache/indv_hybrid_fits.rds'))
ma <- map(all_results_weighted_hybrid, "best_params") |> map("c")


data.frame(id=names(ma),c=as.numeric(ma))

ma = cbind(id=names(all_results_weighted_hybrid),map(all_results_weighted_hybrid, "best_params") |> map_dfr(magrittr::extract,c("c","lr","weight_exam")))

ds |> group_by(id,condit) |> distinct(id,condit) |> left_join(ma,by=join_by(id))


map(all_results_weighted_hybrid,"best_params") |> pluck("c")
all_results_weighted_hybrid[["1"]]$b

 map(~ map(.x$best_params, pluck, "c"))

 map_df(~ map_df(.x$train, pluck, "d"), .id = "density")
#
#
#
#
#
