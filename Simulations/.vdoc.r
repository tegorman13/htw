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
# load and view data
pacman::p_load(tidyverse,lme4,future,furrr,patchwork,here, furrr, future, pander)
purrr::walk(here(c("Functions/Display_Functions.R", "Functions/alm_core.R","Functions/misc_model_funs.R")),source)
select <- dplyr::select; mutate <- dplyr::mutate 

ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))

dsAvg <- ds |> group_by(condit,expMode2,tr, x) |> 
  summarise(y=mean(y),.groups="keep") 

vAvg <- dsAvg |> filter(condit=="Varied")
cAvg <- dsAvg |> filter(condit=="Constant")

head(dsAvg) |> pander::pandoc.table(style = "simple")

i1 <- ds |> filter(id=="1")

input.layer <- c(100,350,600,800,1000,1200)
output.layer <- c(100,350,600,800,1000,1200)


default.layer <- c(100,350,600,800,1000,1200)

output.layer <- adjust_layer(default.layer,k=3)

#
#
#
#
#
#


a_testOnly=list(pred_dat="test_avg",pred_fun="alm.responseOnly",loss_fun="RMSE",loss_data="test_error")
a_trainOnly=list(pred_dat="test_avg",pred_fun="alm.responseOnly",loss_fun="RMSE",loss_data="train_error")
a_testTrain=list(pred_dat="test_avg",pred_fun="alm.responseOnly",loss_fun="RMSE",loss_data="test_error+train_error")


e_testOnly=list(pred_dat="test_avg",pred_fun="exam.response",loss_fun="RMSE",loss_data="test_error")
e_trainOnly=list(pred_dat="test_avg",pred_fun="exam.response",loss_fun="RMSE",loss_data="train_error")
e_testTrain=list(pred_dat="test_avg",pred_fun="exam.response",loss_fun="RMSE",loss_data="test_error+train_error")

h_testOnly=list(pred_dat="test_avg",pred_fun="predict_alm_exam_weighted_hybrid",loss_fun="RMSE",loss_data="test_error")
h_trainOnly=list(pred_dat="test_avg",pred_fun="predict_alm_exam_weighted_hybrid",loss_fun="RMSE",loss_data="train_error")
h_testTrain=list(pred_dat="test_avg",pred_fun="exam.response",loss_fun="RMSE",loss_data="test_error+train_error")

hybrid_te_v <- wrap_grid(vAvg, c_values, lr_values, input.layer, output.layer,predParams=e_testOnly)
hybrid_tetr_v <- wrap_grid(vAvg, c_values, lr_values, input.layer, output.layer,predParams=e_testTrain)
hybrid_tr_v<- wrap_grid(vAvg, c_values, lr_values, input.layer, output.layer,predParams=e_trainOnly)


#
#
#
plan(multisession)

c_values <- seq(0.000001, 1.0, length.out=120)
lr_values <- seq(0.0000001, 4.0, length.out=200)

list2env(readRDS(here::here("data/model_cache/con_group_exam_fits.rds")), envir = .GlobalEnv)
list2env(readRDS(here::here("data/model_cache/var_group_exam_fits.rds")), envir = .GlobalEnv)




#
#
#




#
#
#
#


vte <- pluck(a_te_v, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex_te_v, "test") |> pull(pred)) |>  
  pivot_longer(Observed:EXAM, names_to="Model", values_to = "vx") |> 
  ggplot(aes(x,vx,fill=Model,col=Model,shape=Model)) +geom_point() + 
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(ds$x)))+ylim(0,1500) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Optimize for Test Only")

vtetr <- pluck(a_tetr_v, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex_tetr_v, "test") |> pull(pred)) |>  
  pivot_longer(Observed:EXAM, names_to="Model", values_to = "vx") |> 
  ggplot(aes(x,vx,fill=Model,col=Model,shape=Model)) +geom_point() + 
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(ds$x)))+ylim(0,1500) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Optimize for Test and Train")

vtr <- pluck(a_tr_v, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex_tr_v, "test") |> pull(pred)) |>  
  pivot_longer(Observed:EXAM, names_to="Model", values_to = "vx") |> 
  ggplot(aes(x,vx,fill=Model,col=Model,shape=Model)) +geom_point(alpha=.5) + 
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(ds$x)))+ ylim(0,1500) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Optimize for Train Only")


vte/vtetr/vtr


tvtr<- pluck(a_tr_v, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex_tr_v, "test") |> pull(pred)) 

tvtetr<- pluck(a_tetr_v, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex_tetr_v, "test") |> pull(pred)) 

tvte<- pluck(a_te_v, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex_te_v, "test") |> pull(pred)) 

pander("Varied-fit to test only") + pander(tvte) +
  pander("Varied - fit to train and test")+ pander(tvtetr)+
  pander("Varired-fit to train only")+ pander(tvtr)
#
#
#
#
#


cte <- pluck(a_te_c, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex0_te_c, "test") |> pull(pred)) |>  
  pivot_longer(Observed:EXAM, names_to="Model", values_to = "vx") |> 
  ggplot(aes(x,vx,fill=Model,col=Model,shape=Model)) +geom_point() + 
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(ds$x)))+ylim(0,1500) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Optimize for Test Only")

ctetr <- pluck(a_tetr_c, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex0_tetr_c, "test") |> pull(pred)) |>  
  pivot_longer(Observed:EXAM, names_to="Model", values_to = "vx") |> 
  ggplot(aes(x,vx,fill=Model,col=Model,shape=Model)) +geom_point() + 
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(ds$x)))+ylim(0,1500) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Optimize for Test and Train")

ctr <- pluck(a_tr_c, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex0_tr_c, "test") |> pull(pred)) |>  
  pivot_longer(Observed:EXAM, names_to="Model", values_to = "vx") |> 
  ggplot(aes(x,vx,fill=Model,col=Model,shape=Model)) +geom_point() + 
  scale_x_continuous(breaks=sort(unique(ds$x)), labels=sort(unique(ds$x)))+ ylim(0,1500) +
  theme(legend.title = element_blank(), legend.position="top") +ggtitle("Optimize for Train Only")



cte/ctetr/ctr


tctr<- pluck(a_tr_c, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex0_tr_c, "test") |> pull(pred)) 

tctetr<- pluck(a_tetr_c, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex0_tetr_c, "test") |> pull(pred)) 

tcte<- pluck(a_te_c, "test") |> rename(ALM=pred,Observed=y) %>% 
  cbind(.,EXAM=pluck(ex0_te_c, "test") |> pull(pred)) 

pander("Constant-fit to test only") + pander(tcte) +
  pander("Constant - fit to train and test")+ pander(tctetr)+
  pander("Constant-fit to train only")+ pander(tctr)

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
