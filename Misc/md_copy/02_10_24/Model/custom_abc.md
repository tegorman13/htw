
pacman::p_load(tidyverse,data.table,abc,future,furrr,here,patchwork, conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
options(scipen = 999)
walk(c("fun_alm","fun_model"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
dsv <- ds |> filter(condit=="Varied")  
dsc <- ds |> filter(condit=="Constant") 
tMax=84
avg_dsv <- ds |> filter(condit=="Varied",expMode2=="Train") |> group_by(tr) %>%
  mutate(bandInt2 = sample(rep(c(800, 1000, 1200), each = tMax / 3), tMax, replace = FALSE)[tr]) %>%
  filter(bandInt2 == x) |> select(-bandInt2) |> group_by(tr,condit,x,expMode2) |> summarise(y=mean(y),.groups="keep") |>
  rbind(dsv |> filter(expMode2=="Test") |> group_by(condit,x,expMode2) |> summarise(y=mean(y),tr=1,.groups="keep") ) |> setDT()

avg_dsc <- ds |> filter(condit=="Constant",expMode2=="Train",tr<=tMax) |> group_by(tr, condit,x,expMode2) |> 
  summarise(y=mean(y),.groups="keep") |> rbind(dsc |> filter(expMode2=="Test") |> group_by(condit,x,expMode2) |> summarise(y=mean(y),tr=1,.groups="keep") ) |> setDT()



 #sd <- readRDS(here::here("data/sim_data/sim_data_10k.rds"))
#sd <- readRDS(here::here("data/sim_data/sim_data_300k.rds"))
#sd <- readRDS(here::here("data/sim_data/sim_data_1M_13_38_28.rds"))
sd <- readRDS(here::here("/Users/thomasgorman/Documents/Backup/htw_local/data/sim_data/sim_data_2M_12_21.rds"))
# sd <- readRDS(here::here("/Users/thomasgorman/Documents/Backup/htw_local/data/sim_data/sim_data_1p5M_10_59_41.rds"))

sd$sim_dataAll <- map(sd$sim_dataAll, setDT)
names(sd)

#map(sd$sim_dataAll, class)


dist_rmse <- function(simulated, observed) {
  return(sqrt(mean((simulated - observed)^2))) #RMSE
}


run_abc_fits <- function(data,sim_data, prior_samples, dist_fun="dist_rmse",Model,Group,pct_keep) {
  tol = 500
  rev=FALSE
  target_data_train_test <- data[expMode2 %in% c("Test", "Train"), ]$y
  target_data_test <- data[expMode2 == "Test", ]$y
  target_data_train <- data[expMode2 == "Train", ]$y
  
  # correct order of train & test
  if (rev) { sim_data <- sim_data[7:90, ] |> bind_rows(sim_data[1:6, ]) }

  fn_name=dist_fun
  dist_fun=get(dist_fun)
  
  te_distances <- purrr::map_dbl(sim_data[85:90, ], dist_fun, observed = target_data_test)
  tr_distances <- purrr::map_dbl(sim_data[1:84, ], dist_fun, observed = target_data_train)
  teter_distances <- 0.5 * te_distances + 0.5 * tr_distances
  
  teter_results <- tibble(distance = teter_distances, c = prior_samples$c, 
    lr = prior_samples$lr, sim_index= seq_along(teter_distances)) |> 
    arrange(distance) |> 
    filter(if (!is.null(pct_keep)) row_number() <= n() * pct_keep else distance <= tol) |>
    mutate(rank=row_number(),Model,Group,Fit_Method="Test & Train") |>
    rowwise() |>
    mutate(sim_dat = list(tibble(Model,Fit_Method,c,lr,distance,rank,data,setNames(sim_data[, .SD, .SDcols = sim_index], "pred")) |> 
    mutate(resid=y-pred)))
    
 # teter_pp <- sim_data[,.SD,.SDcols=teter_results$sim_index] |> bind_rows()
#  k = teter_results[1,]$sim_dat

  te_results <- tibble(distance = te_distances, c = prior_samples$c, 
  lr = prior_samples$lr,  sim_index= seq_along(te_distances)) |> 
    arrange(distance) |> 
    filter(if (!is.null(pct_keep)) row_number() <= n() * pct_keep else distance <= tol) |>
    mutate(rank=row_number(),Model,Group,Fit_Method="Test Only") |>
    rowwise() |>
    mutate(sim_dat = list(tibble(Model,Fit_Method,c,lr,distance,rank,data,setNames(sim_data[, .SD, .SDcols = sim_index], "pred")) |> 
    mutate(resid=y-pred)))
  
  tr_results <- tibble(distance = tr_distances, c = prior_samples$c, 
    lr = prior_samples$lr,  sim_index= seq_along(tr_distances)) |>
    arrange(distance) |> 
    filter(if (!is.null(pct_keep)) row_number() <= n() * pct_keep else distance <= tol) |>
    mutate(rank=row_number(),Model,Group,Fit_Method="Train Only") |>
    rowwise() |>
    mutate(sim_dat = list(tibble(Model,Fit_Method,c,lr,distance,rank,data,setNames(sim_data[, .SD, .SDcols = sim_index], "pred")) |> 
    mutate(resid=y-pred)))
  
  tibble::lst(teter_results = teter_results, te_results = te_results, tr_results = tr_results, Model, Group, pct_keep, fn_name)
}


pct_keep=.001
abc_ev <- run_abc_fits(avg_dsv, sim_data= sd$sim_dataAll$Exam_Varied ,sd$prior_samples, dist_fun="dist_rmse",Model="EXAM", Group="Varied",pct_keep)
abc_almv <- run_abc_fits(avg_dsv, sim_data=sd$sim_dataAll$ALM_Varied,sd$prior_samples, dist_fun="dist_rmse", Model="ALM", Group="Varied", pct_keep)
abc_altv <- run_abc_fits(avg_dsv, sim_data=sd$sim_dataAll$Alt_Varied,sd$prior_samples, dist_fun="dist_rmse",Model="Alt_EXAM", Group="Varied",pct_keep)

abc_ec <- run_abc_fits(avg_dsc, sim_data=sd$sim_dataAll$Exam_Constant,sd$prior_samples, dist_fun="dist_rmse",Model="EXAM", Group="Constant",pct_keep)
abc_almc <- run_abc_fits(avg_dsc, sim_data=sd$sim_dataAll$ALM_Constant,sd$prior_samples, dist_fun="dist_rmse",Model="ALM", Group="Constant",pct_keep)
abc_altc <- run_abc_fits(avg_dsc, sim_data=sd$sim_dataAll$Alt_Constant,sd$prior_samples, dist_fun="dist_rmse",Model="Alt_EXAM", Group="Constant",pct_keep)

saveRDS(tibble::lst(abc_ev,abc_almv,abc_altv,abc_ec,abc_almc,abc_altc),
  here::here(paste0("data/abc_2M_rmse_",sub("^0", "", gsub("\\.", "p", pct_keep)),".rds")))








tO <- system.time(abc_ev <- run_abc_fits(avg_dsv, sim_data= sd$sim_dataAll$Exam_Varied ,sd$prior_samples, dist_fun="dist_rmse")
)
  #  user  system elapsed 
  # 7.706   0.073   7.778 

tdt <- system.time(abc_ev <- run_abc_fits(avg_dsv, sim_data= setDT(sd$sim_dataAll$Exam_Varied) ,sd$prior_samples, dist_fun="dist_rmse")
)
#    user  system elapsed 
#  11.588   9.773  51.324 

# 1st vectorized
#    user  system elapsed 
# 133.136   1.472 134.624 




teter_results <- data.frame(distance = teter_distances, c = prior_samples$c, 
                        lr = prior_samples$lr, 
                        sim_index= seq_along(teter_distances)) 

abc_1M <- tibble::lst(abc_ev,abc_almv,abc_altv,abc_ec,abc_almc,abc_altc)
teter_combined <- map_dfr(abc_1M, "teter_results", .id = "Model") %>% update_columns()
te_combined <- map_dfr(abc_1M, "te_results", .id = "Model") %>% update_columns()
tr_combined <- map_dfr(abc_1M, "tr_results", .id = "Model") %>% update_columns()





# return_dat="test_data"
# return_dat="test_data,train_data"
# te_combined <- te_combined |> arrange(distance)
# prior_samples=te_combined[1,]
# args_list <- list(
#   Exam_Varied=list(avg_dsv, input_layer, output_layer, full_sim_exam, prior_samples, return_dat),
#   ALM_Varied=list(avg_dsv, input_layer, output_layer, full_sim_alm, prior_samples, return_dat),
#   Alt_Varied=list(avg_dsv, input_layer, output_layer, full_sim_alt_exam, prior_samples, return_dat)
# )

# sim_dataAll <- map(args_list, sim_data_wrapper)

# run_abc_fits(avg_dsv, sim_data=sim_dataAll$Exam_Varied,prior_samples, dist_fun="dist_rmse",pct_keep = NULL)
# run_abc_fits(avg_dsv, sim_data=sim_dataAll$ALM_Varied,prior_samples, dist_fun="dist_rmse",pct_keep = NULL)

sim_dataAll$ALM_Varied == sim_dataAll$Exam_Varied


sapply(sim_dataAll$Exam_Varied, dist_fun, observed = target_data_test)
sapply(sim_dataAll$ALM_Varied, dist_fun, observed = target_data_test)
sapply(sim_dataAll$Alt_Varied, dist_fun, observed = target_data_test)




# A tibble: 6 × 6
# # Groups:   Model, Group [6]
#   Model    distance         c    lr dist_fun  Group   
#   <chr>       <dbl>     <dbl> <dbl> <chr>     <chr>   
# 1 EXAM         188. 0.0000217 0.912 dist_rmse Varied  
# 2 ALM          187. 0.0000217 0.912 dist_rmse Varied  
# 3 Alt_Exam     188. 0.0000217 0.912 dist_rmse Varied  
# 4 EXAM         101. 0.0000217 0.912 dist_rmse Constant
# 5 ALM          100. 0.0000217 0.912 dist_rmse Constant
# 6 Alt_Exam     100. 0.0000217 0.912 dist_rmse Constant








te_combined |> group_by(Model, Group) |> filter(row_number() <= 1)
tr_combined |> group_by(Model, Group) |> filter(row_number() <= 1)











teter_combined |> ggplot(aes(x=Group,y=distance,fill=Model)) + 
  stat_summary(fun=mean, geom="bar", position=position_dodge()) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge()) 

te_combined |> ggplot(aes(x=Group,y=distance,fill=Model)) + 
  stat_summary(fun=mean, geom="bar", position=position_dodge()) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge()) 

tr_combined |> ggplot(aes(x=Group,y=distance,fill=Model)) + 
  stat_summary(fun=mean, geom="bar", position=position_dodge()) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge()) 



# summary stats foro te_combined - table
te_combined |> group_by(Model) |> summarise(mean_distance=mean(distance),sd_distance=sd(distance),
  mean_c=mean(c),sd_c=sd(c),mean_lr=mean(lr),sd_lr=sd(lr)) 

tr_combined |> group_by(Model) |> summarise(mean_distance=mean(distance),sd_distance=sd(distance),
  mean_c=mean(c),sd_c=sd(c),mean_lr=mean(lr),sd_lr=sd(lr)) 



abc_ev$teter_results |> ggplot(aes(x=c)) + geom_density()
abc_ev$te_results |> ggplot(aes(x=c)) + geom_density()




# abc_ev <- run_abc_fits(avg_dsv, sim_data=sd$exam_v_500k,sd$prior_samples, input_layer, output_layer)
# abc_almv <- run_abc_fits(avg_dsv, sim_data=sd$alm_v_500k,sd$prior_samples, input_layer, output_layer)
# abc_altv <- run_abc_fits(avg_dsv, sim_data=sd$alt_v_500k,sd$prior_samples, input_layer, output_layer)
# 
# abc_ec <- run_abc_fits(avg_dsc, sim_data=sd$exam_c_500k,sd$prior_samples, input_layer, output_layer)
# abc_almc <- run_abc_fits(avg_dsc, sim_data=sd$alm_c_500k,sd$prior_samples, input_layer, output_layer)
# abc_altc <- run_abc_fits(avg_dsc, sim_data=sd$alt_c_500k,sd$prior_samples, input_layer, output_layer)
# 
# 
# saveRDS(tibble::lst(abc_ev,abc_almv,abc_altv,abc_ec,abc_almc,abc_altc),here::here("data/abc_results_500k.rds"))


#abc_1M <- readRDS(here::here("data/abc_results_1M.rds"))
#abc_500k <- readRDS(here::here("data/abc_results_500k.rds"))


