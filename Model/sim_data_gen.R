pacman::p_load(tidyverse,data.table,abc,future,furrr,here,patchwork, conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
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


generate_prior_c_lr <- function(n) {
  prior_samples <- tibble(
    c = runif(n, 0.000001, 5),
    lr = runif(n, 0.000001, 5),
  )
  return(prior_samples)
}

full_sim_exam <- function(data, c, lr,pred_fun=exam.response, input_layer, output_layer,return_dat="test_data",mode="sim") {
  train_data <- data[expMode2=="Train", c("condit","tr","expMode2", "x","y")] 
  test_data <- data[expMode2=="Test", c("condit","tr","expMode2", "x","y")] 
  trainVec=sort(unique(train_data$x))
  if (train_data$condit[1] != "Varied") {
    trainVec <- c(0, trainVec)
  }
  
  train_results <- alm.sim(train_data, c, lr, input_layer, output_layer)
  
  test_prediction <- map_dbl(test_data$x, ~ pred_fun(.x, c, input_layer, 
                                                     output_layer, train_results$wm,  trainVec=trainVec))
  
  train_data$pred <- train_results$d$almResp
  test_data$pred <- test_prediction
  
  fd = eval(parse(text=paste0("rbind(",return_dat,")")))
  if(mode=="sim"){return(fd$pred)
  }else {return(fd)}
  
}

full_sim_alm <- function(data, c, lr,pred_fun=alm.responseOnly, input_layer, output_layer,return_dat="test_data",mode="sim") {
  train_data <- data[expMode2=="Train", c("condit","tr","expMode2", "x","y")] 
  test_data <- data[expMode2=="Test", c("condit","tr","expMode2", "x","y")] 
  trainVec=sort(unique(train_data$x))
  
  train_results <- alm.sim(train_data, c, lr, input_layer, output_layer)
  test_prediction <- map_dbl(test_data$x, ~ alm.responseOnly(.x, c, input_layer, output_layer, train_results$wm,  trainVec=trainVec))
  
  train_data$pred <- train_results$d$almResp
  test_data$pred <- test_prediction
  
  #fd <- rbind(train_data,test_data)
  fd = eval(parse(text=paste0("rbind(",return_dat,")")))
  if(mode=="sim"){return(fd$pred)
  }else {return(fd)}
  
}

full_sim_alt_exam <- function(data, c, lr, pred_fun=alt_exam, input_layer, output_layer,return_dat="test_data",mode="sim") {
  train_data <- data[expMode2=="Train", c("condit","tr","expMode2", "x","y")] 
  test_data <- data[expMode2=="Test", c("condit","tr","expMode2", "x","y")] 
  trainVec=sort(unique(train_data$x))
  trainVecY=train_data$y
  
  train_results <- alm.sim(train_data, c, lr, input_layer, output_layer)
  test_prediction <- map_dbl(test_data$x, ~ alt_exam(.x, c, input_layer, output_layer, train_results$wm,  trainVecX=trainVec, trainVecY=trainVecY))
  
  train_data$pred <- train_results$d$almResp
  test_data$pred <- test_prediction
  
  fd = eval(parse(text=paste0("rbind(",return_dat,")")))
  if(mode=="sim"){return(fd$pred)
  }else {return(fd)}
  
}


sim_data_gen <- function(data, input_layer, output_layer, simulation_function, prior_samples, return_dat) {
  simulation_function <- match.fun(simulation_function)
  plan(multisession)
  suppressMessages(simulation_results <- future_map_dfc(seq_len(nrow(prior_samples)), function(idx) {
    params <- prior_samples[idx, ]
    simulation_function(data=as.data.table(data), c=params$c, lr=params$lr, input_layer=input_layer, output_layer=output_layer, return_dat = return_dat)
  }, .options = furrr_options(seed = TRUE)))
}




input_layer =  c(100,350,600,800,1000,1200)
output_layer = input_layer
n_prior_samples=1000000
prior_samples <- generate_prior_c_lr(n_prior_samples)
return_dat="test_data,train_data"



# Create a list of arguments for each sim_data_gen call
args_list <- list(
  Exam_Varied=list(avg_dsv, input_layer, output_layer, full_sim_exam, prior_samples, return_dat),
  Exam_Constant=list(avg_dsc, input_layer, output_layer, full_sim_exam, prior_samples, return_dat),
  ALM_Varied=list(avg_dsv, input_layer, output_layer, full_sim_alm, prior_samples, return_dat),
  ALM_Constant=list(avg_dsc, input_layer, output_layer, full_sim_alm, prior_samples, return_dat),
  Alt_Varied=list(avg_dsv, input_layer, output_layer, full_sim_alt_exam, prior_samples, return_dat),
  Alt_Constant=list(avg_dsc, input_layer, output_layer, full_sim_alt_exam, prior_samples, return_dat)
)

sim_data_wrapper <- function(args) {
  sim_data_gen(args[[1]], args[[2]], args[[3]], simulation_function = args[[4]], args[[5]], return_dat = args[[6]])
}





t1=system.time({
plan(multisession)
sim_dataAll <- future_map(args_list, sim_data_wrapper, .progress = TRUE,.options = furrr_options(seed = TRUE))
})
t1


saveRDS(tibble::lst(sim_dataAll,prior_samples,args_list,time=t1[3]),
        file = here::here(paste0("data/sim_data/","sim_data_", n_prior_samples,"_",format(Sys.time(),"%H_%M_%OS"), ".rds")))


#r=readRDS(here::here(paste0("data/sim_data/","sim_data_3010_40_32.rds")))


#exam_v_500k <- results[[1]]






# exam_v_500k <- sim_data_gen(avg_dsv, input_layer, output_layer,simulation_function= full_sim_exam, prior_samples, return_dat = return_dat)
# exam_c_500k <- sim_data_gen(avg_dsc, input_layer, output_layer,simulation_function= full_sim_exam, prior_samples, return_dat = return_dat)
# alm_v_500k <- sim_data_gen(avg_dsv, input_layer, output_layer,simulation_function= full_sim_alm, prior_samples, return_dat = return_dat)
# alm_c_500k <- sim_data_gen(avg_dsc, input_layer, output_layer,simulation_function= full_sim_alm, prior_samples, return_dat = return_dat)
# alt_v_500k <- sim_data_gen(avg_dsv, input_layer, output_layer,simulation_function= full_sim_alt_exam, prior_samples, return_dat = return_dat)
# alt_c_500k <- sim_data_gen(avg_dsc, input_layer, output_layer,simulation_function= full_sim_alt_exam, prior_samples, return_dat = return_dat)
# 
# 
# saveRDS(tibble::lst(exam_v_500k, exam_c_500k, alm_v_500k, alm_c_500k, alt_v_500k, alt_c_500k, prior_samples), 
#         file = paste0(here::here("data/sim_data/"),"sim_data_", n_prior_samples, ".rds"))
# 






# 500 secs for 90k on tg_m1 - with future_map
# 300k took about 6 hours on tg_m1
# 8741s for 1M on tg_m1








# 
# abc <- run_abc_fits(avg_dsv, exam_v_500k, input_layer, output_layer)
# 
# 
# calculate_distance <- function(simulated, observed) {
#   return(mean((simulated - observed)^2)) #MSE
# }
# 
# run_abc_fits <- function(data,sim_data, input_layer, output_layer, tol = 100000) {
#   
#   target_data_train_test <- data[expMode2 %in% c("Test", "Train"), ]$y
#   target_data_test <- data[expMode2 == "Test", ]$y
#   target_data_train <- data[expMode2 == "Train", ]$y
#   
#   teter_distances <- sapply(sim_data, calculate_distance, observed = target_data_train_test)
#   te_distances <- sapply(sim_data[85:90, ], calculate_distance, observed = target_data_test)
#   tr_distances <- sapply(sim_data[1:84, ], calculate_distance, observed = target_data_train)
#   
#   teter_results <- tibble(distance = teter_distances, c = prior_samples$c, lr = prior_samples$lr) %>% filter(distance <= tol) %>% arrange(distance)
#   te_results <- tibble(distance = te_distances, c = prior_samples$c, lr = prior_samples$lr) %>% filter(distance <= tol) %>% arrange(distance)
#   tr_results <- tibble(distance = tr_distances, c = prior_samples$c, lr = prior_samples$lr) %>% filter(distance <= tol) %>% arrange(distance)
#   
#   list(teter_results = teter_results, te_results = te_results, tr_results = tr_results)
# }
