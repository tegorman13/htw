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



input_layer =  c(100,350,600,800,1000,1200)
output_layer = input_layer



generate_prior_c_lr <- function(n) {
  prior_samples <- tibble(
    c = runif(n, 0.000001, 7),
    lr = runif(n, 0.000001, 7),
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
  
  #fd <- rbind(train_data,test_data)
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
  
  #fd <- rbind(train_data,test_data)
  fd = eval(parse(text=paste0("rbind(",return_dat,")")))
  if(mode=="sim"){return(fd$pred)
  }else {return(fd)}
  
}



# General function for ABC fits with model flexibility
run_abc_fits <- function(data, input_layer, output_layer, simulation_function, n_prior_samples = 1000, tol = 0.01, return_dat = "train_data,test_data") {
  
  input_layer =  c(100,350,600,800,1000,1200)
  output_layer = input_layer
  
  # Generate prior samples
  prior_samples <- generate_prior_c_lr(n_prior_samples)
  
  # Convert the string to a function
  if (is.character(simulation_function)) {
    simulation_function <- get(simulation_function, mode = "function")
  } else {
    simulation_function <- match.fun(simulation_function)
  }
  
  plan(multisession)
  # Run simulations for each prior sample
  simulation_results <- future_map_dfc(seq_len(nrow(prior_samples)), function(idx) {
    params <- prior_samples[idx, ]
    simulation_function(data=as.data.table(data), c=params$c, lr=params$lr, input_layer=input_layer, output_layer=output_layer, return_dat = return_dat)
  }, .options = furrr_options(seed = TRUE))
  
  # Extract target data
  target_data_train_test <- data[expMode2 %in% c("Test", "Train"), ]$y
  target_data_test <- data[expMode2 == "Test", ]$y
  target_data_train <- data[expMode2 == "Train", ]$y
  
  # ABC for train and test data
  abc_train_test <- abc(
    target = target_data_train_test,
    param = prior_samples,
    sumstat = do.call(rbind, simulation_results),
    tol = tol,
    method = "rejection"
  )
  
  # ABC for test data only
  abc_test <- abc(
    target = target_data_test,
    param = prior_samples,
    sumstat = do.call(rbind, simulation_results[85:90, ]),
    tol = tol,
    method = "rejection"
  )
  
  # ABC for train data only
  abc_train <- abc(
    target = target_data_train,
    param = prior_samples,
    sumstat = do.call(rbind, simulation_results[1:84, ]),
    tol = tol,
    method = "rejection"
  )
  
  # Return results
  tibble::lst(abc_train_test = abc_train_test, abc_test = abc_test, abc_train,data=data,n_prior_samples,prior_samples,tol, simulation_function, targets=tibble::lst(target_data_train_test, target_data_test))
}



n_priors=500000
args_list <- tibble::lst(
  abc_v_exam=list(data = avg_dsv, input_layer = input_layer, output_layer = output_layer, simulation_function = full_sim_exam, n_prior_samples = n_priors),
  abc_v_alt_exam=list(data = avg_dsv, input_layer = input_layer, output_layer = output_layer, simulation_function = full_sim_alt_exam, n_prior_samples = n_priors),
  abc_v_alm=list(data = avg_dsv, input_layer = input_layer, output_layer = output_layer, simulation_function = full_sim_alm, n_prior_samples = n_priors),
  abc_c_alm=list(data = avg_dsc, input_layer = input_layer, output_layer = output_layer, simulation_function = full_sim_alm, n_prior_samples = n_priors),
  abc_c_exam=list(data = avg_dsc, input_layer = input_layer, output_layer = output_layer, simulation_function = full_sim_exam, n_prior_samples = n_priors),
  abc_c_alt_exam=list(data = avg_dsc, input_layer = input_layer, output_layer = output_layer, simulation_function = full_sim_alt_exam, n_prior_samples = n_priors)
)

plan(multisession)
abc_list <- future_map(args_list, ~do.call(run_abc_fits, .x))


saveRDS(abc_list,here::here(paste0("data/model_cache/abc_group500k_",format(Sys.time(), "%H_%M_%OS"),".rds")))


