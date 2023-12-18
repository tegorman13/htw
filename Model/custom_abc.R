
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



sd <- readRDS(here::here("data/sim_data/sim_data_10k.rds"))
sd <- readRDS(here::here("data/sim_data/sim_data_300k.rds"))

calculate_distance <- function(simulated, observed) {
  return(mean((simulated - observed)^2)) #MSE
}


run_abc_fits <- function(data,sim_data, prior_samples, input_layer, output_layer, tol = 100000) {
  
  target_data_train_test <- data[expMode2 %in% c("Test", "Train"), ]$y
  target_data_test <- data[expMode2 == "Test", ]$y
  target_data_train <- data[expMode2 == "Train", ]$y
  
  teter_distances <- sapply(sim_data, calculate_distance, observed = target_data_train_test)
  te_distances <- sapply(sim_data[85:90, ], calculate_distance, observed = target_data_test)
  tr_distances <- sapply(sim_data[1:84, ], calculate_distance, observed = target_data_train)
  
  teter_results <- tibble(distance = teter_distances, c = prior_samples$c, lr = prior_samples$lr) |> 
    filter(distance <= tol) %>% arrange(distance)
  
  te_results <- tibble(distance = te_distances, c = prior_samples$c, lr = prior_samples$lr) |> 
    filter(distance <= tol) %>% arrange(distance)
  
  tr_results <- tibble(distance = tr_distances, c = prior_samples$c, lr = prior_samples$lr) |>
    filter(distance <= tol) %>% arrange(distance)
  
  list(teter_results = teter_results, te_results = te_results, tr_results = tr_results)
}


abc_ev <- run_abc_fits(avg_dsv, sim_data=sd$exam_v_500k,sd$prior_samples, input_layer, output_layer)
abc_almv <- run_abc_fits(avg_dsv, sim_data=sd$alm_v_500k,sd$prior_samples, input_layer, output_layer)
abc_altv <- run_abc_fits(avg_dsv, sim_data=sd$alt_v_500k,sd$prior_samples, input_layer, output_layer)

abc_ec <- run_abc_fits(avg_dsc, sim_data=sd$exam_c_500k,sd$prior_samples, input_layer, output_layer)
abc_almc <- run_abc_fits(avg_dsc, sim_data=sd$alm_c_500k,sd$prior_samples, input_layer, output_layer)
abc_altc <- run_abc_fits(avg_dsc, sim_data=sd$alt_c_500k,sd$prior_samples, input_layer, output_layer)


saveRDS(tibble::lst(abc_ev,abc_almv,abc_altv,abc_ec,abc_almc,abc_altc),here::here("data/abc_results_500k.rds"))



abc_500k <- readRDS(here::here("data/abc_results_500k.rds"))









# Function to calculate distance between simulated and observed data
calculate_distance <- function(simulated, observed) {
  return(mean((simulated - observed)^2)) #MSE
}


# General function for ABC fits with model flexibility
run_abc_fits <- function(data, input_layer, output_layer, simulation_function, n_prior_samples = 1000, tol = 0.01, return_dat = "train_data,test_data") {
  
  input_layer =  c(100,350,600,800,1000,1200)
  output_layer = input_layer
  
  prior_samples <- generate_prior_c_lr(n_prior_samples)
  
  # Convert the string to a function
  if (is.character(simulation_function)) {
    simulation_function <- get(simulation_function, mode = "function")
  } else {
    simulation_function <- match.fun(simulation_function)
  }
  
  
  plan(multisession)
  simulation_results <- future_map_dfc(seq_len(nrow(prior_samples)), function(idx) {
    params <- prior_samples[idx, ]
    simulation_function(data=as.data.table(data), c=params$c, lr=params$lr, input_layer=input_layer, output_layer=output_layer, return_dat = return_dat)
  }, .options = furrr_options(seed = TRUE))
  
  
  
  # Extract target data
  target_data_train_test <- data[expMode2 %in% c("Test", "Train"), ]$y
  target_data_test <- data[expMode2 == "Test", ]$y
  target_data_train <- data[expMode2 == "Train", ]$y
  
  teter_distances <- sapply(simulation_results, calculate_distance, observed = target_data_train_test)
  te_distances <- sapply(simulation_results[85:90, ], calculate_distance, observed = target_data_test)
  tr_distances <- sapply(simulation_results[1:84, ], calculate_distance, observed = target_data_train)
  
  teter_samples <- prior_samples[which(teter_distances <= tol), ]
  te_samples <- prior_samples[which(te_distances <= tol), ]
  tr_samples <- prior_samples[which(tr_distances <= tol), ]


  

  
  tibble::lst(abc_train_test = abc_train_test, abc_test = abc_test, abc_train,data=data,n_prior_samples,prior_samples,tol, simulation_function, targets=tibble::lst(target_data_train_test, target_data_test))
}



# Usage with existing data and functions
n_prior_samples <- 10000
prior_samples <- generate_prior_c_lr(n_prior_samples)
tolerance <- 180000
posterior_samples <- custom_abc_fit(data = avg_dsv, simulation_function = full_sim_exam, prior_samples, tol = tolerance)
