
generate_prior_c_lr <- function(n) {
  prior_samples <- tibble(
    c = runif(n, 0.000001, 7),
    lr = runif(n, 0.000001, 7),
  )
  return(prior_samples)
}

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
