library(dplyr)
library(purrr)
pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted, future, furrr, tictoc,extraDistr)
conflict_prefer_all("dplyr", quiet = TRUE)
walk(c("fun_alm","fun_model"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
set.seed(123)
ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()


lg_generate_prior_c_lr <- function(n,cMean=-5,cSig=2,lrSig=1) {
  prior_samples <- tibble(
    # c = extraDistr::rhnorm(n,sigma=cSig),
    c = rlnorm(n,cMean,sdlog=cSig),
    lr = extraDistr::rhnorm(n,sigma=lrSig),
  )
  return(prior_samples)
}

prior_function <- function(theta) {
  c_prior <- dlnorm(theta$c, -5,2) 
  lr_prior <- extraDistr::rhnorm(theta$lr,sigma=lrSig) 
  return(c_prior * lr_prior) 
}

kernel_density_estimate <- function(T_star, t, h) {
  K <- function(u) dnorm(u) # Gaussian kernel
  mean(sapply((t - T_star) / h, K)) / (length(T_star) * h)
}

# Silverman's Rule for Bandwidth
compute_bandwidth <- function(T_star) {
  n <- length(T_star)
  min(sd(T_star), IQR(T_star) / 1.34) * 0.9 * n^(-1/5)
}


metropolis_hastings <- function(current_theta, proposed_theta, current_T_star, proposed_T_star, t) {
  h_current <- compute_bandwidth(current_T_star)
  h_proposed <- compute_bandwidth(proposed_T_star)
  
  f_hat_current <- kernel_density_estimate(current_T_star, t, h_current)
  f_hat_proposed <- kernel_density_estimate(proposed_T_star, t, h_proposed)

  if (f_hat_current <= 0 || f_hat_proposed <= 0 || is.nan(f_hat_current) || is.nan(f_hat_proposed)) {
    return(current_theta)
  }
  
  a <- min(1, (prior_function(proposed_theta) * f_hat_proposed) / (prior_function(current_theta) * f_hat_current))

  if (is.nan(a) || is.infinite(a)) {
    return(current_theta)
  }

  if (runif(1) < a) return(proposed_theta) else return(current_theta)
}


pda_abc <- function(simulation_function, prior_samples, data, num_iterations = 5000, num_chains = 4, return_dat="test_data") {
  input_layer =  c(100,350,600,800,1000,1200)
  output_layer = input_layer
  #return_dat = "test_data"
  test_idx <- which(data$expMode2 == "Test")
  target_data_test <- data[data$expMode2 == "Test", ]$y
  target_data_train <- data[data$expMode2 == "Train", ]$y
  target_data_train_test <- data[expMode2 %in% c("Test", "Train"), ]$y

 target_data <- case_when(
  return_dat == "test_data" ~ list(target_data_test), 
  return_dat == "train_data" ~ list(target_data_train),
  return_dat == "train_data, test_data" ~ list(target_data_train_test)
) |> unlist() |> as.numeric()


  data <- data |> as.data.table()

#   print(data$id[1])
#   print(head(data))
  chains <- vector("list", num_chains)
  for (chain_idx in 1:num_chains) {

    # current_theta <- c(mean(prior_samples$c), mean(prior_samples$lr)) # Initialize with mean of priors
    current_theta <- prior_samples[sample(1:nrow(prior_samples), 1), ]
    chain <- vector("list", num_iterations)

    for (i in 1:num_iterations) {
      proposed_theta <- prior_samples[sample(1:nrow(prior_samples), 1), ] # Sample from priors
      current_T_star <- simulation_function(data, current_theta$c, current_theta$c, input_layer = input_layer, output_layer = output_layer, return_dat = return_dat) # Simulate data with current theta
      proposed_T_star <- simulation_function(data, proposed_theta$c, proposed_theta$lr, input_layer = input_layer, output_layer = output_layer, return_dat = return_dat) # Simulate data with proposed theta
      t <- target_data  # Your observed data
      current_theta <- metropolis_hastings(current_theta, proposed_theta, current_T_star, proposed_T_star, t)
      # if i is divisble by 100, print current_theta
     # if (i %% 1000 == 0) print(current_theta)
      chain[[i]] <- current_theta
    }
    chains[[chain_idx]] <- chain
  }
  #chains
  chain_dfs <- lapply(chains, function(chain) {
  chain_df <- do.call(rbind, chain) |> as.data.frame()
  colnames(chain_df) <- c("c", "lr")
  chain_df
})

chains_df <- imap_dfr(chain_dfs, ~mutate(.x, chain = .y)) |> mutate(id=data$id[1])
}

####################################
####################################


# ids1 <- 1
ids1 <- c(1,33,66)
#ids1 <- as.numeric(levels(ds$id))

cMean <<- -5; cSig <<- 4; lrSig <<- 3
prior_samples <- lg_generate_prior_c_lr(n=5000, cMean=cMean, cSig=cSig, lrSig=lrSig) 
subjects_data <-  ds |> filter(id %in% ids1)  %>% split(f =c(.$id), drop=TRUE)

num_iterations = 400
num_chains = 2

save_folder <- paste0("n_iter_",num_iterations,"_nc_",num_chains,"_",format(Sys.time(),"%H%M%OS"))
dir.create(paste0("data/abc_pda/",save_folder))




### EXAM Test 
(nc <- future::availableCores())
future::plan(multicore, workers = nc-1)
t1=system.time({
exam_test <- future_map(subjects_data, ~pda_abc(simulation_function = full_sim_exam, 
                                                    prior_samples = prior_samples, 
                                                    data = .x, 
                                                    num_iterations = num_iterations, 
                                                    num_chains = num_chains,
                                                    return_dat="test_data"), .options = furrr_options(seed = TRUE)) %>% 
                                                    setNames(ids1)
                                                  
})
print(t1[3])

plan(sequential)
ri=ri_pda_indv() %>% append(.,t1[3])
run_save <- tibble::lst(exam_test,Model="EXAM",Fit_Method="Test",prior_samples,cMean,cSig,lrSig,ri)
file_name <- paste0("data/abc_pda/",save_folder,"/","pda_",run_save$Model,"_",run_save$Fit_Method,"_",
                    num_iterations,"_",num_chains,"_",format(Sys.time(),"%H%M%OS"), ".rds")
saveRDS(run_save,file=here::here(file_name))



#### ALM Test

(nc <- future::availableCores())
future::plan(multicore, workers = nc-1)
t1=system.time({
  alm_test <- future_map(subjects_data, ~pda_abc(simulation_function = full_sim_alm, 
                                                  prior_samples = prior_samples, 
                                                  data = .x, 
                                                  num_iterations = num_iterations, 
                                                  num_chains = num_chains,
                                                  return_dat="test_data"), .options = furrr_options(seed = TRUE)) %>% 
    setNames(ids1)
  
})
print(t1[3])

plan(sequential)
ri=ri_pda_indv() %>% append(.,t1[3])
run_save <- tibble::lst(alm_test,Model="ALM",Fit_Method="Test",prior_samples,cMean,cSig,lrSig,ri)
file_name <- paste0("data/abc_pda/",save_folder,"/","pda_",run_save$Model,"_",run_save$Fit_Method,"_",
                    num_iterations,"_",num_chains,"_",format(Sys.time(),"%H%M%OS"), ".rds")
saveRDS(run_save,file=here::here(file_name))



#### EXAM Test & Train

(nc <- future::availableCores())
future::plan(multicore, workers = nc-1)
t1=system.time({
  exam_test_train <- future_map(subjects_data, ~pda_abc(simulation_function = full_sim_exam, 
                                                  prior_samples = prior_samples, 
                                                  data = .x, 
                                                  num_iterations = num_iterations, 
                                                  num_chains = num_chains,
                                                  return_dat="train_data, test_data"), .options = furrr_options(seed = TRUE)) %>% 
    setNames(ids1)
  
})
print(t1[3])

plan(sequential)
ri=ri_pda_indv() %>% append(.,t1[3])
run_save <- tibble::lst(exam_test_train,Model="EXAM",Fit_Method="Test_Train",prior_samples,cMean,cSig,lrSig,ri)
file_name <- paste0("data/abc_pda/",save_folder,"/","pda_",run_save$Model,"_",run_save$Fit_Method,"_",
                    num_iterations,"_",num_chains,"_",format(Sys.time(),"%H%M%OS"), ".rds")
saveRDS(run_save,file=here::here(file_name))




#### ALM Test & Train

(nc <- future::availableCores())
future::plan(multicore, workers = nc-1)
t1=system.time({
  alm_test_train <- future_map(subjects_data, ~pda_abc(simulation_function = full_sim_alm, 
                                                  prior_samples = prior_samples, 
                                                  data = .x, 
                                                  num_iterations = num_iterations, 
                                                  num_chains = num_chains,
                                                  return_dat="train_data, test_data"), .options = furrr_options(seed = TRUE)) %>% 
    setNames(ids1)
  
})
print(t1[3])

plan(sequential)
ri=ri_pda_indv() %>% append(.,t1[3])
run_save <- tibble::lst(alm_test_train,Model="ALM",Fit_Method="Test_Train",prior_samples,cMean,cSig,lrSig,ri)
file_name <- paste0("data/abc_pda/",save_folder,"/","pda_",run_save$Model,"_",run_save$Fit_Method,"_",
                    num_iterations,"_",num_chains,"_",format(Sys.time(),"%H%M%OS"), ".rds")
saveRDS(run_save,file=here::here(file_name))




#### EXAM Train

(nc <- future::availableCores())
future::plan(multicore, workers = nc-1)
t1=system.time({
  exam_train <- future_map(subjects_data, ~pda_abc(simulation_function = full_sim_exam, 
                                                  prior_samples = prior_samples, 
                                                  data = .x, 
                                                  num_iterations = num_iterations, 
                                                  num_chains = num_chains,
                                                  return_dat="train_data"), .options = furrr_options(seed = TRUE)) %>% 
    setNames(ids1)
  
})
print(t1[3])

plan(sequential)
ri=ri_pda_indv() %>% append(.,t1[3])
run_save <- tibble::lst(exam_train,Model="EXAM",Fit_Method="Train",prior_samples,cMean,cSig,lrSig,ri)
file_name <- paste0("data/abc_pda/",save_folder,"/","pda_",run_save$Model,"_",run_save$Fit_Method,"_",
                    num_iterations,"_",num_chains,"_",format(Sys.time(),"%H%M%OS"), ".rds")
saveRDS(run_save,file=here::here(file_name))





#### ALM Train

(nc <- future::availableCores())
future::plan(multicore, workers = nc-1)
t1=system.time({
  alm_train <- future_map(subjects_data, ~pda_abc(simulation_function = full_sim_alm, 
                                                       prior_samples = prior_samples, 
                                                       data = .x, 
                                                       num_iterations = num_iterations, 
                                                       num_chains = num_chains,
                                                       return_dat="train_data"), .options = furrr_options(seed = TRUE)) %>% 
    setNames(ids1)
  
})
print(t1[3])

plan(sequential)
ri=ri_pda_indv() %>% append(.,t1[3])
run_save <- tibble::lst(alm_train,Model="ALM",Fit_Method="Train",prior_samples,cMean,cSig,lrSig,ri)
file_name <- paste0("data/abc_pda/",save_folder,"/","pda_",run_save$Model,"_",run_save$Fit_Method,"_",
                    num_iterations,"_",num_chains,"_",format(Sys.time(),"%H%M%OS"), ".rds")
saveRDS(run_save,file=here::here(file_name))

