

library(dplyr)
library(purrr)
pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
walk(c("fun_alm","fun_model"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
set.seed(123)

ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
dsv <- ds |> filter(condit=="Varied")  
dsc <- ds |> filter(condit=="Constant") 

ids <- c(1,2,4,5,6,7,8, 10,11,12,13)



generate_prior_c_lr <- function(n) {
  prior_samples <- tibble(
    c = runif(n, 0.000001, 2),
    lr = runif(n, 0.000001, 3),
  )
  return(prior_samples)
}
lg_generate_prior_c_lr <- function(n) {
  prior_samples <- tibble(
    c = rlnorm(n, -5, 2),
    lr = rlnorm(n, -5, 2),
  )
  return(prior_samples)
}

n=5000
prior_samples <- generate_prior_c_lr(n) # Your existing function to generate priors

prior_samples <- lg_generate_prior_c_lr(n) # Your existing function to generate priors

mean(prior_samples$c)
median(prior_samples$c)
min(prior_samples$c)
max(prior_samples$c)

# KDE Function
kernel_density_estimate <- function(T_star, t, h) {
  K <- function(u) dnorm(u) # Gaussian kernel
  mean(sapply((t - T_star) / h, K)) / (length(T_star) * h)
}

# Silverman's Rule for Bandwidth
compute_bandwidth <- function(T_star) {
  n <- length(T_star)
  min(sd(T_star), IQR(T_star) / 1.34) * 0.9 * n^(-1/5)
}

# Update the prior_function to handle vector input
prior_function <- function(theta) {
  c_prior <- dnorm(theta[1]) # Assuming normal distribution for 'c'
  lr_prior <- dnorm(theta[2]) # Assuming normal distribution for 'lr'
  return(c_prior * lr_prior) # Combine the priors (you might need to adjust this based on your actual prior distributions)
}


metropolis_hastings <- function(current_theta, proposed_theta, current_T_star, proposed_T_star, t) {
  h_current <- compute_bandwidth(current_T_star)
  h_proposed <- compute_bandwidth(proposed_T_star)
  
  f_hat_current <- kernel_density_estimate(current_T_star, t, h_current)
  f_hat_proposed <- kernel_density_estimate(proposed_T_star, t, h_proposed)

  # Ensure non-zero and non-NaN values for f_hat_current and f_hat_proposed
  if (f_hat_current <= 0 || f_hat_proposed <= 0 || is.nan(f_hat_current) || is.nan(f_hat_proposed)) {
    return(current_theta)
  }
  
  a <- min(1, (prior_function(proposed_theta) * f_hat_proposed) / (prior_function(current_theta) * f_hat_current))
  if (is.nan(a) || is.infinite(a)) {
    return(current_theta)
  }

  if (runif(1) < a) return(proposed_theta) else return(current_theta)
}




pda_abc <- function(simulation_function, prior_samples, data, num_iterations = 5000, num_chains = 4) {
  
  input_layer =  c(100,350,600,800,1000,1200)
  output_layer = input_layer
  return_dat = "test_data"
  test_idx <- which(data$expMode2 == "Test")
  target_data_test <- data[expMode2 == "Test", ]$y

  chains <- vector("list", num_chains)
  for (chain_idx in 1:num_chains) {
    current_theta <- c(mean(prior_samples$c), mean(prior_samples$lr)) # Initialize with mean of priors
    chain <- vector("list", num_iterations)
    for (i in 1:num_iterations) {
      proposed_theta <- prior_samples[sample(1:nrow(prior_samples), 1), ] # Sample from priors
      current_T_star <- simulation_function(data, current_theta[1], current_theta[2], input_layer = input_layer, output_layer = output_layer, return_dat = return_dat) # Simulate data with current theta
      proposed_T_star <- simulation_function(data, proposed_theta$c, proposed_theta$lr, input_layer = input_layer, output_layer = output_layer, return_dat = return_dat) # Simulate data with proposed theta
      t <- target_data_test  # Your observed data
      current_theta <- metropolis_hastings(current_theta, as.numeric(proposed_theta), current_T_star, proposed_T_star, t)
      chain[[i]] <- current_theta
    }
    chains[[chain_idx]] <- chain
  }
  chains
}


data <- dsv |> filter(id == 1)
chain_list <- pda_abc(full_sim_exam, prior_samples, data, num_iterations = 5000, num_chains = 4)

# Processing the chains
chain_dfs <- lapply(chain_list, function(chain) {
  chain_df <- do.call(rbind, chain) |> as.data.frame()
  colnames(chain_df) <- c("c", "lr")
  chain_df[1000:nrow(chain_df), ] # trimming the first 1000 iterations
})
combined_chain_df <- do.call(rbind, chain_dfs) |> as.data.frame()


library(coda)
mcmc_chains <- lapply(chain_list, function(chain) {
  chain_df <- do.call(rbind, chain) |> as.data.frame()
  colnames(chain_df) <- c("c", "lr")
  mcmc(as.matrix(chain_df))
})

# Create an mcmc.list object for diagnostics
mcmc_list <- mcmc.list(mcmc_chains)

# Gelman-Rubin Diagnostic
gelman_diag <- gelman.diag(mcmc_list)
print(gelman_diag)

plot(mcmc_list)

for (param in colnames(chain_df)) {
    print(param)
 g1<- ggplot(data = bind_rows(chain_dfs, .id = 'Chain'), aes_string(x = param, group = 'Chain', color = 'Chain')) +
    geom_density() +
    labs(title = paste("Density Plot for", param))
print(g1)
}



summary_stats <- combined_chain_df %>% 
  summarise(across(everything(), list(mean = mean, sd = sd, quantile25 = ~quantile(., probs = 0.25), median = median, quantile75 = ~quantile(., probs = 0.75))))
print(summary_stats)


summary_stats_c <- summary(combined_chain_df$c)
summary_stats_lr <- summary(combined_chain_df$lr)

# Printing summary statistics
print(summary_stats_c)
print(summary_stats_lr)



generate_posterior_predictive <- function(simulation_function, combined_chain_df, data, num_samples = 5) {

  input_layer =  c(100,350,600,800,1000,1200)
  output_layer = input_layer
  return_dat = "test_data"
  test_idx <- which(data$expMode2 == "Test")
  target_data_test <- data[expMode2 == "Test", ]$y
  # Randomly select parameter sets from the posterior
  sampled_params <- combined_chain_df[sample(1:nrow(combined_chain_df), num_samples), ]
  
  # Generate simulated data for each sampled parameter set
  simulated_data_list <- lapply(1:nrow(sampled_params), function(i) {
    params <- sampled_params[i, ]
    sim_data <- simulation_function(data, params$c, params$lr, input_layer = input_layer, output_layer = output_layer, return_dat = "test_data") |> as_tibble() |>
        mutate(y=target_data_test,x=data[expMode2 == "Test", ]$x,c=params$c, lr=params$lr, sim=i ) |> rename(pred="value")
  })

  # Combine simulated data
  sim_trial_data <- bind_rows(simulated_data_list)

  sim_vx_agg <- sim_trial_data |> group_by(sim,x,c,lr) |> summarise(pred=mean(pred), y=mean(y), error=abs(y-pred)) |> ungroup() |> group_by(sim,c,lr) |>
    mutate(meanError=mean(error)) |> ungroup() |>
    mutate(rank = dense_rank(meanError)) |>
    arrange(rank)
   

  sim_rank <- sim_trial_data |> group_by(sim,c,lr) |> mutate(error=abs(y-pred)) |> summarise(mae=mean(error)) |> arrange(mae)



}

# Generate posterior predictive data
posterior_predictive_data <- generate_posterior_predictive(full_sim_exam, combined_chain_df, data, num_samples = 100)


# PDA Algorithm
# pda_abc <- function(simulation_function, prior_samples, data, num_iterations = 5000) {

#   input_layer =  c(100,350,600,800,1000,1200)
#   output_layer = input_layer
#     return_dat="test_data"
#   train_idx <- which(data$expMode2=="Train")
#   test_idx <- which(data$expMode2=="Test")
  
#   target_data_train_test <- data[expMode2 %in% c("Test", "Train"), ]$y
#   target_data_test <- data[expMode2 == "Test", ]$y
#   target_data_train <- data[expMode2 == "Train", ]$y


#   current_theta <- c(mean(prior_samples$c), mean(prior_samples$lr)) # Initialize with mean of priors
#   chain <- vector("list", num_iterations)
#   for(i in 1:num_iterations) {
#     proposed_theta <- prior_samples[sample(1:nrow(prior_samples), 1), ] # Sample from priors
#     current_T_star <- simulation_function(data, current_theta[1], current_theta[2], input_layer=input_layer, output_layer=output_layer, return_dat = return_dat) # Simulate data with current theta
#     proposed_T_star <- simulation_function(data, proposed_theta$c, proposed_theta$lr, input_layer=input_layer, output_layer=output_layer, return_dat = return_dat) # Simulate data with proposed theta
#     t <- target_data_test  # Your observed data
#     current_theta <- metropolis_hastings(current_theta, as.numeric(proposed_theta), current_T_star, proposed_T_star, t)
#     chain[[i]] <- current_theta
#   }
#   chain
# }

# data <- dsv |> filter(id==1)
# # full_sim_exam is simulation function
# chain <- pda_abc(full_sim_exam, prior_samples, data)

# chain_df <- do.call(rbind, chain) |> as.data.frame()
# colnames(chain_df) <- c("c", "lr")
# # trim first 1000 iterations of chain
# chain_df <- chain_df[1000:nrow(chain_df), ]
# summary_stats <- chain_df %>% summarise(across(everything(), list(mean = mean, median = median, sd = sd)))
# print(summary_stats)



#simulation_function=full_sim_exam
