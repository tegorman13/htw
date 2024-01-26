library(dplyr)
library(purrr)
pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted, future, furrr, tictoc)
conflict_prefer_all("dplyr", quiet = TRUE)
walk(c("fun_alm","fun_model"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
set.seed(123)

ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
dsv <- ds |> filter(condit=="Varied")  
dsc <- ds |> filter(condit=="Constant") 

dsId <- ds |> select(id,condit) |> unique()

ids <- c(1,2,4,5,6,7,8, 10,11,12,13)


# generate_prior_c_lr <- function(n) {
#   prior_samples <- tibble(
#     c = runif(n, 0.000001, 2),
#     lr = runif(n, 0.000001, 3),
#   )
#   return(prior_samples)
# }

lg_generate_prior_c_lr <- function(n,cSig=2,lrSig=1) {
  prior_samples <- tibble(
   # c = extraDistr::rhnorm(n,sigma=cSig),
    c = rlnorm(n,-5,sdlog=cSig),
    lr = extraDistr::rhnorm(n,sigma=lrSig),
  )
  return(prior_samples)
}

n=5000
# prior_samples <- generate_prior_c_lr(n) 

prior_samples <- lg_generate_prior_c_lr(n,cSig=2.0,lrSig=2) 

mean(prior_samples$c)
median(prior_samples$c)
min(prior_samples$c)
max(prior_samples$c)
quantile(prior_samples$c)
plot(density(prior_samples$c))

mean(prior_samples$lr)
median(prior_samples$lr)
min(prior_samples$lr)
max(prior_samples$lr)
plot(density(prior_samples$lr))


# Update the prior_function to handle vector input
prior_function <- function(theta) {
  c_prior <- dlnorm(theta$c, -5,3) 
  lr_prior <- dnorm(theta$lr, 2,1) 
  return(c_prior * lr_prior) 
}

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
  target_data_test <- data[data$expMode2 == "Test", ]$y
  data <- data |> as.data.table()

  print(data$id[1])
  print(head(data))
  chains <- vector("list", num_chains)
  for (chain_idx in 1:num_chains) {

    # current_theta <- c(mean(prior_samples$c), mean(prior_samples$lr)) # Initialize with mean of priors
    current_theta <- prior_samples[sample(1:nrow(prior_samples), 1), ]
    chain <- vector("list", num_iterations)

    for (i in 1:num_iterations) {
      proposed_theta <- prior_samples[sample(1:nrow(prior_samples), 1), ] # Sample from priors
      current_T_star <- simulation_function(data, current_theta$c, current_theta$c, input_layer = input_layer, output_layer = output_layer, return_dat = return_dat) # Simulate data with current theta
      proposed_T_star <- simulation_function(data, proposed_theta$c, proposed_theta$lr, input_layer = input_layer, output_layer = output_layer, return_dat = return_dat) # Simulate data with proposed theta
      t <- target_data_test  # Your observed data
      current_theta <- metropolis_hastings(current_theta, proposed_theta, current_T_star, proposed_T_star, t)
      # if i is divisble by 100, print current_theta
      if (i %% 200 == 0) print(current_theta)
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


# data <- dsv |> filter(id == 1)
# chain_list <- pda_abc(full_sim_exam, prior_samples, data, num_iterations = 4000, num_chains = 2)


# chains_df <- imap_dfr(chain_list, ~mutate(.x, chain = .y))
# chain_dfs2 <- split(chains_df |> select(-chain), f = chains_df$chain)
# mcmc_list <- lapply(chain_dfs2, as.mcmc)
# (gelman_result <- gelman.diag(mcmc_list))
# plot(mcmc.list(mcmc_list))

# simulation_function <- full_sim_exam
# prior_samples <- lg_generate_prior_c_lr(n=1000) 
# num_iterations = 100
# num_chains = 4

# posterior_predictive_data <- generate_posterior_predictive(full_sim_exam, chains_by_sbj[[2]], data, num_samples = 100)

# ids1 <- c(1,36,66, 76,101,192)


ids1 <- as.numeric(levels(ds$id))
prior_samples <- lg_generate_prior_c_lr(n=5000) 
subjects_data <-  ds |> filter(id %in% ids1)  %>% split(f =c(.$id), drop=FALSE)

# chains_by_sbj <- map(subjects_data, ~pda_abc(simulation_function = full_sim_exam, 
#                                                  prior_samples = prior_samples, 
#                                                  data = .x, 
#                                                  num_iterations = 10, 
#                                                  num_chains = 1)) %>% # 
#                                                  setNames(ids1)
                                                

(nc <- future::availableCores())
future::plan(future::cluster, workers = nc-1)

tic()
chains_by_sbj <- future_map(subjects_data, ~pda_abc(simulation_function = full_sim_exam, 
                                                    prior_samples = prior_samples, 
                                                    data = .x, 
                                                    num_iterations = 2000, 
                                                    num_chains = 4), .options = furrr_options(seed = TRUE)) %>% setNames(ids1)
                                                     

toc()


plan(sequential)

# 3857s (1 hour 4 min) for 5000 iter; 4 chains on windows_3070. 
#2238 for 3000 iter, 4 chain on 3070
#887 secs for 2K iterations, 4 chains on tg_m1



tic()
k=rnorm(1000)
Sys.sleep(1)
t1=toc(log=TRUE)



map(chains_by_sbj, ~{
  chain_dfs2 <- split(.x |> select(c,lr), f = .x$chain)
  mcmc_list <- lapply(chain_dfs2, as.mcmc)
  plot(mcmc.list(mcmc_list))
} )


map(chains_by_sbj, possibly(~{
  nr = nrow(.x)
  trim = .x |> group_by(chain) |> filter(row_number() < nr/2) |> ungroup()  
  chain_dfs2 <- split(trim |> select(c,lr), f = trim$chain)
  mcmc_list <- lapply(chain_dfs2, as.mcmc, start=nr/2,thin=4)
  (gelman_result <- gelman.diag(mcmc_list))
} ))


map(chains_by_sbj, ~{
  .x |> select(c,lr) |> head(2)
} )


map(chains_by_sbj, ~{
  nr = nrow(.x)
  trim=.x[nr/2:nr, ]
  paste0("c: ", mean(trim$c), " lr: ", mean(trim$lr))
} )


map(chains_by_sbj, ~{
  unique(.x$chain) %>% print()
  nr = nrow(.x)
  trim = .x |> group_by(chain) |> filter(row_number() < nr/2) |> ungroup()
  
  unique(trim$chain) %>% print()
  { trim |> ggplot(aes(x=c,fill=as.factor(chain))) +geom_density() } +
  trim |> ggplot(aes(x=lr,fill=as.factor(chain))) +geom_density() 
  
} )



map(chains_by_sbj, ~{
  unique(.x$chain) %>% print()
  nr = nrow(.x)
  trim = .x |> group_by(chain) |> filter(row_number() < nr/2) |> ungroup()  
  posterior_predictive_data <- generate_posterior_predictive(full_sim_exam, .x, data, num_samples = 100)
  head(posterior_predictive_data$sim_rank)
  head(posterior_predictive_data$sim_vx_agg)
} )

chains_by_sbj = exam_test
chains_by_sbj = alm_test
chains_by_sbj |> tidyselect:::select("1") %>% map(., ~{
  trim = .x |> group_by(chain) |> filter(row_number() < (nrow(.x))/1.5) |> ungroup() 
  data <- ds |> filter(id == trim$id[1])
  posterior_predictive_data <- generate_posterior_predictive(full_sim_exam, trim, data, num_samples = 100)
   print({ trim |> ggplot(aes(x=c,fill=as.factor(chain))) +geom_density() } +
  trim |> ggplot(aes(x=lr,fill=as.factor(chain))) +geom_density(alpha=.5) )
  print(head(posterior_predictive_data$sim_rank))
  head(posterior_predictive_data$sim_vx_agg)
} )




k = rbindlist(chains_by_sbj)|> left_join(dsId, join_by(id))


k |> filter(c<.1) |> ggplot(aes(x=c,fill=condit)) + geom_density()
k  |> ggplot(aes(x=lr,fill=condit)) + geom_density()

k  |> filter(rank<10) |>  ggplot(aes(x=condit,y=lr)) + geom_col()

posterior_predictive_data <- generate_posterior_predictive(full_sim_exam, chains_by_sbj[[2]], data, num_samples = 100)




mean(chains_by_sbj[[1]]$c)
mean(chains_by_sbj[[2]]$c)


map(chains_by_sbj, ~{
  head(.x,5)

} )


# Processing the chains
chain_dfs <- lapply(chain_list, function(chain) {
  chain_df <- do.call(rbind, chain) |> as.data.frame()
  colnames(chain_df) <- c("c", "lr")
  chain_df[100:nrow(chain_df), ] # trimming the first 1000 iterations
})
combined_chain_df <- do.call(rbind, chain_dfs) |> as.data.frame()

#do.call(rbind, chain_list[[1]]) |> as.data.frame()


library(coda)
mcmc_chains <- lapply(chain_list, function(chain) {
  chain_df <- do.call(rbind, chain) |> as.data.frame()
  colnames(chain_df) <- c("c", "lr")
  chain_df[1000:nrow(chain_df), ] # trimming the first 1000 iterations
  mcmc(as.matrix(chain_df))
})

# Create an mcmc.list object for diagnostics
mcmc_list <- mcmc.list(mcmc_chains)

# Gelman-Rubin Diagnostic
gelman_diag <- gelman.diag(mcmc_list)
print(gelman_diag)

plot(mcmc_list)

for (param in colnames(combined_chain_df)) {
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
  target_data_test <- data[data$expMode2 == "Test", ]$y
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

  tibble::lst(sim_vx_agg,sim_rank)

}

# Generate posterior predictive data
posterior_predictive_data <- generate_posterior_predictive(full_sim_exam, combined_chain_df, data, num_samples = 100)


