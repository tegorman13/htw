pacman::p_load(dplyr,purrr,tidyr,ggplot2, data.table, here, patchwork, conflicted, future, furrr, tictoc,extraDistr)
conflict_prefer_all("dplyr", quiet = TRUE)
walk(c("fun_alm","fun_model"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
set.seed(123)
ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()


lg_generate_prior_c_lr <- function(n,cMean=-5,cSig=2,lrSig=1) {
  prior_samples <- tibble(
    #c = extraDistr::rhnorm(n,sigma=cSig),
     #c = runif(n,.000001,1),
     #lr=runif(n, .001, 6)
    c = rlnorm(n,cMean,sdlog=cSig),
    lr = extraDistr::rhnorm(n,sigma=lrSig),
  )
  return(prior_samples)
}

prior_function <- function(theta) {
  log_c_prior <- (dlnorm(theta$c, cMean, cSig))
 # log_c_prior = extraDistr::dhnorm(theta$c,sigma=cSig)
 # log_c_prior <- (dunif(theta$c, .000001, 1))
  #log_lr_prior <-(dunif(theta$lr, .001, 6))
  log_lr_prior <- (extraDistr::dhnorm(theta$lr, sigma=lrSig))
  return(log_c_prior + log_lr_prior)
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


# metropolis_hastings <- function(current_theta, proposed_theta, current_T_star, proposed_T_star, t) {
#   h_current <- compute_bandwidth(current_T_star)
#   h_proposed <- compute_bandwidth(proposed_T_star)
#   
#   f_hat_current <- kernel_density_estimate(current_T_star, t, h_current)
#   f_hat_proposed <- kernel_density_estimate(proposed_T_star, t, h_proposed)
# 
#   if (f_hat_current <= 0 || f_hat_proposed <= 0 || is.nan(f_hat_current) || is.nan(f_hat_proposed)) {
#     warning("Non-positive density estimate encountered.")
#     return(current_theta)
#   }
#   
#   a <- min(1, (prior_function(proposed_theta) * f_hat_proposed) / (prior_function(current_theta) * f_hat_current))
# 
#   if (is.nan(a) || is.infinite(a)) {
#     return(current_theta)
#   }
# 
#   if (runif(1) < a) return(proposed_theta) else return(current_theta)
# }

# metropolis_hastings <- function(current_theta, proposed_theta, current_T_star, proposed_T_star, t, acceptance_factor = 2.5, acceptance_threshold = 0.2) {
#   h_current <- compute_bandwidth(current_T_star)
#   h_proposed <- compute_bandwidth(proposed_T_star)
#   
#   f_hat_current <- kernel_density_estimate(current_T_star, t, h_current)
#   f_hat_proposed <- kernel_density_estimate(proposed_T_star, t, h_proposed)
#   
#   if (f_hat_current <= 0 || f_hat_proposed <= 0 || is.nan(f_hat_current) || is.nan(f_hat_proposed)) {
#     return(current_theta)
#   }
#   
#   ratio <- (prior_function(proposed_theta) * f_hat_proposed) / (prior_function(current_theta) * f_hat_current)
#   a <- min(acceptance_factor, ratio)
#   
#   if (is.nan(a) || is.infinite(a)) {
#     return(current_theta)
#   }
#   
#   if (runif(1) < a || ratio > acceptance_threshold) {
#     return(proposed_theta)
#   } else {
#     return(current_theta)
#   }
# }



metropolis_hastings <- function(current_theta, proposed_theta, current_T_star, proposed_T_star, t, iteration) {
  h_current <- compute_bandwidth(current_T_star)
  h_proposed <- compute_bandwidth(proposed_T_star)
  
  f_hat_current <- kernel_density_estimate(current_T_star, t, h_current)
  f_hat_proposed <- kernel_density_estimate(proposed_T_star, t, h_proposed)
  
  if (f_hat_current <= 0 || f_hat_proposed <= 0 || is.nan(f_hat_current) || is.nan(f_hat_proposed)) {
    return(list(theta = current_theta, accepted = FALSE))
  }
  
  prior_current <- prior_function(current_theta)
  prior_proposed <- prior_function(proposed_theta)
  a <- exp(prior_proposed + log(f_hat_proposed) - prior_current - log(f_hat_current))
  
  if (is.nan(a) || is.infinite(a)) {
    return(list(theta = current_theta, accepted = FALSE))
  }
  
  if (runif(1) < a) {
    acceptance_rate[iteration] <<- acceptance_rate[iteration] + 1
    return(list(theta = proposed_theta, accepted = TRUE))
  } else {
    return(list(theta = current_theta, accepted = FALSE))
  }
}

# Now you can adjust the proposal distribution based on the acceptance rate
adjust_proposal <- function(iteration) {
  # Only adjust every 100 iterations, for example
  if (iteration %% 100 == 0) {
    acceptance <- mean(acceptance_rate[(iteration - 99):iteration])
    if (acceptance < 0.2) {
      proposal_width_c <<- proposal_width_c * 0.9
      proposal_width_lr <<- proposal_width_lr * 0.9
    } else if (acceptance > 0.5) {
      proposal_width_c <<- proposal_width_c * 1.1
      proposal_width_lr <<- proposal_width_lr * 1.1
    }
  }
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
      adjust_proposal(i)
      if (!is.numeric(proposal_width_c) || proposal_width_c <= 0) {
        stop("proposal_width_c is not a positive number")
      }
      if (!is.numeric(proposal_width_lr) || proposal_width_lr <= 0) {
        stop("proposal_width_lr is not a positive number")
      }
      
      print(paste("proposal_width_c:", proposal_width_c))
      print(paste("proposal_width_lr:", proposal_width_lr))
      cat("Calling rnorm with mean:", current_theta$c, "sd:", proposal_width_c, "\n")
      cat("Calling rnorm with mean:", current_theta$lr, "sd:", proposal_width_lr, "\n")
      
      proposed_theta <- list(
        c = rnorm(1, mean = current_theta$c, sd = proposal_width_c),
        lr = rnorm(1, mean = current_theta$lr, sd = proposal_width_lr)
      )
      #proposed_theta <- prior_samples[sample(1:nrow(prior_samples), 1), ] # Sample from priors
      current_T_star <- simulation_function(data, current_theta$c, current_theta$c, input_layer = input_layer, output_layer = output_layer, return_dat = return_dat) # Simulate data with current theta
      proposed_T_star <- simulation_function(data, proposed_theta$c, proposed_theta$lr, input_layer = input_layer, output_layer = output_layer, return_dat = return_dat) # Simulate data with proposed theta
      t <- target_data  # Your observed data
      current_theta <- metropolis_hastings(current_theta, proposed_theta, current_T_star, proposed_T_star, t)
      # if i is divisble by 100, print current_theta
      #if (i %% 1000 == 0) print(current_theta)
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
#ids1 <- c(49)

cMean <<- -5.5; cSig <<- 2.0; lrSig <<- 2.0

# You may want to keep track of the acceptance to tune the proposal distribution
acceptance_rate <- numeric(num_iterations)
proposal_width_c <- cSig
proposal_width_lr <- lrSig

prior_samples <- lg_generate_prior_c_lr(n=6000, cMean=cMean, cSig=cSig, lrSig=lrSig) 
subjects_data <-  ds |> filter(id %in% ids1)  %>% split(f =c(.$id), drop=TRUE)

num_iterations = 1000
num_chains = 2

# save_folder <- paste0("n_iter_",num_iterations,"_nc_",num_chains,"_",format(Sys.time(),"%H%M%OS"))
# dir.create(paste0("data/abc_pda/",save_folder))


### EXAM Test 
(nc <- future::availableCores())
future::plan(multisession, workers = nc)
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





