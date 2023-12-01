


abc_500k <- readRDS('data/model_cache/abc_group500k_06_17_54.rds')


all_models_df_500k <- purrr::imap_dfr(abc_500k, ~ {
  model_name <- gsub("abc_", "", .y)  # Prune the model name here
  relevant_lists <- .x[names(.x) %in% c("abc_train_test", "abc_test", "abc_train")]  # Filter relevant lists
  
  purrr::map2_dfr(
    names(relevant_lists), 
    relevant_lists, 
    ~ {
      data_frame <- cbind(.y$unadj.values, dist = .y$dist[.y$region]) |> as.data.frame() 
      mutate(data_frame,
             Group = sub("_.*", "", model_name),  # Extracts everything before the underscore
             Model = sub("^[^_]*_", "", model_name),  # Extracts everything after the underscore
             Fit_Method = gsub("abc_", "", .x)) |> 
        mutate(
          n_param = .y$numparam,
          numstat = .y$numstat, 
          nsamples = length(.y$dist),
          method = .y$method
        ) |> group_by(Group, Model, Fit_Method) |>
        # find min dist, and the c and lr values at that min dist
        mutate(min_dist = min(dist, na.rm = TRUE), 
               c_at_min_dist = c[which.min(dist)],
               lr_at_min_dist = lr[which.min(dist)], .before = "method")
    }
  )
})



purrr::imap_dfr(abc_500k, ~ {
  model_name <- gsub("abc_", "", .y)  # Prune the model name here
  relevant_lists <- .x[names(.x) %in% c("abc_train_test", "abc_test", "abc_train")]  # Filter relevant lists
  data=.x$data
  model_name
  purrr::map2_dfr()
})




#############

abc_loss_plot <- function(abc_result,prior_samples)
{
  pgrid <- cbind(prior_samples,dist=abc_result$dist) |> mutate(dist=sqrt(dist))
  upper_limit <- quantile(pgrid$dist, 0.95)  # 95th percentile
  pgrid <- pgrid |> filter(dist <= upper_limit)
  gg_prior <- ggplot(pgrid, aes(x = c, y = lr, color = dist)) +
    geom_point() + 
    #scale_color_viridis(limits = c(NA, upper_limit)) +
    scale_color_gradient(low = "blue", high = "red",limits = c(NA, upper_limit)) +  
    geom_contour(aes(x = c, y = lr, z = dist), color = 'white') +
    labs(title = "Point Plot of pgrid", x = "c", y = "lr") +
    theme_minimal()+ ggtitle("Prior Loss")
  
  
  post_grid <- cbind(prior_samples[abc_result$region,],dist=abc_result$dist[abc_result$region],k=abc_result$unadj.values)
  gg_post <- ggplot(post_grid, aes(x = c, y = lr, color = dist)) +
    geom_point() + 
    scale_color_gradient(low = "blue", high = "red",limits = c(NA, upper_limit)) +  
    geom_contour(aes(x = c, y = lr, z = dist), color = 'white') +
    labs(title = "Point Plot of pgrid", x = "c", y = "lr") +
    theme_minimal() + ggtitle("Posterior Loss")
  gg_prior+gg_post
  
}



# Define the function
abc_plot <- function(abc_result, data) {
  

  # Extracting posterior samples
  posterior_samples <- abc_result$unadj.values
  colnames(posterior_samples) <- c("c", "lr")
  posterior_samples_long <- tidyr::pivot_longer(as.data.frame(posterior_samples), everything())
  
  # Generate posterior density plots
  postV <- ggplot(posterior_samples_long, aes(x=value)) +
    geom_density() +
    facet_wrap(~name, scales="free") +
    theme_minimal() +
    labs(x="Value", y="Density", title="Posterior Density Plots")
  
  
  # Compute summary statistics
  sum_stats <- data.frame(
    mean = apply(posterior_samples, 2, mean),
    median = apply(posterior_samples, 2, median),
    mode = apply(posterior_samples, 2, Mode) # Ensure the 'Mode' function is defined
  )
  sum_stats <- tibble::rownames_to_column(sum_stats, var="parameter")
  print(sum_stats)
  
  # Calculate mode for parameters c and lr
  density_c <- density(posterior_samples[, "c"])
  mode_c <- density_c$x[which.max(density_c$y)]
  density_lr <- density(posterior_samples[, "lr"])
  mode_lr <- density_lr$x[which.max(density_lr$y)]
  
  print(paste0("Mode of c: ", mode_c, ", Mode of lr: ", mode_lr))
  
  # Run simulation with full_sim_exam function
  fsv <- full_sim_exam(as.data.table(data), c=mode_c, lr=mode_lr, exam.response, input_layer, output_layer, mode="ret")
  
  # Generate plot for simulation results
  plot_fsv <- fsv |> 
    tidyr::pivot_longer(y:pred, names_to="Resp", values_to="vx") |> 
    ggplot(aes(x, vx, fill=Resp, group=Resp)) +
    stat_summary(geom="bar", fun="mean", position=position_dodge()) +
    scale_fill_manual(values=col_themes$wes2) +
    scale_x_continuous(breaks=sort(unique(data$x)), labels=sort(unique(fsv$x))) +
    theme(legend.title = element_blank(), legend.position="top") +
    ggtitle("Fit to Test Only")
  
  postV + plot_fsv
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
  
  # ABC for test data only
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
