
abc_500k <- readRDS('data/model_cache/abc_group100k_06_39_17.rds')
l1 = abc_500k$abc_v_exam 
k=abc_500k$abc_v_exam$abc_test
names(k)
k$ss[k$region]

post_grid <- tibble(l1$prior_samples[k$region,],dist=k$dist[k$region],k2=k$unadj.values,k$ss) 

post_grid <- tibble(l1$prior_samples[k$region,],dist=k$dist[k$region]) 


post_params <- l1$prior_samples[k$region,][1:10,]
post_params <- post_grid[1:5,]
simulation_results <- pmap_dfr(post_params, ~ {
  params <- tibble::lst(..1,..2)
  l1$simulation_function(l1$data, c=params[[1]], lr=params[[2]], 
                         input_layer=input_layer, output_layer=output_layer, return_dat = "test_data",mode="f") |> 
    mutate(c = params[[1]], lr = params[[2]], dev=y-pred, dev2=dev^2, md=sqrt(sum(dev2)),abc_dist=..3)
}) 
simulation_results



##### Test distance 

sim_data=full_sim_exam(avg_dsv, .03, 2,exam.response,input_layer, output_layer,return_dat = "test_data")


normalise <- function(x,y){
  if(mad(y) == 0)
    return (x)
  else
    return (x/mad(y))
}

calculate_distances <- function(target, sumstat) {
  # Normalize target and sumstat
  nss <- ncol(sumstat)
  for (j in 1:nss) {
    sumstat[, j] <- normalise(sumstat[, j], sumstat[, j])
    target[j] <- normalise(target[j], sumstat[, j])
  }
  
  # Calculate Euclidean distances
  dist <- apply(sumstat, 1, function(row) {
    sqrt(sum((row - target)^2))
  })
  
  return(dist)
}


library(abc)
prior_samples <- generate_prior_c_lr(200)
sim_data <- pmap_dfc(prior_samples, ~ {
  full_sim_exam(avg_dsv, c=..1, lr=..2, input_layer=input_layer,output_layer=output_layer)})
  
target_data_test <- avg_dsv[expMode2 == "Test", ]$y

abc_test <- abc(
  target = target_data_test,
  param = prior_samples,
  sumstat = do.call(rbind, sim_data),
  tol = .1,
  method = "rejection"
)
abc_test$ss
abc_test$dist
post_grid <- tibble(prior_samples[abc_test$region,],dist=abc_test$dist[abc_test$region]) 
post_grid

simulation_results <- pmap_dfr(post_grid, ~ {
  full_sim_exam(avg_dsv, c=..1, lr=..2, input_layer=input_layer,output_layer=output_layer,mode="f") |>
    mutate(c = ..1, lr = ..2, dev=y-pred, dev2=dev^2, abc_dist=..3) |> group_by(c,lr, abc_dist) |>
          summarize(manual_dist=sqrt(sum(dev2)), .groups="keep") 
}) 
simulation_results



target_data <- target_data_test  # Your target data
simulated_data <- do.call(rbind, sim_data)  # Your simulated data
distances <- calculate_distances(target_data, simulated_data)
distances[abc_test$region]




#################



prior_samples <- generate_prior_c_lr(200)
sim_data <- pmap_dfc(prior_samples, ~ {
  full_sim_exam(avg_dsv, c=..1, lr=..2, input_layer=input_layer, output_layer=output_layer)
})
target_data_test <- avg_dsv[expMode2 == "Test", ]$y

# Normalizing target and simulated data
normalized_target <- normalise(target_data_test, target_data_test)
normalized_sim_data <- apply(do.call(rbind, sim_data), 2, function(x) normalise(x, x))

# Running the abc function
abc_test <- abc(
  target = normalized_target,
  param = prior_samples,
  sumstat = normalized_sim_data,
  tol = .1,
  method = "rejection"
)

# Extracting distances from abc_test
abc_test$dist

# Manual computation of Euclidean distances
post_grid <- tibble(prior_samples[abc_test$region,], dist=abc_test$dist[abc_test$region]) 
simulation_results <- pmap_dfr(post_grid, ~ {
  sim <- full_sim_exam(avg_dsv, c=..1, lr=..2, input_layer=input_layer, output_layer=output_layer)
  sim_normalized <- normalise(avg_dsv$y,sim)
  dev2 <- (sim_normalized - normalized_target)^2
  manual_dist <- sqrt(sum(dev2))
  tibble(c = ..1, lr = ..2, manual_dist, abc_dist = ..3)
})

simulation_results
