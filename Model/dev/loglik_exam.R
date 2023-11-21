



# Log likelihood loss function
log_likelihood_loss <- function(prediction, target, noise) {
  # Assuming a normal distribution of predictions with a given noise (standard deviation)
  ll <- sum(dnorm(target, mean = prediction, sd = noise, log = TRUE))
  return(-ll)  # We return negative log likelihood to minimize it during optimization
}

compute_loss <- function(train_results, train_data, output_layer) {
  # Extract the output probabilities from train_results
  output_probs <- train_results$d[, grep("almResp", colnames(train_results$d))]
  
  # Compute the loss for each observation
  losses <- mapply(function(probs, y, output_layer) {
    # Find the index of the output node closest to the observed value y
    closest_node_index <- which.min(abs(output_layer - y))
    # Compute the negative log likelihood for the observation
    # Ensure probabilities are not zero to avoid -Inf values
    safe_prob <- pmax(probs[closest_node_index], 1e-10)
    -log(safe_prob)  # Apply the negative sign here
  }, probs = as.data.frame(t(output_probs)), y = train_data$y, MoreArgs = list(output_layer = output_layer))
  
  # Calculate the average loss to get the negative log likelihood per observation
  avg_loss <- mean(losses, na.rm = TRUE)
  return(avg_loss)
}


compute_test_log_likelihood <- function(pred_dat_y, test_prediction, output_layer=output.layer) {
  # Find the index of the output node closest to each empirical mean
  closest_node_indices <- sapply(pred_dat_y, function(y) which.min(abs(output_layer - y)))
  
  predicted_probs <- mapply(function(row_index, col_index) test_prediction[row_index, col_index],
                            row_index = seq_along(closest_node_indices),
                            col_index = closest_node_indices)

  # Compute the log likelihood
  log_likelihood <- sum(log(pmax(predicted_probs, 1e-10)))

  # Calculate the average log likelihood to get the negative log likelihood per observation
  avg_log_likelihood <- log_likelihood / length(pred_dat_y)

  return(avg_log_likelihood)
}

# Updated ALM response function to return probabilities instead of mean response
alm.responseOnly <- function(input = 1, c, input.layer, output.layer, weight.mat, noise, trainVec=NULL) {
  input.activation <- exp(-c * (input.layer - input)^2)
  output.activation <- weight.mat %*% input.activation
  output.probability <- output.activation / sum(output.activation)
  return(output.probability)
}


# ALM response function returning output probabilities
alm.response <- function(input = 1, c, input.layer, output.layer, weight.mat, noise) {
  input.activation <- exp(-c * (input.layer - input)^2)
  output_activation <- weight.mat %*% input.activation
  
  # Add Gaussian noise to the output activation
  noise_vector <- rnorm(length(output_activation), mean = 0, sd = noise)
  output_activation <- output_activation + noise_vector
  
  output_activation <- pmax(output_activation, 1e-10)
  # Normalize the activations to probabilities
  output.probability <- output_activation / sum(output_activation)
  
  list(output.probability = output.probability, input.activation = input.activation, output.activation = output_activation)
}

# ALM update function remains the same as it updates the weights based on the feedback signal
alm.update <- function(corResp, c, lr, output.layer, input.activation, output.activation, weight.mat) {
  fz <- exp(-c * (output.layer - corResp)^2)
  teacherSignal <- (fz - output.activation) * lr
  wChange <- teacherSignal %*% t(input.activation)
  weight.mat <- weight.mat + wChange
  return(weight.mat)
}

# ALM trial function now returns output probabilities instead of mean response
alm.trial <- function(input, corResp, c, lr, input.layer, output.layer, weight.mat, noise) {
  alm_resp <- alm.response(input, c, input.layer, output.layer, weight.mat, noise)
  updated_weight.mat <- alm.update(corResp, c, lr, output.layer, alm_resp$input.activation, alm_resp$output.activation, weight.mat)
  return(list(output.probability = alm_resp$output.probability, weight.mat = updated_weight.mat))
}

# ALM simulation function now returns output probabilities instead of mean responses
alm.sim <- function(dat, c, lr, input.layer, output.layer, noise) {
  weight.mat <- matrix(0.000001, nrow = length(output.layer), ncol = length(input.layer))
  n <- nrow(dat)
  output_probs <- matrix(NA, nrow = n, ncol = length(output.layer)) # Matrix to store output probabilities
  for(i in 1:n) {
    trial <- alm.trial(dat$x[i], dat$y[i], c, lr, input.layer, output.layer, weight.mat, noise)
    weight.mat <- trial$weight.mat
    output_probs[i, ] <- trial$output.probability
  }
  dat <- cbind(dat, almResp=output_probs)
  return(list(d = dat, wm = weight.mat, c = c, lr = lr))
}

# EXAM response function now returns output probabilities instead of mean response
exam.response <- function(input, c, input.layer, output.layer, weight.mat, trainVec, noise) {
  nearestTrain <- trainVec[which.min(abs(input - trainVec))]
  output_probs <- alm.response(nearestTrain, c, input.layer, output.layer, weight.mat, noise)$output.probability
  
  # Compute the slope for EXAM generalization
  xUnder <- ifelse(min(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) - 1])
  xOver <- ifelse(max(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) + 1])
  
  probsUnder <- alm.response(xUnder, c, input.layer, output.layer, weight.mat, noise)$output.probability
  probsOver <- alm.response(xOver, c, input.layer, output.layer, weight.mat, noise)$output.probability
  
  # Linear interpolation for the probabilities
  slope <- (probsOver - probsUnder) / (xOver - xUnder)
  exam_probs <- output_probs + slope * (input - nearestTrain)
  
  # Ensure probabilities are non-negative and sum to 1
  exam_probs <- pmax(exam_probs, 0)
  exam_probs <- exam_probs / sum(exam_probs)
  
  return(exam_probs)
}




input.layer <- c(100,350,600,800,1000,1200)
output.layer <- c(100,350,600,800,1000,1200)

noise_values <- seq(0.0000001, .01, length.out=30)  # Adjust the range as needed
c_values <- seq(0.0000001, 1.0, length.out=30)
lr_values <- seq(0.0000001, 4.0, length.out=20)

# Update the grid to include the noise parameter
grid <- expand.grid(c = c_values, lr = lr_values, noise = noise_values, Value = NA_real_)



a_testOnly=list(pred_dat="test_avg",pred_fun="alm.responseOnly",loss_fun="RMSE",loss_data="test_error")
a_trainOnly=list(pred_dat="test_avg",pred_fun="alm.responseOnly",loss_fun="RMSE",loss_data="train_error")
a_testTrain=list(pred_dat="test_avg",pred_fun="alm.responseOnly",loss_fun="RMSE",loss_data="test_error+train_error")

e_testOnly=list(pred_dat="test_avg",pred_fun="exam.response",loss_fun="RMSE",loss_data="test_error")
e_trainOnly=list(pred_dat="test_avg",pred_fun="exam.response",loss_fun="RMSE",loss_data="train_error")
e_testTrain=list(pred_dat="test_avg",pred_fun="exam.response",loss_fun="RMSE",loss_data="test_error+train_error")


plan(multisession)

ex_te_v <- wrap_grid.ll(vAvg, c_values, lr_values,noise_values, input.layer, output.layer,predParams=e_testOnly)
ex_tetr_v <- wrap_grid.ll(vAvg, c_values, lr_values,noise_values, input.layer, output.layer,predParams=e_testTrain)
ex_tr_v <- wrap_grid.ll(vAvg, c_values, lr_values,noise_values, input.layer, output.layer,predParams=e_trainOnly)

a_te_v <- wrap_grid.ll(vAvg, c_values, lr_values,noise_values, input.layer, output.layer,predParams=a_testOnly)
a_tetr_v <- wrap_grid.ll(vAvg, c_values, lr_values,noise_values, input.layer, output.layer,predParams=a_testTrain)
a_tr_v <- wrap_grid.ll(vAvg, c_values, lr_values,noise_values, input.layer, output.layer,predParams=a_trainOnly)

ex_te_c <- wrap_grid.ll(cAvg, c_values, lr_values,noise_values, input.layer, output.layer,predParams=e_testOnly)
ex_tetr_c <- wrap_grid.ll(cAvg, c_values, lr_values,noise_values, input.layer, output.layer,predParams=e_testTrain)
ex_tr_c <- wrap_grid.ll(cAvg, c_values, lr_values,noise_values, input.layer, output.layer,predParams=e_trainOnly)

a_te_c <- wrap_grid.ll(cAvg, c_values, lr_values,noise_values, input.layer, output.layer,predParams=a_testOnly)
a_tetr_c <- wrap_grid.ll(cAvg, c_values, lr_values,noise_values, input.layer, output.layer,predParams=a_testTrain)
a_tr_c <- wrap_grid.ll(cAvg, c_values, lr_values,noise_values, input.layer, output.layer,predParams=a_trainOnly)



wrap_grid.ll <- function(dat, c_values, lr_values, noise_values, input.layer, output.layer,
                         predParams=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_data="test_error")) {
  
  # Precompute as much as possible outside the loop
  train_data <- dat[dat$expMode2 == "Train", ]
  test_data <- dat[dat$expMode2 == "Test", ] 
  test_avg <- test_data |> group_by(x) |> summarise(y=mean(y))
  
    trainVec <- sort(unique(train_data$x))  

  pred_dat <- get(predParams$pred_dat)
  pred_fun <- match.fun(predParams$pred_fun)
  
  # Create a grid dataframe with preallocated space
  grid <- expand.grid(c = c_values, lr = lr_values, noise = noise_values, Value = NA_real_)
  
  # Use future_map (parallel version of map) to iterate over rows in the grid
  loop_time <-  system.time({ 
    results <- future_map_dfr(seq_len(nrow(grid)), function(i) {
      c <- grid$c[i]
      lr <- grid$lr[i]
      noise <- grid$noise[i]
      train_results <- alm.sim(train_data, c, lr, input.layer, output.layer, noise)
      weight.mat <- train_results$wm
      
      # Compute log likelihood loss for train data

      train_error <- compute_loss(train_results, train_data, output.layer)
      
      # Compute log likelihood loss for test data
      test_prediction <- sapply(pred_dat$x, function(x) pred_fun(x, c, input.layer, output.layer, weight.mat, trainVec,noise))

    test_error <- compute_test_log_likelihood(pred_dat$y, test_prediction, output_layer)
      
      # Combine train and test errors if needed
      error <- eval(parse(text=predParams$loss_data))
      
      tibble(c = c, lr = lr, noise = noise, Value = error)  # Return a tibble with the results
    })
  })# stop timer
  
  min_value <- min(results$Value, na.rm = TRUE)
  
  # Check if the minimum value occurs more than once
  if(sum(results$Value == min_value, na.rm = TRUE) > 1 && !is.nan(min_value)) {
    warning("The minimum value in the grid occurs more than once.")
  }
  
  # Extract the best fit parameters
  bestFit <- results[which.min(results$Value), ] |>  mutate(across(where(is.numeric), \(x) round(x, 7)))
  
  bf_train_result <- alm.sim(train_data, bestFit$c, bestFit$lr, input.layer, output.layer,bestFit$noise)
  weight.mat <- bf_train_result$wm
  bf_train <- bf_train_result$d
  best_preds <- sapply(pred_dat$x, function(x) pred_fun(x, bestFit$c, input.layer, output_layer, weight.mat,trainVec, bestFit$noise))
  
  (bf_test <- cbind(pred_dat, pred = best_preds ))
  print(paste0("Value:", bestFit$Value, " c:", bestFit$c, " lr:", bestFit$lr, " noise:", bestFit$noise, " time:", round(loop_time[3],3)))
  
  return(list(train=bf_train,test=bf_test,errorGrid=results, Fit = bestFit, c = bestFit$c, lr = bestFit$lr, noise = bestFit$noise, Value = bestFit$Value, time=loop_time[3]))
  
}




# _____________________________________________

# Negative log likelihood loss function
neg_log_likelihood <- function(prediction, target) {
  # Assuming a normal distribution with mean equal to prediction and some standard deviation (noise)
  # We want to minimize this value, hence the negative log likelihood
  -sum(dnorm(target, mean = prediction, sd = noise, log = TRUE))
}

# Updated ALM response function to return probabilities instead of mean response
alm.responseOnly <- function(input = 1, c, input.layer, output.layer, weight.mat, noise, trainVec=NULL) {
  input.activation <- exp(-c * (input.layer - input)^2)
  output.activation <- weight.mat %*% input.activation
  output.probability <- output.activation / sum(output.activation)
  output.probability
}

# Similarly update alm.response if needed


# ... within wrap_grid function

# Define noise values range as another grid dimension
noise_values <- seq(0.1, 10.0, length.out=100)  # Adjust this range as needed

# Update the grid to include noise
grid <- expand.grid(c = c_values, lr = lr_values, noise = noise_values, Value = NA_real_)

# ... within loop iterating over grid rows
  # Extract noise from grid
  noise <- grid$noise[i]
  
  # Use the new neg_log_likelihood loss function
  train_error <- neg_log_likelihood(train_results$d$almResp, train_results$d$y)
  test_prediction <- map_dbl(pred_dat$x, ~ pred_fun(.x, c, input.layer, output.layer, weight.mat, noise, trainVec=trainVec))
  test_error <- neg_log_likelihood(test_prediction, pred_dat$y)




# _____________________________________________



alm.response <- function(input = 1, c, input.layer, output.layer, weight.mat, trainVec=NULL) {
    input.activation <- exp(-c * (input.layer - input)^2) / (sum(exp(-c * (input.layer - input)^2)) +.01)
    output.activation <- (weight.mat %*% input.activation) + .01
    output.probability <- output.activation / sum(output.activation)
    return(output.probability)
}

alm.sim <- function(dat, c, lr, input.layer = INPUT_LAYER_DEFAULT, output.layer = OUTPUT_LAYER_DEFAULT) {
    weight.mat <- matrix(0.000001, nrow = length(output.layer), ncol = length(input.layer))
    n <- nrow(dat)
    prob_distributions <- matrix(nrow = n, ncol = length(output.layer)) # Store probability distributions

    for(i in 1:n) {
        trial <- alm.trial(dat$x[i], dat$y[i], c, lr, input.layer, output.layer, weight.mat)
        weight.mat <- trial$weight.mat
        prob_distributions[i, ] <- trial$output.probability
    }
    dat <- cbind(dat, almResp=prob_distributions)
    return(list(d = dat, wm = weight.mat, c = c, lr = lr))
}


log_likelihood <- function(observed, predicted_probs) {
    ll <- sum(log(predicted_probs[cbind(1:length(observed), observed)]))
    return(ll)
}

exam.response <- function(input, c, input.layer = INPUT_LAYER_DEFAULT, output.layer = OUTPUT_LAYER_DEFAULT, weight.mat, trainVec) {
  nearestTrain <- trainVec[which.min(abs(input - trainVec))]
  nearestTrainProb <- alm.response(nearestTrain, c, input.layer, output.layer, weight.mat)$output.probability

  # Modify the following logic to interpolate between probability distributions
  # rather than mean responses.
  # ...
}

# _____________________________________________

alm.responseOnly <- function(input = 1, c, input.layer, output.layer,weight.mat, trainVec=NULL) {
  input.activation <- exp(-c * (input.layer - input)^2) / (sum(exp(-c * (input.layer - input)^2)) +.01)
  output.activation <- (weight.mat %*% input.activation) + .01
  output.probability <- output.activation / sum(output.activation)
  output.layer * output.probability # Returns probabilities over output layer
}

exam.response <- function(input, c, input.layer = INPUT_LAYER_DEFAULT,output.layer = OUTPUT_LAYER_DEFAULT, weight.mat, trainVec) {
  nearestTrain <- trainVec[which.min(abs(input - trainVec))]
  aresp <- alm.response(nearestTrain, c, input.layer = input.layer,output.layer,weight.mat)$mean.response
  
  xUnder <- ifelse(min(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) - 1])
  xOver <- ifelse(max(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) + 1])
  
  mUnder <- alm.response(xUnder, c, input.layer = input.layer, output.layer, weight.mat)$mean.response
  mOver <- alm.response(xOver, c, input.layer = input.layer,output.layer, weight.mat)$mean.response
  
  output <- aresp + ((mOver - mUnder) / (xOver - xUnder)) * (input - nearestTrain)
  output.layer * output # Returns probabilities over output layer
}

# _____________________________________________




ll <- function(true, pred) {
  -sum(dnorm(true, mean = pred, sd = 1, log = TRUE)) 
}

# Update wrap_grid to optimize log likelihood
wrap_grid <- function(dat, c_values, lr_values, input.layer, output.layer,
                      predParams=list(pred_dat="test_avg", pred_fun="exam.response", 
                                      loss_fun="ll", loss_data="test_ll")) {

  # Rest of function unchanged until results assignment

  results <- future_map_dfr(seq_len(nrow(grid)), function(i) {
    
    # Unchanged
    
    train_ll <- ll(train_results$d$x, train_results$d$almResp) 
    test_ll <- ll(test_prediction, pred_dat$y)
    
    ll <- eval(parse(text=predParams$loss_data))
    
    tibble(c = c, lr = lr, Value = ll)
  })

  # Rest unchanged
  
}

# Update alm.responseOnly to return log probabilities
alm.responseOnly <- function(input = 1, c, input.layer, output.layer,weight.mat, trainVec=NULL) {

  # Unchanged until output.probability
  
  output.logprob <- log(output.probability)
  
  output.logprob
  
}




# _____________________________________________



# New loss function for log likelihood
logLik <- function(pred, obs) {
  -sum(dnorm(obs, pred, sd=sigma, log=TRUE)) 
}

# sigma is another parameter to optimize
sigma <- 1 

# Modify exam response function to return probability 
exam.response <- function(input, c, input.layer = INPUT_LAYER_DEFAULT, output.layer = OUTPUT_LAYER_DEFAULT, weight.mat, trainVec, sigma) {

  # Get nearest neighbor
  nearestTrain <- trainVec[which.min(abs(input - trainVec))]
  
  # Get neighbor responses
  aresp <- alm.response(nearestTrain, c, input.layer, output.layer, weight.mat)$mean.response
  xUnder <- ifelse(min(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) - 1])
  xOver <- ifelse(max(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) + 1])
  mUnder <- alm.response(xUnder, c, input.layer, output.layer, weight.mat)$mean.response
  mOver <- alm.response(xOver, c, input.layer, output.layer, weight.mat)$mean.response

  # Calculate exam prediction
  exam.output <- aresp + ((mOver - mUnder) / (xOver - xUnder)) * (input - nearestTrain)
  
  # Return probability density
  dnorm(exam.output, exam.output, sigma)
}

# Update wrap_grid to optimize log likelihood
wrap_grid <- function(dat, c_values, lr_values, input.layer, output.layer,
                      predParams=list(pred_dat="test_avg", pred_fun="exam.response", 
                                      loss_fun="logLik", loss_data="test_error")) {

  # Rest of function...

  # Add sigma to grid
  sigma_values <- seq(0.1, 5, length.out=100) 
  grid <- expand.grid(c = c_values, lr = lr_values, sigma = sigma_values, Value = NA_real_)

  # Pass sigma to prediction function
  test_prediction <- map_dbl(pred_dat$x, ~ pred_fun(.x, c, input.layer, output.layer, weight.mat,  
                              trainVec, sigma))

  # Rest of function...

}

