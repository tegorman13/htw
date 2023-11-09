

wrap_grid <- function(dat, c_values, lr_values, input.layer, output.layer,
                            predParams=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="RMSE", loss_data="test_error")) {
  
  # Precompute as much as possible outside the loop
  train_data <- dat[dat$expMode2 == "Train", ]
  test_data <- dat[dat$expMode2 == "Test", ] 
  test_avg <- test_data |> group_by(x) |> summarise(y=mean(y))
  trainVec <- sort(c(0,unique(train_data$x)))
  
  pred_dat <- get(predParams$pred_dat)
  pred_fun <- match.fun(predParams$pred_fun)
  loss_fun <- match.fun(predParams$loss_fun)
  
  
  
  # Create a grid dataframe with preallocated space
  grid <- expand.grid(c = c_values, lr = lr_values, Value = NA_real_)
  
  # Use future_map (parallel version of map) to iterate over rows in the grid
  loop_time <-  system.time({ 
    results <- future_map_dfr(seq_len(nrow(grid)), function(i) {
      c <- grid$c[i]
      lr <- grid$lr[i]
      train_results <- alm.sim(train_data, c, lr, input.layer, output.layer)
      weight.mat <- train_results$wm
      
      train_error <- loss_fun(train_results$d$x, train_results$d$almResp)
      test_prediction <- map_dbl(pred_dat$x, ~ pred_fun(.x, c, input.layer, output.layer, weight.mat,  trainVec=trainVec))
      test_error <- loss_fun(test_prediction, pred_dat$y)
      
      error <- eval(parse(text=predParams$loss_data)) 
      
      tibble(c = c, lr = lr, Value = error)  # Return a tibble with the results
    })
  })# stop timer
  
  min_value <- min(results$Value, na.rm = TRUE)
  
  # Check if the minimum value occurs more than once
  if(sum(results$Value == min_value, na.rm = TRUE) > 1 && !is.nan(min_value)) {
    warning("The minimum value in the grid occurs more than once.")
  }
  
  # Extract the best fit parameters
  bestFit <- results[which.min(results$Value), ] |>  mutate(across(where(is.numeric), \(x) round(x, 5)))
  
  bf_train_result <- alm.sim(train_data, bestFit$c, bestFit$lr, input.layer, output.layer)
  weight.mat <- bf_train_result$wm
  bf_train <- bf_train_result$d
  best_preds <- map_dbl(pred_dat$x, ~ pred_fun(.x, bestFit$c, input.layer, output.layer, weight.mat,trainVec))
  
  (bf_test <- cbind(pred_dat, pred = best_preds ))
  print(paste0("Value:", bestFit$Value, " c:", bestFit$c, " lr:", bestFit$lr, " time:", round(loop_time[3],3)))
  
  return(list(train=bf_train,test=bf_test,errorGrid=results, Fit = bestFit, c = bestFit$c, lr = bestFit$lr, Value = bestFit$Value, time=loop_time[3]))
  
}

predict_alm_exam_weighted_hybrid <- function(input, c, weight_exam, trainData, input.layer, output.layer, weight.mat) {
  alm_pred <- alm.responseOnly(input, c, input.layer, output.layer,weight.mat, trainVec=NULL)
  exam_pred <- predict_alm_exam(input, c, input.layer, output.layer,weight.mat, trainVec=NULL)
  (1 - weight_exam) * alm_pred + weight_exam * exam_pred
}


alm.responseOnly <- function(input = 1, c, input.layer, output.layer,weight.mat, trainVec=NULL) {
  input.activation <- exp(-c * (input.layer - input)^2) / (sum(exp(-c * (input.layer - input)^2)) +.01)
  output.activation <- (weight.mat %*% input.activation) + .01
  output.probability <- output.activation / sum(output.activation)
  sum(output.layer * output.probability)
}


alm.response <- function(input = 1, c, input.layer, output.layer,weight.mat, trainVec=NULL) {
  input.activation <- exp(-c * (input.layer - input)^2) / (sum(exp(-c * (input.layer - input)^2)) +.01)
  output.activation <- (weight.mat %*% input.activation) + .01
  output.probability <- output.activation / sum(output.activation)
  mean.response <- sum(output.layer * output.probability)
  list(mean.response = mean.response, input.activation = input.activation, output.activation = output.activation)
}

alm.update <- function(corResp, c, lr, output.layer, input.activation, output.activation, weight.mat) {
  fz <- exp(-c * (output.layer - corResp)^2)
  teacherSignal <- (fz - output.activation) * lr
  wChange <- teacherSignal %*% t(input.activation)
  weight.mat <- weight.mat + wChange
  #weight.mat[weight.mat < 0] = 0
  return(weight.mat)
}

alm.trial <- function(input, corResp, c, lr, input.layer, output.layer, weight.mat) {
  alm_resp <- alm.response(input, c, input.layer,output.layer, weight.mat)
  updated_weight.mat <- alm.update(corResp, c, lr, output.layer, alm_resp$input.activation, alm_resp$output.activation, weight.mat)
  return(list(mean.response = alm_resp$mean.response, weight.mat = updated_weight.mat))
}

alm.sim <- function(dat, c, lr, input.layer = INPUT_LAYER_DEFAULT, output.layer = OUTPUT_LAYER_DEFAULT) {
  weight.mat <- matrix(0.000001, nrow = length(output.layer), ncol = length(input.layer))
  xt <- dat$x
  n <- nrow(dat)
  st <- numeric(n) # Initialize the vector to store mean responses
  for(i in 1:n) {
    trial <- alm.trial(dat$x[i], dat$y[i], c, lr, input.layer, output.layer, weight.mat)
    weight.mat <- trial$weight.mat
    st[i] <- trial$mean.response
  }
   dat <- cbind(dat,almResp=st)
  return(list(d = dat, wm = weight.mat, c = c, lr = lr))
}


exam.response <- function(input, c, input.layer = INPUT_LAYER_DEFAULT,output.layer = OUTPUT_LAYER_DEFAULT, weight.mat, trainVec) {
  nearestTrain <- trainVec[which.min(abs(input - trainVec))]
  aresp <- alm.response(nearestTrain, c, input.layer = input.layer,output.layer,weight.mat)$mean.response
  
  xUnder <- ifelse(min(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) - 1])
  xOver <- ifelse(max(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) + 1])
  
  mUnder <- alm.response(xUnder, c, input.layer = input.layer, output.layer, weight.mat)$mean.response
  mOver <- alm.response(xOver, c, input.layer = input.layer,output.layer, weight.mat)$mean.response
  
  exam.output <- round(aresp + ((mOver - mUnder) / (xOver - xUnder)) * (input - nearestTrain), 3)
  exam.output
}



predict_alm_exam_weighted_hybrid <- function(input, c, weight_exam, input.layer, output.layer, weight.mat,trainVec) {
  
  alm_pred <- alm.responseOnly(input, c, input.layer, output.layer,weight.mat, trainVec=NULL)
  exam_pred <- exam.response(input, c, input.layer, output.layer,weight.mat, trainVec=trainVec)
  (1 - weight_exam) * alm_pred + weight_exam * exam_pred
}


exam_kwantes <- function(input, c, trainData, input.layer = INPUT_LAYER_DEFAULT, output.layer = OUTPUT_LAYER_DEFAULT, weight.mat) {
  # Initialize a vector to store the slopes
  slopes <- numeric(length(trainData))
  
  # Loop over the training data
  for (i in 2:length(trainData)) {
    # Compute the change in Y from the previous trial
    deltaY <- trainData[i] - trainData[i - 1]
    
    # Compute the change in X from the previous trial
    deltaX <- trainData[i] - trainData[i - 1]
    
    # Update the slope estimate
    slopes[i] <- deltaY / deltaX
  }
  # Compute the ALM response for the current input
  almResp <- alm.response(input, c, input.layer, output.layer, weight.mat)$response
  
  # Find the nearest training x-value
  nearestTrainX <- trainData[which.min(abs(input - trainData))]
  
  # Find the slope associated with the nearest training x-value
  slope <- slopes[which.min(abs(input - trainData))]
  
  # Compute the EXAM output based on the estimated slope
  exam.output <- round(almResp + slope * (input - nearestTrainX), 3)
  
  exam.output
}


wrap_grid_serial <- function(dat,c_values,lr_values, input.layer, output.layer,
                             predParams=list(pred_dat="test_avg",pred_fun="exam.response",loss_fun="RMSE",loss_data="test_error")){
  
  cRange <- c_values
  lrRange <-lr_values
  grid <- expand.grid(c = cRange, lr = lrRange)
  
  train_data <- dat[dat$expMode2 == "Train", ]
  test_data <- dat[dat$expMode2 == "Test", ] 
  test_avg<- test_data |> group_by(x) |> summarise(y=mean(y))
  
  trainVec <- sort(unique(train_data$x))
  
  pred_dat <- get(predParams$pred_dat)
  pred_fun <- match.fun(predParams$pred_fun)
  loss_fun <- get(predParams$loss_fun)
  
  grid$Value <- NA
  loop_time <- system.time({ 
    for (i in 1:nrow(grid)) {
      
      c=grid[i,c("c") ]; lr = grid[i,c("lr") ]; 
      train_results <- alm.sim(train_data, c, lr, input.layer, output.layer)
      weight.mat <- train_results$wm
      
      train_error = loss_fun(train_results$d$x,train_results$d$almResp)
      
      test_prediction <- map_dbl(pred_dat$x, ~ pred_fun(.x, c, input.layer, output.layer, weight.mat, trainVec=trainVec))
      test_error <- loss_fun(test_prediction, pred_dat$y)
      
      error <- eval(parse(text=predParams$loss_data))
      
      grid[i, "Value"] <- error
      
    }
  })# stop timer
  
  if(sum(grid$Value == min(grid$Value)) > 1) warning("The minimum value in the grid occurs more than once.")
  bestFit <- grid[which.min(grid$Value), ] |>  mutate(across(where(is.numeric), \(x) round(x, 5)))
  
  bf_train_result <- alm.sim(train_data, bestFit$c, bestFit$lr, input.layer, output.layer)
  weight.mat <- bf_train_result$wm
  bf_train <- bf_train_result$d
  
  best_preds <- map_dbl(pred_dat$x, ~ pred_fun(.x, bestFit$c, input.layer, output.layer, weight.mat,trainVec))
  
  (bf_test <- cbind(pred_dat, pred = best_preds ))
  print(paste0("Value:", bestFit$Value, " c:", bestFit$c, " lr:", bestFit$lr, " time:", round(loop_time[3],3)))
  
  return(list(train=bf_train,test=bf_test,errorGrid=grid,Fit = bestFit, c = bestFit$c, lr = bestFit$lr, Value = bestFit$Value, time=loop_time[3]))
  
}