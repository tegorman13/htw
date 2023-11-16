

fit_lr_c <- function(params_list,pred_dat,pred_fun,loss_fun,loss_dat,model_params,test_data,train_data)
{
  c = params_list[1] |> as.numeric()
  lr = params_list[2] |> as.numeric()
  list2env(model_params,envir = environment())
  train_results <- alm.sim(train_data, c, lr, input.layer, output.layer)
  weight.mat <- train_results$wm
  
    if(train_data$condit[[1]]=="Varied") {trainVec=c(800,1000,1200)}

  tryCatch({
    # Your existing function logic...
    # Ensure at the end of your computations, you check for NA
    test_prediction <- map_dbl(pred_dat$x, ~ pred_fun(.x, c, input.layer, output.layer, weight.mat,  trainVec=trainVec))
   #test_prediction[is.nan(test_prediction)] <- 20000
 #nan_penalty <- ifelse(any(is.nan(test_prediction)), 1e6, 0) # Large penalty for NaN
   #invalid_value_penalty <- ifelse(any(is.nan(test_prediction) | is.na(test_prediction)), 1e2, 0)

  if(length(params_list) >= 3 && !is.null(params_list[3])){
  #train_error <- loss_fun(train_results$d$y, train_results$d$almResp,params_list[3] |> as.numeric())
  test_error <- loss_fun(pred_dat$y,test_prediction,params_list[3] |> as.numeric())
  } else{
    train_error <- loss_fun(train_results$d$y, train_results$d$almResp)
    test_error <- loss_fun(pred_dat$y,test_prediction)
  }
  
  error <- eval(parse(text=loss_dat)) + invalid_value_penalty
    if (is.na(error)) {
      error <- 1e2 # Large penalty
    }
    return(error)
  }, error = function(e) {
    return(1e6) # Large penalty in case of any error
  })


  
}


fit_lr_c1 <- function(params_list,pred_dat,pred_fun,loss_fun,loss_dat,model_params,test_data,train_data)
{
  c = params_list[1] |> as.numeric()
  lr = params_list[2] |> as.numeric()
  list2env(model_params,envir = environment())
  train_results <- alm.sim(train_data, c, lr, input.layer, output.layer)
  weight.mat <- train_results$wm
  test_prediction <- map_dbl(pred_dat$x, ~ pred_fun(.x, c, input.layer, output.layer, weight.mat,  trainVec=trainVec))
   test_prediction[is.nan(test_prediction)] <- 2000000

  if(length(params_list) >= 3 && !is.null(params_list[3])){
  train_error <- loss_fun(train_results$d$y, train_results$d$almResp,params_list[3] |> as.numeric())
  test_error <- loss_fun(pred_dat$y,test_prediction,params_list[3] |> as.numeric())
  } else{
    train_error <- loss_fun(train_results$d$y, train_results$d$almResp)
    test_error <- loss_fun(pred_dat$y,test_prediction)
  }
  
  error <- eval(parse(text=loss_dat))
  
}


fit_lr_c0 <- function(params_list,pred_dat,pred_fun,loss_fun,loss_dat,model_params,test_data,train_data)
{
  c = params_list[1] |> as.numeric()
  lr = params_list[2] |> as.numeric()
  list2env(model_params,envir = environment())
  train_results <- alm.sim(train_data, c, lr, input.layer, output.layer)
  weight.mat <- train_results$wm
  test_prediction <- map_dbl(pred_dat$x, ~ pred_fun(.x, c, input.layer, output.layer, weight.mat,  trainVec=trainVec))
   #test_prediction[is.nan(test_prediction)] <- 2000000

  if(length(params_list) >= 3 && !is.null(params_list[3])){
  train_error <- loss_fun(train_results$d$y, train_results$d$almResp,params_list[3] |> as.numeric())
  test_error <- loss_fun(pred_dat$y,test_prediction,params_list[3] |> as.numeric())
  } else{
    train_error <- loss_fun(train_results$d$y, train_results$d$almResp)
    test_error <- loss_fun(pred_dat$y,test_prediction)
  }
  
  error <- eval(parse(text=loss_dat))
  
}


grid_fit <- function(grid,func, fit_params, model_params,dat_params) {
  
  pred_dat <- get(fit_params$pred_dat)
  loss_dat <- fit_params$loss_data
  pred_fun <- match.fun(fit_params$pred_fun)
  loss_fun <- match.fun(fit_params$loss_fun)
  train_data=dat_params$train_data
  test_data=dat_params$test_data
  
  list2env(model_params,envir = environment())
  plan(multisession)
  loop_time <- system.time({
    results <- future_map_dfr(seq_len(nrow(grid)), function(i) {
      params_list <- tibble(grid[i,])
      error <- func(params_list, pred_dat,pred_fun,loss_fun,loss_dat,model_params,test_data,train_data)
      tibble(params_list, Value = error)
    })
  })
  min_value <- min(results$Value, na.rm = TRUE)
  if(sum(results$Value == min_value, na.rm = TRUE) > 1 && !is.nan(min_value)) {
    warning("The minimum value in the grid occurs more than once.")
  }
  
  bestFit <- results[which.min(results$Value), ] |>  mutate(across(where(is.numeric), \(x) round(x, 8)))
  
  bf_train_result <- alm.sim(train_data, bestFit$c, bestFit$lr, input.layer, output.layer)
  weight.mat <- bf_train_result$wm
  bf_train <- bf_train_result$d
  best_preds <- map_dbl(pred_dat$x, ~ pred_fun(.x, bestFit$c,input.layer, output.layer, weight.mat,trainVec))
  
  (bf_test <- cbind(pred_dat, pred = best_preds ))
  print(paste0("Value:", bestFit$Value, " c:", bestFit$c, " lr:", bestFit$lr, " time:",  round(loop_time[3],3)))
  return(list(train=bf_train,test=bf_test,errorGrid=results, Fit = bestFit, c = bestFit$c, lr = bestFit$lr, Value = bestFit$Value, time=loop_time[3]))
  
}





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
