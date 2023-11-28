





de_optim_fit <- function(func, fit_params, model_params, dat_params, initial_params) {
  
  purrr::walk(here::here(c("Functions/fun_model.R","Functions/alm_core.R","Functions/fit_funs.R")),source)
  model_fun = get(func)
  train_data <- dat_params$train_data
  test_data <- dat_params$test_data
  test_avg <- dat_params$test_avg
  pred_dat <- get(fit_params$pred_dat)
  loss_dat <- fit_params$loss_data
  pred_fun <- match.fun(fit_params$pred_fun)
  loss_fun <- match.fun(fit_params$loss_fun)
  initial_params <- initial_params
  list2env(model_params,envir = environment())
  # Define the objective function for DEoptim
  objective_function <- function(params) {
    model_fun(params, pred_dat, pred_fun, loss_fun, loss_dat, model_params, test_data, train_data)
  }
  
  # Define bounds for parameters
  lower_bounds <- c(.0000001,.000001,100)
  upper_bounds <- c(10,10,2000)
  initialpop=initial_params
  
  # Define number of particles (NP) and iterations (maxiter)
  NP <- 10 * length(initial_params) + 1
  maxiter <- 300
  strategy=6
  
  optim_time <- system.time({
    # Using DEoptim instead of JDEoptim
    de_optim_res <- DEoptim(fn = objective_function, lower = lower_bounds, 
                            upper = upper_bounds, control = list(NP = NP, itermax = maxiter,strategy=strategy))
  })
  
  bestFit <- as.data.frame(as.list(de_optim_res$optim$bestmem)) |> mutate(Value=round(de_optim_res$optim$bestval,3))
  names(bestFit) <- c("c","lr","sigma","Value")
  
  
  bf_train_result <- alm.sim(train_data, bestFit$c, bestFit$lr, input.layer, output.layer)
  weight.mat <- bf_train_result$wm
  bf_train <- bf_train_result$d
  best_preds <- map_dbl(pred_dat$x, ~ pred_fun(.x, bestFit$c,input.layer, output.layer, weight.mat,trainVec))
  (bf_test <- cbind(pred_dat, pred = best_preds ))
  
  if (any(is.nan(best_preds))) {
    warning("Warning: test_prediction contains NaN values.")
  }
  
  if(length(is.null(bestFit$sigma))){
    train_error <- loss_fun(bf_train$y, bf_train$almResp,bestFit$sigma) |> as.numeric()
    test_error <- loss_fun(pred_dat$y, best_preds,bestFit$sigma |> as.numeric())
  } else{
    train_error <- loss_fun(bf_train$y, bf_train$almResp)
    test_error <- loss_fun(pred_dat$y, best_preds)
  }
  
  
  print_b <-paste0(strip_list_notation(tibble::lst(round_tibble(bestFit,3))))
  print(paste0(print_b, ",  time:",  round(optim_time[3],1)))
  
  optim_info <- list(id=as.numeric(train_data$sbj[1]),condit=as.character(train_data$condit[1]),
                     Value = bestFit$Value, c = bestFit$c, lr = bestFit$lr, 
                     print_b,
                     NAN_Test=any(is.nan(best_preds)),
                     model_fun=func, 
                     pred_fun=fit_params$pred_fun,
                     loss_dat=loss_dat,loss_fun=fit_params$loss_fun, 
                     train=bf_train,test=bf_test, Fit = bestFit, 
                     train_error=train_error,test_error=test_error,
                     run_args=tibble::lst(func,fit_params,model_params,dat_params,trainVec,input.layer,output.layer),
                     opt_info=tibble::lst(de_optim_res,initial_params, NP, maxiter, lower_bounds,upper_bounds),
                     info=tibble::lst(time=Sys.time(),Sys.info(), path=getwd()),
                     runtime=optim_time[3])
  
  return(optim_info) 
}



de_optim_fit_id <- function(id_data, func, fit_params, model_params, initial_params) {
  train_data <- id_data %>% filter(expMode2=="Train") # Adjust condition as needed
  test_data <- id_data %>% filter(expMode2=="Test") # Adjust condition as needed
  test_avg <- test_data |> ungroup() |>  group_by(sbj,x,condit) |> summarise(y=mean(y), .groups="drop")
  
  dat_params <- tibble::lst(test_data = test_data, train_data = train_data, test_avg=test_avg)
  
  result <- tryCatch({
    de_optim_fit(func, fit_params, model_params, dat_params, initial_params)
  }, error = function(e) {
    message("Error during optimization for ID ", unique(id_data$id), ": ", e$message)
    return(NULL) 
  })
  
}











optim_fit <- function(func, fit_params, model_params,dat_params,opt.m="BFGS",initial_params) {
  purrr::walk(here::here(c("Functions/fun_model.R","Functions/alm_core.R","Functions/fit_funs.R")),source)
  model_fun = get(func)
  train_data <- dat_params$train_data
  test_data <- dat_params$test_data
  test_avg <- dat_params$test_avg
  pred_dat <- get(fit_params$pred_dat)
  loss_dat <- fit_params$loss_data
  pred_fun <- match.fun(fit_params$pred_fun)
  loss_fun <- match.fun(fit_params$loss_fun)
  initial_params <- initial_params
  list2env(model_params,envir = environment())
  
  
  optim_time <- system.time({
    optim_res <- optim(initial_params, model_fun, pred_dat = pred_dat, pred_fun = pred_fun, 
                       loss_fun = loss_fun, loss_dat = loss_dat, 
                       model_params = model_params, test_data = test_data, train_data = train_data,
                       method = opt.m) # Choose an appropriate method
  })
  
  bestFit <- as.data.frame(as.list(optim_res$par)) |> mutate(Value=round(optim_res$value,3))
  
  bf_train_result <- alm.sim(train_data, bestFit$c, bestFit$lr, input.layer, output.layer)
  weight.mat <- bf_train_result$wm
  bf_train <- bf_train_result$d
  best_preds <- map_dbl(pred_dat$x, ~ pred_fun(.x, bestFit$c,input.layer, output.layer, weight.mat,trainVec))
  (bf_test <- cbind(pred_dat, pred = best_preds ))
  
  if (any(is.nan(best_preds))) {
    warning("Warning: test_prediction contains NaN values.")
  }
  
  if(length(is.null(bestFit$sigma))){
    train_error <- loss_fun(bf_train$y, bf_train$almResp,bestFit$sigma) |> as.numeric()
    test_error <- loss_fun(pred_dat$y, best_preds,bestFit$sigma |> as.numeric())
  } else{
    train_error <- loss_fun(bf_train$y, bf_train$almResp)
    test_error <- loss_fun(pred_dat$y, best_preds)
  }
  
  
  print_b <-paste0(strip_list_notation(tibble::lst(round_tibble(bestFit,3))))
  print(paste0(print_b, ",  time:",  round(optim_time[3],1)))
  
  optim_info <- list(id=as.numeric(train_data$sbj[1]),condit=as.character(train_data$condit[1]),
                     print_b,
                     Value = bestFit$Value, c = bestFit$c, lr = bestFit$lr, 
                     NAN_Test=any(is.nan(best_preds)),
                     opt.m=opt.m,
                     model_fun=func, 
                     pred_fun=fit_params$pred_fun,
                     loss_dat=loss_dat,loss_fun=fit_params$loss_fun, 
                     train=bf_train,test=bf_test, Fit = bestFit, c = bestFit$c, lr = bestFit$lr, Value = bestFit$Value, 
                     train_error=train_error,test_error=test_error,
                     #functions= tibble::lst(pred_fun, loss_fun), 
                     run_args=tibble::lst(func,fit_params,model_params,dat_params,opt.m,initial_params, trainVec,input.layer,output.layer),
                     opt_info=optim_res,
                     info=tibble::lst(time=Sys.time(),Sys.info(), path=getwd()),
                     runtime=optim_time[3])
  
  
  
  return(optim_info) # Return desired results
}



optim_fit_id <- function(id_data, func, fit_params, model_params, opt.m = "BFGS", initial_params) {
  train_data <- id_data %>% filter(expMode2=="Train") 
  test_data <- id_data %>% filter(expMode2=="Test") 
  test_avg <- test_data |> ungroup() |>  group_by(sbj,x,condit) |> summarise(y=mean(y), .groups="drop")
  
  dat_params <- tibble::lst(test_data = test_data, train_data = train_data, test_avg=test_avg)
  
  result <- tryCatch({
    optim_fit(func, fit_params, model_params, dat_params, opt.m = opt.m, initial_params)
  }, error = function(e) {
    message("Error during optimization for ID ", unique(id_data$id), ": ", e$message)
    return(NULL) 
  })
  
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




optim_fit_l_BFGS_b <- function(func, fit_params, model_params, dat_params, initial_params, lower_bounds, upper_bounds, opt.m="L-BFGS-B") {
  purrr::walk(here::here(c("Functions/fun_model.R","Functions/alm_core.R","Functions/fit_funs.R")),source)
  model_fun = get(func)
  train_data <- dat_params$train_data
  test_data <- dat_params$test_data
  test_avg <- dat_params$test_avg
  pred_dat <- get(fit_params$pred_dat)
  loss_dat <- fit_params$loss_data
  pred_fun <- match.fun(fit_params$pred_fun)
  loss_fun <- match.fun(fit_params$loss_fun)
  initial_params <- initial_params
  list2env(model_params,envir = environment())
  
  optim_time <- system.time({
    optim_res <- optim(initial_params, model_fun, pred_dat = pred_dat, pred_fun = pred_fun, 
                       loss_fun = loss_fun, loss_dat = loss_dat, 
                       model_params = model_params, test_data = test_data, train_data = train_data,
                       method = opt.m, lower = lower_bounds, upper = upper_bounds) # Add bounds and set method
  })
  
  bestFit <- as.data.frame(as.list(optim_res$par)) |> mutate(Value=round(optim_res$value,3))
  
  bf_train_result <- alm.sim(train_data, bestFit$c, bestFit$lr, input.layer, output.layer)
  weight.mat <- bf_train_result$wm
  bf_train <- bf_train_result$d
  best_preds <- map_dbl(pred_dat$x, ~ pred_fun(.x, bestFit$c,input.layer, output.layer, weight.mat,trainVec))
  (bf_test <- cbind(pred_dat, pred = best_preds ))
  
  if (any(is.nan(best_preds))) {
    warning("Warning: test_prediction contains NaN values.")
  }
  
  if(length(is.null(bestFit$sigma))){
    train_error <- loss_fun(bf_train$y, bf_train$almResp,bestFit$sigma) |> as.numeric()
    test_error <- loss_fun(pred_dat$y, best_preds,bestFit$sigma |> as.numeric())
  } else{
    train_error <- loss_fun(bf_train$y, bf_train$almResp)
    test_error <- loss_fun(pred_dat$y, best_preds)
  }
  
  
  print_b <-paste0(strip_list_notation(tibble::lst(round_tibble(bestFit,3))))
  print(paste0(print_b, ",  time:",  round(optim_time[3],1)))
  
  optim_info <- list(id=as.numeric(train_data$sbj[1]),condit=as.character(train_data$condit[1]),
                     print_b,
                     Value = bestFit$Value, c = bestFit$c, lr = bestFit$lr, 
                     NAN_Test=any(is.nan(best_preds)),
                     opt.m=opt.m,
                     model_fun=func, 
                     pred_fun=fit_params$pred_fun,
                     loss_dat=loss_dat,loss_fun=fit_params$loss_fun, 
                     train=bf_train,test=bf_test, Fit = bestFit, c = bestFit$c, lr = bestFit$lr, Value = bestFit$Value, 
                     train_error=train_error,test_error=test_error,
                     #functions= tibble::lst(pred_fun, loss_fun), 
                     run_args=tibble::lst(func,fit_params,model_params,dat_params,opt.m,initial_params, trainVec,input.layer,output.layer),
                     opt_info=optim_res,
                     info=tibble::lst(time=Sys.time(),Sys.info(), path=getwd()),
                     runtime=optim_time[3])
  
  
  
  return(optim_info) # Return desired results
}




optim_fit_id_l_BFGS_b <- function(id_data, func, fit_params, model_params, initial_params, lower_bounds, upper_bounds, opt.m = "L-BFGS-B") {
  train_data <- id_data %>% filter(expMode2=="Train") 
  test_data <- id_data %>% filter(expMode2=="Test") 
  test_avg <- test_data |> ungroup() |>  group_by(sbj,x,condit) |> summarise(y=mean(y), .groups="drop")
  
  dat_params <- tibble::lst(test_data = test_data, train_data = train_data, test_avg=test_avg)
  
  result <- tryCatch({
    optim_fit_l_BFGS_b(func, fit_params, model_params, dat_params, initial_params, lower_bounds, upper_bounds, opt.m=opt.m)
  }, error = function(e) {
    message("Error during optimization for ID ", unique(id_data$id), ": ", e$message)
    return(NULL) 
  })
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
