# nll2 used as loss_fun
nll2 <- function(obsv,pred,sigma)
{
  nll= -sum(dnorm(obsv,mean=pred,sd=sigma,log=TRUE)) +.001
  if (is.nan(nll)) {
    nll <- 1e4 # Large penalty
  }
  return(nll)
}

# fit_lr_c0 used as model_fun
fit_lr_c0 <- function(params_list,pred_dat,pred_fun,loss_fun,loss_dat,model_params,test_data,train_data)
{
  c = params_list[1] |> as.numeric()
  lr = params_list[2] |> as.numeric()
  list2env(model_params,envir = environment())
  train_results <- alm.sim(train_data, c, lr, input.layer, output.layer)
  weight.mat <- train_results$wm
  test_prediction <- map_dbl(pred_dat$x, ~ pred_fun(.x, c, input.layer, output.layer, weight.mat,  trainVec=trainVec))

  if(length(params_list) >= 3 && !is.null(params_list[3])){
  train_error <- loss_fun(train_results$d$y, train_results$d$almResp,params_list[3] |> as.numeric())
  test_error <- loss_fun(pred_dat$y,test_prediction,params_list[3] |> as.numeric())
  } else{
    train_error <- loss_fun(train_results$d$y, train_results$d$almResp)
    test_error <- loss_fun(pred_dat$y,test_prediction)
  }
  
  error <- eval(parse(text=loss_dat))
  
}

library(DEoptim)

de_optim_fit <- function(func, fit_params, model_params, dat_params, initial_params) {
 
 purrr::walk(here::here(c("Functions/misc_model_funs.R","Functions/alm_core.R","Functions/fit_funs.R")),source)
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
  lower_bounds <- c(.01,.1,10)
  upper_bounds <- c(5,5,2000)
  initialpop=initial_params

  # Define number of particles (NP) and iterations (maxiter)
  NP <- 10 * length(initial_params) + 1
  maxiter <- 200
  strategy=6

  optim_time <- system.time({
    # Using DEoptim instead of JDEoptim
    de_optim_res <- DEoptim(fn = objective_function, lower = lower_bounds, 
                            upper = upper_bounds, control = list(NP = NP, itermax = maxiter,strategy=strategy))
  })

  bestFit <- as.data.frame(as.list(de_optim_res$optim$bestmem)) |> mutate(Value=round(de_optim_res$optim$bestval,3))
  names(bestFit) <- c("c","lr","sigma","value")
  
  
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
    return(NULL) # Return NULL or another error indicator
  })

}

# code to fit models with DEoptim

pacman::p_load(dplyr,tidyr,purrr, furrr,future,here)
ds <- ds |> mutate(sbj=id)

tMax=84
train_dataV <- ds |> filter(condit=="Varied",expMode2=="Train") |> group_by(tr) %>%
  mutate(bandInt2 = sample(rep(c(800, 1000, 1200), each = tMax / 3), tMax, replace = FALSE)[tr]) %>%
  filter(bandInt2 == x) |> select(-bandInt2)

# fit params specifies which data is used for optimizing, and which function generates predictions. 
fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll2", loss_data="test_error")
model_params = tibble::lst(input.layer = c(100,350,600,800,1000,1200), 
                           output.layer = c(100,350,600,800,1000,1200), 
                           trainVec=c(0,800,1000,1200))

initial_params <- c(c = 0.01, lr = .5, sigma = 200) 

split_data <- split(ds, ds$id)

plan(multisession)

de_ex_te <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))


# ALM and EXAM implementation
alm.responseOnly <- function(input = 1, c, input.layer, output.layer,weight.mat, trainVec=NULL) {
  input.activation <- exp(-c * (input.layer - input)^2) / (sum(exp(-c * (input.layer - input)^2)) +.01)
  output.activation <- (weight.mat %*% input.activation) 
  output.probability <- output.activation / sum(output.activation) +.01
  sum(output.layer * output.probability)
}


alm.response <- function(input = 1, c, input.layer, output.layer,weight.mat, trainVec=NULL) {
  input.activation <- exp(-c * (input.layer - input)^2) / (sum(exp(-c * (input.layer - input)^2)) +.01)
  output.activation <- (weight.mat %*% input.activation) 
  output.probability <- output.activation / sum(output.activation) +.01
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
  st <- numeric(n) 
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