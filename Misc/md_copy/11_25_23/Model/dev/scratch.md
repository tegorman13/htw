pacman::p_load(dplyr,tidyr,purrr,ggplot2)


### Response Functions

exam.response <- function(input, pred_args) {


  c <- pred_args$c
  input.layer <- pred_args$input.layer
  output.layer <- pred_args$output.layer
  weight.mat <- pred_args$weight.mat
  trainVec <- pred_args$trainVec

  nearestTrain <- trainVec[which.min(abs(input - trainVec))]
  aresp <- alm.response(nearestTrain, c, input.layer = input.layer,output.layer,weight.mat)$mean.response
  
  xUnder <- ifelse(min(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) - 1])
  xOver <- ifelse(max(trainVec) == nearestTrain, nearestTrain, trainVec[which(trainVec == nearestTrain) + 1])
  
  mUnder <- alm.response(xUnder, c, input.layer = input.layer, output.layer, weight.mat)$mean.response
  mOver <- alm.response(xOver, c, input.layer = input.layer,output.layer, weight.mat)$mean.response
  
  exam.output <- round(aresp + ((mOver - mUnder) / (xOver - xUnder)) * (input - nearestTrain), 3)
  exam.output
}



alt_exam <- function(input, pred_args) {

      c <- pred_args$c
    input.layer <- pred_args$input.layer
    output.layer <- pred_args$output.layer
    weight.mat <- pred_args$weight.mat
    trainVecX <- pred_args$trainVecX
    trainVecY <- pred_args$trainVecY

       if (length(input) > 1) {
        stop("Input must be scalar")
    }
    if (length(trainVecY) == 0) {
        stop("No training data")
    }

    nearestTrain <- trainVecX[which.min(abs(input - trainVecX))]
    aresp <- alm.response(nearestTrain, c, input.layer, output.layer, weight.mat)$mean.response
    
    if (nearestTrain>input) {
        upper_output <- trainVecY[which(trainVecY > aresp)]
        if (length(upper_output) == 0) {
            return(alm.response(input, c, input.layer, output.layer, weight.mat)$mean.response)
        }
        upper_output <- upper_output[which.max(abs(upper_output - nearestTrain))]
        infer_input <- nearestTrain + (upper_output - aresp)
        slope <- (alm.response(infer_input, c, input.layer, output.layer, weight.mat)$mean.response - aresp) / (infer_input - nearestTrain)
        return(round(aresp + slope * (input - nearestTrain), 3))
    }

    lower_output <- trainVecY[which(trainVecY < aresp)]
    if (length(lower_output) == 0) {
        return(alm.response(input, c, input.layer, output.layer, weight.mat)$mean.response)
    }
    lower_output <- lower_output[which.min(abs(lower_output - nearestTrain))]
    infer_input <- nearestTrain + (aresp - lower_output)
    slope <- (alm.response(infer_input, c, input.layer, output.layer, weight.mat)$mean.response - aresp) / (infer_input - nearestTrain)
    return(round(aresp + slope * (input - nearestTrain), 3))
}


alm.responseOnly <- function(input = 1, pred_args) {

   c <- pred_args$c
  input.layer <- pred_args$input.layer
  output.layer <- pred_args$output.layer
  weight.mat <- pred_args$weight.mat
  trainVec <- pred_args$trainVec

  input.activation <- normalize_vector(activation_function(input, input.layer, c))
  output.activation <- weight.mat %*% input.activation
  output.probability <- normalize_vector(output.activation)
  mean.response <- sum(output.layer * output.probability)
  return(mean.response)
}



alm.response <- function(input = 1, c = 1, input.layer, output.layer, weight.mat, trainVec=NULL) {
  input.activation <- normalize_vector(activation_function(input, input.layer, c))
  output.activation <- (weight.mat %*% input.activation)
  output.probability <- normalize_vector(output.activation)
  mean.response <- sum(output.layer * output.probability)
  list(mean.response = mean.response, input.activation = input.activation, output.activation = output.activation)
}




c_values <- c(.0005, 1, 5)
lr_values <- c(.05,3)
trainVec_values <- list(c(4, 5, 6), seq(1:6))
param_combinations <- expand.grid(c = c_values, lr = lr_values, trainVec = trainVec_values, stringsAsFactors = FALSE)

all_fit_results <- purrr::pmap_dfr(param_combinations, ~simulate_and_fit(..1, ..2, ..3, loss_fun=nll2))

ggplot(all_fit_results,aes(x=true_c,y=fit_c,col=fit_model))+geom_point()+facet_wrap(~trainVec)



simulate_and_fit <- function(c, lr, trainVec,loss_fun) {
  params <- list(
    c = c,
    lr = lr,
    input_layer = c(1, 2, 3, 4, 5, 6),
    output_layer = c(1, 2, 3, 4, 5, 6) * 2,
    models = c('ALM', 'EXAM', 'exam0', 'exam2')
  )
  
  simulated_data_list <- generate_simulated_data(params, n_agents = 1, n_train_trials = 90, 
  trainVec = trainVec, testVec = seq(1,10,by=.25))

  # fit_results <- map_dfr(simulated_data, ~fit_model(..1,"EXAM", loss_fun)) %>%
  # rename_with(~paste0("fit_",.x)) |> mutate(agent=1,true_c=params$c, true_lr=params$lr, trainVec=toString(trainVec))

  fit_results <- map_dfr(params$models, function(modelType) {
    map_dfr(simulated_data_list, function(sim_data) {
      fit_model(sim_data, modelType, loss_fun) |> mutate(model=modelType,id=sim_data$agent[1])
    })
  }) %>%
  rename_with(~paste0("fit_", .x)) |> 
  mutate(true_c = params$c, true_lr = params$lr, trainVec = toString(trainVec))

  return(fit_results)
}

 
## Fitting functions
fit_model <- function(simulated_data,modelType,loss_fun) {

  initial_params <- c(c = 0.1, lr = 1, sigma = 1)
  lower_bounds <- c(0.0001, 0.001, .05)      
  upper_bounds <- c(5, 10, 15)     
  model_params = tibble::lst(input.layer = c(1:6), output.layer = c(1:6) * 2)
  train_data <- simulated_data[simulated_data$expMode == 'train', ]
  test_data <- simulated_data[simulated_data$expMode == 'test' & simulated_data$modelType == modelType, ]

  prediction_functions <- list(
    ALM = alm.responseOnly,
    EXAM = exam.response,
    exam2 = alt_exam,
    exam0 = function(x, pred_args) {
      pred_args$trainVec <- c(0, pred_args$trainVec)
      exam.response(x, pred_args)
    }
  )
  pred_fun <- prediction_functions[[modelType]]

    trainVec <- sort(unique(train_data$x))
  if (modelType == "exam0") {
    trainVec <- c(0, trainVec)
  }

  pred_args <- list(
    input = NULL,  
    c = NULL,
    input.layer = model_params$input.layer,
    output.layer = model_params$output.layer,
    weight.mat = NULL,
    trainVec = trainVec, 
    trainVecX = trainVec,
    trainVecY = NULL
  )

  #print(loss_fun)
  optim_res <- optim_fit_l_BFGS_b(
    func = "fit_lr_c1", 
    fit_params = list(
      pred_dat = test_data, 
      pred_fun = pred_fun, 
      loss_fun = loss_fun, 
      loss_data = "test_error+train_error"
    ),
    model_params = model_params,
    dat_params = list(
      train_data = train_data, 
      test_data = test_data,
      trainVec = sort(unique(train_data$x))
    ),
    pred_args = pred_args,
    initial_params = initial_params, 
    lower_bounds = lower_bounds, 
    upper_bounds = upper_bounds
  )



}



fit_lr_c1 <- function(params_list, pred_dat, pred_fun, pred_args,loss_fun, loss_dat, model_params, test_data, train_data) {
  
  c <- params_list[1]
  lr <- params_list[2]
  train_results <- alm.sim(train_data, c, lr, model_params$input.layer, model_params$output.layer)

  pred_args$c <- c
  pred_args$weight.mat <- train_results$wm
  pred_args$trainVecY <- train_results$d$almResp

  test_prediction <- map_dbl(pred_dat$x, ~ pred_fun(.x, pred_args))  
  
  test_prediction[is.nan(test_prediction)] <- 2000000

  params_sigma <- ifelse(length(params_list) >= 3 && !is.null(params_list[3]), params_list[3], NULL)

  train_error <- loss_fun(train_results$d$y, train_results$d$almResp, params_sigma)
  test_error <- loss_fun(pred_dat$y, test_prediction, params_sigma)
  error <- eval(parse(text = loss_dat))
}


optim_fit_l_BFGS_b <- function(func, fit_params, model_params, dat_params, pred_args, initial_params, lower_bounds, upper_bounds, opt.m="L-BFGS-B") {
  model_fun = get(func)
  
  optim_res <- optim(
    initial_params, 
    model_fun, 
    pred_dat = fit_params$pred_dat, 
    pred_fun = fit_params$pred_fun, 
    pred_args = pred_args,
    loss_fun = fit_params$loss_fun, #match.fun(fit_params$loss_fun), 
    loss_dat = fit_params$loss_data, 
    model_params = model_params, 
    test_data = dat_params$test_data, 
    train_data = dat_params$train_data,
    method = opt.m, 
    lower = lower_bounds, 
    upper = upper_bounds
  )
  bestFit <- as.data.frame(as.list(optim_res$par)) %>% mutate(Value = round(optim_res$value, 3))
  
  return(bestFit)
}


# predict_model <- function(modelType) {
# test_inputs <- agent_data$x[agent_data$modelType == modelType & agent_data$expMode == 'test']
# if (modelType == 'ALM') {
# return(sapply(test_inputs, function(x) alm.responseOnly(x, params$c, params$input_layer, params$output_layer, weight.mat)))
# } else if (modelType == 'EXAM') {
# return(sapply(test_inputs, function(x) exam.response(x, params$c, params$input_layer, params$output_layer, weight.mat, trainVec)))
# } else if (modelType == 'exam2') {
# return(sapply(test_inputs, function(x) alt_exam(x, params$c, params$input_layer, params$output_layer, weight.mat, trainVecX = trainVec, trainVecY)))
# } else if (modelType == 'exam0') {
# return(sapply(test_inputs, function(x) exam.response(x, params$c, params$input_layer, params$output_layer, weight.mat, trainVec = c(0, trainVec))))
# }
# }
# predict_model <- function(modelType, params, agent_data, weight.mat, trainVec,trainVecY) {
# test_inputs <- agent_data$x[agent_data$modelType == modelType & agent_data$expMode == 'test']
# if (modelType == 'ALM') {
# return(sapply(test_inputs, function(x) alm.responseOnly(x, params$c, params$input_layer, params$output_layer, weight.mat)))
# } else if (modelType == 'EXAM') {
# return(sapply(test_inputs, function(x) exam.response(x, params$c, params$input_layer, params$output_layer, weight.mat, trainVec)))
# } else if (modelType == 'exam2') {
# return(sapply(test_inputs, function(x) alt_exam(x, params$c, params$input_layer, params$output_layer, weight.mat, trainVecX = trainVec, trainVecY)))
# } else if (modelType == 'exam0') {
# return(sapply(test_inputs, function(x) exam.response(x, params$c, params$input_layer, params$output_layer, weight.mat, trainVec = c(0, trainVec))))
# }
# }


predict_model <- function(modelType, params, agent_data, weight.mat, trainVec, trainVecY) {
  test_inputs <- agent_data$x[agent_data$modelType == modelType & agent_data$expMode == 'test']
  
  pred_args <- list(
    c = params$c,
    input.layer = params$input_layer,
    output.layer = params$output_layer,
    weight.mat = weight.mat,
    trainVec = trainVec,
    trainVecY = trainVecY
  )
  
  if (modelType == 'ALM') {
    return(sapply(test_inputs, function(x) alm.responseOnly(x, pred_args)))
  } else if (modelType == 'EXAM') {
    return(sapply(test_inputs, function(x) exam.response(x, pred_args)))
  } else if (modelType == 'exam2') {
    pred_args$trainVecX <- trainVec  # Specific to alt_exam function
    return(sapply(test_inputs, function(x) alt_exam(x, pred_args)))
  } else if (modelType == 'exam0') {
    pred_args$trainVec <- c(0, trainVec)  # Specific to exam0 function
    return(sapply(test_inputs, function(x) exam.response(x, pred_args)))
  }
}












simulate_and_fit_alm <- function(params) {
  # Code to simulate ALM data
  # Code to fit ALM model
}

simulate_and_fit_exam <- function(params) {
  # Code to simulate EXAM data
  # Code to fit EXAM model 
}






