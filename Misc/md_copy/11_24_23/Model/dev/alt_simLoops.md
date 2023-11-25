
source(here::here("Model/dev/funs_cm.r"))
pacman::p_load(dplyr,tidyr,purrr)

# Constants
nz <- .0001

normalize_vector <- function(x) x / (sum(x) + nz)

activation_function <- function(input, layer, c) exp(-c * (layer - input)^2)

alm.response <- function(input = 1, c = 1, input.layer, output.layer, weight.mat, trainVec=NULL) {
  input.activation <- normalize_vector(activation_function(input, input.layer, c))
  output.activation <- (weight.mat %*% input.activation)
  output.probability <- normalize_vector(output.activation)
  mean.response <- sum(output.layer * output.probability)
  list(mean.response = mean.response, input.activation = input.activation, output.activation = output.activation)
}

alm.responseOnly <- function(input = 1, c = 1, input.layer, output.layer, weight.mat, trainVec=NULL) {
  input.activation <- normalize_vector(activation_function(input, input.layer, c))
  output.activation <- weight.mat %*% input.activation
  output.probability <- normalize_vector(output.activation)
  mean.response <- sum(output.layer * output.probability)
  return(mean.response)
}
 

alm.update <- function(corResp, c = 1, lr, output.layer, input.activation, output_act, weight.mat) {
  teacherSignal <- (activation_function(corResp,output.layer, c) - output_act) * lr
  wChange <- teacherSignal %*% t(input.activation)
  weight.mat <- weight.mat + wChange
  return(weight.mat)
}

alm.trial <- function(input, corResp, c = 1, lr, input.layer, output.layer, weight.mat) {
  alm_resp <- alm.response(input, c, input.layer, output.layer, weight.mat)
  updated_weight.mat <- alm.update(corResp, c, lr, output.layer, alm_resp$input.activation, alm_resp$output.activation, weight.mat)
  return(list(mean.response = alm_resp$mean.response, weight.mat = updated_weight.mat))
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

alt_exam <- function(input, c, input.layer = input_layer, output.layer = output_layer, weight.mat, trainVecX, trainVecY) {
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


predict_model <- function(modelType, params, agent_data, weight.mat, trainVec,trainVecY) {
test_inputs <- agent_data$x[agent_data$modelType == modelType & agent_data$expMode == 'test']
if (modelType == 'ALM') {
return(sapply(test_inputs, function(x) alm.responseOnly(x, params$c, params$input_layer, params$output_layer, weight.mat)))
} else if (modelType == 'EXAM') {
return(sapply(test_inputs, function(x) exam.response(x, params$c, params$input_layer, params$output_layer, weight.mat, trainVec)))
} else if (modelType == 'exam2') {
return(sapply(test_inputs, function(x) alt_exam(x, params$c, params$input_layer, params$output_layer, weight.mat, trainVecX = trainVec, trainVecY)))
} else if (modelType == 'exam0') {
return(sapply(test_inputs, function(x) exam.response(x, params$c, params$input_layer, params$output_layer, weight.mat, trainVec = c(0, trainVec))))
}
}

generate_simulated_data <- function(params, n_agents = 100, n_train_trials = 151, trainVec = c(4, 5, 6), testVec = params$input_layer) {
  training_inputs <- rep(trainVec, each = n_train_trials / length(trainVec))
  correct_outputs <- training_inputs * 2
 
  
  train_data <- data.frame(
    tr=rep(seq(1:n_train_trials),n_agents),
    x = rep(NA, n_agents * length(training_inputs)), 
    y = rep(NA, n_agents * length(correct_outputs)), 
    expMode = 'train', 
    modelType = 'ALM', 
    model_resp = NA_real_,
    agent = rep(1:n_agents, each = n_train_trials)) |> 
    relocate(agent,tr, expMode,modelType,x,y,model_resp) 

  test_data <- expand.grid(x = testVec, 
                           modelType = params$models, 
                           agent=1:n_agents,
                           expMode = 'test',
                           stringsAsFactors = FALSE) |> 
                mutate(tr=rep(seq(1:(length(testVec)*length(params$models))),n_agents),
                y= x * 2, model_resp = NA_real_) |> 
                relocate(agent,tr,expMode,modelType,x,y,model_resp)
                            
  
  simulated_data <- purrr::map(1:n_agents, function(i) {
    set.seed(123 + i)  # Different seed for each agent
    random_indices <- sample(length(training_inputs))
    # Assign randomized inputs and outputs for each agent
    train_data$x[((i - 1) * length(training_inputs) + 1) : (i * length(training_inputs))] <- training_inputs[random_indices]
    train_data$y[((i - 1) * length(correct_outputs) + 1) : (i * length(correct_outputs))] <- correct_outputs[random_indices]

    weight.mat <- matrix(0.000001, nrow = length(params$output_layer), 
                         ncol = length(params$input_layer))
  
    agent_data <- rbind(
      dplyr::filter(train_data, agent == i),
      dplyr::filter(test_data, agent == i)
    ) %>% mutate(c = params$c, lr = params$lr)
  
    for (t in 1:n_train_trials) {
      trial_result <- alm.trial(agent_data$x[t], agent_data$y[t], params$c, params$lr, params$input_layer, params$output_layer, weight.mat)
      agent_data$model_resp[t] <- trial_result$mean.response
      weight.mat <- trial_result$weight.mat
    }
  
    trainVecY <- agent_data$model_resp[agent_data$expMode == 'train']
   
    for (modelType in params$models) {
      predictions <- predict_model(modelType, params, agent_data, weight.mat, trainVec, trainVecY)
      indices <- which(agent_data$modelType == modelType & agent_data$expMode == 'test')  # Define indices here
      if(length(indices) == length(predictions)) {
        agent_data$model_resp[indices] <- predictions
      } else {
        stop(paste("Mismatch between number of predictions and number of indices for model", modelType))
      }
    }
    return(agent_data |> mutate(trainVec = toString(trainVec)))
  })
  
  return(simulated_data)
}


nll2 <- function(obsv,pred,sigma)
{
  nll= -sum(dnorm(obsv,mean=pred,sd=sigma,log=TRUE)) 
  if (is.nan(nll)) {  nll <- 1e4 }
  return(nll)
}


optim_fit_l_BFGS_b <- function(func, fit_params, model_params, dat_params, initial_params, lower_bounds, upper_bounds, opt.m="L-BFGS-B") {
  
  model_fun = get(func)
  
  optim_res <- optim(
    initial_params, 
    model_fun, 
    pred_dat = fit_params$pred_dat, 
    pred_fun = fit_params$pred_fun, 
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

  bf_train_result <- alm.sim(dat_params$train_data, bestFit$c, bestFit$lr, input.layer = model_params$input.layer, output.layer = model_params$output.layer)
  best_preds <- map_dbl(fit_params$pred_dat$x, ~  match.fun(fit_params$pred_fun)(.x, bestFit$c, model_params$input.layer, output.layer = model_params$output.layer, bf_train_result$wm, dat_params$trainVec))
  test_preds_df <- cbind(fit_params$pred_dat, pred = best_preds)
  
  if (any(is.nan(best_preds))) {
    warning("Warning: test_prediction contains NaN values.")
  }

  train_error <- fit_params$loss_fun(bf_train_result$d$y, bf_train_result$d$almResp, bestFit$sigma) |> as.numeric()
  test_error <- fit_params$loss_fun(fit_params$pred_dat$y, best_preds, bestFit$sigma) |> as.numeric()

  optim_info <- list(
    id = as.numeric(dat_params$train_data$sbj[1]),
    condit = as.character(dat_params$train_data$condit[1]),
    Value = bestFit$Value, 
    c = bestFit$c, 
    lr = bestFit$lr, 
    NAN_Test = any(is.nan(best_preds)),
    opt.m = opt.m,
    # model_fun = func, 
    # pred_fun = fit_params$pred_fun,
    # loss_dat = fit_params$loss_data, 
    # loss_fun = fit_params$loss_fun, 
    train = bf_train_result$d, 
    test = test_preds_df, 
    Fit = bestFit, 
    c = bestFit$c, 
    lr = bestFit$lr, 
    Value = bestFit$Value, 
    train_error = train_error,
    test_error = test_error,
    run_args = tibble::lst(func, fit_params, model_params, dat_params, opt.m, initial_params, 
        trainVec = dat_params$trainVec, model_params$input.layer, model_params$output.layer),
    opt_info = optim_res,
    info = tibble::lst(time = Sys.time(), Sys.info(), path = getwd()),
    runtime = optim_res$time
  )
  
  return(optim_info)
}

fit_lr_c1 <- function(params_list, pred_dat, pred_fun, loss_fun, loss_dat, model_params, test_data, train_data) {
  
  c <- params_list[1]
  lr <- params_list[2]
  train_results <- alm.sim(train_data, c, lr, model_params$input.layer, model_params$output.layer)
  test_prediction <- map_dbl(pred_dat$x, ~ pred_fun(.x, c, model_params$input.layer, model_params$output.layer, train_results$wm, trainVec = sort(unique(train_data$x))))
  test_prediction[is.nan(test_prediction)] <- 2000000

  params_sigma <- ifelse(length(params_list) >= 3 && !is.null(params_list[3]), params_list[3], NULL)

  train_error <- loss_fun(train_results$d$y, train_results$d$almResp, params_sigma)
  test_error <- loss_fun(pred_dat$y, test_prediction, params_sigma)
  error <- eval(parse(text = loss_dat))
}




fit_model <- function(simulated_data,modelType,loss_fun) {

  initial_params <- c(c = 0.1, lr = 1, sigma = 1)
  lower_bounds <- c(0.0001, 0.001, .05)      
  upper_bounds <- c(5, 10, 15)     

  train_data <- simulated_data[simulated_data$expMode == 'train', ]
  test_data <- simulated_data[simulated_data$expMode == 'test' & simulated_data$modelType == modelType, ]

  prediction_functions <- list(
    ALM = alm.responseOnly,
    EXAM = exam.response,
    exam2 = alt_exam,
    exam0 = function(x, c, input_layer, output_layer, weight.mat, trainVec) {
      exam.response(x, c, input_layer, output_layer, weight.mat, trainVec = c(0, trainVec))
    }
  )
       
  #print(loss_fun)
  optim_res <- optim_fit_l_BFGS_b(
    func = "fit_lr_c1", 
    fit_params = list(
      pred_dat = test_data, 
      pred_fun = exam.response, 
      loss_fun = loss_fun, 
      loss_data = "test_error"
    ),
    model_params = list(input.layer = c(1:6), output.layer = c(1:6) * 2),
    dat_params = list(
      train_data = train_data, 
      test_data = test_data,
      trainVec = sort(unique(train_data$x))
    ),
    initial_params = initial_params, 
    lower_bounds = lower_bounds, 
    upper_bounds = upper_bounds
  )

  optim_res$Fit
}

simulate_and_fit <- function(c, lr, trainVec,loss_fun) {
  params <- list(
    c = c,
    lr = lr,
    input_layer = c(1, 2, 3, 4, 5, 6),
    output_layer = c(1, 2, 3, 4, 5, 6) * 2,
    models = c('ALM', 'EXAM', 'exam0', 'exam2')
  )
  simulated_data_list <- generate_simulated_data(params, n_agents = 1, n_train_trials = 90, trainVec = c(4, 5, 6), testVec = params$input_layer)

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

c_values <- c(.05, .1, .5)
lr_values <- c(.5, 1)
trainVec_values <- list(c(4, 5, 6), c(4))

param_combinations <- expand.grid(c = c_values, lr = lr_values, trainVec = trainVec_values, stringsAsFactors = FALSE)

all_fit_results <- purrr::pmap_dfr(param_combinations, ~simulate_and_fit(..1, ..2, ..3, loss_fun=nll2))





# c_values <- c(.05, .1, .5)
# lr_values <- c(.5, 1)
# trainVec_values <- list(c(4,5,6), c(4))

# param_combinations <- expand.grid(c = c_values, lr = lr_values, trainVec = trainVec_values, stringsAsFactors = FALSE)

# all_simulated_data <- purrr::pmap_dfr(param_combinations, function(c, lr, trainVec) {
#   params <- list(
#     c = c,
#     lr = lr,
#     input_layer = c(1,2,3,4,5,6),
#     output_layer = c(1,2,3,4,5,6)*2,
#     models = c('ALM', 'EXAM', 'exam0', 'exam2')
#   )
  
#   generate_simulated_data(params, n_agents = 10, n_train_trials = 90, trainVec = trainVec, testVec=params$input_layer)
# })
# head(all_simulated_data)