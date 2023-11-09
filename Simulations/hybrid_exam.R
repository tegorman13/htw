
pacman::p_load(tidyverse,lme4,future,furrr,patchwork,here, furrr, future, pander)
purrr::walk(here(c("Functions/Display_Functions.R", "Functions/alm_core.R","Functions/misc_model_funs.R")),source)
select <- dplyr::select; mutate <- dplyr::mutate 

ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))

dsAvg <- ds |> group_by(condit,expMode2,tr, x) |> 
  summarise(y=mean(y),.groups="keep") 

vAvg <- dsAvg |> filter(condit=="Varied")
cAvg <- dsAvg |> filter(condit=="Constant")

input.layer <- c(0,100,350,600,800,1000,1200)
output.layer <- c(100,350,600,800,1000,1200)

h_testOnly=list(pred_dat="test_avg",pred_fun="predict_alm_exam_weighted_hybrid",loss_fun="RMSE",loss_data="test_error")
h_trainOnly=list(pred_dat="test_avg",pred_fun="predict_alm_exam_weighted_hybrid",loss_fun="RMSE",loss_data="train_error")
h_testTrain=list(pred_dat="test_avg",pred_fun="predict_alm_exam_weighted_hybrid",loss_fun="RMSE",loss_data="test_error+train_error")

c_values <- seq(0.000001, 1.0, length.out=120)
lr_values <- seq(0.000001, 4.0, length.out=120)
weight_values <- seq(0,1,length.out=15)


hybrid_te_v <- wrap_grid_hybrid(vAvg, c_values, lr_values,weight_values, input.layer, output.layer,predParams=h_testOnly)
hybrid_tetr_v <- wrap_grid_hybrid(vAvg, c_values, lr_values,weight_values, input.layer, output.layer,predParams=h_testTrain)
hybrid_tr_v<- wrap_grid_hybrid(vAvg, c_values, lr_values,weight_values, input.layer, output.layer,predParams=h_trainOnly)

hybrid_te_c <- wrap_grid_hybrid(cAvg, c_values, lr_values,weight_values, input.layer, output.layer,predParams=h_testOnly)
hybrid_tetr_c <- wrap_grid_hybrid(cAvg, c_values, lr_values,weight_values, input.layer, output.layer,predParams=h_testTrain)
hybrid_tr_c<- wrap_grid_hybrid(cAvg, c_values, lr_values,weight_values, input.layer, output.layer,predParams=h_trainOnly)


saveRDS(tibble::lst(hybrid_te_v,hybrid_tetr_v,hybrid_tr_v, hybrid_te_c,hybrid_tetr_c,hybrid_tr_c), here::here("data/model_cache/hybrid_group_exam_fits.rds"))


plan(multisession)

wrap_grid_hybrid <- function(dat, c_values, lr_values,weight_values, input.layer, output.layer,
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
  grid <- expand.grid(c = c_values, lr = lr_values, w=weight_values,Value = NA_real_)
  
  # Use future_map (parallel version of map) to iterate over rows in the grid
  loop_time <-  system.time({ 
    results <- future_map_dfr(seq_len(nrow(grid)), function(i) {
      c <- grid$c[i]
      lr <- grid$lr[i]
      w <- grid$w[i]
      train_results <- alm.sim(train_data, c, lr, input.layer, output.layer)
      weight.mat <- train_results$wm
      
      train_error <- loss_fun(train_results$d$x, train_results$d$almResp)

      test_prediction <- map_dbl(pred_dat$x, ~ pred_fun(.x, c, w, input.layer, output.layer, weight.mat,  trainVec=trainVec))
  
      
      test_error <- loss_fun(test_prediction, pred_dat$y)
      
      error <- eval(parse(text=predParams$loss_data)) 
      
      tibble(c = c, lr = lr, w=w, Value = error)  # Return a tibble with the results
    })
  })# stop timer
  
    min_value <- min(results$Value, na.rm = TRUE)

    # Check if the minimum value occurs more than once
    if(sum(results$Value == min_value, na.rm = TRUE) > 1 && !is.nan(min_value)) {
      warning("The minimum value in the grid occurs more than once.")
    }


  # Extract the best fit parameters
  bestFit <- results[which.min(results$Value), ] |>  mutate(across(where(is.numeric), \(x) round(x, 8)))
  
  bf_train_result <- alm.sim(train_data, bestFit$c, bestFit$lr, input.layer, output.layer)
  weight.mat <- bf_train_result$wm
  bf_train <- bf_train_result$d
  best_preds <- map_dbl(pred_dat$x, ~ pred_fun(.x, bestFit$c, bestFit$w,input.layer, output.layer, weight.mat,trainVec))
  
  (bf_test <- cbind(pred_dat, pred = best_preds ))
  print(paste0("Value:", bestFit$Value, " c:", bestFit$c, " lr:", bestFit$lr, "w:",bestFit$w, " time:",  round(loop_time[3],3)))
  
  return(list(train=bf_train,test=bf_test,errorGrid=results, Fit = bestFit, c = bestFit$c, lr = bestFit$lr,w=bestFit$w, Value = bestFit$Value, time=loop_time[3]))
  
}

predict_alm_exam_weighted_hybrid <- function(input, c, weight_exam, input.layer, output.layer, weight.mat,trainVec) {
  
  alm_pred <- alm.responseOnly(input, c, input.layer, output.layer,weight.mat, trainVec=NULL)
  
  exam_pred <- exam.response(input, c, input.layer, output.layer,weight.mat, trainVec=trainVec)
  (1 - weight_exam) * alm_pred + weight_exam * exam_pred
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


