




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


fit_w_c <- function(params_list,pred_dat,pred_fun,loss_fun,loss_dat,model_params,test_data,train_data)
{
  c = params_list[1] |> as.numeric()
  lr = params_list[2] |> as.numeric()
  list2env(model_params,envir = environment())
  train_results <- alm.sim(train_data, c, lr, input.layer, output.layer)
  weight.mat <- train_results$wm
  test_prediction <- map_dbl(pred_dat$x, ~ pred_fun(.x, c, input.layer, output.layer, weight.mat,  trainVec=trainVec))
 
  if(length(params_list) >= 3 && !is.null(params_list[3])){
  train_error <- loss_fun(train_results$d$x, train_results$d$almResp,params_list[3] |> as.numeric())
  test_error <- loss_fun(test_prediction, pred_dat$y,params_list[3] |> as.numeric())
  } else{
    train_error <- loss_fun(train_results$d$x, train_results$d$almResp)
    test_error <- loss_fun(test_prediction, pred_dat$y)
  }
  
  error <- eval(parse(text=loss_dat))
  
}


optim_fit <- function(func, fit_params, model_params,dat_params,opt.m="BFGS",initial_params) {
  pred_dat <- get(fit_params$pred_dat)
  loss_dat <- fit_params$loss_data
  pred_fun <- match.fun(fit_params$pred_fun)
  loss_fun <- match.fun(fit_params$loss_fun)
  train_data <- dat_params$train_data
  test_data <- dat_params$test_data
  initial_params <- initial_params

  
  
  optim_res <- optim(initial_params, func, pred_dat = pred_dat, pred_fun = pred_fun, 
                     loss_fun = loss_fun, loss_dat = loss_dat, 
                     model_params = model_params, test_data = test_data, train_data = train_data,
                     method = opt.m) # Choose an appropriate method
  
  return(list(optimal_params = optim_res$par, value = optim_res$value)) # Return desired results
}




initial_params <- c(c = 0.5, lr = 2.0, sigma = 175) 


dat_params <- tibble::lst(test_data=test_dataV,train_data=train_dataV)
fit_params=list(pred_dat="test_dataV", pred_fun="exam.response", loss_fun="nll", loss_data="test_error+train_error")

opt_results <- optim_fit(fit_w_c, fit_params,model_params,dat_params,opt.m="BFGS",initial_params)










c_values <- seq(0.000001, 1.0, length.out=25)
lr_values <- seq(0.000001, 4.0, length.out=25)
s_values=seq(50,300,length.out=20)
grid <- expand.grid(c = c_values, lr = lr_values)
model_params = tibble::lst(input.layer = c(100,350,600,800,1000,1200), 
                           output.layer = c(100,350,600,800,1000,1200), 
                           trainVec=c(800,1000,1200))
dat_params <- tibble::lst(test_data=test_dataV,train_data=train_dataV)

fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="test_error+train_error")
grid <- expand.grid(c = c_values, lr = lr_values,sigma=s_values)
ll1 <- grid_fit(grid,fit_w_c,fit_params,model_params,dat_params)
ll1$Fit


fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="RMSE", loss_data="test_error+train_error")
grid <- expand.grid(c = c_values, lr = lr_values)
gf1 <- grid_fit(grid,fit_w_c,fit_params,model_params,dat_params)
gf1$Fit



nll(ll1$test$y,ll1$test$pred,ll1$Fit$sigma)










#############

pacman::p_load(dplyr,tidyr,purrr, furrr,future,here)
purrr::walk(here::here(c("Functions/Display_Functions.R","Functions/misc_model_funs.R","Functions/alm_core.R")),source)

ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))
dsAvg <- ds |> group_by(condit,expMode2,tr, x) |> 
  summarise(y=mean(y),.groups="keep") 

vAvg <- dsAvg |> filter(condit=="Varied")
cAvg <- dsAvg |> filter(condit=="Constant")

dat <- vAvg


trainVec <- sort(c(0,unique(train_data$x)))


tMax=84
train_dataV <- vAvg |> filter(condit=="Varied",expMode2=="Train") |> group_by(tr) %>%
  mutate(bandInt2 = sample(rep(c(800, 1000, 1200), each = tMax / 3), tMax, replace = FALSE)[tr]) %>%
  filter(bandInt2 == x) |> select(-bandInt2)

test_dataV <- vAvg |> filter(condit=="Varied",expMode2=="Test")
test_avg <- test_dataV |> group_by(x) |> summarise(y=mean(y))


e_testOnly=list(pred_dat="test_avg",pred_fun=exam.response,loss_fun="RMSE",loss_data="test_error")



ggplot(train_results$d,aes(x=tr,y=almResp)) + geom_smooth(aes(group=x, col=factor(x))) 


sigma=100
c=.047; lr=1.03
train_results <- alm.sim(train_dataV, c, lr, input.layer, output.layer)


test_prediction <- map_dbl(test_avg$x, ~ pred_fun(.x, c, input.layer, output.layer, train_results$wm,  trainVec=trainVec))
densities <- dnorm(test_avg$y, mean = test_prediction, sd = sigma)
log_likelihoods <- log(densities)

nll(test_avg$y, test_prediction, sigma)

test_avg$y; 
test_prediction
dnorm(test_avg$y, test_prediction,sd=200 )

dnorm(663.0304,715.625,sigma)
dnorm(817,764,sigma)
dnorm(934,883,sigma)
dnorm(1195.5169,1000.006,sigma)
dnorm(test_avg$y, mean = test_prediction, sd = sigma)


avgTrain <- ds |> filter(expMode2=="Train",tr<=82) %>% group_by(id,condit,tr,x,y) %>% 
  summarise(y=mean(y),dist=abs(y-x)) %>% ungroup() %>% group_by(condit,tr,x) %>% summarise(y=mean(y),dist=mean(dist)) %>% ungroup()
ggplot(avgTrain,aes(x=tr,y=y)) + geom_smooth(aes(group=x, col=factor(x))) +facet_grid(~condit) 
ggplot(avgTrain,aes(x=tr,y=dist)) + geom_smooth(aes(group=x, col=factor(x))) +facet_grid(~condit) 


dnorm(663,715,sd=200)



ggplot(train_dataV,aes(x=tr,y=y)) + geom_smooth(aes(group=x, col=factor(x))) 


input.layer <- c(100,350,600,800,1000,1200)
output.layer <- c(100,350,600,800,1000,1200)



