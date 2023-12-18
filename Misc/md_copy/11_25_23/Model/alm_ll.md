


```{r}
runif(10)
```
```




```{r}
pacman::p_load(dplyr,tidyr,purrr, furrr,future,here)

purrr::walk(here::here(c("Functions/Display_Functions.R","Functions/misc_model_funs.R",
                         "Functions/alm_core.R","Functions/fit_funs.R")),source)
ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))

ds <- ds |> mutate(sbj=id)
dsAvg <- ds |> group_by(condit,expMode2,tr, x) |> 
  summarise(y=mean(y),.groups="drop") 

vAvg <- dsAvg |> filter(condit=="Varied")
cAvg <- dsAvg |> filter(condit=="Constant")

tMax=84
train_dataV <- ds |> filter(condit=="Varied",expMode2=="Train") |> group_by(tr) %>%
  mutate(bandInt2 = sample(rep(c(800, 1000, 1200), each = tMax / 3), tMax, replace = FALSE)[tr]) %>%
  filter(bandInt2 == x) |> select(-bandInt2)

test_dataV <- ds |> filter(condit=="Varied",expMode2=="Test")
test_avgV <- ds |> group_by(x) |> summarise(y=mean(y))

train_dataC <- ds |> filter(condit=="Constant",expMode2=="Train") 
test_dataC <- ds |> filter(condit=="Constant",expMode2=="Test")
test_avgC <- test_dataC |> group_by(x) |> summarise(y=mean(y))





alm_result_ql <- function(result){
  
  nllTest <- nll(result$test$y,result$test$pred,result$Fit$sigma)
  nllTrain <- nll(result$train$y,result$train$almResp,result$Fit$sigma)
  
  
  
  
 print(paste0("test nll:",nllTest))
 print(paste0("train nll:",nllTrain))

}
#map(results,alm_result_ql)

```

[1] "c = 0.073, lr = 1.011, sigma = 1249.683, Value = 48.354,  time:0.4"


```{r}

library(DEoptimR)

de_optim_fit <- function(func, fit_params, model_params, dat_params, initial_params) {
   purrr::walk(here::here(c("Functions/misc_model_funs.R",
                         "Functions/alm_core.R","Functions/fit_funs.R")),source)
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


  # Define the objective function for DEoptimR
  objective_function <- function(params) {
    model_fun(params, pred_dat, pred_fun, loss_fun, loss_dat, model_params, test_data, train_data)
  }

  # Define bounds for parameters
  lower_bounds <- c(.0001, .001, 250) 
  upper_bounds <- c(15, 15, 2000)  

  # Define number of particles (NP) and iterations (maxiter)
  NP <- 2 * length(initial_params) + 1
  maxiter <- 3000

  optim_time <- system.time({
  de_optim_res <- JDEoptim(lower = lower_bounds, upper = upper_bounds, fn = objective_function, 
                        NP = NP, maxiter = maxiter)
  })
  
  
  bestFit <- as.data.frame(as.list(de_optim_res$par)) |> mutate(Value=round(de_optim_res$value,3))
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
  
  
  # Wrap the optimization call in tryCatch
  result <- tryCatch({
    de_optim_fit(func, fit_params, model_params, dat_params, initial_params)
  }, error = function(e) {
    message("Error during optimization for ID ", unique(id_data$id), ": ", e$message)
    return(NULL) # Return NULL or another error indicator
  })

}



```

id 1 - test_error
           c        lr    sigma  value
1 0.05944775 0.9984982 999.6437 47.047

```{r}


fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="test_error")
model_params = tibble::lst(input.layer = c(100,350,600,800,1000,1200), 
                           output.layer = c(100,350,600,800,1000,1200), 
                           trainVec=c(0,800,1000,1200))

initial_params <- c(c = 0.01, lr = .5, sigma = 200) 
any(!is.finite(initial_params))


split_data <- split(ds, ds$id)
split_data <- split_data[c(1,2,3)]





plan(multisession)
RNGkind("L'Ecuyer-CMRG")
results <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, opt.m = "Nelder-Mead", initial_params))
saveRDS(results,here::here(paste0("data/model_cache/indv_exam_tetr_bfgs",format(Sys.time(), "%H_%M_%OS"),".rds")))


plan(multisession)
RNGkind("L'Ecuyer-CMRG")
set.seed(123) 

fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="test_error")
de_results <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c1", fit_params, model_params, initial_params))
saveRDS(results,here::here(paste0("data/model_cache/indv_exam_tetr_bfgs",format(Sys.time(), "%H_%M_%OS"),".rds")))








de_results[[1]]$test










opt.m="BFGS"







params_list = initial_params
id_data=split_data[[1]]
func="fit_lr_c"

train_data <- train_data |> mutate(y=ifelse(y>1600,1600,y))
test_data <- test_data |> mutate(y=ifelse(y>1600,1600,y))
```

[1] "c = 0.02, lr = 0.75, sigma = 866.564, value = 46.221,  time:1.4"
[1] "c = 0.069, lr = 1.447, sigma = 999.98, value = 47.363,  time:1.3"

> de_results[[1]]$test
  sbj    x condit         y     pred
1   1  100 Varied  564.8623  717.994
2   1  350 Varied  791.2064  825.478
3   1  600 Varied  984.5210  934.672
4   1  800 Varied 1060.7720 1018.950
5   1 1000 Varied 1187.2458 1096.392
6   1 1200 Varied 1378.7755  999.988


[1] "c = 0.006, lr = 4.653, sigma = 901.095, value = 46.596,  time:1.4"
[1] "c = 10.548, lr = 0.485, sigma = 1000, value = 7988933.002,  time:1.4"
> de_results[[1]]$test
  sbj    x condit         y     pred
1   1  100 Varied  564.8623  733.242
2   1  350 Varied  791.2064  878.848
3   1  600 Varied  984.5210 1156.855
4   1  800 Varied 1060.7720 1140.938
5   1 1000 Varied 1187.2458  595.414
6   1 1200 Varied 1378.7755 1348.224




```{r}



plan(multisession)
RNGkind("L'Ecuyer-CMRG")
set.seed(123) 
split_data <- split(ds, ds$id)


fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="test_error")
de_ex_te <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c1", fit_params, model_params, initial_params))
saveRDS(de_ex_te,here::here(paste0("data/model_cache/indv_nll_de_ex_te",format(Sys.time(), "%H_%M_%OS"),".rds")))

fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="test_error+train_error")
de_ex_tetr <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c1", fit_params, model_params, initial_params))
saveRDS(de_ex_tetr,here::here(paste0("data/model_cache/indv_nll_de_ex_tetr",format(Sys.time(), "%H_%M_%OS"),".rds")))

fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="train_error")
de_ex_tr <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c1", fit_params, model_params, initial_params))
saveRDS(de_ex_tr,here::here(paste0("data/model_cache/indv_exam_de_ex_tr",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
RNGkind("L'Ecuyer-CMRG")
set.seed(123) 
split_data <- split(ds, ds$id)

fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll", loss_data="test_error")
de_alm_te <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c1", fit_params, model_params, initial_params))
saveRDS(de_alm_te,here::here(paste0("data/model_cache/indv_nll_de_alm_te",format(Sys.time(), "%H_%M_%OS"),".rds")))

fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll", loss_data="test_error+train_error")
de_alm_tetr <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c1", fit_params, model_params, initial_params))
saveRDS(de_alm_tetr,here::here(paste0("data/model_cache/indv_nll_de_alm_tetr",format(Sys.time(), "%H_%M_%OS"),".rds")))

fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll", loss_data="train_error")
de_alm_tr <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c1", fit_params, model_params, initial_params))
saveRDS(de_alm_tr,here::here(paste0("data/model_cache/indv_nll_de_alm_tr",format(Sys.time(), "%H_%M_%OS"),".rds")))



plan(multisession)
RNGkind("L'Ecuyer-CMRG")
set.seed(123) 
split_data <- split(ds, ds$id)




fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="RMSE", loss_data="test_error")
de_ex_te_RMSE <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c1", fit_params, model_params, initial_params))
saveRDS(de_ex_te_RMSE,here::here(paste0("data/model_cache/indv_RMSE_de_ex_te",format(Sys.time(), "%H_%M_%OS"),".rds")))

fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="RMSE", loss_data="test_error+train_error")
de_ex_tetr_RMSE <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c1", fit_params, model_params, initial_params))
saveRDS(de_ex_tetr_RMSE,here::here(paste0("data/model_cache/indv_RMSE_de_ex_tetr",format(Sys.time(), "%H_%M_%OS"),".rds")))

fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="RMSE", loss_data="train_error")
de_ex_tr_RMSE <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c1", fit_params, model_params, initial_params))
saveRDS(de_ex_tr,here::here(paste0("data/model_cache/indv_RMSE_de_ex_tr",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
RNGkind("L'Ecuyer-CMRG")
set.seed(123) 
split_data <- split(ds, ds$id)


fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="RMSE", loss_data="test_error")
de_alm_te_RMSE <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c1", fit_params, model_params, initial_params))
saveRDS(de_alm_te_RMSE,here::here(paste0("data/model_cache/indv_RMSE_de_alm_te",format(Sys.time(), "%H_%M_%OS"),".rds")))

fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="RMSE", loss_data="test_error+train_error")
de_alm_tetr_RMSE <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c1", fit_params, model_params, initial_params))
saveRDS(de_alm_tetr_RMSE,here::here(paste0("data/model_cache/indv_RMSE_de_alm_tetr",format(Sys.time(), "%H_%M_%OS"),".rds")))

fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="RMSE", loss_data="train_error")
de_alm_tr_RMSE <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c1", fit_params, model_params, initial_params))
saveRDS(de_alm_tr_RMSE,here::here(paste0("data/model_cache/indv_RMSE_de_alm_tr",format(Sys.time(), "%H_%M_%OS"),".rds")))


```







```{r}


c_values <- seq(0.000001, 1.0, length.out=15)
lr_values <- seq(0.000001, 4.0, length.out=15)
s_values=seq(50,300,length.out=10)
grid <- expand.grid(c = c_values, lr = lr_values)
grid <- expand.grid(c = c_values, lr = lr_values, s_values)
model_params = tibble::lst(input.layer = c(100,350,600,800,1000,1200), 
                           output.layer = c(100,350,600,800,1000,1200), 
                           trainVec=c(800,1000,1200))
dat_params <- tibble::lst(test_data=test_dataV,train_data=train_dataV)

initial_params <- c(c = 0.1, lr = 1.0, sigma = 100) 


dat_params <- tibble::lst(test_data=test_dataV,train_data=train_dataV)
opt_results <- optim_fit(fit_w_c, fit_params,model_params,dat_params,opt.m="BFGS",initial_params)





```

