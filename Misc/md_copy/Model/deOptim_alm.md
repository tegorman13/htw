



```{r}
pacman::p_load(dplyr,tidyr,purrr, furrr,future,here)

purrr::walk(here::here(c("Functions/Display_Functions.R","Functions/misc_model_funs.R",
                         "Functions/alm_core.R","Functions/fit_funs.R")),source)
ds <- readRDS(here::here("data/e1_md_11-06-23.rds")) |> mutate(sbj=id)

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








```

de_results[[1]]$Fit
           c        lr    sigma  value
1 0.05578472 0.9892148 162.0428 39.038

de_results[[1]]$test
  sbj    x condit         y     pred
1   1  100 Varied  564.8623  756.010
2   1  350 Varied  791.2064  857.286
3   1  600 Varied  984.5210  948.443
4   1  800 Varied 1060.7720 1039.582
5   1 1000 Varied 1187.2458 1171.196
6   1 1200 Varied 1378.7755 1040.500



```{r}
# 11/16/23 - fairly successful fits

fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="test_error")
model_params = tibble::lst(input.layer = c(100,350,600,800,1000,1200), 
                           output.layer = c(100,350,600,800,1000,1200), 
                           trainVec=c(0,800,1000,1200))

initial_params <- c(c = 0.01, lr = .5, sigma = 200) 
any(!is.finite(initial_params))


split_data <- split(ds, ds$id)
split_data <- split_data[c(1,2)]



# plan(multisession)
# fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll2", loss_data="test_error")
# de_results <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))


plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll2", loss_data="test_error")
de_ex_te <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_ex_te,here::here(paste0("data/model_cache/indv_nll_de2_ex_te",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll2", loss_data="test_error+train_error")
de_ex_tetr <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_ex_tetr,here::here(paste0("data/model_cache/indv_nll_de2_ex_tetr",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll2", loss_data="train_error")
de_ex_tr <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_ex_tr,here::here(paste0("data/model_cache/indv_exam_de2_ex_tr",format(Sys.time(), "%H_%M_%OS"),".rds")))


plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll2", loss_data="test_error")
de_alm_te <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_alm_te,here::here(paste0("data/model_cache/indv_nll_de2_alm_te",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll2", loss_data="test_error+train_error")
de_alm_tetr <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,initial_params))
saveRDS(de_alm_tetr,here::here(paste0("data/model_cache/indv_nll_de2_alm_tetr",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll2", loss_data="train_error")
de_alm_tr <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_alm_tr,here::here(paste0("data/model_cache/indv_nll_de2_alm_tr",format(Sys.time(), "%H_%M_%OS"),".rds")))



  
  




params_list = initial_params
id_data=split_data[[1]]
func="fit_lr_c"
```




```{r}


parVar=list("params_list","pred_dat","pred_fun","loss_fun","loss_dat","model_params","test_data","train_data")


de_optim_res <- DEoptim(fn = objective_function, lower = c(.01,.1,10), 
                            upper = c(30,30,2000), control = list(NP = 30, itermax = 100,strategy=6))

de_optim_res <- DEoptim(fn = objective_function, lower = c(.001,.01,10), 
                            upper = c(5,5,2000), control = list(NP = 30, itermax = 100,strategy=6))


  objective_function <- function(params) {
    model_fun(params, pred_dat, pred_fun, loss_fun, loss_dat, model_params, test_data, train_data)
  }


paramEvolve <- DEoptim(fit_lr_c0,
                       lower = c(0, 0, 0, 0),
                       upper = c(10, 10, 20, 0.2),
                       control = list(parallelType = 1, 
                                      itermax = 100,
                                      parVar=parVar))

paramEvolve <- DEoptim(fit_lr_c0,
                       lower = c(0, 0, 0, 0),
                       upper = c(10, 10, 20, 0.2),
                       control = list(parallelType = 1, 
                                      itermax = 100,
                                      params_list=params_list,
                                      pred_dat= pred_dat,
                                      pred_fun=alm.responseOnly,
                                      loss_fun=nll2,
                                      loss_dat=loss_dat,
                                      model_params=model_params,
                                      test_data=test_data,
                                      train_data=train_data))
```

