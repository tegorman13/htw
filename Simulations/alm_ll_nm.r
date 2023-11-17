

pacman::p_load(tidyverse,here,future,furrr)
purrr::walk(here::here(c("Functions/fit_funs.R", "Functions/alm_core.R","Functions/misc_model_funs.R")),source)

ds <- readRDS(here::here("data/e1_md_11-06-23.rds")) |> mutate(sbj=id) |> relocate(sbj,.after=id)


model_params = tibble::lst(input.layer = c(100,350,600,800,1000,1200), 
                           output.layer = c(100,350,600,800,1000,1200), 
                           trainVec=c(0,800,1000,1200))

initial_params <- c(c = 0.01, lr = .5, sigma = 200) 
#any(!is.finite(initial_params))


split_data <- split(ds, ds$id)
split_data <- split_data[c(1)]

params_list = initial_params
id_data=split_data[[1]]
func="fit_lr_c0"


fit_params = list(pred_dat="test_data", pred_fun="exam.response", loss_fun="nll2", loss_data="test_error")
lower_bounds = c(.0000001, .000001, 100) # example lower bounds
upper_bounds = c(3,10, 2000) # example upper bounds
initial_params <- c(c = 0.01, lr = .5, sigma = 200) 





###########################

library(DEoptim)

plan(multisession)
fit_params=list(pred_dat="test_data", pred_fun="exam.response", loss_fun="nll2", loss_data="test_error")
de_ex_te_raw <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_ex_te_raw,here::here(paste0("data/model_cache/indv_nll_raw_de2_ex_te_",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_data", pred_fun="exam.response", loss_fun="nll2", loss_data="test_error+train_error")
de_ex_tetr_raw <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_ex_tetr_raw,here::here(paste0("data/model_cache/indv_nll_raw_de2_ex_tetr_",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_data", pred_fun="exam.response", loss_fun="nll2", loss_data="train_error")
de_ex_tr_raw <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_ex_tr_raw,here::here(paste0("data/model_cache/indv_nll_raw_de2_ex_tr_",format(Sys.time(), "%H_%M_%OS"),".rds")))


plan(multisession)
fit_params=list(pred_dat="test_data", pred_fun="alm.responseOnly", loss_fun="nll2", loss_data="test_error")
de_alm_te_raw <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_alm_te_raw,here::here(paste0("data/model_cache/indv_nll_raw_de2_alm_te_",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_data", pred_fun="alm.responseOnly", loss_fun="nll2", loss_data="test_error+train_error")
de_alm_tetr_raw <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,initial_params))
saveRDS(de_alm_tetr_raw,here::here(paste0("data/model_cache/indv_nll_raw_de2_alm_tetr_",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_data", pred_fun="alm.responseOnly", loss_fun="nll2", loss_data="train_error")
de_alm_tr_raw <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_alm_tr_raw,here::here(paste0("data/model_cache/indv_nll_raw_de2_alm_tr_",format(Sys.time(), "%H_%M_%OS"),".rds")))












##############################

plan(multisession)
fit_params=list(pred_dat="test_data", pred_fun="exam.response", loss_fun="nll2", loss_data="test_error")
nm_ex_te_raw <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(nm_ex_te_raw,here::here(paste0("data/model_cache/indv_nll_raw_nm_ex_te_",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_data", pred_fun="exam.response", loss_fun="nll2", loss_data="test_error+train_error")
nm_ex_tetr_raw <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(nm_ex_tetr_raw,here::here(paste0("data/model_cache/indv_nll_raw_nm_ex_tetr_",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_data", pred_fun="exam.response", loss_fun="nll2", loss_data="train_error")
nm_ex_tr_raw <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(nm_ex_tr_raw,here::here(paste0("data/model_cache/indv_exam_nm_raw_ex_tr_",format(Sys.time(), "%H_%M_%OS"),".rds")))


plan(multisession)
fit_params=list(pred_dat="test_data", pred_fun="alm.responseOnly", loss_fun="nll2", loss_data="test_error")
nm_alm_te_raw <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(nm_alm_te_raw,here::here(paste0("data/model_cache/indv_nll_nm_raw_alm_te_",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_data", pred_fun="alm.responseOnly", loss_fun="nll2", loss_data="test_error+train_error")
nm_alm_tetr_raw <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(nm_alm_tetr_raw,here::here(paste0("data/model_cache/indv_nl_raw_nm_alm_tetr_",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_data", pred_fun="alm.responseOnly", loss_fun="nll2", loss_data="train_error")
nm_alm_tr_raw <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, opt.m, initial_params))
saveRDS(nm_alm_tr_raw,here::here(paste0("data/model_cache/indv_nll_raw_nm_alm_tr_",format(Sys.time(), "%H_%M_%OS"),".rds")))










##################################


plan(multisession)
fit_params = list(pred_dat="test_data", pred_fun="exam.response", loss_fun="nll2", loss_data="test_error")
indv_lbfg_ex_te_raw_nll <- future_map(split_data, 
                            ~ optim_fit_id_l_BFGS_b(.x, func="fit_lr_c0", fit_params, model_params, initial_params, lower_bounds, upper_bounds))
saveRDS(indv_lbfg_ex_te_raw_nll,here::here(paste0("data/model_cache/indv_lbfg_ex_te_raw_nll",format(Sys.time(), "%H_%M_%OS"),".rds")))


plan(multisession)
fit_params = list(pred_dat="test_data", pred_fun="exam.response", loss_fun="nll2", loss_data="test_error+train_error")
indv_lbfg_ex_tetr_raw_nll <- future_map(split_data, 
                                      ~ optim_fit_id_l_BFGS_b(.x, func="fit_lr_c0", fit_params, model_params, initial_params, lower_bounds, upper_bounds))
saveRDS(indv_lbfg_ex_tetr_raw_nll,here::here(paste0("data/model_cache/indv_lbfg_ex_tetr_raw_nll",format(Sys.time(), "%H_%M_%OS"),".rds")))



plan(multisession)
fit_params = list(pred_dat="test_data", pred_fun="exam.response", loss_fun="nll2", loss_data="train_error")
indv_lbfg_ex_tr_raw_nll <- future_map(split_data, 
                                      ~ optim_fit_id_l_BFGS_b(.x, func="fit_lr_c0", fit_params, model_params, initial_params, lower_bounds, upper_bounds))
saveRDS(indv_lbfg_ex_tr_raw_nll,here::here(paste0("data/model_cache/indv_lbfg_ex_tr_raw_nll",format(Sys.time(), "%H_%M_%OS"),".rds")))




plan(multisession)
fit_params = list(pred_dat="test_data", pred_fun="alm.responseOnly", loss_fun="nll2", loss_data="test_error")
indv_lbfg_alm_te_raw_nll <- future_map(split_data, 
                                      ~ optim_fit_id_l_BFGS_b(.x, func="fit_lr_c0", fit_params, model_params, initial_params, lower_bounds, upper_bounds))
saveRDS(indv_lbfg_alm_te_raw_nll,here::here(paste0("data/model_cache/indv_lbfg_alm_te_raw_nll",format(Sys.time(), "%H_%M_%OS"),".rds")))


plan(multisession)
fit_params = list(pred_dat="test_data", pred_fun="alm.responseOnly", loss_fun="nll2", loss_data="test_error+train_error")
indv_lbfg_alm_tetr_raw_nll <- future_map(split_data, 
                                        ~ optim_fit_id_l_BFGS_b(.x, func="fit_lr_c0", fit_params, model_params, initial_params, lower_bounds, upper_bounds))
saveRDS(indv_lbfg_alm_tetr_raw_nll,here::here(paste0("data/model_cache/indv_lbfg_alm_tetr_raw_nll",format(Sys.time(), "%H_%M_%OS"),".rds")))


plan(multisession)
fit_params = list(pred_dat="test_data", pred_fun="alm.responseOnly", loss_fun="nll2", loss_data="train_error")
indv_lbfg_alm_tr_raw_nll <- future_map(split_data, 
                                      ~ optim_fit_id_l_BFGS_b(.x, func="fit_lr_c0", fit_params, model_params, initial_params, lower_bounds, upper_bounds))
saveRDS(indv_lbfg_alm_tr_raw_nll,here::here(paste0("data/model_cache/indv_lbfg_alm_tr_raw_nll",format(Sys.time(), "%H_%M_%OS"),".rds")))



##################################


initial_params <- c(c = 0.01, lr = .5, sigma = 200) 
opt.m="Nelder-Mead"
fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll2", loss_data="test_error")
nm_ex_te_avg<- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
nm_ex_te_avg[[1]]$test
nm_ex_te_avg[[1]]$Fit








saveRDS(nm_ex_te,here::here(paste0("data/model_cache/indv_nll_nm_ex_te",format(Sys.time(), "%H_%M_%OS"),".rds")))










#### 11/15/23 
opt.m="Nelder-Mead"
plan(multisession)
RNGkind("L'Ecuyer-CMRG")

fit_params=list(pred_dat="test_dat", pred_fun="exam.response", loss_fun="nll", loss_data="test_error")
nm_ex_te <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(nm_ex_te,here::here(paste0("data/model_cache/indv_nll_nm_ex_te",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="test_error+train_error")
nm_ex_tetr <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(nm_ex_tetr,here::here(paste0("data/model_cache/indv_nll_nm_ex_tetr",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="train_error")
nm_ex_tr <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(nm_ex_tr,here::here(paste0("data/model_cache/indv_exam_nm_ex_tr",format(Sys.time(), "%H_%M_%OS"),".rds")))


plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll", loss_data="test_error")
nm_alm_te <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(nm_alm_te,here::here(paste0("data/model_cache/indv_nll_nm_alm_te",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll", loss_data="test_error+train_error")
nm_alm_tetr <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(nm_alm_tetr,here::here(paste0("data/model_cache/indv_nll_nm_alm_tetr",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll", loss_data="train_error")
nm_alm_tr <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, opt.m, initial_params))
saveRDS(nm_alm_tr,here::here(paste0("data/model_cache/indv_nll_nm_alm_tr",format(Sys.time(), "%H_%M_%OS"),".rds")))



####


opt.m="BFGS"


plan(multisession)
RNGkind("L'Ecuyer-CMRG")

fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="test_error")
bfgs_ex_te <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(bfgs_ex_te,here::here(paste0("data/model_cache/indv_nll_bfgs_ex_te",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="test_error+train_error")
bfgs_ex_tetr <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(bfgs_ex_tetr,here::here(paste0("data/model_cache/indv_nll_bfgs_ex_tetr",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="train_error")
bfgs_ex_tr <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(bfgs_ex_tr,here::here(paste0("data/model_cache/indv_exam_bfgs_ex_tr",format(Sys.time(), "%H_%M_%OS"),".rds")))


plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll", loss_data="test_error")
bfgs_alm_te <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(bfgs_alm_te,here::here(paste0("data/model_cache/indv_nll_bfgs_alm_te",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll", loss_data="test_error+train_error")
bfgs_alm_tetr <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(bfgs_alm_tetr,here::here(paste0("data/model_cache/indv_nll_bfgs_alm_tetr",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll", loss_data="train_error")
bfgs_alm_tr <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, opt.m, initial_params))
saveRDS(bfgs_alm_tr,here::here(paste0("data/model_cache/indv_nll_bfgs_alm_tr",format(Sys.time(), "%H_%M_%OS"),".rds")))




####### SANN


opt.m="SANN"


plan(multisession)
RNGkind("L'Ecuyer-CMRG")

fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="test_error")
sann_ex_te <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(sann_ex_te,here::here(paste0("data/model_cache/indv_nll_sann_ex_te",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="test_error+train_error")
sann_ex_tetr <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(sann_ex_tetr,here::here(paste0("data/model_cache/indv_nll_sann_ex_tetr",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="train_error")
sann_ex_tr <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(sann_ex_tr,here::here(paste0("data/model_cache/indv_exam_sann_ex_tr",format(Sys.time(), "%H_%M_%OS"),".rds")))


plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll", loss_data="test_error")
sann_alm_te <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(sann_alm_te,here::here(paste0("data/model_cache/indv_nll_sann_alm_te",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll", loss_data="test_error+train_error")
sann_alm_tetr <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,opt.m, initial_params))
saveRDS(sann_alm_tetr,here::here(paste0("data/model_cache/indv_nll_sann_alm_tetr",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll", loss_data="train_error")
sann_alm_tr <- future_map(split_data, ~ optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, opt.m, initial_params))
saveRDS(sann_alm_tr,here::here(paste0("data/model_cache/indv_nll_sann_alm_tr",format(Sys.time(), "%H_%M_%OS"),".rds")))


### de - no adjustment for error - 2K iterations

library(DEoptimR)

plan(multisession)
RNGkind("L'Ecuyer-CMRG")

fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="test_error")
de_ex_te <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_ex_te,here::here(paste0("data/model_cache/indv_nll_de2_ex_te",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="test_error+train_error")
de_ex_tetr <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_ex_tetr,here::here(paste0("data/model_cache/indv_nll_de2_ex_tetr",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="train_error")
de_ex_tr <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_ex_tr,here::here(paste0("data/model_cache/indv_exam_de2_ex_tr",format(Sys.time(), "%H_%M_%OS"),".rds")))


plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll", loss_data="test_error")
de_alm_te <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_alm_te,here::here(paste0("data/model_cache/indv_nll_de2_alm_te",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll", loss_data="test_error+train_error")
de_alm_tetr <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,initial_params))
saveRDS(de_alm_tetr,here::here(paste0("data/model_cache/indv_nll_de2_alm_tetr",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll", loss_data="train_error")
de_alm_tr <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_alm_tr,here::here(paste0("data/model_cache/indv_nll_de2_alm_tr",format(Sys.time(), "%H_%M_%OS"),".rds")))




### de - no adjustment for error - 3K iterations - small minimum for all 3 params, raised maximum

library(DEoptimR)

plan(multisession)

fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="test_error")
de_ex_te <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_ex_te,here::here(paste0("data/model_cache/indv_nll_de3_ex_te",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="test_error+train_error")
de_ex_tetr <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_ex_tetr,here::here(paste0("data/model_cache/indv_nll_de3_ex_tetr",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="exam.response", loss_fun="nll", loss_data="train_error")
de_ex_tr <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_ex_tr,here::here(paste0("data/model_cache/indv_exam_de3_ex_tr",format(Sys.time(), "%H_%M_%OS"),".rds")))


plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll", loss_data="test_error")
de_alm_te <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_alm_te,here::here(paste0("data/model_cache/indv_nll_de3_alm_te",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll", loss_data="test_error+train_error")
de_alm_tetr <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params,initial_params))
saveRDS(de_alm_tetr,here::here(paste0("data/model_cache/indv_nll_de3_alm_tetr",format(Sys.time(), "%H_%M_%OS"),".rds")))

plan(multisession)
fit_params=list(pred_dat="test_avg", pred_fun="alm.responseOnly", loss_fun="nll", loss_data="train_error")
de_alm_tr <- future_map(split_data, ~ de_optim_fit_id(.x, func="fit_lr_c0", fit_params, model_params, initial_params))
saveRDS(de_alm_tr,here::here(paste0("data/model_cache/indv_nll_de3_alm_tr",format(Sys.time(), "%H_%M_%OS"),".rds")))
