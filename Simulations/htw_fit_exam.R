pacman::p_load(dplyr,purrr,ggplot2,patchwork,here)
purrr::walk(here(c("Functions/Display_Functions.R", "Functions/alm_core.R","Functions/misc_model_funs.R")),source)
select <- dplyr::select; mutate <- dplyr::mutate 

ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))

dsAvg <- ds |> group_by(condit,expMode2,tr, x) |> 
  summarise(y=mean(y),.groups="keep") 

#head(dsAvg) |> pander::pandoc.table(style = "simple")

i1 <- ds |> filter(id=="1")

input.layer <- c(100,350,600,800,1000,1200)
output.layer <- c(100,350,600,800,1000,1200)

trainVec0 <- c(0,sort(unique(train_data$x)))
trainVec <- sort(unique(train_data$x))
testVec <- c(unique(test_data$x))


c_values <- seq(0.000001, 1.0, length.out=150)
lr_values <- seq(0.0000001, 4.0, length.out=200)


# default.layer <- c(100,350,600,800,1000,1200)
# output.layer <- adjust_layer(default.layer,k=3)


a_testOnly=list(pred_dat="test_avg",pred_fun="alm.responseOnly",loss_fun="RMSE",loss_data="test_error")
a_trainOnly=list(pred_dat="test_avg",pred_fun="alm.responseOnly",loss_fun="RMSE",loss_data="train_error")
a_testTrain=list(pred_dat="test_avg",pred_fun="alm.responseOnly",loss_fun="RMSE",loss_data="test_error+train_error")

e_testOnly=list(pred_dat="test_avg",pred_fun="exam.response",loss_fun="RMSE",loss_data="test_error")
e_trainOnly=list(pred_dat="test_avg",pred_fun="exam.response",loss_fun="RMSE",loss_data="train_error")
e_testTrain=list(pred_dat="test_avg",pred_fun="exam.response",loss_fun="RMSE",loss_data="test_error+train_error")


plan(multisession)

ex_te_v <- wrap_grid(vAvg, c_values, lr_values, input.layer, output.layer,predParams=e_testOnly)
ex_tetr_v <- wrap_grid(vAvg, c_values, lr_values, input.layer, output.layer,predParams=e_testTrain)
ex_tr_v<- wrap_grid(vAvg, c_values, lr_values, input.layer, output.layer,predParams=e_trainOnly)

a_te_v <- wrap_grid(vAvg, c_values, lr_values, input.layer, output.layer,predParams=a_testOnly)
a_tetr_v <- wrap_grid(vAvg, c_values, lr_values, input.layer, output.layer,predParams=a_testTrain)
a_tr_v <- wrap_grid(vAvg, c_values, lr_values, input.layer, output.layer,predParams=a_trainOnly)

saveRDS(tibble::lst(ex_te_v,ex_tetr_v,ex_tr_v, a_te_v,a_tetr_v,a_tr_v), here::here("data/model_cache/var_group_exam_fits.rds"))


ex0_te_c <- wrap_grid(cAvg, c_values, lr_values, input.layer, output.layer,predParams=e_testOnly)
ex0_tetr_c <- wrap_grid(cAvg, c_values, lr_values, input.layer, output.layer,predParams=e_testTrain)
ex0_tr_c <- wrap_grid(cAvg, c_values, lr_values, input.layer, output.layer,predParams=e_trainOnly)

a_te_c <- wrap_grid(cAvg, c_values, lr_values, input.layer, output.layer,predParams=a_testOnly)
a_tetr_c <- wrap_grid(cAvg, c_values, lr_values, input.layer, output.layer,predParams=a_testTrain)
a_tr_c <- wrap_grid(cAvg, c_values, lr_values, input.layer, output.layer,predParams=a_trainOnly)


saveRDS(tibble::lst(ex0_te_c,ex0_tetr_c,ex0_tr_c, a_te_c,a_tetr_c,a_tr_c), here::here("data/model_cache/con_group_exam_fits.rds"))












