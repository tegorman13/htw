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














c_values <- seq(0.000001, 0.01, length.out=50)


lr_values <- seq(0.01, 1.5, length.out=30)





predParams=list(pred_dat="test_avg",pred_fun="exam.response",loss_fun=RMSE,loss_data="

test_error")



predParams=list(pred_dat="test_avg",pred_fun="exam.response",


loss_fun=RMSE,loss_data="test_error+train_error")





# predParams=list(pred_dat="rbind(test_avg,train_results$d[,c('x','y')])",pred_fun="ex

am.response",lossFun=RMSE)



# ex_te_v <- wrap_grid(vAvg, c_values, lr_values, input.layer, output.layer,predParams=e_testOnly)
# ex_tetr_v <- wrap_grid(vAvg, c_values, lr_values, input.layer, output.layer,predParams=e_testTrain)
# ex_tr_v<- wrap_grid(vAvg, c_values, lr_values, input.layer, output.layer,predParams=e_trainOnly)
# 
# a_te_v <- wrap_grid(vAvg, c_values, lr_values, input.layer, output.layer,predParams=a_testOnly)
# a_tetr_v <- wrap_grid(vAvg, c_values, lr_values, input.layer, output.layer,predParams=a_testTrain)
# a_tr_v <- wrap_grid(vAvg, c_values, lr_values, input.layer, output.layer,predParams=a_trainOnly)

#saveRDS(tibble::lst(ex_te_v,ex_tetr_v,ex_tr_v, a_te_v,a_tetr_v,a_tr_v), here::here("data/model_cache/var_group_exam_fits.rds"))


# ex0_te_c <- wrap_grid(cAvg, c_values, lr_values, input.layer, output.layer,predParams=e_testOnly)
# ex0_tetr_c <- wrap_grid(cAvg, c_values, lr_values, input.layer, output.layer,predParams=e_testTrain)
# ex0_tr_c <- wrap_grid(cAvg, c_values, lr_values, input.layer, output.layer,predParams=e_trainOnly)
# 
# a_te_c <- wrap_grid(cAvg, c_values, lr_values, input.layer, output.layer,predParams=a_testOnly)
# a_tetr_c <- wrap_grid(cAvg, c_values, lr_values, input.layer, output.layer,predParams=a_testTrain)
# a_tr_c <- wrap_grid(cAvg, c_values, lr_values, input.layer, output.layer,predParams=a_trainOnly)


# saveRDS(tibble::lst(ex0_te_c,ex0_tetr_c,ex0_tr_c, a_te_c,a_tetr_c,a_tr_c), here::here("data/model_cache/con_group_exam_fits.rds"))




best_paramsV <- wrap_grid(dsAvg |> filter(condit=="Varied"), c_values, lr_values, 


  input.layer, output.layer,


  predParams)





eval(parse(text="train_data[,c('x','y')]"))










