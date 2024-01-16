#
#
#

pacman::p_load(tidyverse, data.table, here, patchwork, conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)

walk(c("fun_alm","fun_model", "Display_Functions"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
ids <- c(1,2,4,5,6,7,8, 10,11,12,13)
ids2 <- c(1,66,141,142,36,57,81)
ids3 <- c(20,71,101,4,76,192)


#
#
#
#
#
ind_fits <- readRDS("~/Library/CloudStorage/GoogleDrive-tegorman13@gmail.com/My Drive/HTW/gl/data/indv_sim/ind_abc_ss_1000_ng100_bufp55.rds")
ind_fits <- ind_fits |> list_assign(runInfo = zap())

# names(ind_fits)
# names(ind_fits$id_fits_alm_varied)
# names(ind_fits$id_fits_alm_varied[[1]])

teter_results_df <- imap_dfr(ind_fits, ~imap_dfr(.x, ~.x[["teter_results"]] ,.id="id"))


combined_results_df <- imap_dfr(ind_fits, ~imap_dfr(.x, ~{
  teter_data <- .x[["teter_results"]] %>% mutate(result_type = "teter_results")
  te_data <- .x[["te_results"]] %>% mutate(result_type = "te_results")
  tr_data <- .x[["tr_results"]] %>% mutate(result_type = "tr_results")
  combined_data <- bind_rows(teter_data, te_data, tr_data)
}, .id = "id"))


id1 <- combined_results_df |> filter(id %in% ids2) 
names(id1)
typeof(id1$sim_dat)
str(id1$sim_dat[[1]])
# sim_dat is a list dataframe. Add the value of rank into sim_dat
id1 <- id1 |> rowwise() |> mutate(sim_dat = map2(sim_dat, rank, ~mutate(.x, rank = .y)))



id1 <- id1 %>%
  rowwise() %>%
  mutate(
    sim_dat = list(map2(
      sim_dat, rank, 
      ~ if (is.data.frame(.x)) {
        mutate(.x, rank = .y)
      } else {
        .x
      }
    ))
  ) %>%
  ungroup()

#
#
#
#
#

plt_c_gi <- combined_results_df |> 
  ggplot(aes(x=c)) + geom_density(aes(fill=Group), alpha=.5) + 
  facet_wrap(Fit_Method~Model,scales="free",ncol=2) + 
  theme_bw() + theme(legend.position = "bottom") + labs(title = "Posterior distribution of c") 

plt_lr_gi <- combined_results_df |> 
  ggplot(aes(x=lr)) + geom_density(aes(fill=Group), alpha=.5) + 
  facet_wrap(Fit_Method~Model,scales="free",ncol=2) + 
  theme_bw() + theme(legend.position = "bottom") + labs(title = "Posterior distribution of lr") 

plt_c_gi + plt_lr_gi

#
#
#
#
#

post_dat <- combined_results_df |> filter(id %in% (c(ids2,ids3)), rank<=1) |> 
  select(sim_dat) |> unnest(sim_dat) |> filter(expMode2=="Test") |>
  mutate(facet_label = paste0("Model: ", Model, "\n", "Fit: ", Fit_Method, "\n", "c: ", round(c, 4), "\n", "lr: ", round(lr, 4)))

post_dat |> filter(id %in% ids2) |> ggplot(aes(x = x, y = y, fill=condit)) + 
  stat_summary(fun=mean, geom="bar", position=position_dodge(), alpha=.75) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge()) +
  stat_halfeye(aes(x=x,y=pred,color=condit),position=position_dodge()) +
  ggh4x::facet_nested_wrap(id~facet_label, ncol=6, scales="free") + labs(title = "Average of posterior predictive distributions")


post_dat |> filter(id %in% ids3) |> ggplot(aes(x = x, y = y, fill=condit)) + 
  stat_summary(fun=mean, geom="bar", position=position_dodge(), alpha=.75) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge()) +
  stat_halfeye(aes(x=x,y=pred,color=condit),position=position_dodge()) +
  ggh4x::facet_nested_wrap(id~facet_label, ncol=6, scales="free") + labs(title = "Average of posterior predictive distributions")


#
#
#
#
#




post_ev <- ind_fits$id_fits_exam_varied |> 
  map_dfr(~tibble(pluck(.x$teter_results)),.id = "id") |> 
  filter(id %in% c(13)) |> 
  filter(rank==1) |> 
  select(sim_dat) |> unnest(sim_dat) 


post_ev |> filter(expMode2=="Test") |> 
  ggploat(aes(x = x, y = pred, fill=condit)) + 
  stat_summary(fun=mean, geom="bar", position=position_dodge(), alpha=.75) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge()) +
  #stat_summary(fun=mean, geom="point", aes(x=x,y=y,color=condit), position=position_dodge()) +
  stat_halfeye(aes(x=x,y=y,color=condit),position=position_dodge()) +
  ggh4x::facet_nested_wrap(~id) + labs(title = "Average of posterior predictive distributions")  



post_alm <- ind_fits$id_fits_alm_varied |> 
  map_dfr(~tibble(pluck(.x$teter_results)),.id = "id") |> 
  filter(id %in% c(1,2,4,5,6,7,8, 10,11,12,13)) |> 
  filter(rank<10) |> 
  select(sim_dat) |> unnest(sim_dat) 

post_alm |> filter(expMode2=="Test") |> 
  ggplot(aes(x = x, y = pred, fill=condit)) + 
  stat_summary(fun=mean, geom="bar", position=position_dodge(), alpha=.75) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge()) +
  #stat_summary(fun=mean, geom="point", aes(x=x,y=y,color=condit), position=position_dodge()) +
  stat_halfeye(aes(x=x,y=y,color=condit),position=position_dodge()) +
  ggh4x::facet_nested_wrap(~id) + labs(title = "Average of posterior predictive distributions")  



#
#
#
#
#
#
#
#| eval: false
# ind_abc_s500_ng300_buf5 <- readRDS("~/Library/CloudStorage/GoogleDrive-tegorman13@gmail.com/My Drive/HTW/gl/data/indv_sim/ind_abc_s500_ng300_buf5.rds")

ind_abc_ss_1000_13_56_52 <- readRDS("~/Library/CloudStorage/GoogleDrive-tegorman13@gmail.com/My Drive/HTW/gl/data/indv_sim/ind_abc_ss_1000_13_56_52.rds")


abc_ev <- ind_abc_ss_1000_13_56_52$id_fits_exam_varied |>
  map_dfr(~tibble(pluck(.x$teter_results)),.id = "id") |> filter(id %in% c(1,2)) |> select(sim_dat) |> unnest(sim_dat) |> filter(rank==1)



abc_ev <- ind_abc_s500_ng300_buf5$id_fits_exam_varied |>
  map_dfr(~tibble(pluck(.x$teter_results)),.id = "id") |> filter(id %in% c(1,2)) |>
  map(sim_dat)


abc_ev <- ind_abc_s500_ng300_buf5$id_fits_exam_varied |>
  map_dfr(~tibble(pluck(.x$teter_results)),.id = "id") |> filter(id %in% c(1,2)) |> select(sim_dat) |>
  pluck("sim_dat")


abc_ev <- ind_abc_s500_ng300_buf5$id_fits_exam_varied |>
  map_dfr(~tibble(pluck(.x$teter_results)),.id = "id") |> filter(id %in% c(1,2)) |> select(sim_dat) |> unnest(sim_dat) |> filter(rank==1) |>
  mutate(pred=pred[[1]])
  #mutate(pred=unnest(pred)$...195,resid=unnest(resid)$...195 )

abc_ev |> filter(expMode2=="Test") |> 
  ggplot(aes(x = x, y = pred, fill=condit)) + 
  stat_summary(fun=mean, geom="bar", position=position_dodge(), alpha=.75) +
  stat_summary(fun.data=mean_se, geom="errorbar", position=position_dodge()) +
  #stat_summary(fun=mean, geom="point", aes(x=x,y=y,color=condit), position=position_dodge()) +
  stat_halfeye(aes(x=x,y=y,color=condit),position=position_dodge()) +
  ggh4x::facet_nested_wrap(~id) + labs(title = "Average of posterior predictive distributions")  

colnames(abc_ev)
k=abc_ev["pred"]
names(k)


abc_ev <- id_fits_exam_varied |>
  map_dfr(~tibble(pluck(.x$teter_results)),.id = "id") |> filter(id %in% c(1,2)) |> select(sim_dat) |> unnest(sim_dat)


#
#
#
#
#
#
#
#
#
#
#
#
#
#
#


  sim_data_gen_s <- function(data, input_layer, output_layer, simulation_function, prior_samples, return_dat) {
    simulation_function <- match.fun(simulation_function)
    suppressMessages(simulation_results <- map_dfc(seq_len(nrow(prior_samples)), function(idx) {
      params <- prior_samples[idx, ]
      simulation_function(data=as.data.table(data), c=params$c, lr=params$lr, input_layer=input_layer, output_layer=output_layer, return_dat = return_dat)    
      }))
  }




input_layer =  c(100,350,600,800,1000,1200)
output_layer = input_layer
return_dat="train_data, test_data"
data <- ds |> filter(id==1) 
train_idx <- which(data$expMode2=="Train")
test_idx <- which(data$expMode2=="Test")

kde_results <- purrr::map(group_posterior_all, ~purrr::map_if(.x, is.list, ~compute_kde(.x$c, .x$lr, ngrid = 100, nsamples=100)))

sim_data <- sim_data_gen_s(data,input_layer,output_layer,full_sim_exam, kde_results$abc_ev$te_results$kde_samples, return_dat)

simulated = as_vector(sim_data[,1])
observed=data 


dist_mean_sd <- function(simulated, observed) {
  # Calculate the means and standard deviations for simulated and observed data
  
  d <- tibble(observed,pred=simulated) 
  nbins <- 4
  dtrain <- d |> filter(expMode2=="Train") |> 
  mutate(Block=cut(tr,breaks=seq(1,max(tr), length.out=nbins+1),include.lowest=TRUE,labels=FALSE))
 
  sum_stat_train <- dtrain |> group_by(Block,x) |> summarise(mean_sim=mean(pred),sd_sim=sd(pred),mean_obs=mean(y),sd_obs=sd(y), mean_dif = mean_obs-mean_sim, sd_dif = sd_obs-sd_sim) 

  sum_stat_test <- d |> filter(expMode2=="Test") |> group_by(x) |> summarise(mean_sim=mean(pred),sd_sim=sd(pred),mean_obs=mean(y),sd_obs=sd(y), mean_dif = mean_obs-mean_sim, sd_dif = sd_obs-sd_sim)

  train_mean_rmse <- sqrt(mean((sum_stat$mean_dif^2)))
  train_sd <- mean(sum_stat_train$sd_obs)

  test_mean_rmse <- sqrt(mean((sum_stat_test$mean_dif^2)))
  test_sd <- mean(sum_stat_test$sd_obs)
  
  return(tibble::lst(train_mean_rmse,test_mean_rmse, train_sd, test_sd))
  
}

#
#
#
#
#
#
#
#
id1 <- ds |> filter(id==1) 

input_layer =  c(100,350,600,800,1000,1200)
output_layer = input_layer
n_prior_samples=500
prior_samples <- kde_results$abc_ev$teter_results$kde_samples
return_dat="train_data, test_data"
Model="EXAM"; Group="Varied"
data=id1

train_idx <- which(id1$expMode2=="Train")
test_idx <- which(id1$expMode2=="Test")



sim_data_wrapper <- function(args) {
  sim_data_gen_s(args[[1]], args[[2]], args[[3]], simulation_function = args[[4]], args[[5]], return_dat = args[[6]])
}

sim_data_gen_s <- function(data, input_layer, output_layer, simulation_function, prior_samples, return_dat) {
  simulation_function <- match.fun(simulation_function)
  suppressMessages(simulation_results <- map_dfc(seq_len(nrow(prior_samples)), function(idx) {
    params <- prior_samples[idx, ]
    simulation_function(data=as.data.table(data), c=params$c, lr=params$lr, input_layer=input_layer, output_layer=output_layer, return_dat = return_dat)
  }))
}

full_sim_exam <- function(data, c, lr,pred_fun=exam.response, input_layer, output_layer,return_dat="test_data",mode="sim") {
  train_data <- data[expMode2=="Train", c("condit","tr","expMode2", "x","y")] 
  test_data <- data[expMode2=="Test", c("condit","tr","expMode2", "x","y")] 
  trainVec=sort(unique(train_data$x))
  if (train_data$condit[1] != "Varied") {
    trainVec <- c(0, trainVec)
  }
  train_results <- alm.sim(train_data, c, lr, input_layer, output_layer)
  test_prediction <- map_dbl(test_data$x, ~ pred_fun(.x, c, input_layer, 
                                                     output_layer, train_results$wm,  trainVec=trainVec))
  
  train_data$pred <- train_results$d$almResp
  test_data$pred <- test_prediction
  fd = eval(parse(text=paste0("rbind(",return_dat,")")))
  if(mode=="sim"){return(fd$pred)
  }else {return(fd)}
}

args_list <- list(
  Exam_Varied=list(id1, input_layer, output_layer, full_sim_exam, prior_samples, return_dat)
)

sim1 <- full_sim_exam(data=as.data.table(id1), c=.02, lr=2.9, input_layer=input_layer, output_layer=output_layer, return_dat = return_dat)
id1p<- id1 |> mutate(pred=sim1, resid=y-pred)



t1=system.time({
sim_dataAll <- map(args_list, sim_data_wrapper) |> as.data.table()
})
t1





target_data_train_test <- id1[expMode2 %in% c("Test", "Train"), ]$y
target_data_test <- id1[expMode2 == "Test", ]$y
target_data_train <- id1[expMode2 == "Train", ]$y

te_distances <- purrr::map_dbl(sim_dataAll[test_idx, ], dist_rmse, observed = target_data_test)


pct_keep=.1
te_results <- tibble(distance = te_distances, c = prior_samples$c, 
  lr = prior_samples$lr,  sim_index= seq_along(te_distances)) |> 
    arrange(distance) |> 
    filter(if (!is.null(pct_keep)) row_number() <= n() * pct_keep else distance <= tol) |>
    mutate(rank=row_number(),Model,Group,Fit_Method="Test Only") |>
    rowwise() |>
    mutate(sim_dat = list(tibble(Model,Fit_Method,c,lr,distance,rank,data,pred=(sim_dataAll[, .SD, .SDcols = sim_index])) |> 
    mutate(resid=y-pred)))

head(te_results)

id1_test <- sim_dataAll[test_idx, ]
id1_test[,1]


#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
