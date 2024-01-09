#
#
#
#
#
#

pacman::p_load(tidyverse, data.table, here, patchwork, conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)


walk(c("fun_alm","fun_model", "Display_Functions"), ~ source(here::here(paste0("Functions/", .x, ".R"))))
ds <- readRDS(here::here("data/e1_md_11-06-23.rds"))  |> as.data.table()
dsv <- ds |> filter(condit=="Varied")  
dsc <- ds |> filter(condit=="Constant") 
ind_ds <- ds |> 
  filter(expMode2 == "Test") |> 
  group_by(id, condit, x, expMode2) |> 
  summarise(y = mean(y), .groups = "keep") 



tMax=84
avg_dsv <- ds |> filter(condit=="Varied",expMode2=="Train") |> group_by(tr) %>%
  mutate(bandInt2 = sample(rep(c(800, 1000, 1200), each = tMax / 3), tMax, replace = FALSE)[tr]) %>%
  filter(bandInt2 == x) |> select(-bandInt2) |> group_by(tr,condit,x,expMode2) |> summarise(y=mean(y),.groups="keep") |>
  rbind(dsv |> filter(expMode2=="Test") |> group_by(condit,x,expMode2) |> summarise(y=mean(y),tr=1,.groups="keep") ) |> setDT()

avg_dsc <- ds |> filter(condit=="Constant",expMode2=="Train",tr<=tMax) |> group_by(tr, condit,x,expMode2) |> 
  summarise(y=mean(y),.groups="keep") |> 
  rbind(dsc |> filter(expMode2=="Test") |> 
          group_by(condit,x,expMode2) |> 
          summarise(y=mean(y),tr=1,.groups="keep") ) |> 
  setDT()

# load group level abc posterior - posterior parameters and predicted data for each combo of condit, Fit_Method, and Model
group_posterior_all <- readRDS(here::here("data/abc_2M_rmse_p001.rds"))  

Model_Names <- list("EXAM", "Alt_EXAM", "ALM")
Fit_Methods <- list("Test & Train", "Test Only","Train Only")
condits <- list("Constant", "Varied")


# posteriors for each Fit_Method 
teter <- group_posterior_all |> map_dfr(~tibble(pluck(.x$teter_results))) 
te <- group_posterior_all |> map_dfr(~tibble(pluck(.x$te_results))) 
tr <- group_posterior_all |> map_dfr(~tibble(pluck(.x$tr_results))) 



#
#
#
#
#
#

sample_from_kde <- function(kde_result, n_samples) {
  # Flatten the z matrix and sample indices based on these probabilities
  sampled_indices <- sample(length(kde_result$z), size = n_samples, replace = TRUE, prob = c(kde_result$z))
  
  # Convert indices to matrix row and column numbers
  rows <- ((sampled_indices - 1) %% nrow(kde_result$z)) + 1
  cols <- ((sampled_indices - 1) %/% nrow(kde_result$z)) + 1
  
  # Map back to x and y values and create a tibble
  samples <- tibble(c = kde_result$x[rows], lr = kde_result$y[cols])
  return(samples)
}

compute_kde <- function(c, lr, ngrid = 100, nsamples = 500) {
  kde=MASS::kde2d(c, lr, n = ngrid)
  kde_samples=sample_from_kde(kde, nsamples)
  return(list(kde=kde, kde_samples=kde_samples))
}

kde_results <- purrr::map(group_posterior_all, ~purrr::map_if(.x, is.list, ~compute_kde(.x$c, .x$lr, ngrid = 100, nsamples=1000)))





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
